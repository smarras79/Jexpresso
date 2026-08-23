#=============================================================================
 ark.jl -- linearly-implicit additive Runge-Kutta (IMEX) integrators,
 registered as OrdinaryDiffEq algorithms.

 WHY NOT KenCarp4 / ARS343 OUT OF OrdinaryDiffEq
 -----------------------------------------------
 OrdinaryDiffEq's IMEX solvers assume the implicit part is nonlinear: they run
 a Newton iteration, form W = I - γΔt J, and hand W to a LinearSolve
 algorithm. Getting our column structure through that stack means teaching it
 about a Jacobian that is block-banded per column and distributed over MPI --
 fighting three layers of abstraction to arrive somewhere we can reach
 directly. Our implicit part is *linear by construction* (it is the acoustic
 operator linearised about the reference state), so each stage is one direct
 banded solve. No Newton, no convergence tuning, no W assembly.

 What we keep from OrdinaryDiffEq is everything else, and it is not small:
 `solve`, `init`/`step!`/`solve!`, the CallbackSet, tstops, saveat, the
 pre-compilation pass in TimeIntegrators.jl, and the restart/diagnostics/LES
 statistics callbacks all keep working untouched, because a HEVI run is still
 an ordinary `solve(prob, alg; ...)` -- only `alg` changes.

 THE CONTRACT WITH `p`
 ---------------------
 The algorithm reaches the implicit operator through the ODE parameter
 object, which in Jexpresso is `params`. `ark_hooks(alg, p)` resolves which
 field carries them -- `p.hevi` for HEVI_ARK, `p.imex` for IMEX_ARK -- once at
 `init`, and the stepper then only ever sees:

     hooks.fimp!(du, u, p, t)          f_imp evaluated at u
     hooks.solve!(dst, src, p, gdt)    dst = (I - gdt*A)^-1 (src - qref) + qref

 That is the whole interface, which is why the same stepper drives both the
 vertically-implicit split (whose solve is a banded LU per column) and the
 fully 3D one (whose solve is a preconditioned Krylov iteration). Anything
 supplying those two functions can be integrated by this file.

 `f_exp` is never written down: it is `f(u) - f_imp(u)`, evaluated as a
 subtraction. That is one extra (cheap) operator application per stage, in
 exchange for the guarantee that the two halves sum to the explicit RHS
 exactly -- no source term, boundary condition or closure can go missing in
 the split, whatever the deck turns on.

 CHOOSING A TABLEAU
 ------------------
 The instinct is to judge the explicit half on its own, by the imaginary-axis
 radius of its stability polynomial divided by the RHS evaluations per step.
 THAT IS THE WRONG NUMBER, and it is wrong in a way that is expensive to find
 out by running:

     scheme    order  RHS/step  imag. radius  radius/RHS   joint dt_max
     ARS111      1       1          0.00         0.00       unstable
     ARS121      1       2          1.00         0.50       0.038 s
     ARS222      2       2          0.00         0.00       unstable
     ARS232      2       3          1.73         0.58       0.052 s   <- default
     ARS343      3       4          2.83         0.71       0.0004 s
     ARS443      3       4          1.57         0.39       0.059 s

 The two right-hand columns rank the family DIFFERENTLY, and the second one is
 the one that decides whether a run survives. A HEVI acoustic mode with
 horizontal wavenumber kx and vertical wavenumber kz loads both halves at the
 same time -- zE = dt*c*kx on the explicit half, zI = dt*c*kz on the implicit
 one -- so what has to sit inside the unit disc is the joint amplification
 factor over the whole (zE, zI) rectangle the mesh can support, not the value
 on either axis. ark_joint_amplification measures it; ark_imaginary_radius
 measures the zI = 0 edge, which for HEVI is precisely the mode that does not
 matter.

 ARS343 has the best explicit radius in the family and is unusable here: once
 zI > 0, |R| > 1 for EVERY zE > 0. It was the default in this file until the
 joint region was measured, and the symptom was not obviously a tableau
 problem -- max|R| = 1.024 at dt = 0.03 is a 2.4% growth per step, slow enough
 that the run blew up at a fixed MODEL time (~15 s) almost independently of
 dt, and marginal enough that changing the rank count moved the blow-up by ten
 seconds. Both of those read as an MPI or formulation bug, not as a tableau.

 ARS222 is a second trap for a second reason: its explicit half is Heun's
 method, unstable everywhere on the imaginary axis. Fine for a
 diffusion-dominated test, wrong for advection.

 The joint numbers above are measured on the rtb_hevi mesh (explicit half
 26.6 /s, implicit half 69.5 /s) by ark_joint_dt_max, not quoted. Redo them
 for a new mesh: the ranking depends on the ratio of the two rates. The setup
 report does exactly this and refuses to run a combination that amplifies.

 For reference, CarpenterKennedy2N54 has radius 3.34 over 5 RHS and is stable
 at dt = 0.02 on that mesh; ARS232 at 0.052 s over 3 RHS is 3.5x its
 wall-clock.
=============================================================================#

# OrdinaryDiffEq v6.60+ split its internals into OrdinaryDiffEqCore. Resolve
# the module once so this file works against either layout.
const _ODEC = isdefined(OrdinaryDiffEq, :OrdinaryDiffEqCore) ?
              OrdinaryDiffEq.OrdinaryDiffEqCore : OrdinaryDiffEq

"""
    ARKTableau

An additive Runge-Kutta pair. `aE`/`bE`/`cE` drive the explicit half,
`aI`/`bI`/`cI` the implicit half; `aI` is lower triangular with a constant
diagonal `γ`, which is what lets one banded factorisation serve every stage
of every step.

`stiffly_accurate` records that both halves finish with `b == a[s, :]`, so
the last stage value *is* the step's answer. That saves one full RHS
evaluation per step, because the last stage's evaluation doubles as the next
step's first-same-as-last.
"""
struct ARKTableau
    name::Symbol
    s::Int
    order::Int
    γ::Float64
    aE::Matrix{Float64}
    aI::Matrix{Float64}
    bE::Vector{Float64}
    bI::Vector{Float64}
    cE::Vector{Float64}
    cI::Vector{Float64}
    stiffly_accurate::Bool
end

function _mk_tableau(name, order, aE, aI, bE, bI)
    s  = size(aE, 1)
    cE = vec(sum(aE, dims = 2))
    cI = vec(sum(aI, dims = 2))
    γ  = aI[s, s]
    sa = isapprox(bE, aE[s, :]; atol = 1e-14) &&
         isapprox(bI, aI[s, :]; atol = 1e-14) &&
         isapprox(cE[s], 1.0; atol = 1e-14)
    return ARKTableau(name, s, order, γ, aE, aI, bE, bI, cE, cI, sa)
end

"""
    ark_tableau(name::Symbol) -> ARKTableau

Ascher-Ruuth-Spiteri IMEX pairs (Appl. Numer. Math. 25, 1997). Naming follows
their ARS(s_explicit, s_implicit, order) convention.
"""
function ark_tableau(name::Symbol)
    if name === :ARS111
        aE = [0.0 0.0; 1.0 0.0]
        aI = [0.0 0.0; 0.0 1.0]
        return _mk_tableau(:ARS111, 1, aE, aI, [1.0, 0.0], [0.0, 1.0])

    elseif name === :ARS121
        aE = [0.0 0.0; 1.0 0.0]
        aI = [0.0 0.0; 0.0 1.0]
        return _mk_tableau(:ARS121, 1, aE, aI, [0.0, 1.0], [0.0, 1.0])

    elseif name === :ARS122
        aE = [0.0 0.0; 0.5 0.0]
        aI = [0.0 0.0; 0.0 0.5]
        return _mk_tableau(:ARS122, 2, aE, aI, [0.0, 1.0], [0.0, 1.0])

    elseif name === :ARS222
        γ = (2.0 - sqrt(2.0)) / 2.0
        δ = 1.0 - 1.0 / (2.0 * γ)
        aE = [0.0 0.0 0.0; γ 0.0 0.0; δ (1.0-δ) 0.0]
        aI = [0.0 0.0 0.0; 0.0 γ 0.0; 0.0 (1.0-γ) γ]
        return _mk_tableau(:ARS222, 2, aE, aI, [δ, 1.0-δ, 0.0], [0.0, 1.0-γ, γ])

    elseif name === :ARS232
        γ = (2.0 - sqrt(2.0)) / 2.0
        δ = -2.0 * sqrt(2.0) / 3.0
        aE = [0.0 0.0 0.0; γ 0.0 0.0; δ (1.0-δ) 0.0]
        aI = [0.0 0.0 0.0; 0.0 γ 0.0; 0.0 (1.0-γ) γ]
        return _mk_tableau(:ARS232, 2, aE, aI, [0.0, 1.0-γ, γ], [0.0, 1.0-γ, γ])

    elseif name === :ARS343
        γ  = 0.4358665215084590
        b1 = -1.5*γ*γ + 4.0*γ - 0.25
        b2 =  1.5*γ*γ - 5.0*γ + 1.25
        a31, a32 = 0.3212788860, 0.3966543747
        a41, a42, a43 = -0.105858296, 0.5529291479, 0.5529291479
        aE = [0.0 0.0 0.0 0.0
              γ   0.0 0.0 0.0
              a31 a32 0.0 0.0
              a41 a42 a43 0.0]
        aI = [0.0 0.0        0.0 0.0
              0.0 γ          0.0 0.0
              0.0 (1.0-γ)/2  γ   0.0
              0.0 b1         b2  γ]
        return _mk_tableau(:ARS343, 3, aE, aI, [0.0, b1, b2, γ], [0.0, b1, b2, γ])

    elseif name === :ARS443
        aE = [0.0     0.0     0.0  0.0  0.0
              0.5     0.0     0.0  0.0  0.0
              11/18   1/18    0.0  0.0  0.0
              5/6    -5/6     0.5  0.0  0.0
              0.25    1.75    0.75 -1.75 0.0]
        aI = [0.0  0.0   0.0   0.0  0.0
              0.0  0.5   0.0   0.0  0.0
              0.0  1/6   0.5   0.0  0.0
              0.0 -0.5   0.5   0.5  0.0
              0.0  1.5  -1.5   0.5  0.5]
        return _mk_tableau(:ARS443, 3, aE, aI,
                           [0.25, 1.75, 0.75, -1.75, 0.0], [0.0, 1.5, -1.5, 0.5, 0.5])
    end
    error("unknown ARK tableau :$name -- available: :ARS111 :ARS121 :ARS122 ",
          ":ARS222 :ARS232 :ARS343 :ARS443")
end

"""
    ark_stability_polynomial(tab) -> Vector{Float64}

Coefficients of the explicit half's stability function, `R(z) = Σ p[k] z^(k-1)`.
Used by `ark_imaginary_radius` and by the setup report, so the deck's choice
of tableau can be judged against the Δt it is being run at rather than
against a remembered rule of thumb.
"""
function ark_stability_polynomial(tab::ARKTableau)
    s = tab.s
    p = [1.0]
    v = ones(Float64, s)
    for _ = 1:s
        push!(p, dot(tab.bE, v))
        v = tab.aE * v
    end
    return p
end

"""
    ark_joint_amplification(tab, zEmax, zImax; n = 321) -> Float64

The IMEX figure of merit: `max |R(i·zE, i·zI)|` over the rectangle
`[0, zEmax] x [0, zImax]` of the imaginary axes.

WHY `ark_imaginary_radius` IS THE WRONG NUMBER FOR A SPLIT SCHEME
----------------------------------------------------------------
The explicit half's imaginary radius is what an explicit scheme's step size is
judged by, and it is silently misleading for an ARK pair. A HEVI acoustic mode
with horizontal wavenumber kx and vertical kz puts `zE = Δt·c·kx` on the
explicit half and `zI = Δt·c·kz` on the implicit one AT THE SAME TIME, and the
joint amplification factor is not the explicit one:

    ARS343   zI = 0:  zE stable to 2.828      <- the number that reads "best"
             zI > 0:  |R| > 1 for EVERY zE > 0

ARS343 is the best of the ARS family on the explicit axis and unusable for
HEVI. The mesh supports every (kx, kz) pair, so the whole rectangle has to be
inside the unit disc, not just its two edges.

The failure this catches is mild and slow, which is what makes it worth
automating: on the rtb_hevi case at Δt = 0.03, ARS343 gives max|R| = 1.024,
a 2.4% growth PER STEP that takes ~17 s of model time to reach amplitudes the
equation of state chokes on. It presents as a blow-up at a fixed physical
time, nearly independent of Δt, and it is marginal enough that round-off --
the rank count -- moves the blow-up by ten seconds. Nothing about that reads
as "wrong tableau" from the output.

The rectangle is sampled rather than solved because the stable set is not
simply connected in `zE` (ARS343 is stable, then unstable, then stable again),
so bisecting for a boundary reports a number that is not a limit.
"""
function ark_joint_amplification(tab::ARKTableau, zEmax::Real, zImax::Real; n::Int = 321)
    s = tab.s
    r = 0.0
    for e in range(0.0, Float64(zEmax); length = n), i in range(0.0, Float64(zImax); length = n)
        zE = im * e
        zI = im * i
        Y  = (I(s) - zE * tab.aE - zI * tab.aI) \ ones(ComplexF64, s)
        r  = max(r, abs(1 + zE * sum(tab.bE .* Y) + zI * sum(tab.bI .* Y)))
    end
    return r
end

"""
    ark_joint_dt_max(tab, rate_exp, rate_imp; rtol = 1e-4, hi = nothing) -> Float64

Largest Δt for which `ark_joint_amplification(tab, Δt·rate_exp, Δt·rate_imp)`
stays at 1. `rate_exp` and `rate_imp` are the largest eigenvalue magnitudes
(1/s) the explicit and implicit halves see on this mesh.
"""
function ark_joint_dt_max(tab::ARKTableau, rate_exp::Real, rate_imp::Real;
                          rtol::Real = 1.0e-4, hi::Union{Nothing,Real} = nothing,
                          n::Int = 161)
    ok(dt) = ark_joint_amplification(tab, dt * rate_exp, dt * rate_imp; n = n) <= 1 + 1.0e-9
    lo = 1.0e-8
    ok(lo) || return 0.0

    # BRACKET FIRST, then bisect. This used to take a fixed `hi = 0.5` and
    # bisect inside it, which silently returns the CAP rather than the limit
    # whenever the true limit is above it -- and 0.5 is not a large number
    # here: on hevi_10x1x50 the answer is 0.357, well inside it, but the same
    # tableau at rate_exp = 1 returns 0.4999 for both ARS232 and ARS443, which
    # is the cap and not a property of either method. A capped value is
    # indistinguishable from a real one in the setup report, and it errs
    # towards saying a Δt is safe.
    h = hi === nothing ? 1.0e-7 : Float64(hi)
    if hi === nothing
        while ok(h)
            h *= 2.0
            h > 1.0e8 && return h    # unconditionally stable in this direction
        end
        lo = 0.5 * h
    end

    # Relative tolerance: an absolute one is meaningless across limits that
    # span 1e-4 to 1e2 in the tableaux this is called on.
    while h - lo > rtol * max(lo, eps())
        m = 0.5 * (lo + h)
        ok(m) ? (lo = m) : (h = m)
    end
    return lo
end

"""
    ark_wedge_amplification(tab, zEmax, zImax, mach; n = 241) -> (amp, zE, zI)

The IMEX figure of merit for a split whose implicit half takes **all** of the
acoustics: `max |R(i·zE, i·zI)|` over the WEDGE

    { (zE, zI) : 0 <= zI <= zImax,  0 <= zE <= min(mach·zI, zEmax) }
    ∪  the zI = 0 edge, 0 <= zE <= zEmax

rather than over the full rectangle. Returns the amplification and where it
was attained.

WHY THE RECTANGLE IS THE WRONG SET ONCE THE SPLIT IS FULLY THREE-DIMENSIONAL
---------------------------------------------------------------------------
`ark_joint_amplification` scans the whole rectangle `[0,zEmax] x [0,zImax]`,
and for HEVI that is exactly right: there `zE` is the HORIZONTAL acoustic
rate `Δt·c·kx` and `zI` the VERTICAL one `Δt·c·kz`, the mesh supports every
`(kx, kz)` pair independently, so every corner of the rectangle is a mode the
run actually has.

Move all three directions of the acoustics to the implicit half and that stops
being true. What is left explicit is advection, `zE = Δt·|v|·k`, and the SAME
wavenumber `k` puts `zI = Δt·c·k` on the implicit half at the same time. The
two are locked together by the Mach number:

    zE / zI = |v| / c

so the reachable set is a wedge of opening `mach`, not a rectangle. Its one
detached piece is the `zI = 0` edge -- the entropy and vorticity modes, which
are advected and carry no acoustic partner.

This is not a technicality; it reverses the ranking of the family. Measured
here (see `test/imex3d/test_wedge_stability.jl`, which recomputes them):

    tableau   RHS/step   rectangle, sup over zI      wedge at mach <= 0.3
    ARS232        3      neutral                     zE <= 1.15
    ARS443        4      neutral                     zE <= 1.57
    ARS343        4      amplifies for zE > 0.07     zE <= 2.83

Per RHS evaluation -- which is what converts into wall-clock -- that is 0.38,
0.39 and 0.71. ARS343 is not marginally better here, it is nearly twice the
step per unit work, and its wedge limit is its full explicit imaginary radius:
at these Mach numbers the acoustic coupling costs it nothing at all.

ARS343 is unusable for HEVI and the best of the family here, for the same
reason in both cases: its unstable region sits at `zI` of order 1 with `zE`
of order 1 too. HEVI puts modes there. A fully 3D split cannot -- reaching
`zE = 1` costs `zI = 1/mach >= 3`, and by then the L-stable implicit half has
damped the mode by an order of magnitude.

`mach` should be an upper bound on `|v|/c` over the run, not its initial
value: it is the ONE input here that grows as the flow develops.
"""
function ark_wedge_amplification(tab::ARKTableau, zEmax::Real, zImax::Real,
                                 mach::Real; n::Int = 241)
    s   = tab.s
    zEm = Float64(zEmax); zIm = Float64(zImax); M = Float64(mach)
    amp = 0.0; aE = 0.0; aI = 0.0

    @inline function _amp(e, i)
        zE = im * e
        zI = im * i
        Y  = (I(s) - zE * tab.aE - zI * tab.aI) \ ones(ComplexF64, s)
        return abs(1 + zE * sum(tab.bE .* Y) + zI * sum(tab.bI .* Y))
    end

    for i in range(0.0, zIm; length = n)
        emax = min(M * i, zEm)
        # `emax` shrinks to zero at the wedge apex, so a fixed sample count
        # there would spend the whole budget resolving nothing. Scale it, but
        # keep at least two points so the edge itself is always sampled.
        ne = max(2, ceil(Int, n * emax / max(zEm, eps())))
        for e in range(0.0, emax; length = ne)
            r = _amp(e, i)
            r > amp && (amp = r; aE = e; aI = i)
        end
    end
    # The advected, non-acoustic modes: zI = 0 is not in the wedge (it is its
    # apex) but the modes are real and this is where the explicit half stands
    # alone, so it is the edge that sets zEmax in the first place.
    for e in range(0.0, zEm; length = n)
        r = _amp(e, 0.0)
        r > amp && (amp = r; aE = e; aI = 0.0)
    end
    return (amp, aE, aI)
end

"""
    ark_wedge_dt_max(tab, rate_exp, rate_imp, mach; rtol = 1e-4) -> Float64

Largest Δt for which `ark_wedge_amplification(tab, Δt·rate_exp, Δt·rate_imp,
mach)` stays at 1. The wedge analogue of `ark_joint_dt_max`, and the number
the 3D IMEX setup report quotes.

Brackets before bisecting, for the reason spelled out in `ark_joint_dt_max`:
a fixed upper cap returns the cap rather than the limit whenever the limit is
above it, and a capped value is indistinguishable from a real one in a report
that a deck's Δt is then chosen against.
"""
function ark_wedge_dt_max(tab::ARKTableau, rate_exp::Real, rate_imp::Real, mach::Real;
                          rtol::Real = 1.0e-4, n::Int = 121)
    ok(dt) = ark_wedge_amplification(tab, dt * rate_exp, dt * rate_imp, mach; n = n)[1] <=
             1 + 1.0e-9
    ok(1.0e-8) || return 0.0

    h = 1.0e-7
    while ok(h)
        h *= 2.0
        h > 1.0e8 && return h        # unconditionally stable in this direction
    end
    lo = 0.5 * h

    while h - lo > rtol * max(lo, eps())
        m = 0.5 * (lo + h)
        ok(m) ? (lo = m) : (h = m)
    end
    return lo
end

"""
    ark_imaginary_radius(tab; tol = 1e-4) -> Float64

Largest `y` for which `|R(iy)| <= 1`: the pure-advection stability limit of
the explicit half, in units of `Δt·λ`.

This is the limit only when the implicit half sees nothing (`zI = 0`), which
for HEVI is exactly the mode that does not matter. Use
`ark_joint_amplification` to judge a tableau; this one is reported alongside
it because the difference between the two is the point.
"""
function ark_imaginary_radius(tab::ARKTableau; tol::Real = 1.0e-4)
    p = ark_stability_polynomial(tab)
    R(y) = abs(sum(c * (im*y)^(k-1) for (k, c) in enumerate(p)))
    R(tol) > 1.0 + 1e-12 && return 0.0
    lo, hi = 0.0, 1.0
    while R(hi) <= 1.0 && hi < 64.0
        lo = hi; hi *= 2.0
    end
    while hi - lo > tol
        mid = 0.5 * (lo + hi)
        R(mid) <= 1.0 ? (lo = mid) : (hi = mid)
    end
    return lo
end

#-----------------------------------------------------------------------------
# The OrdinaryDiffEq algorithm
#-----------------------------------------------------------------------------

"""
    HEVI_ARK(name::Symbol = :ARS232)

Horizontally-explicit / vertically-implicit additive Runge-Kutta.

Drop it into a deck exactly where an explicit integrator went:

    :ode_solver => HEVI_ARK(:ARS232),
    :Δt         => 0.05,

Fixed step only -- the tableaux carry no embedded error estimate, and an
adaptive controller on a turbulence run would chase the SGS noise anyway.
`:ode_adaptive_solver` already defaults to `false`.

WHICH TABLEAU
-------------
Judge one by `ark_joint_amplification`, never by the explicit half's imaginary
radius: a HEVI acoustic mode loads BOTH halves at once, and the two figures of
merit rank the family differently. Measured on the rtb_hevi mesh, where the
explicit half sees 26.6 /s and the implicit half 69.5 /s:

    tableau   RHS/step   explicit radius   joint Δt_max   wall-clock to t=700
    ARS232        3           1.732           0.052 s        1848 s  (measured)
    ARS443        4           1.570           0.059 s       ~2200 s  (projected)
    ARS343        4           2.828           0.0004 s       unusable
    (CarpenterKennedy2N54, fully explicit: 5 RHS, 0.02 s, 6475 s measured)

The two measured figures are single-rank runs of problems/CompEuler/rtb_hevi on
one machine: 0.185 s/step at Δt = 0.02 explicit against 0.132 s/step at
Δt = 0.05 under ARS232, so 35000 steps against 14000 -- 3.50x. ARS443's is
scaled from its RHS count and was not run.

ARS343 has the largest explicit radius in the family and is unusable here:
`|R| > 1` for every `zE > 0` once `zI > 0`. It was this package's default until
that was measured, and the symptom was a blow-up at a fixed model time (~15 s)
that barely moved with Δt and shifted with the rank count -- a 2.4% per-step
growth is slow enough to look like anything but a tableau choice.

ARS232 is the default: second order, three RHS evaluations, and the best
wall-clock of the family on an anisotropic mesh. Use :ARS443 when third-order
accuracy in time is wanted; it costs about 20% more for the same physics.
"""
struct HEVI_ARK <: _ODEC.OrdinaryDiffEqAlgorithm
    tab::ARKTableau
end
# ONE constructor, deliberately. A keyword form alongside this --
#
#     HEVI_ARK(; tableau::Symbol = :ARS232) = ...
#
# -- looks like a harmless convenience and is not: keyword arguments do not
# participate in dispatch, so both it and the default-argument form below
# define the SAME zero-argument positional method `HEVI_ARK()`, and the second
# silently replaces the first. Julia does not warn about that by default
# (--warn-overwrite=no), so it survives every include-based test and then
# fails the package outright:
#
#     ERROR: Method overwriting is not permitted during Module precompilation.
#
# Same trap as the duplicate assemble_mpi! documented in
# src/kernel/mpi/mpi_communications.jl. test/hevi/test_load.jl guards it.
HEVI_ARK(name::Symbol = :ARS232) = HEVI_ARK(ark_tableau(name))

Base.show(io::IO, alg::HEVI_ARK) = print(io, "HEVI_ARK(:", alg.tab.name, ")")

"""
    IMEX_ARK(name::Symbol = :ARS343)

Fully three-dimensional implicit-explicit additive Runge-Kutta: the **whole**
linear acoustic-gravity subsystem is implicit, in all three directions, and
what is left explicit is advection, diffusion, the sources and the closures.

Same place in a deck as any other integrator:

    :ode_solver => IMEX_ARK(:ARS343),
    :Δt         => 1.5,                 # limited by ADVECTION now, not sound

It runs the same stepper as `HEVI_ARK` -- same tableaux, same
`fimp!` / `solve!` contract, same FSAL and stiff-accuracy handling. The one
difference is where the stage solve goes: HEVI's implicit operator is
column-local and is inverted by a banded LU per column, while this one couples
the whole domain and is inverted by a preconditioned Krylov iteration, with
that very column solve as the preconditioner. See `imex3d.jl`.

WHICH TABLEAU, AND WHY IT IS THE OPPOSITE ANSWER TO HEVI'S
----------------------------------------------------------
`:ARS343`, and the reason is `ark_wedge_amplification`. HEVI leaves the
horizontal acoustics explicit, so the explicit and implicit halves are loaded
by two INDEPENDENT wavenumbers and the whole `(zE, zI)` rectangle is
reachable; ARS343 amplifies there for any `zE > 0.07` and is unusable. Here
the acoustics are entirely implicit, the explicit half sees only advection,
and the same wavenumber sets both -- so the reachable set is the wedge
`zE <= mach·zI` and ARS343 is neutral on it out to `zE = 2.83` -- its full
explicit imaginary radius, i.e. the acoustic coupling costs it nothing at all
-- against ARS443's 1.57 and ARS232's 1.15. Over four RHS evaluations rather
than three, that is still nearly twice the step per unit work.

Fixed step only, as `HEVI_ARK`: the tableaux carry no embedded error estimate,
and `rhs!` contains MPI collectives that an adaptive controller would reach a
different number of times on each rank.
"""
struct IMEX_ARK <: _ODEC.OrdinaryDiffEqAlgorithm
    tab::ARKTableau
end
# ONE constructor here too -- see the comment above HEVI_ARK's.
IMEX_ARK(name::Symbol = :ARS343) = IMEX_ARK(ark_tableau(name))

Base.show(io::IO, alg::IMEX_ARK) = print(io, "IMEX_ARK(:", alg.tab.name, ")")

"""
    ARK_ALG

The two algorithms that share this file's stepper. They differ only in which
object on `p` carries the `fimp!` / `solve!` hooks, which `ark_hooks` resolves
once per run rather than per step.
"""
const ARK_ALG = Union{HEVI_ARK, IMEX_ARK}

"""
    ark_hooks(alg, p)

The implicit-half hook object for this algorithm, off the ODE parameter
object. Resolved in `alg_cache` and stored, so `perform_step!` neither
searches `p` nor hard-codes a field name -- and so a deck that names the
integrator without `params_setup` having built the matching cache fails at
`init`, with a message, instead of on the first stage.
"""
ark_hooks(::HEVI_ARK, p) = hasproperty(p, :hevi) ? p.hevi :
    error("HEVI_ARK was selected but `params` carries no `hevi` cache. ",
          "params_setup builds it only when hevi_enabled(inputs) is true, so put ",
          ":ode_solver => HEVI_ARK(...) in user_inputs.jl rather than swapping the ",
          "integrator in afterwards.")
ark_hooks(::IMEX_ARK, p) = hasproperty(p, :imex) ? p.imex :
    error("IMEX_ARK was selected but `params` carries no `imex` cache. ",
          "params_setup builds it only when imex3d_enabled(inputs) is true, so put ",
          ":ode_solver => IMEX_ARK(...) in user_inputs.jl rather than swapping the ",
          "integrator in afterwards.")

"""
    ark_relinearize!(params, hooks, u)

Refresh the implicit operator's coefficients from the state `u` and refactorise
whatever the scheme factorises. The `:PS` half of the linearisation switch;
`:RS` never calls it. One method per hook type, in the file that defines it.
"""
function ark_relinearize! end

_ODEC.alg_order(alg::ARK_ALG)   = alg.tab.order
_ODEC.isfsal(::ARK_ALG)         = true
_ODEC.isadaptive(::ARK_ALG)     = false
_ODEC.alg_extrapolates(::ARK_ALG) = false

mutable struct HEVI_ARK_Cache{uType, rateType, TabType, HookType} <: _ODEC.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    tab::TabType
    fE::Vector{rateType}
    fI::Vector{rateType}
    U::uType
    tmp::uType
    ftmp::rateType
    fsalfirst::rateType
    k::rateType
    # The implicit half's hooks (`fimp!`, `solve!`, the linearisation policy),
    # resolved from `p` once at `init` by `ark_hooks`. Held here rather than
    # looked up per step so the stepper never names a field of `params` and
    # the same `perform_step!` serves HEVI and the 3D IMEX.
    hv::HookType
end

_ODEC.get_fsalfirstlast(cache::HEVI_ARK_Cache, u) = (cache.fsalfirst, cache.k)
_ODEC.full_cache(c::HEVI_ARK_Cache) = (c.u, c.uprev, c.U, c.tmp, c.ftmp, c.fsalfirst, c.k,
                                       c.fE..., c.fI...)

# Trailing `args...` absorbs the `verbose` argument that OrdinaryDiffEq v7
# added to this signature; v6 calls it without.
function _ODEC.alg_cache(alg::ARK_ALG, u, rate_prototype, ::Type{uEltypeNoUnits},
                         ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                         uprev, uprev2, f, t, dt, reltol, p, calck,
                         ::Val{true}, args...) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    s = alg.tab.s
    HEVI_ARK_Cache(u, uprev, alg.tab,
                   [zero(rate_prototype) for _ = 1:s],
                   [zero(rate_prototype) for _ = 1:s],
                   zero(u), zero(u),
                   zero(rate_prototype), zero(rate_prototype), zero(rate_prototype),
                   ark_hooks(alg, p))
end

function _ODEC.alg_cache(alg::ARK_ALG, u, rate_prototype, ::Type{uEltypeNoUnits},
                         ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                         uprev, uprev2, f, t, dt, reltol, p, calck,
                         ::Val{false}, args...) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    error(typeof(alg), " is in-place only. The state here is a large distributed ",
          "field; an out-of-place formulation would allocate a full copy per stage.")
end

function _ODEC.initialize!(integrator, cache::HEVI_ARK_Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.fsalfirst
    integrator.k[2] = cache.k
    hevi_trace("initialize!: first rhs! evaluation (this is where the JIT of the ",
               "whole RHS chain happens)")
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    hevi_trace("initialize!: first rhs! done")
    integrator.stats.nf += 1
    return nothing
end

function _ODEC.perform_step!(integrator, cache::HEVI_ARK_Cache, repeat_step = false)

    t     = integrator.t
    dt    = integrator.dt
    uprev = integrator.uprev
    u     = integrator.u
    f     = integrator.f
    p     = integrator.p
    tab   = cache.tab
    s     = tab.s
    hv    = cache.hv

    HEVI_STEP_COUNT[] += 1

    # LHEVI-PS: refresh the operator's coefficients from the previous solution
    # and refactorise, every `update_freq` steps. Done HERE, before the stage
    # loop, so the whole step sees one consistent operator -- refreshing
    # mid-step would put different stages on different implicit operators and
    # break the ARK order conditions.
    #
    # Step 1 is skipped: the operator was just built from qe at setup, and the
    # state has not moved yet.
    if hv.linearization === :PS && HEVI_STEP_COUNT[] > 1 &&
       (HEVI_STEP_COUNT[] - 1) % hv.update_freq == 0
        ark_relinearize!(p, hv, uprev)
    end

    tracing = HEVI_TRACE[] && HEVI_STEP_COUNT[] <= HEVI_TRACE_STEPS[]
    tracing && hevi_trace("step ", HEVI_STEP_COUNT[], " begin: t=", t, " dt=", dt)
    _tstep = hevi_tic()

    #--- stage 1: U1 = uprev, and f(U1) is the FSAL from the previous step ----
    _t = hevi_tic()
    hv.fimp!(cache.fI[1], uprev, p, t + tab.cI[1] * dt)
    if hevi_prof_on()
        HEVI_PROFILE.t_fimp += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_fimp += 1
    end
    @inbounds @. cache.fE[1] = cache.fsalfirst - cache.fI[1]
    tracing && hevi_trace("  stage 1 done (f_imp on uprev)")

    #--- stages 2..s ----------------------------------------------------------
    for i = 2:s
        tmp = cache.tmp
        copyto!(tmp, uprev)
        @inbounds for j = 1:i-1
            aE = tab.aE[i, j]
            aI = tab.aI[i, j]
            aE != 0.0 && (@. tmp += (dt * aE) * cache.fE[j])
            aI != 0.0 && (@. tmp += (dt * aI) * cache.fI[j])
        end

        γ  = tab.aI[i, i]
        Ui = cache.U
        if γ != 0.0
            # (I - γΔt A)(U - qref) = tmp - qref. HEVI solves it column by
            # column with a banded LU; the 3D IMEX solves it with a
            # preconditioned Krylov iteration over the whole domain.
            tracing && hevi_trace("  stage ", i, ": entering implicit solve, gdt=", dt * γ)
            _t = hevi_tic()
            hv.solve!(Ui, tmp, p, dt * γ)
            if hevi_prof_on()
                HEVI_PROFILE.t_solve += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_solve += 1
            end
            tracing && hevi_trace("  stage ", i, ": implicit solve done")
        else
            copyto!(Ui, tmp)
        end

        # The final stage of a stiffly-accurate pair IS the step's answer, so
        # its RHS evaluation is written straight into the FSAL slot instead of
        # being recomputed after the update.
        last_is_answer = tab.stiffly_accurate && i == s
        dest = last_is_answer ? cache.k : cache.ftmp
        tracing && hevi_trace("  stage ", i, ": entering rhs!")
        _t = hevi_tic()
        f(dest, Ui, p, t + tab.cE[i] * dt)
        if hevi_prof_on()
            HEVI_PROFILE.t_rhs += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_rhs += 1
        end
        tracing && hevi_trace("  stage ", i, ": rhs! done")
        integrator.stats.nf += 1
        _t = hevi_tic()
        hv.fimp!(cache.fI[i], Ui, p, t + tab.cI[i] * dt)
        if hevi_prof_on()
            HEVI_PROFILE.t_fimp += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_fimp += 1
        end
        @inbounds @. cache.fE[i] = dest - cache.fI[i]

        last_is_answer && copyto!(u, Ui)
    end

    #--- combine ---------------------------------------------------------------
    if !tab.stiffly_accurate
        copyto!(u, uprev)
        @inbounds for i = 1:s
            bE = tab.bE[i]
            bI = tab.bI[i]
            bE != 0.0 && (@. u += (dt * bE) * cache.fE[i])
            bI != 0.0 && (@. u += (dt * bI) * cache.fI[i])
        end
        tracing && hevi_trace("  final rhs! (FSAL for the next step)")
        _t = hevi_tic()
        f(cache.k, u, p, t + dt)
        if hevi_prof_on()
            HEVI_PROFILE.t_rhs += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_rhs += 1
        end
        integrator.stats.nf += 1
    end

    tracing && hevi_trace("step ", HEVI_STEP_COUNT[], " end")
    HEVI_MONITOR[] && hevi_monitor_step!(u, p)
    if hevi_prof_on()
        HEVI_PROFILE.t_step += (time_ns() - _tstep) * 1e-9; HEVI_PROFILE.n_step += 1
        HEVI_PROFILE.n_step >= HEVI_PROF_EVERY[] && hevi_profile_report!()
    end
    return nothing
end
