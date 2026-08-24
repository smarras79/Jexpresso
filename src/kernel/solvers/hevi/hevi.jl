#=============================================================================
 hevi.jl -- setup, self-check and reporting for the HEVI time integrator.

 Include order inside Jexpresso: columns.jl, operator.jl, factorize.jl,
 ark.jl, then this file.

 WHAT A DECK TURNS ON
 --------------------
     :ode_solver => HEVI_ARK(:ARS232),      # instead of CarpenterKennedy2N54()
     :Δt         => 0.02,                   # measure it; see below

 Optional:
     :hevi_verify        => true    # setup self-check (default true, cheap)
     :hevi_wall_flux     => true    # zero implicit vertical mass flux at
                                    # floor and lid (default true)
     :hevi_vars          => [1,4,5] # override the implicit variable set
     :lcfl_report        => true    # print the stability table at startup

 Everything else in the deck is untouched. The RHS, the SGS closure, the wall
 model, the sponge, the boundary conditions, the diagnostics, restart and LES
 statistics all keep working exactly as they do under an explicit solver,
 because f_exp is computed as rhs! - f_imp rather than re-derived.

 WHAT TO EXPECT FROM IT
 ----------------------
 HEVI removes the *vertical acoustic* term from the explicit eigenvalue
 budget and nothing else, so the gain is bounded by the grid's acoustic
 anisotropy. For the LESICP2 target mesh -- 128x128x120 over
 10240 x 10240 x 5000 m, nop 4 in every direction -- that is

     dx_element = 80 m,  dz_element = 41.7 m,  ratio 1.92

 Removing the vertical acoustic term halves the acoustic rate -- and it is NOT
 the 10-30x the same split buys on a mesh with dx >> dz.

 HALVING THE RATE DOES NOT DOUBLE Δt. The default tableau gives most of it
 back: ARS232's explicit imaginary radius is 1.732 against
 CarpenterKennedy2N54's 3.34, so measured on this mesh Δt goes from 0.0229 s
 explicit to 0.0202 s under HEVI -- about 0.9x, not 2x. (This file used to say
 "about two"; that figure was computed for ARS343, radius 2.83, before the
 joint-stability analysis ruled ARS343 out for HEVI and made ARS232 the
 default.) What is left is the cheaper STEP -- 3 rhs! plus 3 operator
 applications and 2 column solves against 5 rhs! -- which comes to roughly
 1.0-1.35x end to end depending on how expensive rhs! is relative to the
 vertical operator. On a case with a cheap rhs! HEVI is SLOWER than explicit.
 See README.md, "What it is worth on the LESICP2 target mesh", for the
 measured table, and use JEXPRESSO_HEVI_PROFILE=1 to place your own case on it.

 Getting past that means also removing the *horizontal* acoustic term:
 acoustic substepping (substep.jl, implemented but not yet delivering its
 intended outer step) or the fully 3D implicit solve (imex3d.jl, which reuses
 the column solver built here as its preconditioner).

 Run with :lcfl_report => true first. If the report says SGS diffusion rather
 than vertical acoustics is the binding term, HEVI as configured here will
 buy almost nothing and the vertical diffusion is what wants to be implicit.
=============================================================================#

"""
    HEVICache

Everything the HEVI stage solve needs, assembled once at setup and reached
from `params.hevi` by the integrator.

`fimp!` and `solve!` are the two hooks `HEVI_ARK` calls; they are stored as
fields (rather than looked up by name) so the integrator stays independent of
this file and can be tested against any other implicit operator.
"""
struct HEVICache{OPT, CCT, FAT, F1, F2}
    topo::ColumnTopology
    op::OPT
    cc::CCT
    fac::FAT
    fimp!::F1
    solve!::F2
    tableau::Symbol
    gdt::Float64
    # :RS -- linearise once, about the deck's reference state qe (the default,
    #        and what this code did before the switch existed).
    # :PS -- re-linearise about the PREVIOUS SOLUTION every `update_freq`
    #        steps, following Giraldo et al. (arXiv:2311.11425), whose timings
    #        make LHEVI the fastest of the HEVI variants and who refresh every
    #        five steps by default.
    linearization::Symbol
    update_freq::Int
end

"""
    hevi_fimp!(du, u, p, t)

The `f_imp` hook. Thin wrapper so `HEVI_ARK` can call a plain function while
the operator state travels on `p`.
"""
hevi_fimp!(du, u, p, t) = hevi_implicit_rhs!(du, u, p, p.hevi)

"""
    hevi_relinearize!(params, hevi, u) -> hevi

Recompute the implicit operator's coefficients from the state `u` and
refactorise. This is the :PS (previous solution) half of the linearisation
switch; :RS never calls it.

Only `beta = dp/d(rho theta)` and `thetabar` move. They are what make `A` an
approximation of the vertical Jacobian, and they drift as the solution departs
from the reference state -- a stale `beta` leaves a residual vertical acoustic
term in the explicit part, which shows up as a step size that will not go where
the stability analysis says it should, with nothing in the output to explain it.

`u` holds TOTAL variables (build_hevi refuses PERT), so the thermodynamics can
be read straight off it.

The band is re-extracted by the same probing path used at setup, so the matrix
stays exactly the operator rather than drifting into an approximation of it.
"""
function hevi_relinearize!(params, hevi::HEVICache, u)

    op    = hevi.op
    npoin = Int(params.mesh.npoin)
    PhysConst = PhysicalConst{Float64}()
    γ = PhysConst.γ

    ρoff  = 0 * npoin
    ρθoff = 4 * npoin

    @inbounds for ip = 1:npoin
        ρ  = u[ρoff + ip]
        ρθ = u[ρθoff + ip]
        if ρ > 0 && ρθ > 0
            p = perfectGasLaw_ρθtoP(PhysConst, ρ, ρθ / ρ)
            op.beta[ip]     = γ * p / ρθ
            op.thetabar[ip] = ρθ / ρ
        end
        # A non-positive rho or rho*theta means the state is already broken;
        # keeping the previous coefficients there is strictly better than
        # writing a NaN into the operator and losing the error message.
    end

    # The eddy viscosity moves far more than beta and thetabar do -- from its
    # molecular floor to tens of m^2/s as the boundary layer spins up -- so on
    # a diffusion-carrying operator this refresh is the point of :PS, not a
    # refinement of it.
    op.vd === nothing || vdiff_refresh!(op.vd, params, op.vars, u)

    refactorize!(hevi.fac, params, op, hevi.cc, hevi.topo, hevi.gdt)
    return hevi
end

# The stepper in ark.jl calls `ark_relinearize!` so that it never has to know
# which scheme's hooks it is holding; this is HEVI's method.
ark_relinearize!(params, hevi::HEVICache, u) = hevi_relinearize!(params, hevi, u)


"""
    hevi_enabled(inputs) -> Bool

True when the deck asks for a HEVI run, either by naming a `HEVI_ARK`
integrator or by setting `:lhevi`.
"""
function hevi_enabled(inputs)
    get(inputs, :lhevi, false) == true && return true
    return haskey(inputs, :ode_solver) && inputs[:ode_solver] isa HEVI_ARK
end

"""
    _band_matvec!(y, AB, kl, ku, n, x)

Multiply by a matrix in LAPACK general-band storage, before it is factorised.
Used only by the setup self-check.
"""
function _band_matvec!(y, AB, kl, ku, n, x)
    fill!(y, 0.0)
    d = kl + ku + 1
    @inbounds for j = 1:n
        xj = x[j]
        xj == 0.0 && continue
        for i = max(1, j - ku):min(n, j + kl)
            y[i] += AB[d + i - j, j] * xj
        end
    end
    return y
end

"""
    hevi_verify_band(params, op, cc, topo, AB, kl, ku, n, gdt) -> Float64

Check that the assembled band really is the matrix of the operator.

Applies A to a deterministic, DSS-consistent test field and compares the
banded product `W x` against `x - γΔt (A x)` computed by the operator itself.
Returns the largest relative discrepancy, reduced over all ranks.

This single number exercises the whole chain at once: the column topology,
the level catalogue, the MPI gather/scatter election, the probe colouring and
the band index arithmetic. A bug in any of them shows up here, at setup, with
a number -- rather than later as an unexplained step-size limit, which is how
a wrong implicit matrix otherwise presents.

The test field is built from the *global* point id, so it takes the same
value on every rank that holds a shared node. A field of `rand()` would not,
and the check would fail for a reason that has nothing to do with the code
under test.
"""
function hevi_verify_band(params, op::HEVIOperator, cc::ColumnComm, topo::ColumnTopology,
                          AB, kl, ku, n, gdt)

    npoin = Int(params.mesh.npoin)
    nimp  = op.nimp
    gip   = params.mesh.ip2gip

    V = op.V
    @inbounds for q = 1:nimp, ip = 1:npoin
        V[ip, q] = sinpi(1.0e-3 * (Float64(gip[ip]) + 137.0 * q))
    end

    column_gather!(cc, V, 1:nimp)
    Xv = copy(cc.X)

    hevi_apply_A!(op.R, V, params, op)
    column_gather!(cc, op.R, 1:nimp)
    Xr = cc.X

    y   = zeros(Float64, n)
    err = 0.0
    scl = 0.0
    @inbounds for ic = 1:cc.nown
        _band_matvec!(y, AB[ic], kl, ku, n, view(Xv, :, ic))
        for i = 1:n
            expected = Xv[i, ic] - gdt * Xr[i, ic]
            err = max(err, abs(y[i] - expected))
            scl = max(scl, abs(Xv[i, ic]) + abs(gdt * Xr[i, ic]))
        end
    end

    if MPI.Comm_size(cc.comm) > 1
        err = MPI.Allreduce(err, MPI.MAX, cc.comm)
        scl = MPI.Allreduce(scl, MPI.MAX, cc.comm)
    end
    return scl > 0 ? err / scl : err
end

"""
    hevi_verify_solve(params, op, cc, fac) -> Float64

Check the factorised solve, not just the assembled matrix.

`hevi_verify_band` says the band is the right matrix; this says the banded LU
of it actually inverts. They fail for different reasons -- the first for an
indexing or topology mistake, the second for conditioning -- and on the
production mesh the columns are 481 levels deep rather than the 11 of the
unit test, so the second is the one worth having a number for.

Returns the relative residual of `x - γΔt A x - b` after solving for `x`.
"""
function hevi_verify_solve(params, op::HEVIOperator, cc::ColumnComm, fac::ColumnFactorization)

    npoin = Int(params.mesh.npoin)
    nimp  = op.nimp
    gip   = params.mesh.ip2gip
    gdt   = fac.gdt[]

    B = op.V
    @inbounds for q = 1:nimp, ip = 1:npoin
        B[ip, q] = sinpi(1.0e-3 * (Float64(gip[ip]) + 61.0 * q))
    end
    Bsave = copy(B)

    column_gather!(cc, B, 1:nimp)
    @inbounds for ic = 1:cc.nown
        LAPACK.gbtrs!('N', fac.kl, fac.ku, fac.n, fac.AB[ic], fac.ipiv[ic], view(cc.X, :, ic))
    end
    column_scatter!(cc, B, 1:nimp)      # B now holds x

    X = similar(B)
    hevi_apply_A!(X, B, params, op)

    err = 0.0; scl = 0.0
    @inbounds for q = 1:nimp, ip = 1:npoin
        err = max(err, abs(B[ip, q] - gdt * X[ip, q] - Bsave[ip, q]))
        scl = max(scl, abs(Bsave[ip, q]))
    end
    if MPI.Comm_size(cc.comm) > 1
        err = MPI.Allreduce(err, MPI.MAX, cc.comm)
        scl = MPI.Allreduce(scl, MPI.MAX, cc.comm)
    end
    return scl > 0 ? err / scl : err
end

"""
    hevi_explicit_implicit_rates(params, spec) -> (rate_exp, rate_imp)

The two eigenvalue magnitudes (1/s) that the joint IMEX stability check needs:
what the explicit half of the split still sees, and what the implicit half
takes away.

`rate_imp` is MEASURED -- it is `max|Im λ|` of the assembled column operator,
which is the honest vertical acoustic rate including the spectral-element
N²/h growth that a node-spacing estimate misses. On the rtb_hevi mesh it comes
out 1.44x the `c/Δz_min` figure the CFL table prints.

`rate_exp` is the node-spacing estimate of everything left explicit, with the
horizontal acoustic terms scaled by that SAME measured factor. Reusing it is
the point: the correction is a property of the polynomial order, not of the
direction, so calibrating it once against the one operator whose spectrum is
actually available makes the horizontal estimate as good as the vertical one.

WHY THE THIRD RETURN VALUE EXISTS -- READ IT BEFORE BLAMING THE TABLEAU.
`rate_exp` is a SUM over three physically unrelated processes, and a refusal
prints only the sum. That is not enough to act on, because the three want
opposite fixes: horizontal acoustics wants a smaller Δt or the 3D IMEX,
advection wants a smaller Δt, and DIFFUSION wants neither -- it wants the
closure or `:implicit_vdiff`. Diagnosed on this project's LESICP2-coarse mesh
(10 x 1 x 60 over 6400 x 1600 x 5000 m, first z element 40 m):

    rate_exp = 32.3 1/s     of which     acoustic  3.1     viscous  29.2

-- i.e. 90% of the refusal was SGS diffusion, on a mesh whose horizontal
acoustic Courant number at that Δt was 0.6. Nothing in "tableau ARS232
amplifies by 49.99" says that, and every remedy the message used to suggest
(drop Δt, change tableau) was the wrong one: 2 nu / h_z^2 goes as 1/h_z^2, so
buying stability with Δt costs a factor of ten in cost for a factor of ten in
h_z that was not the problem.

WHAT MADE IT 29 1/s THERE, since the same figure on a laminar sounding is
normally ~0.02. `hevi_verify_physics` calls `rhs!` a few lines before this
runs, so `sgs.mu_turb` is POPULATED by the time `cfl_limits` reads it -- these
rates are not taken on a zero-viscosity state, whatever the CFL report's own
"taken at t = 0" caveat says. And nu_t goes as the filter width SQUARED. On
that mesh the elements are 640 x 1600 x 40 m with ONE element across y, so
`:les_filter_width => :max` (the default) hands the whole domain the 1600 m
dummy direction: 1600/40 = 40 on the width, 1600x on nu_t, and a viscous rate
that was negligible becomes the binding term. `compute_element_size_driver`
warns about exactly this at setup; `:geometric` is the fix on a quasi-2D mesh.

The vertical viscous term is dropped from `rate_exp` when `op.vd !== nothing`,
because `:implicit_vdiff` has already put d/dz(mu d/dz) into the column
operator -- charging it to the explicit half would refuse a Δt the scheme can
actually take. It is a real (decaying) eigenvalue, so it does not belong in
`rate_imp` either: the implicit half takes it unconditionally.

Advection goes in unscaled. It is included rather than dropped so an
advection-dominated mesh does not quietly slip past.
"""
function hevi_explicit_implicit_rates(params, spec, op::HEVIOperator)

    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    qe    = params.qp.qe

    u = zeros(Float64, npoin * neqs)
    @inbounds for ieq = 1:neqs, ip = 1:npoin
        u[(ieq - 1) * npoin + ip] = qe[ip, ieq]
    end
    L = cfl_limits(params, u)

    _r(dt) = (isfinite(dt) && dt > 0) ? 1.0 / dt : 0.0

    # spec = (max|Im λ| of the column operator, max|Re|/max|Im|, c/Δz_min)
    rate_imp = spec[1]
    κ = (isfinite(spec[3]) && spec[3] > 0) ? spec[1] / spec[3] : 1.0

    rate_ac  = κ * (_r(L.dt_acoustic_x) + _r(L.dt_acoustic_y))
    rate_adv = _r(L.dt_advective_x) + _r(L.dt_advective_y) + _r(L.dt_advective_z)
    # :implicit_vdiff has already taken d/dz(mu d/dz); only the horizontal and
    # cross terms are still explicit.
    vdz      = op.vd === nothing ? _r(L.dt_viscous_z) : 0.0
    rate_vis = _r(L.dt_viscous_x) + _r(L.dt_viscous_y) + vdz

    rate_exp = rate_ac + rate_adv + rate_vis

    return rate_exp, rate_imp,
           (acoustic = rate_ac, advective = rate_adv, viscous = rate_vis,
            vdiff_implicit = op.vd !== nothing, numax = L.numax,
            hmin_z = L.hmin_z, kappa = κ)
end

"""
    hevi_verify_physics(params, op, topo; eps_rel) -> Vector{NTuple{5,Float64}}

Does `f_imp` actually equal the vertical part of the real `rhs!`?

THE GAP THIS CLOSES
-------------------
`hevi_verify_band` and `hevi_verify_solve` are both self-consistent: the first
compares the assembled matrix against the operator it was probed from, the
second inverts the matrix it was just handed. Neither ever looks at `rhs!`.
So an implicit operator that is a factor of two small, or missing a term, or
carrying the wrong sign on the buoyancy coupling, passes both -- and the split
stays CONSISTENT (f_exp = rhs! - f_imp by construction, so the answer is still
the right answer), it just stops being a SPLIT. Whatever f_imp fails to
capture stays in f_exp and stays explicit, and the vertical acoustic step-size
limit comes back with nothing in the output to say why. That is not a
hypothetical failure mode; it is the one that looks exactly like "HEVI is
mysteriously no faster than explicit".

THE MEASUREMENT
---------------
Perturb the reference state by a smooth profile that is a function of the
LEVEL INDEX only, so it is exactly uniform in x and y. On a vertically
extruded mesh every horizontal derivative of such a field is identically zero
at the discrete level, so

    rhs!(qe + δ) - rhs!(qe)

is purely the vertical response, and to first order in |δ| it must equal
`A·δ` -- which is precisely what f_imp claims to be. The ρw profile is chosen
to vanish at the floor and the lid so the free-slip projection inside `rhs!`
does not alter the state and contaminate the difference.

Two numbers per implicit equation:

  `ratio` = <r_rhs, r_imp> / <r_imp, r_imp>, the amplitude f_imp captures.
            1.0 is right. 0.5 means half the vertical term is still explicit;
            -1.0 means the sign is inverted and the split is making it worse.
  `rerr`  = max|r_rhs - r_imp| / max|r_rhs|, what is left over.

The viscous terms in `rhs!` are in `r_rhs` and not in `r_imp`, but for a
domain-scale profile they are three to four orders below the acoustic terms,
so they set the floor on `rerr` and do not move `ratio`.
"""
function hevi_verify_physics(params, op::HEVIOperator, topo::ColumnTopology;
                             eps_rel::Float64 = 1.0e-6)

    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    nlev  = topo.nlev
    qe    = params.qp.qe
    comm  = params.mesh.parts.comm

    lev = level_of_node(topo, npoin)

    # Scales of the reference state, global so every rank perturbs identically.
    ρs = 0.0; θs = 0.0
    @inbounds for ip = 1:npoin
        ρs = max(ρs, abs(qe[ip, 1]))
        θs = max(θs, abs(qe[ip, 5]))
    end
    if MPI.Comm_size(comm) > 1
        ρs = MPI.Allreduce(ρs, MPI.MAX, comm)
        θs = MPI.Allreduce(θs, MPI.MAX, comm)
    end
    PC = PhysicalConst{Float64}()
    cs = sqrt(PC.γ * PC.Rair * 300.0)                      # ~ speed of sound

    u0 = zeros(Float64, npoin * neqs)
    u1 = zeros(Float64, npoin * neqs)
    @inbounds for ieq = 1:neqs, ip = 1:npoin
        u0[(ieq - 1) * npoin + ip] = qe[ip, ieq]
    end
    copyto!(u1, u0)

    δ = zeros(Float64, npoin, op.nimp)
    @inbounds for ip = 1:npoin
        l = lev[ip]
        l == 0 && continue
        σ = nlev > 1 ? (l - 1) / (nlev - 1) : 0.0
        # ρw vanishes at both ends so the free-slip projection is a no-op on it
        d1 = eps_rel * ρs      * cospi(2σ)
        d4 = eps_rel * ρs * cs * sinpi(2σ)
        d5 = eps_rel * θs      * cospi(4σ)
        for (q, ieq) in enumerate(op.vars)
            d = ieq == 1 ? d1 : ieq == 4 ? d4 : ieq == 5 ? d5 : 0.0
            δ[ip, q] = d
            u1[(ieq - 1) * npoin + ip] += d
        end
    end

    # `rhs!` mutates its state argument (the Dirichlet projection writes back
    # into u), so both evaluations get their own copy.
    du0 = zeros(Float64, npoin * neqs)
    du1 = zeros(Float64, npoin * neqs)
    rhs!(du0, u0, params, 0.0)
    rhs!(du1, u1, params, 0.0)

    R = similar(δ)
    hevi_apply_A!(R, δ, params, op)

    z   = @view params.mesh.coords[3, :]
    lev = level_of_node(topo, npoin)

    # A/B harness: dump A*δ keyed by GLOBAL point id so two runs at different
    # rank counts can be diffed node by node. δ is a function of the level
    # index and of globally-reduced scales only, so it is identical between
    # runs by construction and any difference in the dump is a difference in
    # the operator itself. Off unless asked for; this is a debugging tool.
    if get(ENV, "JEXPRESSO_HEVI_DUMP", "0") == "1"
        gip = params.mesh.ip2gip
        open("hevi_dump_r$(MPI.Comm_rank(comm)).txt", "w") do io
            for ip = 1:npoin
                Base.print(io, gip[ip])
                for q = 1:op.nimp
                    Base.print(io, " ", R[ip, q])
                end
                # ...and the SAME response as computed by Jexpresso's own rhs!.
                # If A and rhs! move together between rank counts, the operator
                # faithfully mirrors the discretisation and the partition
                # dependence lives upstream of HEVI; if only A moves, it is mine.
                for (q, ieq) in enumerate(op.vars)
                    off = (ieq - 1) * npoin
                    Base.print(io, " ", du1[off + ip] - du0[off + ip])
                end
                Base.println(io)
            end
        end
    end

    out = NTuple{5,Float64}[]
    for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        num = 0.0; den = 0.0; err = 0.0; scl = 0.0; ipw = 1
        @inbounds for ip = 1:npoin
            a = du1[off + ip] - du0[off + ip]     # the true vertical response
            b = R[ip, q]                          # what f_imp claims it is
            num += a * b
            den += b * b
            d = abs(a - b)
            d > err && (err = d; ipw = ip)
            scl = max(scl, abs(a))
        end
        # WHERE the worst node is, not just how bad it is. A discrepancy that
        # only appears in parallel has to be localised before it can be
        # explained, and "which level" separates the two candidates outright:
        # the floor/lid rows point at the wall flux treatment, an interior
        # level at the assembly across a rank cut.
        zw = z[ipw]; lw = Float64(lev[ipw])
        if MPI.Comm_size(comm) > 1
            num = MPI.Allreduce(num, MPI.SUM, comm)
            den = MPI.Allreduce(den, MPI.SUM, comm)
            scl = MPI.Allreduce(scl, MPI.MAX, comm)
            # carry z and the level of the argmax with the max itself
            errmax = MPI.Allreduce(err, MPI.MAX, comm)
            if err < errmax
                zw = 0.0; lw = 0.0
            end
            zw = MPI.Allreduce(zw, MPI.MAX, comm)
            lw = MPI.Allreduce(lw, MPI.MAX, comm)
            err = errmax
        end
        push!(out, (Float64(ieq), den > 0 ? num / den : 0.0,
                    scl > 0 ? err / scl : err, zw, lw))
    end
    return out
end

"""
    hevi_column_spectrum(AB, kl, ku, n, gdt, params, topo, op) -> (maxim, reim, cref)

Eigenvalues of one column's implicit operator, and what they should be.

The band self-check proves the matrix IS the operator. It says nothing about
whether the operator is the right PHYSICS -- if the linearisation coefficients
were wrong (a bad β, a θ̄ read off the wrong slot, a qe that is not the state
the case was initialised on) the band would still match the operator faithfully
and the check would pass, while the vertical acoustic mode stayed partly
explicit and the run went unstable at a step size the report said was safe.

A hyperbolic acoustic operator has essentially imaginary spectrum with
|λ| ≈ c/Δz times an O(1) factor from the SEM derivative. Reporting both the
measured radius and c/Δz_min turns "is the implicit operator the right one" into
a number, on the real mesh, at setup.

Returns (max|Im λ|, max|Re λ|/max|Im λ|, c/Δz_min).
"""
function hevi_column_spectrum(AB, kl, ku, n, gdt, params, topo::ColumnTopology,
                              op::HEVIOperator)

    # EVERY RANK REACHES THE REDUCTIONS AT THE BOTTOM.
    #
    # This used to `return (NaN, NaN, NaN)` early when the rank owned no
    # columns, and the caller then branched on that NaN before running further
    # collectives. Under a p4est space-filling-curve partition rank 0 routinely
    # owns nothing, so the ranks took different branches, the collectives after
    # this point paired up wrongly, and the setup report printed 4.6e18 columns
    # while the joint-stability check -- the one that catches an amplifying
    # tableau -- was skipped entirely on a parallel run. A collective behind a
    # rank-local condition is the same trap params_setup.jl warns about where
    # it builds this cache eagerly.
    mi = -Inf; mr = -Inf; cmax = 0.0; dzmin = Inf

    # Sample several owned columns, not just the first. On one rank AB[1] is
    # global column 1; on many ranks it is a different column per rank, so a
    # figure taken from AB[1] alone silently changes meaning with the rank
    # count and cannot be compared across runs.
    if !isempty(AB)
        d    = kl + ku + 1
        W    = zeros(Float64, n, n)
        nsam = min(length(AB), 8)
        for k = 1:nsam
            ic = 1 + ((k - 1) * (length(AB) - 1)) ÷ max(nsam - 1, 1)
            fill!(W, 0.0)
            @inbounds for j = 1:n, i = max(1, j - ku):min(n, j + kl)
                W[i, j] = AB[ic][d + i - j, j]
            end
            λ  = eigvals((Matrix(1.0I, n, n) .- W) ./ gdt)
            mi = max(mi, maximum(abs, imag.(λ)))
            # SIGNED, not |Re|. A discrete operator M^-1 K with M SPD and K
            # skew has a purely imaginary spectrum; a real part means the
            # skewness is broken, and its SIGN says whether the modes decay or
            # GROW. |Re| cannot tell those apart, and the difference is the
            # difference between a harmless discretisation blemish and an
            # exponential instability.
            mr = max(mr, maximum(real, λ))
        end
    end

    # c and Δz for the reference figure, from whatever consecutive level pairs
    # this rank holds. Reducing the two ingredients separately -- max c, min Δz
    # -- rather than the ratio is what makes the answer independent of where
    # the partition cut the columns; requiring a whole column to be local (as
    # this did) yields nothing at all when every column is split, which is the
    # normal case here.
    z = @view params.mesh.coords[3, :]
    @inbounds for ic = 1:topo.ncol, il = 1:topo.nlev-1
        ipa, ipb = topo.node[il, ic], topo.node[il+1, ic]
        (ipa == 0 || ipb == 0) && continue
        dz = abs(z[ipb] - z[ipa])
        dz > 0 && (dzmin = min(dzmin, dz))
        cmax = max(cmax, sqrt(max(op.beta[ipa] * op.thetabar[ipa], 0.0)))
    end

    comm = params.mesh.parts.comm
    if MPI.Comm_size(comm) > 1
        mi    = MPI.Allreduce(mi,    MPI.MAX, comm)
        mr    = MPI.Allreduce(mr,    MPI.MAX, comm)
        cmax  = MPI.Allreduce(cmax,  MPI.MAX, comm)
        dzmin = MPI.Allreduce(dzmin, MPI.MIN, comm)
    end

    isfinite(mi) || return (NaN, NaN, NaN)
    cref = (isfinite(dzmin) && dzmin > 0) ? cmax / dzmin : NaN
    return (mi, mi > 0 ? mr / mi : NaN, cref)
end

"""
    build_hevi(params, inputs) -> HEVICache

Discover the column structure, build the implicit operator, plan the MPI
redistribution, assemble and factorise one banded matrix per column, and
report what it found.
"""
function build_hevi(params, inputs)

    mesh = params.mesh
    comm = mesh.parts.comm
    rank = MPI.Comm_rank(comm)
    t0   = time_ns()

    alg = inputs[:ode_solver]
    alg isa HEVI_ARK ||
        error("HEVI setup was asked for but :ode_solver is $(typeof(alg)). ",
              "Set :ode_solver => HEVI_ARK(:ARS343) (or drop :lhevi).")
    # HEVI linearises about qe and acts on (u - qe), which is the deviation
    # only when u holds the TOTAL variables. Under PERT() u IS the deviation
    # already and qe is the base state, so (u - qe) would be q - 2qe: silently
    # wrong, and wrong in a way no self-check here would catch, because the
    # band would still faithfully match the operator it was probed from.
    #
    # Supporting PERT is a small, well-defined change -- feed u straight in as
    # the deviation, and re-derive the buoyancy term against whatever the
    # case's perturbation-form user_source! uses -- but it has to be verified
    # per case rather than assumed.
    get(inputs, :SOL_VARS_TYPE, TOTAL()) isa TOTAL ||
        error("HEVI currently supports :SOL_VARS_TYPE => TOTAL() only; this case asks ",
              "for $(typeof(get(inputs, :SOL_VARS_TYPE, TOTAL()))). Under PERT() the ",
              "solution vector already holds the deviation from qe, so the implicit ",
              "operator would subtract the reference state a second time.")

    # HEVI on the p4est space-filling-curve partition diverges, and slowly
    # enough to look like anything but a partition problem. The columnar
    # partition has neither problem: the ghost-inflation factors match, and no
    # column is split across ranks (0 of 180 measured at four ranks), so the
    # column solve needs no redistribution at all.
    #
    # The check is NOT just `:lxy_partition == true` -- that flag is a request
    # the mesh code ignores whenever there is refinement. See
    # check_columnar_partition in columns.jl for what goes wrong and why
    # nothing else here catches it.
    check_columnar_partition(inputs, comm, "HEVI")

    params.Minv isa AbstractVector ||
        error("HEVI needs a lumped (diagonal) mass matrix; this case has a dense Minv, ",
              "which comes from :QT => Exact() with :llump => false. Applying the ",
              "vertical operator would then cost O(npoin^2) per stage. Use inexact ",
              "quadrature, or lump the mass matrix.")

    get(params, :laguerre, false) == true &&
        error("HEVI does not support the Laguerre semi-infinite elements: u carries ",
              "points beyond mesh.npoin, which the column topology does not cover.")

    get(inputs, :ladapt, false) == true &&
        error("HEVI does not support :ladapt. Adaptive mesh refinement changes the ",
              "column structure mid-run, and the factorised column matrices, the ",
              "gather/scatter plan and the level catalogue would all be stale after ",
              "the first adaptation.")

    tab = alg.tab
    # The step size the integrator will actually be handed, computed exactly as
    # time_loop! computes it -- Float32, and divided down by the AMR level.
    # Deriving γΔt from inputs[:Δt] instead would be off in the last bits, the
    # first stage solve would find a stale factorisation, and every column
    # would be silently refactorised inside the time loop.
    ad_lvl_max = MPI.Allreduce(maximum(mesh.ad_lvl; init = 0), MPI.MAX, comm)
    Δt  = Float64(Float32(inputs[:Δt] / (2.0^ad_lvl_max)))
    gdt = tab.γ * Δt

    # LINEARISATION MODE.
    #
    # :RS keeps the operator's coefficients frozen at the deck's reference
    # state qe for the whole run. :PS recomputes them from the current solution
    # every :hevi_update_freq steps and refactorises.
    #
    # WHAT :PS CHANGES, AND WHAT IT DELIBERATELY DOES NOT. The operator's
    # stiffness-removing power lives entirely in its COEFFICIENTS (beta and
    # thetabar): f_imp(u) = A(u - qe) has Jacobian A regardless of what is
    # subtracted, so how well the split removes the vertical acoustic term
    # depends on A alone. The offset cancels in f_exp = rhs! - f_imp. So :PS
    # refreshes the coefficients and keeps measuring the deviation from qe.
    #
    # That is a deliberate departure from Giraldo et al., who linearise about
    # the previous solution outright. Keeping qe as the origin costs nothing in
    # stiffness and preserves a property their form gives up: A acts on
    # (u - qe), so f_imp(qe) = 0 EXACTLY, the reference state stays an exact
    # steady state of the implicit half, and whatever discrete hydrostatic
    # imbalance the code has stays entirely inside the explicit part. Their
    # motivation for the previous-solution form -- that a balanced reference
    # state may not exist for real data or for whole-atmosphere runs with large
    # day/night swings -- does not apply to a deck that always defines qe.
    lin_mode = Symbol(get(inputs, :hevi_linearization, :RS))
    lin_mode in (:RS, :PS) ||
        error("HEVI: :hevi_linearization must be :RS (reference state, the default) ",
              "or :PS (previous solution); got $lin_mode.")
    update_freq = Int(get(inputs, :hevi_update_freq, 5))
    lin_mode === :PS && update_freq < 1 &&
        error("HEVI: :hevi_update_freq must be >= 1 when :hevi_linearization => :PS; ",
              "got $update_freq.")

    # Implicit vertical diffusion is opt-in and, under a dynamic closure, only
    # meaningful with :PS -- see vdiff_check_linearisation for why the wrong
    # combination is worse than not asking for it at all.
    vdiff_check_linearisation(inputs, params, "HEVI", lin_mode)

    hevi_trace_init!(comm)
    get(inputs, :hevi_monitor, false) == true && (HEVI_MONITOR[] = true)
    hevi_trace("build_hevi: start")

    rank == 0 && (print(" # HEVI setup: column topology ..."); flush(stdout))
    topo = build_column_topology(mesh, comm)
    hevi_trace("build_hevi: topology done (", topo.ncol, " local columns, ",
               topo.nlev, " levels)")

    vars = haskey(inputs, :hevi_vars) ? collect(Int, inputs[:hevi_vars]) :
                                        hevi_choose_vars(params.metrics, comm)
    # Diffusion acts on u, v, w and theta, so switching it on widens the
    # implicit set even where the acoustic operator alone would not need to.
    lvdiff = vdiff_enabled(inputs)
    vars   = vdiff_vars(params, inputs, vars)
    op = build_hevi_operator(params, topo, vars;
                             lwall_flux = get(inputs, :hevi_wall_flux, true),
                             vdiff = lvdiff)

    owner, own = assign_column_owners(topo, comm)
    cc  = build_column_comm(topo, owner, own, comm, length(vars))

    rank == 0 && (print(" operator ..."); flush(stdout))
    verr = NaN
    spec = (NaN, NaN, NaN)
    hook = get(inputs, :hevi_verify, true) != true ? nothing :
        (AB, kl, ku, n) -> begin
            rank == 0 && (print(" self-check ..."); flush(stdout))
            spec = hevi_column_spectrum(AB, kl, ku, n, gdt, params, topo, op)
            verr = hevi_verify_band(params, op, cc, topo, AB, kl, ku, n, gdt)
            if !(verr < 1.0e-8)
                error("HEVI self-check failed: the assembled column matrix differs from ",
                      "the operator by a relative $(verr). The implicit solve would be ",
                      "wrong. This usually means the mesh is not a clean vertical ",
                      "extrusion -- see build_column_topology in hevi/columns.jl.")
            end
        end
    hevi_trace("build_hevi: probing the operator and factorising")
    fac = build_column_factorization(params, op, cc, topo, gdt; verify_hook = hook)
    hevi_trace("build_hevi: factorisation done")

    serr = NaN
    phys = NTuple{5,Float64}[]
    if hook !== nothing
        phys = hevi_verify_physics(params, op, topo)
        serr = hevi_verify_solve(params, op, cc, fac)
        if !(serr < 1.0e-6)
            error("HEVI self-check failed: the banded solve leaves a relative residual ",
                  "of $(serr). The column matrices are assembled correctly (that check ",
                  "passed) but their factorisation does not invert them, which points at ",
                  "conditioning rather than at indexing.")
        end
    end
    rank == 0 && (println(" done"); flush(stdout))

    # Fast-operator check: build the full acoustic operator and verify it
    # against the vertical one this setup already trusts. Off by default -- it
    # allocates a second operator -- and only meaningful while the substepping
    # scheme that uses it is being developed.
    if get(ENV, "JEXPRESSO_HEVI_FAST_CHECK", "0") == "1"
        opf = build_hevi_fast_operator(params, topo; lwall_flux = op.lwall_flux)
        vf  = hevi_verify_fast(params, opf, op, topo)
        # Cost of one f_fast application against one full rhs!. This is the
        # number that decides whether substepping can pay at all: an outer step
        # costs 3 rhs! + ns f_fast, so the whole scheme is worth building only
        # if f_fast is a small fraction of rhs!.
        ntime = 20
        _np = Int(params.mesh.npoin); _nq = Int(params.neqs)
        du = zeros(Float64, _np * _nq)
        u0 = zeros(Float64, _np * _nq)
        @inbounds for ieq = 1:_nq, ip = 1:_np
            u0[(ieq - 1) * _np + ip] = params.qp.qn[ip, ieq]
        end
        hevi_fast_rhs!(du, u0, params, opf)                      # warm up
        t0 = time_ns()
        for _ = 1:ntime; hevi_fast_rhs!(du, u0, params, opf); end
        t_fast = (time_ns() - t0) / 1e9 / ntime
        rhs!(du, u0, params, 0.0)                                # warm up
        t0 = time_ns()
        for _ = 1:ntime; rhs!(du, u0, params, 0.0); end
        t_rhs = (time_ns() - t0) / 1e9 / ntime
        rank == 0 && println("\n [fast op] full-vs-vertical on a horizontally uniform field: ",
                             vf.rel_vertical,
                             "\n [fast op] horizontal response to a field varying in x: ",
                             vf.horiz_response,
                             "\n [fast op] ok = ", vf.ok,
                             "\n [fast op] cost: f_fast ", round(t_fast*1e3, digits=3),
                             " ms   rhs! ", round(t_rhs*1e3, digits=3),
                             " ms   ratio ", round(t_fast/t_rhs, digits=4), "\n")
        for dtau in (0.02, 0.04, 0.08)
            e = hevi_verify_fast_energy(params, opf, dtau, 40)
            rank == 0 && println("   [fast op] FB energy growth/step at dtau=", dtau,
                                 ": ", round(e.growth_per_step, digits=6))
        end
    end

    # The joint IMEX stability of THIS tableau at THIS Δt on THIS mesh. The
    # explicit half's imaginary radius -- the number a tableau is usually
    # picked on -- is the zI = 0 edge of this rectangle, and for HEVI that is
    # the one mode that cannot go unstable. See ark_joint_amplification.
    jamp = NaN; jdt = NaN; rate_e = NaN; rate_i = NaN
    terms = nothing
    # `spec` is reduced across ranks inside hevi_column_spectrum, so this
    # condition is the SAME on every rank and the collectives inside
    # hevi_explicit_implicit_rates are reached by all of them or by none.
    if !isnan(spec[1])
        rate_e, rate_i, terms = hevi_explicit_implicit_rates(params, spec, op)
        jamp = ark_joint_amplification(tab, Δt * rate_e, Δt * rate_i)
        jdt  = ark_joint_dt_max(tab, rate_e, rate_i)
    end

    hevi_report(params, topo, op, cc, fac, tab, Δt, verr, serr, spec, phys,
                (jamp, jdt, rate_e, rate_i), (time_ns() - t0) / 1e9;
                lin = lin_mode, ufreq = update_freq, terms = terms)

    if !isnan(jamp) && jamp > 1 + 1.0e-6
        rate = log(jamp) / Δt
        # NAME THE TERM. `rate_e` is a sum over three processes with three
        # different remedies, and a message that reports only the sum sends the
        # reader after the tableau even when the tableau is fine. Measured on
        # LESICP2-coarse: 29.2 of 32.3 1/s was SGS diffusion, from an LES filter
        # width taken off a one-element-wide y direction.
        diag = ""
        if terms !== nothing
            diag = string("WHAT IS IN THE EXPLICIT HALF ($(round(rate_e; sigdigits = 3)) 1/s): ",
                          "acoustic $(round(terms.acoustic; sigdigits = 3)) + ",
                          "advective $(round(terms.advective; sigdigits = 3)) + ",
                          "viscous $(round(terms.viscous; sigdigits = 3)) 1/s. ")
            if terms.viscous > 0.5 * rate_e && terms.viscous > 1.0
                diag *= string("DIFFUSION IS THE BINDING TERM, not the acoustics and not the ",
                               "tableau: nu_max = $(round(terms.numax; sigdigits = 3)) m^2/s over ",
                               "h_z = $(round(terms.hmin_z; sigdigits = 3)) m. Neither a smaller Δt ",
                               "nor another tableau addresses that -- 2nu/h_z^2 is a property of the ",
                               "closure and the mesh. ",
                               terms.vdiff_implicit ?
                                 "The vertical part is already implicit, so this is the horizontal " *
                                 "and cross terms; the remaining lever is the closure. " :
                                 "Look FIRST at the LES filter width line in the mesh report: on a " *
                                 "quasi-2D mesh (one element across a direction) the default " *
                                 ":les_filter_width => :max takes the width from that dummy " *
                                 "direction and inflates nu_t by its square -- 1600x on the " *
                                 "10x1x60 LESICP2 mesh. :les_filter_width => :geometric is the fix " *
                                 "there. If the width is right, :implicit_vdiff => true moves " *
                                 "d/dz(mu d/dz) into the column operator for the cost of a wider " *
                                 "band and no new solve. ")
            end
        end
        error("HEVI refuses to start: tableau $(String(tab.name)) amplifies by ",
              "$(round(jamp; sigdigits = 6)) per step at Δt = $Δt on this mesh, a growth ",
              "rate of $(round(rate; sigdigits = 3)) 1/s. ", diag,
              "That is slow enough to look like ",
              "anything but a tableau problem -- it presents as a blow-up at a fixed model ",
              "time, nearly independent of Δt, that moves with the rank count. ",
              "Either drop Δt to $(round(HEVI_DT_SAFETY[] * jdt; sigdigits = 3)) s or below ",
              "(the neutral limit is $(round(jdt; sigdigits = 3)) s; running at it leaves no ",
              "margin for the estimate or for the rank count), or pick a tableau whose ",
              "JOINT stability region covers this mesh: :ARS232 (3 RHS) and :ARS443 (4 RHS, ",
              "third order) both do, :ARS343 does not despite having the largest explicit ",
              "imaginary radius of the family. Set :hevi_verify => false to run anyway.")
    end


    return HEVICache(topo, op, cc, fac, hevi_fimp!, hevi_column_solve!, tab.name, gdt,
                     lin_mode, update_freq)
end

"""
    hevi_report(...)

One block on rank 0 stating what the setup found. Two numbers here decide
whether the run will behave: the fraction of columns that had to be split
across ranks (communication volume per stage) and the memory the banded
factors take (which is per rank, and there are 64 ranks on a node).
"""
function hevi_report(params, topo, op, cc, fac, tab, Δt, verr, serr, spec, phys, joint, secs;
                     lin::Symbol = :RS, ufreq::Int = 0, terms = nothing)

    comm = cc.comm
    rank = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    nown_g  = nranks > 1 ? MPI.Allreduce(cc.nown, MPI.SUM, comm) : cc.nown
    # The split fraction is per LOCAL COLUMN INSTANCE -- how many of the
    # (rank, column) pairs that exist are partial -- so the denominator is the
    # summed local column count, not the number of global columns owned. With
    # nown_g it read "495 of 165 (300.0%)" on three ranks.
    ncol_g  = nranks > 1 ? MPI.Allreduce(topo.ncol, MPI.SUM, comm) : topo.ncol
    bytes   = sum(sizeof, fac.AB; init = 0) + sum(sizeof, fac.ipiv; init = 0)
    bytes_mx = nranks > 1 ? MPI.Allreduce(bytes, MPI.MAX, comm) : bytes
    # COLLECTIVE, so it belongs above the early return -- every rank reaches it
    # or none does. `mu` is rank-local; a max over it alone would print
    # whichever value this rank happened to hold.
    μmax    = op.vd === nothing ? 0.0 : maximum(op.vd.mu)
    μg      = nranks > 1 ? MPI.Allreduce(μmax, MPI.MAX, comm) : μmax
    rank == 0 || return nothing

    radius = ark_imaginary_radius(tab)
    nrhs   = tab.stiffly_accurate ? tab.s - 1 : tab.s
    println()
    println(" ┌── HEVI (horizontally explicit / vertically implicit) ────────────")
    @printf(" │  tableau %-8s  order %d   %d RHS/step   explicit imag. radius %.3f\n",
            String(tab.name), tab.order, nrhs, radius)
    @printf(" │  implicit variables: %s   (nimp = %d)\n",
            string(op.vars), op.nimp)
    @printf(" │  columns: %d global, %d levels each; %d owned across %d ranks\n",
            topo.ncolg, topo.nlev, nown_g, nranks)
    @printf(" │  split across ranks: %d of %d local column-instances (%.1f%%)\n",
            cc.nsplit_global, max(ncol_g, 1), 100 * cc.nsplit_global / max(ncol_g, 1))
    @printf(" │  moved per stage: %d nodes globally (%.2f MB)\n",
            cc.nmoved_global, cc.nmoved_global * op.nimp * 8 / 1024^2)
    @printf(" │  banded factors: kl=ku=%d, n=%d per column, %.1f MB on the busiest rank\n",
            fac.kl, fac.n, bytes_mx / 1024^2)
    @printf(" │  linearisation: %s\n",
            lin === :PS ?
              "LHEVI-PS -- coefficients refreshed from the solution every $(ufreq) steps" :
              "LHEVI-RS -- coefficients frozen at the reference state qe")
    if op.vd !== nothing
        @printf(" │  vertical diffusion: IMPLICIT (:implicit_vdiff) on %s, max mu = %.3g\n",
                string(filter(v -> v != 1, op.vars)), μg)
        println(" │      the horizontal and cross terms stay explicit; they are dz/dx of this")
    end
    @printf(" │  γΔt = %.6g   (γ = %.6g, Δt = %.6g); factorised once\n", fac.gdt[], tab.γ, Δt)
    if isnan(verr)
        println(" │  self-check: SKIPPED (:hevi_verify => false)")
    else
        @printf(" │  self-check: band matches operator to %.2e; solve residual %.2e\n",
                verr, serr)
    end
    if !isnan(spec[1])
        @printf(" │  column spectrum: max|Im λ| = %.3g 1/s   (c/Δz_min = %.3g, ratio %.2f)\n",
                spec[1], spec[3], spec[1] / max(spec[3], eps()))
        @printf(" │                   max(Re)/max|Im| = %+.2e  (hyperbolic: expect ~0; POSITIVE = growth)\n",
                spec[2])
        # The threshold is 0.1, not round-off. A few percent of real part is
        # EXPECTED on more than one rank and is not a HEVI defect: the
        # assembled RHS and the mass matrix pick up different duplication
        # factors at a small set of shared nodes, so the discretisation
        # Jexpresso actually applies is itself slightly partition-dependent
        # there, and this operator mirrors it faithfully -- measured at
        # 4.3e-02 on three ranks against 5.2e-16 on one, with rhs! shifting by
        # the identical 1.003253 at the identical nodes. Warning at 1e-3 would
        # fire on every parallel run and send the reader after a defect that
        # is not in this file. See hevi_verify_physics.
        if !(spec[2] < 0.1) || !(0.2 < spec[1] / max(spec[3], eps()) < 20.0)
            println(" │  WARNING: that spectrum does not look like a vertical acoustic operator.")
            println(" │  The implicit half may not be removing the mode it is supposed to, in")
            println(" │  which case Δt is limited as if HEVI were not on.")
        end
    end
    if !isempty(phys)
        println(" │  f_imp vs the real rhs!, on a vertical-only perturbation:")
        worst = 0.0
        for (ieq, ratio, rerr, zw, lw) in phys
            @printf(" │      eq %d: captures %6.3f of the vertical term   (residual %.2e",
                    Int(ieq), ratio, rerr)
            @printf(" at z=%.1fm, level %d)\n", zw, Int(lw))
            worst = max(worst, abs(ratio - 1.0))
        end
        if worst > 0.05
            println(" │  WARNING: f_imp is NOT the vertical part of rhs!. Whatever it fails to")
            println(" │  capture stays in f_exp and stays EXPLICIT, so Δt is still limited by")
            println(" │  the vertical acoustic term -- the split is consistent but not a split.")
        end
    end
    if !isnan(joint[1])
        jamp, jdt, rate_e, rate_i = joint
        @printf(" │  joint IMEX stability: explicit half %.1f 1/s, implicit half %.1f 1/s\n",
                rate_e, rate_i)
        # WHAT THE EXPLICIT HALF IS MADE OF. The sum alone cannot be acted on:
        # acoustics and advection want a smaller Δt, diffusion wants the
        # closure or :implicit_vdiff, and the remedy for one is wasted on the
        # others. Printed always, not only on a refusal, so the margin can be
        # attributed while the run is still starting successfully.
        if terms !== nothing
            @printf(" │      explicit half = acoustic %.2f + advective %.2f + viscous %.2f 1/s\n",
                    terms.acoustic, terms.advective, terms.viscous)
            if terms.viscous > 0.5 * rate_e && terms.viscous > 1.0
                @printf(" │      DIFFUSION IS THE BINDING TERM (%.0f%% of it): nu_max = %.3g m^2/s ",
                        100 * terms.viscous / max(rate_e, eps()), terms.numax)
                @printf("over h_z = %.3g m\n", terms.hmin_z)
                if terms.vdiff_implicit
                    println(" │      (the vertical part is already implicit; this is the horizontal one)")
                else
                    println(" │      A smaller Δt is the wrong remedy -- 2nu/h_z^2 does not scale with it.")
                    println(" │      Check the LES filter width warning above, then :implicit_vdiff => true.")
                end
            end
        end
        @printf(" │      max|R| over (zE,zI) = (%.3f, %.3f) is %.6f;  Δt_max here = %.4f s\n",
                Δt * rate_e, Δt * rate_i, jamp, jdt)
        frac = jdt > 0 ? Δt / jdt : Inf
        @printf(" │      Δt is %.0f%% of that limit  (recommended: %.4f s, i.e. 70%%)\n",
                100 * frac, HEVI_DT_SAFETY[] * jdt)
        if jamp > 1 + 1.0e-6
            @printf(" │  FATAL: this amplifies %.3f%% per step (%.3g 1/s). See the error below.\n",
                    100 * (jamp - 1), log(jamp) / Δt)
        elseif frac > 0.85
            println(" │  WARNING: that is no margin. Δt_max is computed from a node-spacing")
            println(" │  estimate of the explicit rate calibrated against ONE measured spectrum,")
            println(" │  and the flow speeds it assumes are the ones at t = 0. At the limit the")
            println(" │  scheme is exactly neutral, so a few percent of error in either -- or the")
            println(" │  partition-dependence of the assembly at shared nodes, which grows with")
            println(" │  the rank count -- is the difference between neutral and growing. A run")
            println(" │  can then be stable on 3 ranks and diverge on 10 at the same Δt.")
        end
    end
    @printf(" │  setup took %.2f s\n", secs)
    if cc.nsplit_global > 0
        println(" │")
        println(" │  NOTE: some columns are split across ranks, so each stage carries the")
        println(" │  gather/scatter above. That is expected under a p4est space-filling")
        println(" │  curve partition and is cheap while the percentage stays small. If it")
        println(" │  is large, pick a rank count that divides the tree-column count --")
        println(" │  tools/pick_nranks.jl --columns works this out from the coarse mesh.")
    end
    println(" └", "─"^66)
    flush(stdout)
    return nothing
end
