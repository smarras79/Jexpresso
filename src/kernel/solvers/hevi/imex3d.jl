#=============================================================================
 imex3d.jl -- a fully three-dimensional implicit-explicit (IMEX) time
 integrator: the complete linear acoustic-gravity subsystem is implicit, in
 all three directions, and everything else stays explicit.

 WHERE THIS SITS RELATIVE TO HEVI
 --------------------------------
 HEVI (hevi.jl) takes the VERTICAL acoustic terms implicitly. That removes one
 term from the explicit eigenvalue budget, so its gain is bounded by the grid's
 acoustic anisotropy -- and on the LESICP2 target mesh that works out at about
 0.9x on Δt, recovered to roughly 1.0-1.35x in wall-clock only through the
 cheaper step. Nothing at all on an isotropic mesh. The README in this directory names the two ways past that
 ceiling: acoustic substepping (substep.jl), and

     "a 3D implicit acoustic solve. Removes the constraint entirely, at the
      cost of a global elliptic solve per stage."

 This file is that. The step size stops being limited by sound in ANY
 direction and becomes limited by advection, which on an atmospheric mesh is
 slower by 1/Mach -- so the STEP grows by a factor of 20 to 30 rather than 2.
 Whether the WALL-CLOCK follows is a separate question with a different answer;
 see THE COST below, and do not assume it.

 WHAT IS IMPLICIT
 ----------------
 Exactly the operator acoustic.jl already builds and verifies:

     d(dρ)/dt   = -∇·(dρu, dρv, dρw)
     d(dρu)/dt  = -d/dx [ β · dρθ ]
     d(dρv)/dt  = -d/dy [ β · dρθ ]
     d(dρw)/dt  = -d/dz [ β · dρθ ]  -  g · dρ
     d(dρθ)/dt  = -∇·(θ̄ dρu, θ̄ dρv, θ̄ dρw)

 linearised about the reference state, acting on `u - qe`, with `β = c²/θ`.
 All five equations are implicit, because the horizontal pressure gradient and
 the horizontal mass flux live in the two HEVI leaves out. Advection stays
 explicit in every direction; so do the SGS closure, the sources, the sponge,
 the wall model and the boundary conditions, because `f_exp` is computed as
 `rhs!(u) - f_imp(u)` rather than re-derived.

 THE COST -- WHICH IS REAL, AND MEASURED
 ---------------------------------------
 `I - γΔt A` is no longer block-banded per column, so there is no direct
 factorisation and no "factorise once" -- every stage is a Krylov solve over
 the whole distributed domain. That is the entire price of the scheme, and it
 comes down to one number: the iteration count.

 The preconditioner is what keeps that number small, and it was already
 sitting next to this file. HEVI's `I - γΔt A_vertical` is exactly a
 line-relaxation smoother in the vertical -- the direction in which an
 atmospheric mesh is stiffest -- it is already factorised once per γΔt, and
 applying it is two triangular solves per column. What is left for the Krylov
 iteration is the horizontal acoustic coupling. This is the arrangement the
 HEVI README anticipated: "the column solver built here is exactly the
 line-relaxation smoother such a preconditioner wants in the vertical".

 Measured on two meshes over an 8x range of Δt, the result is

     iterations  ~  25 * CFL_h ,      CFL_h = γΔt * c / h_x

 -- the HORIZONTAL acoustic Courant number, independent of h_z. That is the
 preconditioner working exactly as designed, and it is also the bill. It is
 LINEAR in Δt, because the operator is skew and the spectrum of I - γΔt A is a
 line segment of that length; so the cost of a step grows with Δt and the cost
 per unit SIMULATED time approaches a floor rather than falling. Past the point
 where the four rhs! evaluations stop dominating, a larger Δt buys nothing.

 BE CLEAR-EYED ABOUT WHAT THAT MEANS. One Krylov iteration measures 3.74x one
 HEVI vertical-operator application, which on the rtb_hevi case is about 1.8
 full rhs! evaluations. On a case whose rhs! is cheap -- rtb_imex's is a bare
 inviscid RHS -- this scheme removes the acoustic step-size restriction
 completely and is still SLOWER in wall-clock than either HEVI or the explicit
 scheme. It wins when rhs! is expensive relative to the linear acoustic
 operator, which is a property of the case and not of this file. README_IMEX3D.md
 states the break-even as an inequality; JEXPRESSO_HEVI_PROFILE=1 measures the
 two sides of it for any deck.

 WHY ARS343, WHEN HEVI REFUSES IT
 --------------------------------
 Not an oversight in either place -- see `ark_wedge_amplification`. HEVI
 splits the acoustics between the halves, so two INDEPENDENT wavenumbers load
 `zE` and `zI` and the whole rectangle is reachable; ARS343 amplifies over
 most of it. Here the acoustics are entirely implicit, the explicit half sees
 only advection, and one wavenumber sets both -- so the reachable set is the
 wedge `zE <= mach·zI`, on which ARS343 is neutral out to `zE = 2.83` --
 which is its full explicit imaginary radius, i.e. the acoustics cost it
 nothing -- against ARS443's 1.57 and ARS232's 1.15. The setup report measures
 this for the deck's own mesh and Δt rather than quoting it, and
 test/imex3d/test_wedge_stability.jl recomputes the table above.
=============================================================================#

#-----------------------------------------------------------------------------
# Monitoring
#-----------------------------------------------------------------------------
const IMEX_MONITOR       = Ref(false)
const IMEX_MONITOR_EVERY = Ref(50)

"""
    IMEX3DPrecond

The preconditioner: HEVI's vertically-implicit operator, factorised once per
γΔt into one banded LU per column.

`pvars` IS NOT ALWAYS ALL FIVE EQUATIONS, and that is the point. On a mesh
whose ζ lines are vertical the vertical operator's ρu and ρv rows are exactly
the identity — `dζ/dx` and `dζ/dy` vanish, so no horizontal momentum enters the
ζ flux divergence — and inverting them is inverting `I`. Preconditioning only
`(ρ, ρw, ρθ)` and leaving the horizontal momenta alone is therefore not an
approximation there; it is the same operator with the trivial rows dropped.

It is worth doing because the banded solve is one of the two costs per Krylov
iteration and it scales as `n·(kl+ku) ∝ nvar²`: three equations instead of five
is `(3/5)² = 0.36` of the arithmetic, `0.36` of the factor storage, and `3/5` of
the bytes moved by the gather/scatter for split columns.

On a terrain-following mesh those rows are not the identity, `hevi_choose_vars`
says so, and `pvars` widens to all five.

`opv` is a SEPARATE operator object from the 3D one even though the two share
their coefficients, because each carries its own element-loop scratch and the
stage solve applies both, interleaved, many times per step.
"""
struct IMEX3DPrecond{OPT}
    opv::OPT
    cc::ColumnComm
    fac::ColumnFactorization
    topo::ColumnTopology
    # Which columns of the packed 5-wide field this preconditions. The
    # remaining ones get the identity, which is what they deserve.
    pvars::Vector{Int}
end

"""
    imex3d_precond!(V, pc, params, gdt) -> V

Apply `(I - γΔt A_vertical)^-1` in place to a packed `npoin x 5` deviation
field.

Refactorises if γΔt has moved, for the same reason `hevi_column_solve!` does:
OrdinaryDiffEq shortens the step that lands on a `tstop` and Jexpresso passes
every diagnostic time as one. Solving with a stale factorisation would not be
wrong here -- a preconditioner is allowed to be approximate -- but it would
quietly cost iterations, which is the one currency this scheme runs on.
"""
function imex3d_precond!(V::AbstractMatrix, pc::IMEX3DPrecond, params, gdt::Real)

    fac = pc.fac
    if abs(fac.gdt[] - Float64(gdt)) > 1.0e-14 * max(abs(gdt), 1.0)
        hevi_trace("    IMEX: refactorising the column preconditioner: gdt ",
                   fac.gdt[], " -> ", Float64(gdt))
        _t = hevi_tic()
        refactorize!(fac, params, pc.opv, pc.cc, pc.topo, gdt)
        if hevi_prof_on()
            HEVI_PROFILE.t_refac += (time_ns() - _t) * 1e-9; HEVI_PROFILE.n_refac += 1
        end
    end

    # `pvars` selects which columns of the packed 5-wide field take part;
    # column_gather!/scatter! index `vars[m]` for m = 1:cc.nimp, so a
    # non-contiguous set like (1, 4, 5) works directly and the untouched
    # columns keep their values -- i.e. the identity on ρu and ρv.
    column_gather!(pc.cc, V, pc.pvars)
    @inbounds for ic = 1:pc.cc.nown
        LAPACK.gbtrs!('N', fac.kl, fac.ku, fac.n, fac.AB[ic], fac.ipiv[ic],
                      view(pc.cc.X, :, ic))
    end
    column_scatter!(pc.cc, V, pc.pvars)
    return V
end

"""
    IMEX3DCache

Everything the 3D IMEX stage solve needs, assembled once at setup and reached
from `params.imex` by the integrator.

`fimp!` and `solve!` are the two hooks `IMEX_ARK` calls, and they carry the
same signatures HEVI's do -- which is why one stepper drives both.
"""
struct IMEX3DCache{OPT, PCT, F1, F2}
    topo::ColumnTopology
    op::OPT                     # the full 3D acoustic operator
    pc::PCT                     # IMEX3DPrecond, or nothing
    ws::GMRESWorkspace
    B::Matrix{Float64}          # packed right-hand side  (npoin x 5)
    X::Matrix{Float64}          # packed solution         (npoin x 5)
    fimp!::F1
    solve!::F2
    tableau::Symbol
    gdt::Float64
    linearization::Symbol
    update_freq::Int
    comm::MPI.Comm
end

"""
    imex3d_fimp!(du, u, p, t)

The `f_imp` hook: the full 3D acoustic operator applied to `u - qe`, written
into the flat `du` layout. Zero in every equation beyond the fifth, which is
what leaves moisture and tracers entirely to the explicit half.
"""
imex3d_fimp!(du, u, p, t) = hevi_fast_rhs!(du, u, p, p.imex.op)

"""
    imex3d_solve!(dst, src, params, gdt) -> dst

The stage solve: `dst = qe + (I - γΔt A)^-1 (src - qe)`, by right-
preconditioned GMRES over the whole distributed domain.

EVERY RANK RUNS THE SAME NUMBER OF ITERATIONS. The iteration is driven
entirely by quantities that come out of `MPI.Allreduce`, and MPI guarantees
every process the identical result of a reduction, so the loop conditions are
bit-identical across ranks. That is not a nicety: the operator contains
`assemble_mpi!`, so a rank that decided to stop one iteration early would
leave the others blocked inside a halo exchange forever. Anything added to
this loop that branches on a rank-local quantity reintroduces that deadlock.
"""
function imex3d_solve!(dst, src, params, gdt::Real)

    imex  = params.imex
    op    = imex.op
    ws    = imex.ws
    pc    = imex.pc
    npoin = Int(params.mesh.npoin)
    qe    = params.qp.qe
    g     = Float64(gdt)

    B = imex.B
    X = imex.X
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            B[ip, q] = src[off + ip] - qe[ip, ieq]
        end
    end

    # Closures over concretely-typed locals: the bodies are ordinary calls into
    # functions that specialise on their arguments, so the per-solve closure
    # construction is four allocations a step and nothing inside the Arnoldi
    # loop is dynamic.
    matvec! = let params = params, op = op, g = g
        (W, V) -> begin
            hevi_apply_A!(W, V, params, op)
            @inbounds @. W = V - g * W
            return W
        end
    end
    precon! = pc === nothing ? (Z -> Z) :
              let pc = pc, params = params, g = g
                  Z -> imex3d_precond!(Z, pc, params, g)
              end

    hevi_trace("    IMEX: entering GMRES, gdt=", g)
    iters, rel, converged = gmres_solve!(X, B, ws, matvec!, precon!)
    hevi_trace("    IMEX: GMRES done in ", iters, " iterations, rel. residual ", rel)

    if !converged
        # Not fatal, and deliberately so: one loose stage solve degrades the
        # order of the step, it does not corrupt the state, and a run that is
        # 200 steps from a diagnostic dump is better finished than aborted.
        # But it IS reported every time it happens on rank 0, because a stage
        # solve that stops converging is the one failure mode of this scheme
        # that has no other symptom until the step size looks mysteriously
        # limited.
        if MPI.Comm_rank(imex.comm) == 0
            @warn string("IMEX3D: the stage solve did not reach its tolerance -- ",
                         iters, " iterations left a relative residual of ", rel,
                         " against :imex_rtol = ", ws.rtol, ". Raise :imex_maxiter, ",
                         "raise :imex_restart, or lower Δt.") maxlog = 20
        end
    end

    dst === src || copyto!(dst, src)
    @inbounds for (q, ieq) in enumerate(op.vars)
        off = (ieq - 1) * npoin
        for ip = 1:npoin
            dst[off + ip] = X[ip, q] + qe[ip, ieq]
        end
    end

    if IMEX_MONITOR[] && ws.nsolve % IMEX_MONITOR_EVERY[] == 0 &&
       MPI.Comm_rank(imex.comm) == 0
        @printf(" [imex] %d solves, %.1f Krylov iterations each on average; last %d at %.2e, worst %.2e\n",
                ws.nsolve, ws.niter / max(ws.nsolve, 1), ws.last_iters,
                ws.last_relres, ws.worst_relres)
        flush(stdout)
    end
    return dst
end

"""
    imex3d_enabled(inputs) -> Bool

True when the deck asks for a fully 3D IMEX run, either by naming an
`IMEX_ARK` integrator or by setting `:limex`.
"""
function imex3d_enabled(inputs)
    get(inputs, :limex, false) == true && return true
    return haskey(inputs, :ode_solver) && inputs[:ode_solver] isa IMEX_ARK
end

"""
    imex_lateral_walls(params, mode) -> (wallx, wally)

Per-node free-slip flags for the lateral boundaries, where the implicit
horizontal mass flux has to vanish.

`mode`:

  `:auto`  use the deck's own boundary faces if the mesh carries them, and
           fall back to the domain bounding box if it does not (which is the
           case for the standalone test fixture). The default.
  `:bc`    boundary faces only; an error if they are not there.
  `:box`   bounding box only.
  `:none`  no lateral walls -- what a laterally periodic case wants, and what
           it gets automatically under `:bc`, since a periodic face is not a
           boundary face.

THE `:bc` PATH IS THE PRINCIPLED ONE. It walks exactly the faces that
`apply_boundary_conditions_dirichlet!` walks, skipping the same `periodic*`
tags it skips, and reads the same face normals -- so the implicit operator's
idea of where a wall is cannot drift from the one `rhs!` enforces. Deriving it
from the bounding box instead would call a periodic pair a wall, and mark
nodes on an open lateral boundary as reflecting.
"""
function imex_lateral_walls(params, mode::Symbol)

    mesh  = params.mesh
    npoin = Int(mesh.npoin)
    wx = falses(npoin)
    wy = falses(npoin)
    mode === :none && return (wx, wy, :none)

    mode in (:auto, :bc, :box, :none) ||
        error("IMEX3D: :imex_lateral_walls must be :auto, :bc, :box or :none; got $mode.")

    have_faces = hasproperty(mesh, :poin_in_bdy_face) && hasproperty(mesh, :bdy_face_type) &&
                 hasproperty(mesh, :nfaces_bdy) &&
                 hasproperty(params.metrics, :nx) && hasproperty(params.metrics, :ny)

    if (mode === :bc || mode === :auto) && have_faces
        met = params.metrics
        ngl = Int(mesh.ngl)
        @inbounds for iface = 1:Int(mesh.nfaces_bdy)
            tag = String(mesh.bdy_face_type[iface])
            startswith(tag, "periodic") && continue
            for j = 1:ngl, i = 1:ngl
                ip = Int(mesh.poin_in_bdy_face[iface, i, j])
                (ip < 1 || ip > npoin) && continue
                # 0.5 rather than a tolerance around 1: on a face that is not
                # axis-aligned the normal is shared between components, and
                # the component that dominates is the flux that has to go.
                abs(met.nx[iface, i, j]) > 0.5 && (wx[ip] = true)
                abs(met.ny[iface, i, j]) > 0.5 && (wy[ip] = true)
            end
        end
        return (wx, wy, :bc)
    end

    mode === :bc && error("IMEX3D: :imex_lateral_walls => :bc, but this mesh carries no ",
                          "boundary-face table (poin_in_bdy_face / bdy_face_type) or the ",
                          "metrics carry no face normals. Use :box or :none.")

    # Bounding box. mesh.xmin/xmax are already MPI-reduced over the whole mesh,
    # so this is the same set of nodes on every rank.
    tol = 1.0e-6 * max(mesh.xmax - mesh.xmin, mesh.ymax - mesh.ymin, eps())
    @inbounds for ip = 1:npoin
        x = mesh.coords[1, ip]; y = mesh.coords[2, ip]
        (abs(x - mesh.xmin) < tol || abs(x - mesh.xmax) < tol) && (wx[ip] = true)
        (abs(y - mesh.ymin) < tol || abs(y - mesh.ymax) < tol) && (wy[ip] = true)
    end
    return (wx, wy, :box)
end

"""
    imex_single_valued_error(V, op, comm) -> Float64

How far a packed field is from being single-valued at nodes more than one rank
holds, relative to its own size.

The stage solve is where this can go wrong and nowhere else can see it: every
other self-check here feeds the machinery a field that is single-valued by
construction. A continuous-Galerkin state that has gone double-valued is not a
small error -- the next explicit RHS differentiates across the jump.

Same measurement as `hevi_single_valued_error`, on the packed layout: assemble
a field of ones to get the multiplicity, assemble the field, and compare
against multiplicity times the local value.
"""
function imex_single_valued_error(V::AbstractMatrix, op, comm::MPI.Comm)

    op.dss_cache === nothing && return 0.0
    npoin, nimp = size(V)

    m = ones(Float64, npoin, nimp)
    assemble_mpi!(m, op.dss_cache)
    w = copy(V)
    assemble_mpi!(w, op.dss_cache)

    err = 0.0; scl = 0.0
    @inbounds for q = 1:nimp, ip = 1:npoin
        scl = max(scl, abs(V[ip, q]))
        err = max(err, abs(w[ip, q] - m[ip, q] * V[ip, q]))
    end
    if MPI.Comm_size(comm) > 1
        err = MPI.Allreduce(err, MPI.MAX, comm)
        scl = MPI.Allreduce(scl, MPI.MAX, comm)
    end
    return scl > 0 ? err / scl : err
end

"""
    imex_verify_solve(params, imex, gdt) -> (relres, iters, sv)

Solve `(I - γΔt A) x = b` for a deterministic, DSS-consistent `b` and measure
what came back:

  `relres`  the residual `‖b - (x - γΔt A x)‖ / ‖b‖`, computed by APPLYING the
            operator to the answer rather than by trusting the number GMRES
            reported. The two differ if the Arnoldi recursion has drifted, and
            it is the applied one that the time integrator will feel.
  `iters`   Krylov iterations. This is the cost of the scheme; the setup
            report prints it because on a mesh where the preconditioner is
            weak it is the number that decides whether the run is affordable,
            and it is available here for one solve's worth of work.
  `sv`      how single-valued the answer is across ranks.

`b` is built from the GLOBAL point id, so it takes the same value on every
rank holding a shared node. A field of `rand()` would not, and the check would
fail for a reason that has nothing to do with the code under test.
"""
function imex_verify_solve(params, imex::IMEX3DCache, gdt::Real)

    npoin = Int(params.mesh.npoin)
    op    = imex.op
    nimp  = op.nimp
    gip   = params.mesh.ip2gip
    g     = Float64(gdt)

    B = zeros(Float64, npoin, nimp)
    @inbounds for q = 1:nimp, ip = 1:npoin
        B[ip, q] = sinpi(1.0e-3 * (Float64(gip[ip]) + 53.0 * q))
    end

    X = zeros(Float64, npoin, nimp)
    matvec! = let params = params, op = op, g = g
        (W, V) -> begin
            hevi_apply_A!(W, V, params, op)
            @inbounds @. W = V - g * W
            return W
        end
    end
    precon! = imex.pc === nothing ? (Z -> Z) :
              let pc = imex.pc, params = params, g = g
                  Z -> imex3d_precond!(Z, pc, params, g)
              end

    iters, _, _ = gmres_solve!(X, B, imex.ws, matvec!, precon!)

    R = similar(X)
    matvec!(R, X)
    err = 0.0; scl = 0.0
    @inbounds for q = 1:nimp, ip = 1:npoin
        err = max(err, abs(R[ip, q] - B[ip, q]))
        scl = max(scl, abs(B[ip, q]))
    end
    if MPI.Comm_size(imex.comm) > 1
        err = MPI.Allreduce(err, MPI.MAX, imex.comm)
        scl = MPI.Allreduce(scl, MPI.MAX, imex.comm)
    end

    sv = imex_single_valued_error(X, op, imex.comm)
    return (scl > 0 ? err / scl : err, iters, sv)
end

"""
    imex_verify_wellbalanced(params, imex) -> Float64

`max |f_imp(qe)|`, which must be **exactly** zero.

The operator acts on `u - qe`, so this is zero by construction rather than
small -- which is the point of testing it. The reference state is an exact
steady state of the implicit half, whatever discrete hydrostatic imbalance the
case carries stays entirely inside the explicit part, and no amount of moving
terms between the halves can create or destroy balance. A non-zero here means
the packing has drifted from `qe`, not that the physics is off.
"""
function imex_verify_wellbalanced(params, imex::IMEX3DCache)
    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    qe    = params.qp.qe
    u  = zeros(Float64, npoin * neqs)
    du = zeros(Float64, npoin * neqs)
    @inbounds for ieq = 1:neqs, ip = 1:npoin
        u[(ieq - 1) * npoin + ip] = qe[ip, ieq]
    end
    # Through the operator directly, not through `imex.fimp!`: this runs at
    # setup, and `params.imex` -- which the hook dereferences -- does not exist
    # until params_setup merges this cache in on return.
    hevi_fast_rhs!(du, u, params, imex.op)
    m = maximum(abs, du; init = 0.0)
    return MPI.Comm_size(imex.comm) > 1 ? MPI.Allreduce(m, MPI.MAX, imex.comm) : m
end

"""
    imex3d_rates(params, u, spec, umax_deck) -> NamedTuple

The eigenvalue magnitudes (1/s) the wedge stability analysis needs, on this
mesh and this state.

  `rate_exp`  what stays explicit: advection in all three directions plus the
              diffusive terms.
  `rate_imp`  what the implicit half takes: the acoustics, all three
              directions.
  `mach`      `|v|/c`, the wedge's opening.

THE SPECTRAL CORRECTION. `cfl_limits` reports node-spacing limits, and the
spectral radius of an SEM derivative grows like `N²/h`, not like `1/h_min`.
`spec` carries the ONE place that ratio is measured rather than estimated --
the assembled column operator's spectrum against `c/Δz_min` -- so `κ` is taken
from there and applied to the first-derivative terms (advection, acoustics)
and squared for the second-derivative ones (diffusion). It is a property of
the polynomial order, not of the direction, which is what makes reusing it
legitimate.

THE ONE INPUT THAT IS NOT MEASURABLE AT SETUP. Under HEVI the explicit half is
dominated by horizontal SOUND, which is a property of the mesh and the base
state and is therefore known at `t = 0`. Here it is dominated by ADVECTION,
which at `t = 0` on a bubble or a spin-up case is zero and which grows to
whatever the flow reaches. Estimating Δt from the initial state would say any
step is safe. So `umax_deck` (`:imex_umax`) is the deck's statement of the
largest flow speed it expects, and it is used when it exceeds the measured
one. The report says which of the two it used.
"""
function imex3d_rates(params, u, spec, umax_deck::Float64)

    L = cfl_limits(params, u)
    _r(dt) = (isfinite(dt) && dt > 0 && dt < 1.0e29) ? 1.0 / dt : 0.0

    κ = (isfinite(spec[1]) && isfinite(spec[3]) && spec[3] > 0) ? spec[1] / spec[3] : 1.0

    # Advection at the speed the run will actually reach, not the one it starts
    # at. Rebuilt from the node spacings rather than taken from dt_advective_*,
    # because those are Inf on a state at rest and carry no scale to raise.
    uref = max(L.umax, L.wmax, umax_deck)
    rate_adv = uref * (_r(L.hmin_x) + _r(L.hmin_y) + _r(L.hmin_z))

    rate_visc = _r(L.dt_viscous_x) + _r(L.dt_viscous_y) + _r(L.dt_viscous_z)

    rate_exp = κ * rate_adv + κ * κ * rate_visc
    rate_imp = κ * (_r(L.dt_acoustic_x) + _r(L.dt_acoustic_y) + _r(L.dt_acoustic_z))

    cref = L.cmax > 0 ? L.cmax : sqrt(1.4 * 287.0 * 300.0)
    mach = uref / cref

    return (rate_exp = rate_exp, rate_imp = rate_imp, mach = mach,
            rate_adv = κ * rate_adv, rate_visc = κ * κ * rate_visc,
            uref = uref, umeasured = max(L.umax, L.wmax), cmax = cref, kappa = κ,
            limits = L)
end

"""
    ark_relinearize!(params, imex, u)

The `:PS` linearisation for the 3D IMEX: refresh `β` and `θ̄` from the current
solution on BOTH operators -- the 3D one that is applied and the vertical one
that preconditions it -- and refactorise the preconditioner.

Refreshing only the 3D operator would leave the preconditioner approximating
an operator that no longer exists, which costs iterations rather than
correctness and is exactly the sort of drift that shows up as "the scheme got
slower and nobody changed anything".
"""
function ark_relinearize!(params, imex::IMEX3DCache, u)

    npoin = Int(params.mesh.npoin)
    PhysConst = PhysicalConst{Float64}()
    γ = PhysConst.γ

    op = imex.op
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
        # A non-positive ρ or ρθ means the state is already broken; keeping the
        # previous coefficients there beats writing a NaN into the operator and
        # losing the error message.
    end

    if imex.pc !== nothing
        copyto!(imex.pc.opv.beta,     op.beta)
        copyto!(imex.pc.opv.thetabar, op.thetabar)
        refactorize!(imex.pc.fac, params, imex.pc.opv, imex.pc.cc, imex.pc.topo,
                     imex.pc.fac.gdt[])
    end
    return imex
end

"""
    build_imex3d(params, inputs) -> IMEX3DCache

Discover the column structure, build the full 3D acoustic operator and the
column preconditioner, allocate the Krylov workspace, run the self-checks and
report what it found.

The refusals below are the same set HEVI refuses on, for the same reasons, plus
one of its own (`:imex_precond => :column` needs the column machinery). They
are checks on the DECK, run before anything expensive, because each of them
otherwise presents as a slow divergence rather than as a mistake.
"""
function build_imex3d(params, inputs)

    mesh = params.mesh
    comm = mesh.parts.comm
    rank = MPI.Comm_rank(comm)
    t0   = time_ns()

    alg = inputs[:ode_solver]
    alg isa IMEX_ARK ||
        error("IMEX3D setup was asked for but :ode_solver is $(typeof(alg)). ",
              "Set :ode_solver => IMEX_ARK(:ARS343) (or drop :limex).")

    Int(params.neqs) >= 5 ||
        error("IMEX3D needs the five compressible Euler equations (ρ, ρu, ρv, ρw, ρθ); ",
              "this case has neqs = $(params.neqs). The implicit operator IS the ",
              "acoustic subsystem and has nowhere to read the momenta from.")

    # Same reason as HEVI: the operator acts on u - qe, which is the deviation
    # only when u holds TOTAL variables. Under PERT() the solution vector
    # already IS the deviation, so u - qe would be q - 2qe -- and wrong in a
    # way no self-check here would catch, because the operator would still be
    # applied faithfully.
    get(inputs, :SOL_VARS_TYPE, TOTAL()) isa TOTAL ||
        error("IMEX3D currently supports :SOL_VARS_TYPE => TOTAL() only; this case asks ",
              "for $(typeof(get(inputs, :SOL_VARS_TYPE, TOTAL()))). Under PERT() the ",
              "solution vector already holds the deviation from qe, so the implicit ",
              "operator would subtract the reference state a second time.")

    # The columnar partition is needed twice over here: the preconditioner
    # solves per column, and -- the reason HEVI refuses the alternative -- on
    # the p4est space-filling-curve partition the assembled RHS and the mass
    # matrix pick up different ghost multiplicities at some rank-shared nodes,
    # so M^-1 K stops being skew and the acoustic operator acquires a positive
    # real eigenvalue. A split scheme cannot survive that, and this one less
    # than HEVI: nothing is left explicit to absorb a growing mode.
    check_columnar_partition(inputs, comm, "IMEX3D")

    params.Minv isa AbstractVector ||
        error("IMEX3D needs a lumped (diagonal) mass matrix; this case has a dense Minv, ",
              "which comes from :QT => Exact() with :llump => false. Every Krylov ",
              "iteration applies the operator, so an O(npoin²) mass solve inside it ",
              "would dominate the run. Use inexact quadrature, or lump the mass matrix.")

    get(params, :laguerre, false) == true &&
        error("IMEX3D does not support the Laguerre semi-infinite elements: u carries ",
              "points beyond mesh.npoin, which the column topology does not cover.")

    get(inputs, :ladapt, false) == true &&
        error("IMEX3D does not support :ladapt. Adaptive mesh refinement changes the ",
              "column structure mid-run, and the preconditioner's column matrices, its ",
              "gather/scatter plan and the level catalogue would all be stale after the ",
              "first adaptation.")

    tab = alg.tab
    # The step size the integrator will actually be handed, computed exactly as
    # time_loop! computes it -- Float32, and divided down by the AMR level.
    # Deriving γΔt from inputs[:Δt] instead would be off in the last bits and
    # the first stage solve would refactorise the preconditioner for nothing.
    ad_lvl_max = MPI.Allreduce(maximum(mesh.ad_lvl; init = 0), MPI.MAX, comm)
    Δt  = Float64(Float32(inputs[:Δt] / (2.0^ad_lvl_max)))
    gdt = tab.γ * Δt

    lin_mode = Symbol(get(inputs, :imex_linearization, :RS))
    lin_mode in (:RS, :PS) ||
        error("IMEX3D: :imex_linearization must be :RS (reference state, the default) ",
              "or :PS (previous solution); got $lin_mode.")
    update_freq = Int(get(inputs, :imex_update_freq, 5))
    lin_mode === :PS && update_freq < 1 &&
        error("IMEX3D: :imex_update_freq must be >= 1 when :imex_linearization => :PS; ",
              "got $update_freq.")

    hevi_trace_init!(comm)
    get(inputs, :imex_monitor, false) == true && (IMEX_MONITOR[] = true)
    # max(1, ...): the monitor takes `nsolve % every`, and a deck that set this
    # to 0 would take down the run with a DivideError from inside the stage
    # solve, which is a long way from where the mistake is.
    IMEX_MONITOR_EVERY[] = max(1, Int(get(inputs, :imex_monitor_every, 50)))
    hevi_trace("build_imex3d: start")

    rank == 0 && (print(" # IMEX3D setup: column topology ..."); flush(stdout))
    topo = build_column_topology(mesh, comm)
    hevi_trace("build_imex3d: topology done (", topo.ncol, " local columns, ",
               topo.nlev, " levels)")

    wallmode = Symbol(get(inputs, :imex_lateral_walls, :auto))
    wallx, wally, wallsrc = imex_lateral_walls(params, wallmode)
    nwall = count(wallx) + count(wally)
    nwall_g = MPI.Comm_size(comm) > 1 ? MPI.Allreduce(nwall, MPI.SUM, comm) : nwall

    rank == 0 && (print(" 3D operator ..."); flush(stdout))
    lwall = get(inputs, :imex_wall_flux, true)
    op = build_hevi_fast_operator(params, topo; lwall_flux = lwall,
                                  wallx = wallx, wally = wally)

    #-------------------------------------------------------------------------
    # The preconditioner: HEVI's vertical operator over all five equations.
    #
    # Five rather than the three HEVI uses, so that it maps the same field the
    # 3D operator does and the two compose without a repacking step in the
    # inner loop. The extra two rows are the identity on a mesh whose ζ lines
    # are vertical (their vertical flux divergence vanishes), which costs band
    # width and buys the ability to precondition a terrain-following mesh
    # without a second code path.
    #-------------------------------------------------------------------------
    pcmode = Symbol(get(inputs, :imex_precond, :column))
    pcmode in (:column, :none) ||
        error("IMEX3D: :imex_precond must be :column (the default) or :none; got $pcmode.")

    verify = get(inputs, :imex_verify, true) == true
    spec   = (NaN, NaN, NaN)

    pc = nothing
    opv = nothing
    # The trivial rows of the vertical operator, dropped where they really are
    # trivial. hevi_choose_vars reads the metrics rather than a deck flag, so a
    # mesh warped by any route still gets all five. See IMEX3DPrecond.
    pvars = pcmode === :column ? hevi_choose_vars(params.metrics, comm) : Int[]
    if pcmode === :column
        rank == 0 && (print(" column preconditioner ..."); flush(stdout))
        opv = build_hevi_operator(params, topo, pvars; lwall_flux = lwall,
                                  full = false)
        owner, own = assign_column_owners(topo, comm)
        cc  = build_column_comm(topo, owner, own, comm, length(pvars))
        # The spectrum has to be taken through the factorisation's hook, not
        # afterwards: `gbtrf!` overwrites AB with the LU factors, so reading the
        # band once `fac` exists would diagonalise the factors rather than the
        # matrix. Same reason build_hevi passes its self-check the same way.
        #
        # The hook runs on EVERY rank -- hevi_column_spectrum reduces at the end
        # -- so no collective here sits behind a rank-local branch.
        hook = verify ? ((AB, kl, ku, n) -> begin
                   spec = hevi_column_spectrum(AB, kl, ku, n, gdt, params, topo, opv)
                   nothing
               end) : nothing
        fac = build_column_factorization(params, opv, cc, topo, gdt; verify_hook = hook)
        pc  = IMEX3DPrecond(opv, cc, fac, topo, copy(pvars))
    end

    #-------------------------------------------------------------------------
    # The Krylov workspace.
    #-------------------------------------------------------------------------
    npoin   = Int(mesh.npoin)
    inner   = build_distributed_inner(op.dss_cache, npoin, comm)
    restart = Int(get(inputs, :imex_restart, 20))
    maxiter = Int(get(inputs, :imex_maxiter, 200))
    rtol    = Float64(get(inputs, :imex_rtol, 1.0e-8))
    atol    = Float64(get(inputs, :imex_atol, 1.0e-30))
    ws = GMRESWorkspace(npoin, op.nimp, inner;
                        m = restart, maxiter = maxiter, rtol = rtol, atol = atol)

    imex = IMEX3DCache(topo, op, pc, ws,
                       zeros(Float64, npoin, op.nimp), zeros(Float64, npoin, op.nimp),
                       imex3d_fimp!, imex3d_solve!, tab.name, gdt,
                       lin_mode, update_freq, comm)

    #-------------------------------------------------------------------------
    # Self-checks.
    #
    # `params.imex` does not exist yet -- params_setup merges this cache in on
    # return -- so anything here that needs the hooks goes through `imex`
    # directly rather than through `params.imex`.
    #-------------------------------------------------------------------------
    vfast = nothing; sres = NaN; siter = 0; ssv = NaN; wb = NaN
    if verify
        rank == 0 && (print(" self-check ..."); flush(stdout))
        # 1. the 3D operator reduces to the vertical one on a horizontally
        #    uniform field, and moves the horizontal momenta on a field that
        #    varies in x. Needs a vertical operator; reuse the preconditioner's
        #    when it exists, build a throw-away one when it does not.
        opvchk = opv === nothing ?
                 build_hevi_operator(params, topo,
                                     hevi_choose_vars(params.metrics, comm);
                                     lwall_flux = lwall, full = false) : opv
        vfast = hevi_verify_fast(params, op, opvchk, topo)

        # 2. the stage solve itself.
        sres, siter, ssv = imex_verify_solve(params, imex, gdt)
        wb = imex_verify_wellbalanced(params, imex)

        if !(sres < 1.0e-6)
            error("IMEX3D self-check failed: the stage solve left a relative residual of ",
                  sres, " after ", siter, " Krylov iterations. Either the iteration did ",
                  "not converge (raise :imex_maxiter / :imex_restart) or the operator and ",
                  "the preconditioner disagree about the field they act on.")
        end
        if !(ssv < 1.0e-10)
            error("IMEX3D self-check failed: the stage solve left the state DOUBLE-VALUED ",
                  "across ranks by a relative ", ssv, ". A continuous-Galerkin field with ",
                  "different values in the two copies of a rank-interface node is ",
                  "differentiated across the jump by the next explicit RHS, and the jump ",
                  "grows. This points at the gather/scatter in the preconditioner, not at ",
                  "the Krylov iteration.")
        end
        if !(wb == 0.0)
            error("IMEX3D self-check failed: f_imp(qe) = ", wb, ", which must be EXACTLY ",
                  "zero -- the operator acts on u - qe, so the reference state is a steady ",
                  "state of the implicit half by construction. A non-zero here means the ",
                  "packing has drifted from qe.")
        end
        if vfast !== nothing && !vfast.ok
            error("IMEX3D self-check failed: the 3D acoustic operator does not reduce to ",
                  "the vertical one on a horizontally uniform field (relative ",
                  vfast.rel_vertical, "), or it produces no horizontal momentum response ",
                  "at all (", vfast.horiz_response, "). The first is a metric or index ",
                  "error in the ξ/η sweeps; the second means the operator has silently ",
                  "collapsed to HEVI's.")
        end
    end
    rank == 0 && (println(" done"); flush(stdout))

    #-------------------------------------------------------------------------
    # Stability, on the wedge rather than on the rectangle.
    #-------------------------------------------------------------------------
    u0 = zeros(Float64, npoin * Int(params.neqs))
    @inbounds for ieq = 1:Int(params.neqs), ip = 1:npoin
        u0[(ieq - 1) * npoin + ip] = params.qp.qn[ip, ieq]
    end
    umax_deck = Float64(get(inputs, :imex_umax, 0.0))
    R = imex3d_rates(params, u0, spec, umax_deck)
    mach = Float64(get(inputs, :imex_mach, R.mach))
    # The wedge is remarkably insensitive to `mach` -- the sup is identical for
    # 0.05 and 0.3 on every tableau in the family -- but a degenerate value
    # would collapse it onto the zI axis and report a limit that is not one.
    mach = clamp(mach, 1.0e-3, 1.0)

    wamp = ark_wedge_amplification(tab, Δt * R.rate_exp, Δt * R.rate_imp, mach)
    wdt  = ark_wedge_dt_max(tab, R.rate_exp, R.rate_imp, mach)
    # The rectangle figure too, for contrast: it is what HEVI is judged on and
    # what would (wrongly) rule ARS343 out here.
    jamp = ark_joint_amplification(tab, Δt * R.rate_exp, Δt * R.rate_imp; n = 161)

    imex3d_report(params, topo, op, pc, ws, tab, Δt, gdt,
                  (vfast, sres, siter, ssv, wb, spec),
                  (R, mach, wamp, jamp, wdt),
                  (wallsrc, nwall_g), (time_ns() - t0) / 1e9;
                  lin = lin_mode, ufreq = update_freq, verify = verify)

    guard = Symbol(get(inputs, :imex_stability_guard, :warn))
    guard in (:warn, :error, :off) ||
        error("IMEX3D: :imex_stability_guard must be :warn (the default), :error or :off; ",
              "got $guard.")
    if guard === :error && wamp[1] > 1 + 1.0e-6
        error("IMEX3D refuses to start: tableau $(String(tab.name)) amplifies by ",
              round(wamp[1]; sigdigits = 6), " per step at Δt = $Δt on this mesh, on the ",
              "WEDGE of (zE, zI) the split can actually reach. Drop Δt to ",
              round(HEVI_DT_SAFETY[] * wdt; sigdigits = 3), " s or below (the neutral ",
              "limit is ", round(wdt; sigdigits = 3), " s), lower :imex_umax if it is ",
              "pessimistic, or pick another tableau -- on the wedge the ranking is ",
              ":ARS343 (best), then :ARS443, then :ARS232, which is the reverse of the ",
              "HEVI ranking. Set :imex_stability_guard => :warn to run anyway.")
    end

    return imex
end

"""
    imex3d_report(...)

One block on rank 0 stating what the setup found.

Three numbers here decide whether the run will behave: the Krylov iteration
count (the entire cost of the scheme over an explicit step), the wedge
amplification at the deck's Δt, and the flow speed the second of those was
computed with -- which is the one input that is a guess rather than a
measurement.
"""
function imex3d_report(params, topo, op, pc, ws, tab, Δt, gdt, checks, stab, walls, secs;
                       lin::Symbol = :RS, ufreq::Int = 0, verify::Bool = true)

    comm   = params.mesh.parts.comm
    rank   = MPI.Comm_rank(comm)
    nranks = MPI.Comm_size(comm)

    npoin_g = MPI.Allreduce(Int(params.mesh.npoin), MPI.SUM, comm)
    kbytes  = sum(sizeof, ws.V; init = 0) + 3 * sizeof(ws.W)
    kb_max  = nranks > 1 ? MPI.Allreduce(kbytes, MPI.MAX, comm) : kbytes
    pbytes  = pc === nothing ? 0 :
              sum(sizeof, pc.fac.AB; init = 0) + sum(sizeof, pc.fac.ipiv; init = 0)
    pb_max  = nranks > 1 ? MPI.Allreduce(pbytes, MPI.MAX, comm) : pbytes
    nsplit  = pc === nothing ? 0 : pc.cc.nsplit_global
    ncol_g  = nranks > 1 ? MPI.Allreduce(topo.ncol, MPI.SUM, comm) : topo.ncol
    rank == 0 || return nothing

    vfast, sres, siter, ssv, wb, spec = checks
    R, mach, wamp, jamp, wdt = stab
    wallsrc, nwall_g = walls

    radius = ark_imaginary_radius(tab)
    nrhs   = tab.stiffly_accurate ? tab.s - 1 : tab.s

    println()
    println(" ┌── IMEX3D (fully three-dimensional implicit acoustics) ───────────")
    @printf(" │  tableau %-8s  order %d   %d RHS/step   explicit imag. radius %.3f\n",
            String(tab.name), tab.order, nrhs, radius)
    @printf(" │  implicit: all of %s over %d global nodes; explicit: everything else\n",
            string(op.vars), npoin_g)
    @printf(" │  γΔt = %.6g   (γ = %.6g, Δt = %.6g)\n", gdt, tab.γ, Δt)
    @printf(" │  linearisation: %s\n",
            lin === :PS ?
              "PS -- coefficients refreshed from the solution every $(ufreq) steps" :
              "RS -- coefficients frozen at the reference state qe")
    # Built before the @printf that uses it: a nested @sprintf inside a
    # ternary inside a @printf argument list parses, but only just.
    uwhere = R.uref > R.umeasured + 1.0e-12 ?
             @sprintf("(:imex_umax; the state at t=0 shows only %.1f)", R.umeasured) :
             "(measured on the initial state)"

    @printf(" │  lateral free-slip walls: %d nodes, from %s\n", nwall_g,
            wallsrc === :bc  ? "the deck's boundary faces (periodic faces skipped)" :
            wallsrc === :box ? "the domain bounding box" : "nothing (:none)")

    println(" │")
    println(" │  Stage solve: preconditioned GMRES")
    @printf(" │      preconditioner: %s\n",
            pc === nothing ? "NONE (:imex_precond => :none) -- expect many iterations" :
            string("HEVI column solve on ", pc.pvars, ", banded LU per column",
                   length(pc.pvars) == 3 ?
                     " (the omitted rows are exactly the identity here)" : ""))
    if pc !== nothing
        @printf(" │      columns: %d global x %d levels; %d of %d local instances split\n",
                topo.ncolg, topo.nlev, nsplit, max(ncol_g, 1))
        @printf(" │      banded factors %.1f MB, Krylov basis %.1f MB (busiest rank)\n",
                pb_max / 1024^2, kb_max / 1024^2)
    end
    @printf(" │      restart %d, max %d iterations, rtol %.1e\n",
            ws.m, ws.maxiter, ws.rtol)
    if verify
        @printf(" │      measured: %d iterations to %.1e on a full-spectrum right-hand side\n",
                siter, ws.rtol)
        @printf(" │      applied residual %.2e; state single-valued across ranks to %.2e\n",
                sres, ssv)
        # The iteration count is the entire cost of the scheme, so say what
        # drives it rather than leaving the reader to guess. Measured on two
        # meshes over an 8x range of Δt, it is ~25 x the HORIZONTAL acoustic
        # Courant number that survives the vertical preconditioner -- linear,
        # because the operator is skew and its spectrum is a line segment whose
        # length is that number.
        println(" │      Cost per step is LINEAR in Δt for this reason, so cost per unit")
        println(" │      simulated time approaches a floor: past the point where the 4 rhs!")
        println(" │      evaluations stop dominating the step, a larger Δt buys nothing.")
        if siter > ws.m
            println(" │      (More than one restart cycle, which is fine here: a skew operator")
            println(" │      has no spectral clustering for a longer Krylov polynomial to")
            println(" │      exploit, so raising :imex_restart buys a few percent at most.)")
        end
    else
        println(" │      self-check SKIPPED (:imex_verify => false)")
    end

    if verify
        println(" │")
        println(" │  Operator self-check:")
        @printf(" │      f_imp(qe) = %.1e   (must be exactly 0: the split cannot disturb balance)\n", wb)
        if vfast !== nothing
            @printf(" │      reduces to the vertical operator on a horizontally uniform field: %.2e\n",
                    vfast.rel_vertical)
            @printf(" │      horizontal momentum response to a field varying in x: %.3g (must be > 0)\n",
                    vfast.horiz_response)
        end
        if !isnan(spec[1])
            @printf(" │      vertical spectrum: max|Im λ| = %.3g 1/s (c/Δz_min = %.3g, κ = %.2f)\n",
                    spec[1], spec[3], R.kappa)
            @printf(" │                        max(Re)/max|Im| = %+.2e (hyperbolic: expect ~0)\n",
                    spec[2])
        end
    end

    println(" │")
    println(" │  Stability on the reachable wedge (see ark_wedge_amplification):")
    @printf(" │      explicit half %.3g 1/s  = advection %.3g + diffusion %.3g\n",
            R.rate_exp, R.rate_adv, R.rate_visc)
    @printf(" │      implicit half %.3g 1/s  (acoustics, all three directions)\n", R.rate_imp)
    @printf(" │      |v| used: %.1f m/s %s;  Mach %.3g -> wedge opening\n",
            R.uref, uwhere, mach)
    @printf(" │      (zE, zI) = (%.3f, %.3f);  max|R| on the wedge = %.6f\n",
            Δt * R.rate_exp, Δt * R.rate_imp, wamp[1])
    @printf(" │      neutral up to Δt = %.4g s; this run is at %.0f%% of it\n",
            wdt, wdt > 0 ? 100 * Δt / wdt : Inf)
    @printf(" │      for contrast, over the full RECTANGLE: %.6f\n", jamp)
    if jamp > 1 + 1.0e-6 && wamp[1] <= 1 + 1.0e-6
        println(" │      -- i.e. this tableau WOULD be refused by HEVI's criterion and is")
        println(" │      fine here. That is not a relaxation of the test, it is a different")
        println(" │      set: HEVI's two halves are loaded by independent wavenumbers, this")
        println(" │      scheme's by the same one. See ark_wedge_amplification.")
    end
    if wamp[1] > 1 + 1.0e-6
        @printf(" │      WARNING: this amplifies %.4g%% per step (%.3g 1/s) at the largest\n",
                100 * (wamp[1] - 1), log(wamp[1]) / Δt)
        @printf(" │      point of the wedge, (zE,zI) = (%.3f, %.3f). Growth that slow presents\n",
                wamp[2], wamp[3])
        println(" │      as a blow-up at a fixed model time, nearly independent of Δt, that")
        println(" │      moves with the rank count. Drop Δt, or set :imex_stability_guard")
        println(" │      => :error to make this fatal.")
    elseif wdt > 0 && Δt / wdt > 0.85
        println(" │      WARNING: no margin. The advective rate this is built on assumes the")
        println(" │      flow speed above; if the run exceeds it, the limit moves down.")
    end

    @printf(" │  setup took %.2f s\n", secs)
    println(" └", "─"^66)
    flush(stdout)
    return nothing
end
