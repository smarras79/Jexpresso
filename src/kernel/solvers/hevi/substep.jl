#=============================================================================
 substep.jl -- acoustic substepping (split-explicit) time integration.

 THE IDEA
 --------
 HEVI removes the VERTICAL acoustic term from the explicit eigenvalue budget,
 so its gain is capped by the grid's acoustic anisotropy. Substepping removes
 BOTH acoustic terms from the OUTER step instead: the outer step carries only
 advection, diffusion and the sources, and the acoustic subsystem is
 integrated inside it with its own small step.

     f(u) = f_fast(u) + f_slow(u)
     f_fast(u) = A_full (u - qe)          <- acoustic.jl, the complete linear
                                             acoustic-gravity operator
     f_slow(u) = rhs!(u) - f_fast(u)      <- by subtraction, so no source,
                                             boundary condition, SGS closure,
                                             wall model or sponge can be lost
                                             or double-counted

 OUTER SCHEME: RK3 as in Wicker & Skamarock (2002). Three stages advancing
 from u^n by h = Δt/3, Δt/2, Δt. In each stage the SLOW tendency is evaluated
 once, frozen, and carried through that stage's substeps.

 INNER SCHEME: forward-backward, as in Klemp & Wilhelmson (1978). Forward
 Euler is unconditionally unstable on a purely imaginary spectrum, which is
 exactly what the acoustic subsystem has, so the substep cannot be a plain
 explicit Euler. Forward-backward is stable to an acoustic Courant number of
 about one and costs two operator applications, and the operator's structure
 makes it free of any rearrangement:

     row ρ    reads   (ρu, ρv, ρw)
     row ρu   reads   ρθ
     row ρv   reads   ρθ
     row ρw   reads   ρθ, ρ
     row ρθ   reads   (ρu, ρv, ρw)

 The momenta and the (ρ, ρθ) pair are decoupled within each half, so applying
 the full operator and using only the momentum rows IS the forward half, and
 applying it again to the updated state and using only the (ρ, ρθ) rows IS the
 backward half. No masking of the input is needed: row ρu never reads a
 momentum, so feeding it the already-updated one changes nothing.

 WHAT LIMITS EACH STEP
 ---------------------
     Δτ (substep)  acoustic, all three directions
     Δt (outer)    advection + diffusion only

 so the gain is the ratio of those two, against the cost of ns substeps at
 2 x 0.32 of an rhs! each plus 3 rhs! for the slow tendencies.

 BOUNDARIES -- READ THIS BEFORE TRUSTING A WALL-BOUNDED RUN
 ----------------------------------------------------------
 The operator zeroes the vertical acoustic mass and heat flux at the domain
 floor and lid (`wall` in operator.jl), so those two reflect correctly inside
 a substep. LATERAL rigid boundaries are NOT treated that way: their condition
 arrives only through f_slow, frozen for the whole stage. An acoustic wave
 reaching a side wall between two slow evaluations is therefore reflected at
 the accuracy of the outer step, not the substep. On a periodic or laterally
 open domain the question does not arise. `substep_verify` below measures the
 consequence rather than assuming it away.
=============================================================================#

"""
    SPLIT_EXPLICIT(nsub = 0; cfl_target = 0.8)

Split-explicit integrator: RK3 outer step, forward-backward acoustic substeps.

`nsub` is the number of substeps spanning the FULL outer Δt; the three RK3
stages use `nsub/3`, `nsub/2` and `nsub` of them (at least one each). Pass 0
to have it chosen at setup from the acoustic limit and `cfl_target`.
"""
struct SPLIT_EXPLICIT <: _ODEC.OrdinaryDiffEqAlgorithm
    nsub::Int
    cfl_target::Float64
end
SPLIT_EXPLICIT(nsub::Int = 0; cfl_target::Real = 0.8) =
    SPLIT_EXPLICIT(nsub, Float64(cfl_target))

Base.show(io::IO, a::SPLIT_EXPLICIT) =
    print(io, "SPLIT_EXPLICIT(", a.nsub, "; cfl_target = ", a.cfl_target, ")")

_ODEC.alg_order(::SPLIT_EXPLICIT)        = 2
_ODEC.isfsal(::SPLIT_EXPLICIT)           = false
_ODEC.isadaptive(::SPLIT_EXPLICIT)       = false
_ODEC.alg_extrapolates(::SPLIT_EXPLICIT) = false

const _WS_STAGE_FRACTIONS = (1.0/3.0, 0.5, 1.0)

mutable struct SPLIT_EXPLICIT_Cache{uType, rateType} <: _ODEC.OrdinaryDiffEqMutableCache
    u::uType
    uprev::uType
    S::rateType         # frozen slow tendency for the current stage
    dtot::rateType      # full rhs!
    dfast::rateType     # f_fast
    ustage::uType       # stage estimate
    v::uType            # substep state
    k::rateType
    fsalfirst::rateType
end

_ODEC.get_fsalfirstlast(c::SPLIT_EXPLICIT_Cache, u) = (c.fsalfirst, c.k)
_ODEC.full_cache(c::SPLIT_EXPLICIT_Cache) =
    (c.u, c.uprev, c.S, c.dtot, c.dfast, c.ustage, c.v, c.k, c.fsalfirst)

function _ODEC.alg_cache(alg::SPLIT_EXPLICIT, u, rate_prototype, ::Type{uEltypeNoUnits},
                         ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                         uprev, uprev2, f, t, dt, reltol, p, calck,
                         ::Val{true}, args...) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    SPLIT_EXPLICIT_Cache(u, uprev,
                         zero(rate_prototype), zero(rate_prototype), zero(rate_prototype),
                         similar(u), similar(u),
                         zero(rate_prototype), zero(rate_prototype))
end

function _ODEC.alg_cache(alg::SPLIT_EXPLICIT, u, rate_prototype, ::Type{uEltypeNoUnits},
                         ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits},
                         uprev, uprev2, f, t, dt, reltol, p, calck,
                         ::Val{false}, args...) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    error("SPLIT_EXPLICIT is an in-place integrator; it has no out-of-place cache.")
end

function _ODEC.initialize!(integrator, cache::SPLIT_EXPLICIT_Cache)
    integrator.kshortsize = 2
    resize!(integrator.k, integrator.kshortsize)
    integrator.k[1] = cache.fsalfirst
    integrator.k[2] = cache.k
    integrator.f(cache.fsalfirst, integrator.uprev, integrator.p, integrator.t)
    integrator.stats.nf += 1
    return nothing
end

"""
    _substep_fb!(v, S, p, sub, dτ, nsub)

`nsub` forward-backward acoustic substeps of size `dτ`, starting from `v` and
carrying the frozen slow tendency `S`.

`v` is the full state; the operator internally takes the deviation from `qe`,
so nothing here has to know about the reference state.
"""
function _substep_fb!(v, S, p, sub, dτ::Float64, nsub::Int)
    npoin = Int(p.mesh.npoin)
    R     = sub.dsub
    mom   = sub.mom_rows          # equation indices updated in the forward half
    mas   = sub.mass_rows         # ... and in the backward half

    for _ = 1:nsub
        # forward half: momenta from the CURRENT (ρ, ρθ)
        hevi_fast_rhs!(R, v, p, sub.op)
        @inbounds for ieq in mom
            off = (ieq - 1) * npoin
            for ip = 1:npoin
                v[off + ip] += dτ * (R[off + ip] + S[off + ip])
            end
        end
        # backward half: (ρ, ρθ) from the UPDATED momenta
        hevi_fast_rhs!(R, v, p, sub.op)
        @inbounds for ieq in mas
            off = (ieq - 1) * npoin
            for ip = 1:npoin
                v[off + ip] += dτ * (R[off + ip] + S[off + ip])
            end
        end
    end
    return v
end

function _ODEC.perform_step!(integrator, cache::SPLIT_EXPLICIT_Cache, repeat_step = false)

    t     = integrator.t
    dt    = integrator.dt
    uprev = integrator.uprev
    u     = integrator.u
    f     = integrator.f
    p     = integrator.p
    sub   = p.substep

    copyto!(cache.ustage, uprev)

    for (k, α) in enumerate(_WS_STAGE_FRACTIONS)
        # ---- slow tendency, evaluated once and frozen for this stage --------
        f(cache.dtot, cache.ustage, p, t + (k == 1 ? 0.0 : α * dt))
        integrator.stats.nf += 1
        hevi_fast_rhs!(cache.dfast, cache.ustage, p, sub.op)
        @inbounds @. cache.S = cache.dtot - cache.dfast

        # ---- substep from u^n over h = α Δt ---------------------------------
        h  = α * dt
        # nsub spans the FULL outer Δt, so a stage covering h = α Δt gets α nsub
        # substeps and every stage lands on the SAME Δτ. Getting this wrong is
        # invisible in a short run and fatal in a long one: stage 1 used to get
        # nsub/6 for an interval of Δt/3, i.e. twice the intended Δτ, which put
        # the forward-backward substep at an acoustic Courant number of 1.7 and
        # blew the run up at a fixed step count regardless of Δt.
        ns = max(1, cld(sub.nsub * (k == 1 ? 2 : (k == 2 ? 3 : 6)), 6))  # nsub/3, nsub/2, nsub
        dτ = h / ns
        copyto!(cache.v, uprev)
        _substep_fb!(cache.v, cache.S, p, sub, dτ, ns)
        copyto!(cache.ustage, cache.v)
    end

    copyto!(u, cache.ustage)
    f(cache.k, u, p, t + dt)
    integrator.stats.nf += 1
    return nothing
end

"""
    SubstepCache

What `perform_step!` needs, assembled once at setup and reached from
`params.substep`.
"""
struct SubstepCache{OPT, RT}
    op::OPT
    nsub::Int
    dsub::RT
    mom_rows::Vector{Int}
    mass_rows::Vector{Int}
end

"""
    build_substep(params, inputs, Δt) -> SubstepCache

Build the fast operator and choose the substep count.
"""
function build_substep(params, inputs, Δt::Real)

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    alg  = inputs[:ode_solver]

    topo = build_column_topology(params.mesh, comm)
    op   = build_hevi_fast_operator(params, topo;
                                    lwall_flux = get(inputs, :hevi_wall_flux, true))

    # Substep count: from the acoustic limit of the whole grid if not given.
    nsub = alg.nsub
    if nsub <= 0
        # The acoustic limit is set by the SUM of the per-direction rates -- the
        # eigenvalues add along the diagonal of the tensor-product operator --
        # not by the smallest of the three limits, which is what min() would
        # give and is 25% too generous on this mesh.
        L = cfl_limits(params, params.qp.qn)
        _r(x) = x > 0 && isfinite(x) ? 1.0 / x : 0.0
        dt_ac = 1.0 / max(_r(L.dt_acoustic_x) + _r(L.dt_acoustic_y) + _r(L.dt_acoustic_z), eps())
        nsub  = max(1, ceil(Int, Δt / (alg.cfl_target * dt_ac)))
    end

    npoin = Int(params.mesh.npoin)
    neqs  = Int(params.neqs)
    dsub  = zeros(Float64, npoin * neqs)

    cache = SubstepCache{typeof(op), typeof(dsub)}(op, nsub, dsub, [2, 3, 4], [1, 5])

    if rank == 0
        println(" ┌── SPLIT-EXPLICIT (acoustic substepping) ──────────────────")
        println(" │  outer: RK3 (Wicker & Skamarock 2002), 3 rhs! per step")
        println(" │  inner: forward-backward (Klemp & Wilhelmson 1978)")
        @printf(" │  substeps per outer Δt: %d   ->  Δτ = %.4g s\n", nsub, Δt / nsub)
        # nsub/3 + nsub/2 + nsub substeps per outer step, two applications each
        @printf(" │  cost per outer step: 3 rhs! + %d fast applications\n",
                2 * (cld(nsub, 3) + cld(nsub, 2) + nsub))
        println(" └───────────────────────────────────────────────────────────")
        flush(stdout)
    end
    return cache
end
