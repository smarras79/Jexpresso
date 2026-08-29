abstract type AbstractSGSModel end

#----------------------------------------------------------------------
# Smagorinsky
#----------------------------------------------------------------------
Base.@kwdef mutable struct SGS_SMAG{T <: AbstractFloat, dims1, backend, VT} <: AbstractSGSModel

    # physical constants — no defaults; set by allocator from PhysicalConst
    Pr_t::T
    Sc_t::T
    μ_mol::T
    κ_mol::T
    Ri_crit::T
    g::T
    cp::T
    Lc::T
    Ls::T
    Rvap::T
    Rair::T
    ε_ratio::T
    C_s::T
    C_s2::T
    karman::T

    # run-time configuration flags (set from inputs after allocation)
    lrichardson::Bool = false   # enable Richardson stability correction
    ltheta_eqn::Bool  = true   # true = dry θ equation; false = moist energy equation
    # Near-wall limit ℓ = min(C_sΔ, κz) on the mixing length. Off by default
    # because it needs the distance to the wall, and on a terrain-following
    # (warped) mesh the height above terrain is not available where μ_turb is
    # built — see sgs_mixing_length2 in SGS.jl. Safe for a flat lower boundary.
    lwall_damping::Bool = false

    # Distance to the wall at every node, filled once at setup when
    # lwall_damping is on (left at zero otherwise, which disables the limit
    # node by node). Precomputed rather than derived inside compute_sgs_cache!
    # so that the "which coordinate is the wall normal" question is answered
    # once, at setup, where the mesh and its zmin are in scope.
    zwall::VT  = KernelAbstractions.zeros(backend, T, dims1)

    # BASE-STATE POTENTIAL TEMPERATURE, and the flag that says whether it is
    # needed. Under :SOL_VARS_TYPE => PERT() the deck's user_primitives! hands
    # compute_sgs_cache! the PERTURBATION theta,
    #     uprimitive[5] = (u[5]+qe[5])/rho - qe[5]/qe[1],
    # which is the right thing to DIFFUSE (the base state is in hydrostatic
    # balance and must not be) but the wrong thing to build N2 from. N2 is
    # (g/theta)*dtheta/dz on the TOTAL field: dividing g by an O(0.1 K)
    # perturbation that changes sign node to node, and dropping the base-state
    # stratification from dtheta/dz entirely, turns f_Ri into a grid-scale
    # multiplier swinging between 0 (mixing off) and ~3 -- incoherent exactly
    # where a stratified LES is most sensitive to it.
    #
    # theta_base[ip] = qe[ip,5]/qe[ip,1] restores both halves: the cache adds
    # it to uprimitive[5] for the denominator and folds it into the collocation
    # derivative for dtheta/dz, so the sum IS the derivative of the total.
    # Filled once at setup (qe is static); left zeroed with lpert = false under
    # TOTAL(), where uprimitive[5] is already the total and the branch is a
    # compile-time-predictable no-op.
    lpert::Bool = false
    theta_base::VT = KernelAbstractions.zeros(backend, T, dims1)

    # per-point caches (size npoin) — zeroed at construction, filled by compute_sgs_cache!
    μ_turb::VT = KernelAbstractions.zeros(backend, T, dims1)  # turbulent viscosity
    f_Ri::VT   = KernelAbstractions.zeros(backend, T, dims1)  # Richardson correction factor
    N2::VT     = KernelAbstractions.zeros(backend, T, dims1)  # N² (dry or moist)
    # Sij components — read by les_statistics for Reynolds stress output
    S11::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S22::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S33::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S12::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S13::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S23::VT    = KernelAbstractions.zeros(backend, T, dims1)
end

#----------------------------------------------------------------------
# Vreman
#----------------------------------------------------------------------
Base.@kwdef mutable struct SGS_VREM{T <: AbstractFloat, dims1, backend, VT} <: AbstractSGSModel

    # physical constants — no defaults; set by allocator from PhysicalConst
    Pr_t::T
    Sc_t::T
    μ_mol::T
    κ_mol::T
    Ri_crit::T
    g::T
    cp::T
    Lc::T
    Ls::T
    Rvap::T
    Rair::T
    ε_ratio::T
    C_s::T
    C_s2::T
    karman::T
    C_vrem::T  # = 2.5 * C_s²

    # run-time configuration flags (set from inputs after allocation)
    lrichardson::Bool = false
    ltheta_eqn::Bool  = true
    # Near-wall limit on the mixing length, off by default and for a stronger
    # reason than in SGS_SMAG: Vreman's operator ALREADY vanishes at a wall by
    # construction -- B_β collapses faster than ‖∇u‖² in the near-wall layer --
    # which is the model's main advantage over Smagorinsky and the reason no
    # limiter was applied here originally.
    #
    # It is available anyway because "vanishes at a wall" is an asymptotic
    # statement about a RESOLVED wall layer, and an LES on a 160 x 160 m
    # horizontal cell with a wall model does not resolve one: the first node
    # sits at z = 6.9 m, deep inside the surface layer, where the asymptotics
    # have not taken over and nothing else bounds ℓ. Switching it on there
    # bounds ℓ by κz exactly as in SGS_SMAG.
    #
    # It reduces EXACTLY to the previous behaviour when false, which is the
    # default: sgs_mixing_length2 returns C_vrem*Δ² untouched.
    lwall_damping::Bool = false

    # Distance to the wall at every node, filled once at setup when
    # lwall_damping is on (left at zero otherwise, which disables the limit
    # node by node). Same contract as SGS_SMAG.zwall.
    zwall::VT  = KernelAbstractions.zeros(backend, T, dims1)

    # BASE-STATE POTENTIAL TEMPERATURE, and the flag that says whether it is
    # needed. Under :SOL_VARS_TYPE => PERT() the deck's user_primitives! hands
    # compute_sgs_cache! the PERTURBATION theta,
    #     uprimitive[5] = (u[5]+qe[5])/rho - qe[5]/qe[1],
    # which is the right thing to DIFFUSE (the base state is in hydrostatic
    # balance and must not be) but the wrong thing to build N2 from. N2 is
    # (g/theta)*dtheta/dz on the TOTAL field: dividing g by an O(0.1 K)
    # perturbation that changes sign node to node, and dropping the base-state
    # stratification from dtheta/dz entirely, turns f_Ri into a grid-scale
    # multiplier swinging between 0 (mixing off) and ~3 -- incoherent exactly
    # where a stratified LES is most sensitive to it.
    #
    # theta_base[ip] = qe[ip,5]/qe[ip,1] restores both halves: the cache adds
    # it to uprimitive[5] for the denominator and folds it into the collocation
    # derivative for dtheta/dz, so the sum IS the derivative of the total.
    # Filled once at setup (qe is static); left zeroed with lpert = false under
    # TOTAL(), where uprimitive[5] is already the total and the branch is a
    # compile-time-predictable no-op.
    lpert::Bool = false
    theta_base::VT = KernelAbstractions.zeros(backend, T, dims1)

    # per-point caches (size npoin)
    μ_turb::VT = KernelAbstractions.zeros(backend, T, dims1)
    f_Ri::VT   = KernelAbstractions.zeros(backend, T, dims1)
    N2::VT     = KernelAbstractions.zeros(backend, T, dims1)
    # Sij components. Vreman does not need them for its own eddy viscosity --
    # it works from B_β/‖∇u‖² -- but les_statistics does, and storing them here
    # is what lets the diagnostic read μ_turb straight out of this cache
    # instead of recomputing it. That distinction stopped being cosmetic when
    # :lwall_damping became available for Vreman: a recomputed μ_turb does not
    # see the near-wall limiter, so the reported subfilter stress would have
    # been the one the model did NOT apply, and overstated near the ground.
    S11::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S22::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S33::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S12::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S13::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S23::VT    = KernelAbstractions.zeros(backend, T, dims1)
end

#----------------------------------------------------------------------
# Allocators — dispatched on the existing ::SMAG / ::VREM type tags
#----------------------------------------------------------------------
# EVERY allocator MUST swallow unknown keywords, and the `::Any` fallback below
# does NOT cover them for it. `::SMAG` is MORE SPECIFIC than `::Any`, so a call
# carrying a keyword this method does not name selects this method and then
# dies in the keyword sorter -- the fallback is never reached. That is not
# hypothetical: the dynamic-SGS work added nelem/neqs/SD/C1/C2 to the single
# call site in params_setup.jl unconditionally, which left every deck using
# :visc_model => SMAG() or VREM() failing at setup with "no method matching
# allocate_SGS", while DSGS (which names them) and AV (which reaches the
# fallback) were fine.
#
# `kwargs...` here is the fix and the guard: these models genuinely do not need
# a dynamic model's parameters, so ignoring them is correct, and it makes the
# next model-specific keyword a non-event rather than a breakage.
function allocate_SGS(npoin, T, backend, PhysConst, ::SMAG; C_s = PhysConst.C_s,
                      kwargs...)
    dims1 = (Int64(npoin),)
    VT    = typeof(KernelAbstractions.zeros(backend, T, dims1))
    return SGS_SMAG{T, dims1, backend, VT}(
        Pr_t    = T(PhysConst.Pr_t),
        Sc_t    = T(PhysConst.Sc_t),
        μ_mol   = T(PhysConst.μ_mol),
        κ_mol   = T(PhysConst.κ_mol),
        Ri_crit = T(PhysConst.Ri_crit),
        g       = T(PhysConst.g),
        cp      = T(PhysConst.cp),
        Lc      = T(PhysConst.Lc),
        Ls      = T(PhysConst.Ls),
        Rvap    = T(PhysConst.Rvap),
        Rair    = T(PhysConst.Rair),
        ε_ratio = T(PhysConst.ε_ratio),
        C_s     = T(C_s),
        C_s2    = T(C_s * C_s),
        karman  = T(PhysConst.karman),
    )
end

# Fallback for :visc_model values with no dynamic SGS model (AV, DSGS, ...).
# Type-only slots rather than `_`: once a method takes keyword arguments,
# Julia lowers it into a keyword sorter that FORWARDS the positional
# arguments to an inner body method, so all-underscore names become used
# rather than discarded and the definition no longer lowers.
allocate_SGS(::Any, ::Any, ::Any, ::Any, ::Any; kwargs...) = nothing

# EVERY allocator MUST swallow unknown keywords, and the `::Any` fallback below
# does NOT cover them for it. `::SMAG` is MORE SPECIFIC than `::Any`, so a call
# carrying a keyword this method does not name selects this method and then
# dies in the keyword sorter -- the fallback is never reached. That is not
# hypothetical: the dynamic-SGS work added nelem/neqs/SD/C1/C2 to the single
# call site in params_setup.jl unconditionally, which left every deck using
# :visc_model => SMAG() or VREM() failing at setup with "no method matching
# allocate_SGS", while DSGS (which names them) and AV (which reaches the
# fallback) were fine.
#
# `kwargs...` here is the fix and the guard: these models genuinely do not need
# a dynamic model's parameters, so ignoring them is correct, and it makes the
# next model-specific keyword a non-event rather than a breakage.
function allocate_SGS(npoin, T, backend, PhysConst, ::VREM; C_s = PhysConst.C_s,
                      kwargs...)
    dims1 = (Int64(npoin),)
    VT    = typeof(KernelAbstractions.zeros(backend, T, dims1))
    return SGS_VREM{T, dims1, backend, VT}(
        Pr_t    = T(PhysConst.Pr_t),
        Sc_t    = T(PhysConst.Sc_t),
        μ_mol   = T(PhysConst.μ_mol),
        κ_mol   = T(PhysConst.κ_mol),
        Ri_crit = T(PhysConst.Ri_crit),
        g       = T(PhysConst.g),
        cp      = T(PhysConst.cp),
        Lc      = T(PhysConst.Lc),
        Ls      = T(PhysConst.Ls),
        Rvap    = T(PhysConst.Rvap),
        Rair    = T(PhysConst.Rair),
        ε_ratio = T(PhysConst.ε_ratio),
        C_s     = T(C_s),
        C_s2    = T(C_s * C_s),
        karman  = T(PhysConst.karman),
        C_vrem  = T(2.5 * C_s * C_s),
    )
end

#----------------------------------------------------------------------
# DynSGS -- residual-based dynamic SGS, 3D
#----------------------------------------------------------------------
#
# WHY DynSGS GETS A CLOSURE STRUCT IN 3D WHEN IT DOES NOT IN 1D/2D.
#
# The 1D and 2D DynSGS paths (SGS.jl, rhs.jl) bypass `params.sgs` entirely:
# they compute one coefficient per element, pack it into `visc_coeff_dsgs` and
# hand that to a scalar `visc_coeff[ieq]*grad(primitive)` kernel. That is fine
# for a shock tube and for the rising bubble, and it is wrong for an
# atmospheric LES, because the 3D viscous kernel is not a scalar Laplacian --
# it is the full stress tensor
#
#     tau_ij = mu (du_i/dx_j + du_j/dx_i) - (2/3) mu div(u) delta_ij
#
# and it is reached only through `_expansion_visc!(..., sgs::AbstractSGSModel,
# ..., NSD_3D, ...)`, which reads its coefficient from `SGS_diffusion` on a
# closure struct. Everything else that has to agree with the viscous term
# reads the same struct:
#
#   cfl_limits          (hevi/cfl_diagnostics.jl)  the viscous row of the CFL
#                                                  table, via sgs.mu_turb
#   vdiff_refresh!      (hevi/vdiffusion.jl)       the IMPLICIT vertical
#                                                  diffusion coefficient, via
#                                                  SGS_diffusion on this struct
#   _fill_sgs_cached!   (io/les_statistics.jl)     the reported subfilter
#                                                  stress, via mu_turb and Sij
#
# So making DynSGS an AbstractSGSModel is not packaging: it is what buys the
# correct 3D stress tensor, a truthful CFL report, working :implicit_vdiff and
# correct LES statistics, none of which the 1D/2D DynSGS path has.
#
# WHAT IS DIFFERENT FROM SMAG/VREM. Only where mu_turb comes from. Smagorinsky
# builds it node by node from the local strain rate; DynSGS builds ONE
# coefficient per element from the L-infinity norm of the PDE residual over
# that element and global normalising scales, which cannot be computed inside a
# per-element loop. So it is filled in two stages:
#
#   compute_dsgs_viscosity!(sgs, NSD_3D, ...)   once per RHS call, before the
#                                               element loop -> sgs.nu_el,
#                                               sgs.mu_el (per element)
#   compute_sgs_cache!(sgs, ..., iel, ...)      per element, inside the loop ->
#                                               broadcasts mu_el[iel] onto the
#                                               element's nodes in mu_turb and
#                                               fills Sij for les_statistics
#
# The broadcast is last-writer-wins on DSS-shared nodes, which is exactly the
# convention broadcast_dsgs_to_nodes! already uses for the VTK fields. It is
# not a defect here: compute_sgs_cache! runs immediately before the ieq loop
# for the SAME element, so every value the viscous kernel reads inside that
# element is that element's own.
#
Base.@kwdef mutable struct SGS_DSGS{T <: AbstractFloat, dims1, dimsE, backend, VT} <: AbstractSGSModel

    # physical constants — no defaults; set by allocator from PhysicalConst
    Pr_t::T
    Sc_t::T
    μ_mol::T
    κ_mol::T
    Ri_crit::T
    g::T
    cp::T
    Lc::T
    Ls::T
    Rvap::T
    Rair::T
    ε_ratio::T
    C_s::T
    C_s2::T
    karman::T

    # DynSGS coefficients, Marras/Nazarov eq. (9): mu_res = C1 Delta^2 max_i
    # ‖R_i‖/‖q_i-<q_i>‖ capped at mu_max = C2 Delta (‖v‖+c). Deck-settable as
    # :dsgs_C1 / :dsgs_C2; the paper's values are 1.0 and 0.5.
    C1::T
    C2::T

    # run-time configuration flags (set from inputs after allocation)
    lrichardson::Bool = false
    ltheta_eqn::Bool  = true
    lwall_damping::Bool = false

    # ADD SMAGORINSKY ON TOP OF THE RESIDUAL VISCOSITY (:dsgs_add_smagorinsky).
    #
    # Off by default, because "parameter-free" is the point of the model and
    # because that is what the 1D/2D paths do. It exists because a residual
    # sensor and a surface layer want different things and this deck has both.
    #
    # DynSGS is a STABILISATION: it measures how badly the discrete solution
    # fails to satisfy the PDE and puts viscosity exactly there, which is
    # nearly nothing in a well-resolved region. In the surface layer of a
    # wall-modelled PBL the subfilter stress carries MOST of the momentum flux
    # and has to be there whether or not the solution is locally smooth --
    # that is a physics requirement (the log law), not a stabilisation one, and
    # no residual sensor produces it. With this on, mu_turb = mu_smag + mu_dsgs
    # so the closure keeps its wall behaviour and DynSGS adds dissipation only
    # where the residual says the discretisation needs it.
    #
    # Run it off first: whether the residual term alone is enough is exactly
    # the question this deck exists to answer.
    ladd_smagorinsky::Bool = false

    # Distance to the wall at every node; same contract as SGS_SMAG.zwall, and
    # read only when ladd_smagorinsky is on (the residual viscosity has no
    # mixing length to limit).
    zwall::VT  = KernelAbstractions.zeros(backend, T, dims1)

    # BASE-STATE POTENTIAL TEMPERATURE, and the flag that says whether it is
    # needed. Under :SOL_VARS_TYPE => PERT() the deck's user_primitives! hands
    # compute_sgs_cache! the PERTURBATION theta,
    #     uprimitive[5] = (u[5]+qe[5])/rho - qe[5]/qe[1],
    # which is the right thing to DIFFUSE (the base state is in hydrostatic
    # balance and must not be) but the wrong thing to build N2 from. N2 is
    # (g/theta)*dtheta/dz on the TOTAL field: dividing g by an O(0.1 K)
    # perturbation that changes sign node to node, and dropping the base-state
    # stratification from dtheta/dz entirely, turns f_Ri into a grid-scale
    # multiplier swinging between 0 (mixing off) and ~3 -- incoherent exactly
    # where a stratified LES is most sensitive to it.
    #
    # theta_base[ip] = qe[ip,5]/qe[ip,1] restores both halves: the cache adds
    # it to uprimitive[5] for the denominator and folds it into the collocation
    # derivative for dtheta/dz, so the sum IS the derivative of the total.
    # Filled once at setup (qe is static); left zeroed with lpert = false under
    # TOTAL(), where uprimitive[5] is already the total and the branch is a
    # compile-time-predictable no-op.
    lpert::Bool = false
    theta_base::VT = KernelAbstractions.zeros(backend, T, dims1)

    # per-ELEMENT DynSGS coefficients (size nelem), filled once per RHS call by
    # compute_dsgs_viscosity!(::SGS_DSGS, ::NSD_3D).
    #   ν_el  kinematic  [m²/s]  — the model's own output, what to compare
    #                              against the C2·Δ·(|v|+c) cap
    #   μ_el  dynamic    [Pa·s]  — ρ̄_el·ν_el, what the stress tensor wants
    ν_el::VT = KernelAbstractions.zeros(backend, T, dimsE)
    μ_el::VT = KernelAbstractions.zeros(backend, T, dimsE)

    # Domain-reduction scratch, one entry per equation. Plain host Vectors and
    # not backend arrays on purpose: they are the MPI.Allreduce! buffers and
    # are read scalar-wise by the element loop, so they must be host-resident.
    #   avg    <q'_i>                 domain mean of the PERTURBATION
    #   denom  ‖q'_i - <q'_i>‖_inf    the normalising scale of eq. (9)
    #   scale  <|q_i|>                domain mean of the TOTAL, which is what
    #                                 the denominator FLOORS are built from
    avg::Vector{T}
    denom::Vector{T}
    scale::Vector{T}
    #   raw     the MEASURED ‖q_i - <q_i>‖ before any floor is applied
    #   floors  the floor for each slot
    # Kept apart so the floor can be used as a GATE -- "this field has no
    # resolved scale yet, leave its residual out of the max" -- rather than as
    # a value. See the block that sets denom[ieq] = Inf in SGS.jl.
    raw::Vector{T}
    floors::Vector{T}

    # per-point caches (size npoin) — same names and meanings as SGS_SMAG, so
    # every consumer of a closure struct works unchanged.
    μ_turb::VT = KernelAbstractions.zeros(backend, T, dims1)
    f_Ri::VT   = KernelAbstractions.zeros(backend, T, dims1)
    N2::VT     = KernelAbstractions.zeros(backend, T, dims1)
    S11::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S22::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S33::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S12::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S13::VT    = KernelAbstractions.zeros(backend, T, dims1)
    S23::VT    = KernelAbstractions.zeros(backend, T, dims1)
end

# DynSGS gets a closure struct ONLY in 3D. The 1D and 2D paths run through
# viscous_rhs_el!'s own `params.VT == DSGS()` branch and never look at
# params.sgs; allocating one there would put npoin-sized arrays behind a
# `params.sgs isa AbstractSGSModel` test that several diagnostics use as
# "there is a closure to read", and they would read zeros.
function allocate_SGS(npoin, T, backend, PhysConst, ::DSGS;
                      C_s = PhysConst.C_s, nelem = 0, neqs = 5, SD = nothing,
                      C1 = 1.0, C2 = 0.5)
    SD === NSD_3D() || return nothing
    dims1 = (Int64(npoin),)
    dimsE = (Int64(nelem),)
    VT    = typeof(KernelAbstractions.zeros(backend, T, dims1))
    return SGS_DSGS{T, dims1, dimsE, backend, VT}(
        avg     = zeros(T, Int64(neqs)),
        denom   = zeros(T, Int64(neqs)),
        scale   = zeros(T, Int64(neqs)),
        raw     = zeros(T, Int64(neqs)),
        floors  = zeros(T, Int64(neqs)),
        Pr_t    = T(PhysConst.Pr_t),
        Sc_t    = T(PhysConst.Sc_t),
        μ_mol   = T(PhysConst.μ_mol),
        κ_mol   = T(PhysConst.κ_mol),
        Ri_crit = T(PhysConst.Ri_crit),
        g       = T(PhysConst.g),
        cp      = T(PhysConst.cp),
        Lc      = T(PhysConst.Lc),
        Ls      = T(PhysConst.Ls),
        Rvap    = T(PhysConst.Rvap),
        Rair    = T(PhysConst.Rair),
        ε_ratio = T(PhysConst.ε_ratio),
        C_s     = T(C_s),
        C_s2    = T(C_s * C_s),
        karman  = T(PhysConst.karman),
        C1      = T(C1),
        C2      = T(C2),
    )
end
