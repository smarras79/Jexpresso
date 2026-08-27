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
