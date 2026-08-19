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
    # Vreman's operator already vanishes at a wall by construction, so no
    # near-wall length-scale limit is applied here (unlike SGS_SMAG).

    # per-point caches (size npoin)
    μ_turb::VT = KernelAbstractions.zeros(backend, T, dims1)
    f_Ri::VT   = KernelAbstractions.zeros(backend, T, dims1)
    N2::VT     = KernelAbstractions.zeros(backend, T, dims1)
    # Vreman computes B_β/‖∇u‖² internally — no Sij components stored
end

#----------------------------------------------------------------------
# Allocators — dispatched on the existing ::SMAG / ::VREM type tags
#----------------------------------------------------------------------
function allocate_SGS(npoin, T, backend, PhysConst, ::SMAG; C_s = PhysConst.C_s)
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

allocate_SGS(_, _, _, _, ::Any; kwargs...) = nothing

function allocate_SGS(npoin, T, backend, PhysConst, ::VREM; C_s = PhysConst.C_s)
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
