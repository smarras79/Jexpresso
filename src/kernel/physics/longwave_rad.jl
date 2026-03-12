# ─────────────────────────────────────────────────────────────────────────────
# LONGWAVE
# ─────────────────────────────────────────────────────────────────────────────

"""
    LWParams

Configuration for the longwave gray radiation calculation.

# Fields
- `ε_surface`: surface emissivity ∈ [0,1]. 1.0 → perfect blackbody.
               Typical land: 0.95–0.98, ocean: 0.97–0.99.
- `T_space`: effective brightness temperature of cold space (K). Use 0.0 to
             impose zero inflow from above, or 2.7 for cosmic background.
             The latter is negligible for terrestrial heating rates.
"""

struct LWParams
    ε_surface :: Float64
    T_space   :: Float64
end



"""
    user_rhs_longwave(x, y, z, θ, ϕ, ip, κ, atmos_data) → Float64

Volumetric source term for the longwave RTE.

Thermal emission is isotropic and proportional to the absorption coefficient:

    Q(x, Ω) = κ_abs(x) B(T(x))

where B(T) = σ_SB T⁴ / π is the gray Planck radiance (W m⁻² sr⁻¹).

The scattering source ∫ Φ I' dΩ' is handled by the LHS and does not
appear here.

# Arguments
- `ip`: spatial node index into `κ` and `atmos_data`
- `κ`: absorption coefficient array (m⁻¹) from `atmos_to_rad_longwave`
- `atmos_data`: atmospheric state struct with field `t_lay`

# Notes
κ_abs (absorption only) is used here, not κ_ext (total extinction). These
differ when scattering is non-negligible. Make sure you are passing the
absorption component returned from `atmos_to_rad_longwave`, not the sum
κ + σ.
"""


function user_rhs_longwave(x, y, z, θ, ϕ, ip, κ, atmos_data)

    σ_SB  = 5.670374419e-8   # W m⁻² K⁻⁴
    T     = atmos_data.t_lay[ip]
    κ_abs = κ[ip]

    # Gray Planck radiance (W m⁻² sr⁻¹)
    B_T = σ_SB * T^4 / π

    return κ_abs * B_T
end


"""
    user_rad_bc_longwave(x, y, z, θ, ϕ, nx_n, ny_n, nz_n,
                         ip, atmos_data, bdy, lw) → Float64

Inflow boundary radiance I_bc(x,Ω) for the longwave RTE on a 3D box domain.

Called only at boundary nodes where Ω·n̂ < 0.

# Boundary treatment

**Top face (z = zmax):** Cold space — zero downwelling longwave. If
`lw.T_space > 0` a small cosmic background term is included:

    I_TOA = ε_space × σ_SB T_space⁴ / π

**Bottom face (z = zmin):** Gray body surface emission in all upward directions:

    I_sfc(Ω) = ε_sfc × σ_SB T_sfc⁴ / π    for Ωz > 0

The surface temperature is taken from `atmos_data.t_lay[ip]` at the lowest
model level. Downward-going directions at the bottom boundary (Ωz < 0) are
not inflow and are never called here.

**Lateral faces:** The domain is periodic in x and y. These branches are
included for completeness but should not be reached in a periodic setup.
They return the local blackbody radiance as a neutral non-reflecting wall.

# Arguments
- `ip`: spatial node index
- `atmos_data`: atmospheric state struct with field `t_lay`
- `bdy`: boundary predicates from `make_boundary_predicates`
- `lw`: `LWParams` struct
"""
function user_rad_bc_longwave(x, y, z, θ, ϕ, nx_n, ny_n, nz_n,
                               ip, atmos_data, bdy, lw)

    σ_SB = 5.670374419e-8

    # Route based on which face the best normal belongs to.
    # Do NOT add secondary Ωz checks here — the caller already selected the
    # correct face normal for this direction via best_dot, so we trust it.
    if bdy.is_top(z)
        return σ_SB * lw.T_space^4 / π

    elseif bdy.is_bottom(z)
        T_sfc = atmos_data.t_lay[ip]
        return lw.ε_surface * σ_SB * T_sfc^4 / π

    else
        # Lateral face — periodic in x,y so this should rarely be reached,
        # but emit as a blackbody wall for robustness
        T_wall = atmos_data.t_lay[ip]
        return σ_SB * T_wall^4 / π
    end
end