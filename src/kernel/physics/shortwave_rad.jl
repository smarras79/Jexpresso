# ─────────────────────────────────────────────────────────────────────────────
# SHORTWAVE
# ─────────────────────────────────────────────────────────────────────────────

"""
    SWParams

Configuration for the shortwave gray radiation calculation.

# Fields
- `S₀_flux`: solar constant at TOA (W m⁻²), the total flux integrated over
             the solar disk (nominally 1361 W m⁻² for Earth).
- `μ₀`: cosine of solar zenith angle. μ₀ = 1 → overhead sun, μ₀ = 0 → horizon.
- `φ₀`: solar azimuth angle (rad), measured from the x-axis. φ₀ = 0 → sun in
        the x-z plane.
- `δ_beam`: Gaussian angular half-width of the solar beam (rad). The beam is
            spread over this width to be resolvable on the finite angular mesh.
            Should be ≥ 2–3× the angular element size. Default 0.05 rad ≈ 3°.

# Notes on δ_beam and flux conservation
The Gaussian spread conserves flux: the TOA radiance field integrates to
S₀_flux × μ₀ over the downward hemisphere when δ_beam is small relative to π.
Verify this on your angular mesh before production runs:

    flux_check = ∑_{Ω·ẑ < 0} I_TOA(Ω) |Ωz| dΩ  ≈  S₀_flux × μ₀
"""

struct SWParams
    S₀_flux    :: Float64   # W m⁻²
    μ₀         :: Float64   # cos(SZA)
    φ₀         :: Float64   # solar azimuth, rad
    δ_beam     :: Float64   # angular beam width, rad
    z_prof     :: Vector{Float64}
    τ_from_TOA :: Vector{Float64}
end


"""
    user_rhs_shortwave(x, y, z, θ, ϕ) → Float64

Volumetric source term for the shortwave RTE.

There is no thermal emission in the shortwave. All solar forcing enters through
the top-of-atmosphere inflow boundary condition. Multiple scattering is handled
implicitly by the LHS scattering integral.
"""
@inline function user_rhs_shortwave(x, y, z, θ, ϕ)
    return 0.0
end


# ── Corrected BC normalisation ────────────────────────────────────────────────
# I_beam is specific intensity (W m⁻² sr⁻¹).
# The solar disk subtends solid angle ~ π δ² (Gaussian approximation).
# The flux through a horizontal surface is:
#   F = ∫ I |cosθ| dΩ = I₀ × |cos(θ_sun)| × π δ² = S₀ × μ₀  ✓
# So the peak radiance is S₀ / (π δ²), with NO μ₀ in the numerator.

function user_rad_bc_shortwave(x, y, z, θ, ϕ, nx_n, ny_n, nz_n, bdy, sw)

    Ωx = sin(θ)*cos(ϕ); Ωy = sin(θ)*sin(ϕ); Ωz = cos(θ)
    Ω_dot_n = nx_n*Ωx + ny_n*Ωy + nz_n*Ωz
    Ω_dot_n < -1e-13 || return 0.0

    if bdy.is_top(z)
       

        θ_sun = π - acos(clamp(sw.μ₀, 0.0, 1.0))
        Δθ    = θ - θ_sun
        Δϕ    = ϕ - sw.φ₀
        Δϕ    = Δϕ - 2π * round(Δϕ / (2π))
        ang_dist² = Δθ^2 + (sin(θ_sun) * Δϕ)^2

        # Peak radiance: S₀/(π δ²), no μ₀ factor here
        # The μ₀ projection enters when computing flux via ∫ I |cosθ| dΩ
        I_beam = (sw.S₀_flux / (π * sw.δ_beam^2)) *
                 exp(-ang_dist² / sw.δ_beam^2)

        return I_beam

    elseif !bdy.is_bottom(z)
        # ── Lateral face: attenuated direct beam ──────────────────────────
        # Only directions close to the solar beam direction receive inflow.
        # A direction Ω at a lateral face is inflow (Ω·n < 0) and carries
        # the solar beam intensity attenuated from TOA down to height z
        # along the beam path.
        #
        # The optical path from TOA (z=zmax) to height z along the beam
        # direction with cosine μ₀ = sw.μ₀ is:
        #
        #   τ(z) = ∫_z^{zmax} κ_ext(z') dz' / μ₀
        #
        # The attenuated beam radiance at height z on the lateral face is:
        #
        #   I_lateral(Ω) = I_TOA(Ω) × exp(-τ(z) / μ₀)
        #
        # where τ(z) is the vertical optical depth from z to TOA.
        # κ_ext_col is a precomputed array of cumulative optical depth
        # from TOA as a function of z for this column.

        sw.μ₀ > 0.0 || return 0.0

        θ_sun     = π - acos(clamp(sw.μ₀, 0.0, 1.0))
        Δθ        = θ - θ_sun
        Δϕ        = ϕ - sw.φ₀
        Δϕ        = Δϕ - 2π * round(Δϕ / (2π))
        ang_dist² = Δθ^2 + (sin(θ_sun) * Δϕ)^2

        # Beam radiance at TOA
        I_TOA = (sw.S₀_flux / (π * sw.δ_beam^2)) * exp(-ang_dist² / sw.δ_beam^2)

        # Beer-Lambert attenuation from TOA to height z
        # τ_col(z) = cumulative vertical optical depth from zmax down to z
        # Divide by μ₀ to get the slant path optical depth
        τ_z   = interp_optical_depth(z, sw.z_prof, sw.τ_from_TOA)
        I_lat = I_TOA * exp(-τ_z / sw.μ₀)

        return I_lat

    else

        return 0.0
    end
end


"""
    user_rhs_shortwave_diffuse(x, y, z, θ, ϕ, ip, F_dir, κ_ext_data, σ_data,
                                sw, atmos_data)

RHS source term for the diffuse shortwave RTE.

The direct beam scatters into the diffuse field at each point:

    S_dir(x, Ω) = σ(x) × Φ(Ω, Ω_sun) × F_dir(x) / μ₀

where:
- σ(x) is the scattering coefficient
- Φ(Ω, Ω_sun) is the phase function evaluated between direction Ω
  and the solar beam direction Ω_sun
- F_dir(x)/μ₀ = S₀ exp(-τ/μ₀) is the direct beam radiance
  (flux divided by μ₀ recovers the beam radiance along Ω_sun)
"""
function user_rhs_shortwave_diffuse(x, y, z, θ, ϕ, ip,
                                     F_dir, σ_data, sw, g_eff)

    σ = σ_data[ip]
    σ < 1e-30 && return 0.0   # no scattering — no source

    # Solar beam direction
    θ_sun = π - acos(clamp(sw.μ₀, 0.0, 1.0))
    φ_sun = sw.φ₀

    # Phase function between current direction and solar beam direction
    Φ_val = user_scattering_functions(θ, θ_sun, ϕ, φ_sun, g_eff)

    # Direct beam radiance along Ω_sun = F_dir / μ₀
    # (F_dir is the flux = radiance × μ₀ × solid_angle_of_beam ≈ radiance × μ₀)
    I_dir_beam = F_dir[ip] / sw.μ₀

    # Scattering source: energy scattered from direct beam into direction Ω
    # The 4π factor accounts for the normalisation convention of Φ
    return σ * Φ_val * I_dir_beam
end


"""
    user_rad_bc_shortwave_diffuse(x, y, z, θ, ϕ, nx_n, ny_n, nz_n, bdy, sw)

Boundary conditions for the diffuse shortwave problem.

All inflow BCs are zero — the direct beam is handled analytically via
F_dir, and there is no diffuse inflow from outside the domain.

Surface reflection (if needed) can be added here as:
    I_reflected = albedo × F_dir_surface / π   (Lambertian surface)
but for the first implementation zero is correct.
"""
function user_rad_bc_shortwave_diffuse(x, y, z, θ, ϕ,
                                        bdy, sw, F_dir_ip,
                                        τ_z, sw_ω₀_lateral, g_eff;
                                        surface_albedo = 0.15)
    if bdy.is_bottom(z)
        # Lambertian surface reflection of direct beam
        return surface_albedo * F_dir_ip / π

    elseif bdy.is_top(z)
        # Zero diffuse inflow at TOA
        return 0.0

    else
        # Lateral boundary: approximate diffuse inflow as the single-scatter
        # source integrated from TOA to this height along the beam.
        # This is the diffuse radiance that the surrounding atmosphere would
        # produce at this height — it accounts for the fact that the domain
        # is only ~5 mean free paths wide and lateral fluxes are significant.
        #
        # Single-scatter approximation for lateral diffuse inflow:
        #   I_lat(Ω) = ∫₀^{τ_z} σ Φ(Ω,Ω_sun) I_dir(τ') / κ_ext dτ'
        #
        # For a homogeneous atmosphere this simplifies to:
        #   I_lat(Ω) = σ/κ_ext × Φ(Ω,Ω_sun) × (S₀/μ₀) × 
        #              ∫₀^{τ_z} exp(-τ'/μ₀) dτ'
        #            = σ/κ_ext × Φ(Ω,Ω_sun) × S₀ × (1 - exp(-τ_z/μ₀))
        #
        # This gives a physically consistent lateral inflow that matches
        # the interior scattering source at the boundary.

        θ_sun = π - acos(clamp(sw.μ₀, 0.0, 1.0))
        Φ_val = user_scattering_functions(θ, θ_sun, ϕ, sw.φ₀, g_eff)

        # ω₀ = σ/κ_ext at this boundary node — pass as argument or
        # compute from σ_ip, κ_ext_ip passed in
        ω₀_val = sw_ω₀_lateral   # domain-mean single scatter albedo,
                                   # precomputed as mean(σ_sw ./ κ_ext_sw)

        I_lat = ω₀_val * Φ_val * sw.S₀_flux * (1.0 - exp(-τ_z / sw.μ₀))

        return I_lat
    end
end