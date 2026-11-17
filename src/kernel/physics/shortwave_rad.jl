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

    sw.μ₀ > 0.0 || return 0.0  # no direct beam at horizon/night
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

        sw.μ₀ > 0.0 || return 0.0  # no direct beam at horizon/night

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

"""
    sw_single_scatter_column(prof, S0, μ0, φ0, θ, ϕ; albedo=0.0) -> Vector{Float64}

Single-scattering diffuse radiance of a horizontally homogeneous (plane-parallel) column, for one direction Ω = (θ, ϕ),
at every level of `prof` (from build_sw_column_optics: z, τ from TOA, ω, g; TOA first).

  downward (cosθ < 0):  I(τ) = ∫₀^τ   J(τ′) e^{-τ′/μ₀} e^{-(τ-τ′)/μ} dτ′/μ,      J = ω S₀ Φ(Ω,Ω_sun; g)
  upward   (cosθ ≥ 0):  I(τ) = ∫_τ^τ* J(τ′) e^{-τ′/μ₀} e^{-(τ′-τ)/μ} dτ′/μ + (A/π) μ₀ S₀ e^{-τ*/μ₀} e^{-(τ*-τ)/μ}

μ = |cosθ|; vacuum top; Lambertian surface reflecting the direct beam. J is taken LINEAR in τ inside each layer and the
exponential integrals are done exactly (closed forms with the e^{-Δ/μ} factor folded in: no overflow as μ → 0, where
I → ω S₀ Φ e^{-τ/μ₀}; Taylor limits when |aΔ| is tiny).
"""
function sw_single_scatter_column(prof, S0, μ0, φ0, θ, ϕ; albedo = 0.0)
    n = length(prof.z); I = zeros(n)
    μ0 > 0 || return I
    θs = π - acos(clamp(μ0, 0.0, 1.0))
    Ωz = cos(θ); μ = max(abs(Ωz), 1e-12)
    τ = prof.τ
    J = [prof.ω[k] * S0 * phase_function_hg(θ, θs, ϕ, φ0, prof.g[k]) for k = 1:n]
    if Ωz < 0
        a = 1/μ - 1/μ0
        for k = 1:n-1
            Δ = τ[k+1] - τ[k]; m = Δ > 0 ? (J[k+1] - J[k])/Δ : 0.0
            eΔ = exp(-Δ/μ0)
            if abs(a*Δ) < 1e-6                                        # Taylor limits
                A0 = Δ*eΔ*(1 - a*Δ/2); A1 = Δ^2*eΔ*(0.5 - a*Δ/3)
            else
                A0 = eΔ*(-expm1(-a*Δ))/a
                A1 = (Δ*eΔ - A0)/a
            end
            I[k+1] = I[k]*exp(-Δ/μ) + exp(-τ[k]/μ0)/μ * (J[k]*A0 + m*A1)
        end
    else
        b = 1/μ + 1/μ0
        I[n] = albedo/π * μ0*S0*exp(-τ[n]/μ0)
        for k = n-1:-1:1
            Δ = τ[k+1] - τ[k]; m = Δ > 0 ? (J[k+1] - J[k])/Δ : 0.0
            B0 = -expm1(-b*Δ)/b
            B1 = (B0 - Δ*exp(-b*Δ))/b
            I[k] = I[k+1]*exp(-Δ/μ) + exp(-τ[k]/μ0)/μ * (J[k]*B0 + m*B1)
        end
    end
    return I
end

"""
    SWLateralSingleScatter(prof, S0, μ0, φ0, albedo)

Lateral-boundary inflow radiance for the diffuse SW solve: the plane-parallel single-scattering radiance of the
reference column (sw_single_scatter_column), computed once per direction and interpolated linearly in z.
Call as `lat(z, θ, ϕ)`.
"""
struct SWLateralSingleScatter
    prof
    S0::Float64; μ0::Float64; φ0::Float64; albedo::Float64
    cache::Dict{Tuple{Float64,Float64},Vector{Float64}}
end
SWLateralSingleScatter(prof, S0, μ0, φ0, albedo) =
    SWLateralSingleScatter(prof, S0, μ0, φ0, albedo, Dict{Tuple{Float64,Float64},Vector{Float64}}())

function (lat::SWLateralSingleScatter)(z, θ, ϕ)
    Icol = get!(lat.cache, (θ, ϕ)) do
        sw_single_scatter_column(lat.prof, lat.S0, lat.μ0, lat.φ0, θ, ϕ; albedo = lat.albedo)
    end
    zp = lat.prof.z                                       # descending
    z >= zp[1] && return Icol[1]
    z <= zp[end] && return Icol[end]
    k = clamp(searchsortedfirst(zp, z, rev=true) - 1, 1, length(zp)-1)
    t = (zp[k] - z) / (zp[k] - zp[k+1])
    return Icol[k] + t*(Icol[k+1] - Icol[k])
end


"""
    build_sw_lateral_plane_parallel(prof, S0, μ0, φ0, albedo; N=16, nphi=144) -> SWLateralPlaneParallel

Lateral-boundary inflow for the diffuse SW solve: the full plane-parallel radiance of the reference column,

    I_lat(z, Ω) = I_ss(z, Ω) + ΔI(z, μ),      ΔI(z, μ_i) = I_1D(z, μ_i) - ⟨I_ss⟩_φ(z, μ_i)

I_ss is the exact single-scattering radiance (sw_single_scatter_column: keeps the sharp solar-azimuth dependence);
ΔI is the multiply scattered part, from the 1D reference solver (SWRef1D.solve_column, azimuth-averaged, at its
Gauss μ_i) minus the azimuth mean of I_ss at the same μ_i, taken independent of azimuth. ΔI is interpolated linearly
in z and in μ within each hemisphere (constant beyond the outermost Gauss nodes). μ here is the 1D convention:
μ > 0 = downward travel = -cosθ. NOTE: the 1D solve reflects the TOTAL downward flux at a Lambertian surface; the 3D
bottom condition currently reflects only the direct beam, so the two agree only for albedo = 0.
"""
struct SWLateralPlaneParallel
    ss::SWLateralSingleScatter
    μp::Vector{Float64}             # Gauss μ on (0,1), ascending
    ΔI_dn::Matrix{Float64}          # (N, nlev) downward-travelling, at μp
    ΔI_up::Matrix{Float64}          # (N, nlev) upward-travelling, at μp
    z::Vector{Float64}              # levels, descending
end

function build_sw_lateral_plane_parallel(prof, S0, μ0, φ0, albedo; N::Int = 16, nphi::Int = 144)
    ss = SWLateralSingleScatter(prof, S0, μ0, φ0, albedo)
    κa = (1 .- prof.ω) .* prof.e; σ = prof.ω .* prof.e
    r = SWRef1D.solve_column(prof.z, κa, σ, prof.g; S0 = S0, μ0 = μ0, A = albedo, N = N, deltam = true, dτmax = 0.05)
    nl = length(prof.z)
    length(r.z) == nl || error("SWLateralPlaneParallel: level mismatch")
    ΔI = zeros(2N, nl)
    φs = range(0, 2π, length = nphi + 1)[1:nphi]
    for i = 1:2N
        θ = acos(-r.μ[i])                                  # 1D μ > 0 (down) -> Ωz = cosθ = -μ
        m = zeros(nl)
        for ϕ in φs; m .+= sw_single_scatter_column(prof, S0, μ0, φ0, θ, ϕ; albedo = albedo); end
        ΔI[i, :] .= r.Idiff[i, :] .- m ./ nphi
    end
    return SWLateralPlaneParallel(ss, r.μ[1:N], ΔI[1:N, :], ΔI[N+1:2N, :], collect(prof.z))
end

function (lat::SWLateralPlaneParallel)(z, θ, ϕ)
    μ = -cos(θ)
    T = μ >= 0 ? lat.ΔI_dn : lat.ΔI_up
    a = abs(μ); μp = lat.μp; N = length(μp)
    # μ interpolation weights (constant beyond the end nodes)
    if a <= μp[1];      i1, i2, t = 1, 1, 0.0
    elseif a >= μp[N];  i1, i2, t = N, N, 0.0
    else
        i2 = searchsortedfirst(μp, a); i1 = i2 - 1; t = (a - μp[i1]) / (μp[i2] - μp[i1])
    end
    # z interpolation (levels descending)
    zp = lat.z
    if z >= zp[1];       k1, k2, s = 1, 1, 0.0
    elseif z <= zp[end]; k1, k2, s = length(zp), length(zp), 0.0
    else
        k1 = clamp(searchsortedfirst(zp, z, rev = true) - 1, 1, length(zp) - 1); k2 = k1 + 1
        s = (zp[k1] - z) / (zp[k1] - zp[k2])
    end
    at(k) = T[i1, k] + t * (T[i2, k] - T[i1, k])
    return lat.ss(z, θ, ϕ) + at(k1) + s * (at(k2) - at(k1))
end
