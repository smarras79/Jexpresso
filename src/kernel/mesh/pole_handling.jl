"""
    pole_shifted_theta_range(nelem_θ, ngl; pole_fraction=0.01)

Compute a θ range [θ_min, θ_max] that excludes the poles by a margin
proportional to the angular element size, ensuring LGL endpoint nodes
never land at θ=0 or θ=π where sin(θ)=0.

# Arguments
- `nelem_θ`: number of elements in the θ direction
- `ngl`: number of GL points per element (polynomial order + 1)
- `pole_fraction`: fraction of one element width to exclude at each pole.
                   Default 0.01 (1% of element width). See Notes.

# Returns
- `θ_min`: shifted lower bound (> 0)
- `θ_max`: shifted upper bound (< π)

# Notes
The element width is Δθ = π / nelem_θ. The pole margin is:

    δ = pole_fraction × Δθ

so the mesh covers [δ, π-δ]. The solid angle of the two excluded polar
caps is:

    ΔΩ = 2 × 2π × (1 - cos(δ)) ≈ 2π δ²   (for small δ)

For nelem_θ=4, ngl=5, pole_fraction=0.01:
    Δθ = π/4 ≈ 0.785 rad
    δ  = 0.01 × 0.785 ≈ 0.00785 rad
    ΔΩ / 4π ≈ δ²/2 ≈ 3×10⁻⁵  (0.003% of total solid angle — negligible)
    sin(θ_min) = sin(0.00785) ≈ 0.00785   (Je_min / Je_typical ~ 0.01)

pole_fraction=0.01 is recommended: it keeps the mass matrix condition
number well-behaved while the flux error from missing polar caps is below
any physically meaningful threshold. Do not use pole_fraction < 1e-4 as
sin(θ_min) approaches machine epsilon and the mass matrix becomes
ill-conditioned again.
"""
function pole_shifted_theta_range(nelem_θ, ngl; pole_fraction=0.01)
    Δθ    = π / nelem_θ
    δ     = pole_fraction * Δθ
    θ_min = δ
    θ_max = π - δ

    # Solid angle of excluded caps as a fraction of 4π
    ΔΩ_frac = (1 - cos(δ))   # = ΔΩ_both_caps / (4π), since 2×2π×(1-cosδ)/(4π)

    @info "Pole exclusion:"
    @info "  Element width Δθ    = $(round(Δθ, digits=6)) rad  " *
          "($(round(rad2deg(Δθ), digits=2))°)"
    @info "  Pole margin δ       = $(round(δ, digits=8)) rad  " *
          "($(round(rad2deg(δ), digits=4))°)"
    @info "  θ range             = [$(round(θ_min,digits=8)), $(round(θ_max,digits=8))]"
    @info "  sin(θ_min)          = $(round(sin(θ_min), digits=8))  " *
          "(Je_pole / Je_equator ≈ $(round(sin(θ_min), digits=4)))"
    @info "  Excluded solid angle: $(round(ΔΩ_frac*100, digits=4))% of 4π  " *
          "(both caps combined)"

    if sin(θ_min) < 1e-6
        @warn "sin(θ_min) = $(sin(θ_min)) is dangerously small. " *
              "Increase pole_fraction above $(pole_fraction)."
    end

    return θ_min, θ_max
end


"""
    verify_pole_exclusion(extra_mesh)

Check that no angular node lands at or near the poles after mesh construction.
Call immediately after building the angular mesh to confirm the shift was
applied correctly.

Prints the five smallest and five largest θ values found on the mesh, and
the corresponding sin(θ) values (= Je_ang / Je_ref, the relative Jacobian
magnitude at those nodes).
"""
function verify_pole_exclusion(extra_mesh)
    coords_ang = extra_mesh.extra_coords
    npoin_ang  = size(coords_ang, 2)

    θ_all = [coords_ang[1, ip] for ip = 1:npoin_ang]
    sort!(θ_all)

    θ_min = θ_all[1]
    θ_max = θ_all[end]

    @info "Pole exclusion verification:"
    @info "  Smallest θ values: $(round.(θ_all[1:min(5,end)],   digits=8))"
    @info "  Largest  θ values: $(round.(θ_all[end-min(4,end-1):end], digits=8))"
    @info "  sin(θ_min) = $(round(sin(θ_min), digits=8))"
    @info "  sin(θ_max) = $(round(sin(θ_max), digits=8))"

    if sin(θ_min) < 1e-8 || sin(θ_max) < 1e-8
        @error "Pole node detected: sin(θ) < 1e-8 at θ = $θ_min or $θ_max. " *
               "The θ mesh was not shifted correctly — check that θ_min and " *
               "θ_max from pole_shifted_theta_range are being passed to the " *
               "angular mesh constructor."
    elseif sin(θ_min) < 1e-4
        @warn "sin(θ_min) = $(sin(θ_min)) is small. Mass matrix condition " *
              "number will be O(1/sin(θ_min)) = O($(round(1/sin(θ_min), digits=0))). " *
              "Consider increasing pole_fraction."
    else
        @info "  ✓ No pole nodes detected. Mass matrix is safe to invert."
    end

    return θ_min, θ_max
end


"""
    check_solid_angle_with_pole_exclusion(extra_mesh, θ_min, θ_max)

Verify that the solid angle integral ∫ sin(θ) dθ dϕ recovers the expected
value for the truncated domain [θ_min, θ_max] × [0, 2π].

The exact value for the truncated sphere is:

    ∫_{θ_min}^{θ_max} ∫_0^{2π} sin(θ) dθ dϕ
        = 2π × [-cos(θ)]_{θ_min}^{θ_max}
        = 2π × (cos(θ_min) - cos(θ_max))
        = 2π × (cos(δ) - cos(π-δ))
        = 2π × 2cos(δ)
        = 4π cos(δ)

For δ = 0.00785 rad: 4π cos(δ) = 4π × 0.99997 ≈ 12.566 (differs from 4π
by only 0.003%).
"""
function check_solid_angle_with_pole_exclusion(extra_mesh, θ_min, θ_max)

    connijk_ang = extra_mesh.extra_connijk
    nop_ang     = extra_mesh.extra_nop
    nelem_ang   = extra_mesh.extra_nelem
    coords_ang  = extra_mesh.extra_coords
    Je_ang      = extra_mesh.extra_metrics.Je
    ωθ_weights  = extra_mesh.ωθ
    ωϕ_weights  = extra_mesh.ωθ

    solid_angle  = 0.0
    solid_angle_in  = 0.0

    for e = 1:nelem_ang
        nop = nop_ang[e]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            ip_ang = connijk_ang[e, iθ, iϕ]
            θ      = coords_ang[1, ip_ang]
            Je     = Je_ang[e, iθ, iϕ]
            ω_θ    = ωθ_weights[iθ]
            ω_ϕ    = ωϕ_weights[iϕ]

            dΩ = Je * ω_θ * ω_ϕ
            solid_angle += dΩ
            cos(θ) < 0.0 && (solid_angle_in += dΩ)
        end
    end

    # ϕ coverage correction (periodic endpoint excluded)
    npoin_ang = size(coords_ang, 2)
    ϕ_max   = maximum(coords_ang[2, ip] for ip = 1:npoin_ang)
    ϕ_scale = 1#2π / ϕ_max

    # Exact value for truncated domain
    exact_full = 4π * cos(θ_min)   # = 4π cos(δ) ≈ 4π for small δ
    exact_in   = 2π * cos(θ_min)   # lower hemisphere only

    solid_angle_corr    = solid_angle    * ϕ_scale
    solid_angle_in_corr = solid_angle_in * ϕ_scale

    err_full = abs(solid_angle_corr - exact_full) / exact_full * 100.0
    err_in   = abs(solid_angle_in_corr - exact_in) / exact_in * 100.0

    @info "Solid angle check (pole-excluded mesh):"
    @info "  θ range             : [$(round(θ_min,digits=6)), $(round(θ_max,digits=6))]"
    @info "  ϕ scale correction  : $(round(ϕ_scale, digits=6))"
    @info "  ∫dΩ full (raw)      : $(round(solid_angle,      digits=6))"
    @info "  ∫dΩ full (corrected): $(round(solid_angle_corr, digits=6))  " *
          "(exact = $(round(exact_full, digits=6)))"
    @info "  Error full sphere   : $(round(err_full, digits=4))%"
    @info "  ∫dΩ inflow (corr.)  : $(round(solid_angle_in_corr, digits=6))  " *
          "(exact = $(round(exact_in, digits=6)))"
    @info "  Error inflow hemi   : $(round(err_in, digits=4))%"

    if err_full > 1.0
        @warn "Solid angle error > 1% after pole exclusion and ϕ correction. " *
              "Check that Je_ang includes sin(θ) and that the ϕ mesh is " *
              "correctly periodic."
    else
        @info "  ✓ Solid angle integral correct to $(round(err_full, digits=3))%"
    end

    return (solid_angle_full = solid_angle_corr,
            solid_angle_in   = solid_angle_in_corr,
            exact_full       = exact_full,
            exact_in         = exact_in)
end