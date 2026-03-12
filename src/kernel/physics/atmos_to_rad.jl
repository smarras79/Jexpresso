"""
    atmos_to_rad_longwave(atmos_data, npoin)

Convert atmospheric state to volumetric absorption and scattering coefficients
for a broadband gray longwave (thermal infrared, ~10 μm effective wavelength)
radiative transfer calculation.

# Physics

The thermal IR gray atmosphere is dominated by:
- Water vapor continuum and rotational band absorption
- CO₂ and other well-mixed greenhouse gas absorption (lumped into k_dry)
- Liquid and ice cloud absorption (geometric optics regime)
- Scattering is weak in the thermal IR (single-scatter albedo << 1)

Mass absorption cross sections are Planck-mean values appropriate for
Earth's atmosphere in the 8–12 μm window and surrounding bands.

# References
- Pierrehumbert (2010), "Principles of Planetary Climate", Ch. 4
- Lacis & Hansen (1974), J. Atmos. Sci., 31, 118–133
- Fu & Liou (1993), J. Atmos. Sci., 50, 2013–2025 (cloud optics)

# Arguments
- `atmos_data`: struct with fields t_lay, p_lay, vmr_h2o (mol/mol), q_liq (kg/kg),
                q_ice (kg/kg), rho (kg/m³)
- `npoin`: number of spatial nodes

# Returns
- `κ`: volumetric absorption coefficient (m⁻¹), passed to the RTE as the
       extinction minus scattering term
- `σ`: volumetric scattering coefficient (m⁻¹)

# Notes
The total extinction seen by the RTE solver is κ + σ. In the longwave the
single-scatter albedo ω₀ = σ/(κ+σ) is typically 0.02–0.10 in cloudy regions
and ~0 in clear sky.
"""

function atmos_to_rad_longwave(atmos_data, npoin)

    PhysConst = PhysicalConst{Float64}()

    κ = zeros(npoin)
    σ = zeros(npoin)

    # ── Molecular weight ratio for VMR → mass mixing ratio conversion ─────
    M_h2o = 18.015e-3   # kg/mol
    M_air  = 28.964e-3   # kg/mol
    r_mol  = M_h2o / M_air   # ≈ 0.6220

    # ── Mass absorption cross sections (m²/kg) ────────────────────────────
    # Water vapor: Planck-mean over thermal IR, including continuum.
    # Typical range 0.01–0.10 m²/kg; 0.05 is a standard gray value.
    k_vap = 0.05

    # Dry air / well-mixed gases (CO₂, O₃, CH₄ etc.) — pressure scaled.
    # At surface conditions roughly 0.002 m²/kg for a CO₂-like gray gas.
    k_dry = 0.002

    # Liquid cloud: geometric optics, Q_abs ≈ 2, r_eff ≈ 10 μm, ρ_liq = 1000 kg/m³
    # k = 3 Q_abs / (4 ρ_liq r_eff) ≈ 3×2/(4×1000×10e-6) ≈ 150 m²/kg
    # Broadband Planck-mean is somewhat lower; 100 m²/kg is representative.
    k_liq = 100.0

    # Ice cloud: larger effective radius (~30 μm), lower absorption efficiency.
    # k_ice ≈ k_liq × (r_eff_liq / r_eff_ice) × (Q_abs_ice / Q_abs_liq) ≈ 35 m²/kg
    k_ice = 35.0

    # ── Scattering ────────────────────────────────────────────────────────
    # Thermal IR scattering is negligible in clear sky (Rayleigh ∝ λ⁻⁴ → 0).
    # In clouds, absorption dominates over scattering at 10 μm.
    # Single-scatter albedo ω₀ ≈ 0.02–0.08 for liquid clouds in thermal IR.
    ω₀_clear = 0.0    # clear-sky: no scattering
    ω₀_cloud = 0.05   # in-cloud: small but non-zero

    # ── Reference pressure for collision-broadening scaling ───────────────
    P_ref = 101325.0   # Pa

    for ip = 1:npoin
        T    = atmos_data.t_lay[ip]     # K
        P    = atmos_data.p_lay[ip]     # Pa
        vmr  = atmos_data.vmr_h2o[ip]  # mol/mol
        qliq = atmos_data.q_liq[ip]    # kg/kg
        qice = atmos_data.q_ice[ip]    # kg/kg
        ρ    = atmos_data.rho[ip]      # kg/m³

        # VMR (mol/mol) → mass mixing ratio (kg/kg)
        qv = vmr * r_mol

        # Pressure scaling: collision broadening ∝ P for gas absorption.
        # Cloud particle optics are independent of pressure.
        p_scale = P / P_ref

        # Total mass absorption coefficient (m²/kg)
        k_abs = k_dry * p_scale + k_vap * qv * p_scale + k_liq * qliq + k_ice * qice

        # In-cloud fraction drives scattering albedo
        in_cloud = (qliq + qice) > 1e-8
        ω₀ = in_cloud ? ω₀_cloud : ω₀_clear

        # Volumetric coefficients (m⁻¹)
        κ_total   = ρ * k_abs
        κ[ip] = (1.0 - ω₀) * κ_total
        σ[ip] =        ω₀  * κ_total
    end

    κ_ext = κ .+ σ
    @info "LW extinction   extrema: $(extrema(κ_ext))"
    @info "LW absorption   extrema: $(extrema(κ))"
    @info "LW scattering   extrema: $(extrema(σ))"
    @info "LW SSA          extrema: $(extrema(σ ./ max.(κ_ext, 1e-30)))"

    return κ, σ
end


"""
    atmos_to_rad_shortwave(atmos_data, npoin)

Convert atmospheric state to volumetric absorption and scattering coefficients
for a broadband gray shortwave (visible/near-infrared, ~550 nm effective
wavelength) radiative transfer calculation.

# Physics

The shortwave gray atmosphere is dominated by:
- Rayleigh scattering by air molecules (strong at 550 nm, ω₀ = 1 exactly)
- Water vapor near-IR absorption bands (0.94, 1.14, 1.38, 1.87 μm weighted)
- Ozone Hartley/Huggins band absorption (important in upper atmosphere)
- Liquid cloud Mie scattering (highly forward-peaked, g ≈ 0.85, ω₀ ≈ 0.999)
- Ice cloud scattering (g ≈ 0.80, ω₀ ≈ 0.99)
- Aerosol absorption is neglected here (add a_aer term for urban cases)

Scattering is the dominant process in the shortwave — the single-scatter
albedo is typically 0.9–1.0 in clear sky and >0.99 in liquid clouds.

# References
- Liou (2002), "An Introduction to Atmospheric Radiation", 2nd ed., Ch. 3 & 6
- Hansen & Travis (1974), Space Sci. Rev., 16, 527–610 (Rayleigh scattering)
- Hu & Stamnes (1993), J. Climate, 6, 728–742 (liquid cloud optics)

# Arguments
- `atmos_data`: struct with fields t_lay, p_lay, vmr_h2o (mol/mol),
                vmr_o3 (mol/mol), q_liq (kg/kg), q_ice (kg/kg), rho (kg/m³)
- `npoin`: number of spatial nodes

# Returns
- `κ`: volumetric absorption coefficient (m⁻¹)
- `σ`: volumetric scattering coefficient (m⁻¹)

# Notes
The Henyey-Greenstein asymmetry parameter g used in the RTE scattering kernel
should be set consistently with the dominant scatterer:
- Clear sky (Rayleigh): g = 0.0 (symmetric)
- Cloudy (liquid):      g = 0.85
- Mixed:                g weighted by scattering optical depth contribution

A single g cannot represent both Rayleigh and cloud scattering simultaneously.
For clear-sky problems use g = 0; for overcast scenes use g = 0.85.
"""
function atmos_to_rad_shortwave(atmos_data, npoin)

    PhysConst = PhysicalConst{Float64}()

    κ = zeros(npoin)
    σ = zeros(npoin)
    g_eff = zeros(npoin)

    # ── Molecular weight ratios ───────────────────────────────────────────
    M_h2o = 18.015e-3   # kg/mol
    M_o3  = 47.998e-3   # kg/mol
    M_air  = 28.964e-3  # kg/mol
    r_h2o = M_h2o / M_air   # ≈ 0.6220
    r_o3  = M_o3  / M_air   # ≈ 1.6573

    # ── Mass absorption cross sections (m²/kg) ────────────────────────────
    # Water vapor: solar band (0.7–3 μm) Planck-weighted mean.
    # Clear-sky SW absorption by H₂O is weaker than in LW: ~0.002–0.01 m²/kg.
    k_vap = 0.006

    # Ozone: Hartley-Huggins band (0.2–0.35 μm) dominates SW ozone absorption.
    # Cross section at 550 nm is ~4×10⁻²³ m²/molecule = 5×10⁻⁴ m²/kg.
    # Chappuis band (0.45–0.75 μm) adds ~1×10⁻⁴ m²/kg.
    # Combined gray value: ~6×10⁻⁴ m²/kg (solar-flux weighted).
    k_o3 = 6.0e-4

    # Liquid cloud: Mie scattering at 550 nm, r_eff ≈ 10 μm.
    # Extinction efficiency Q_ext ≈ 2, single-scatter albedo ω₀ ≈ 0.9999.
    # k_ext = 3 Q_ext / (4 ρ_liq r_eff) ≈ 150 m²/kg
    # Absorption contribution: k_abs ≈ (1 - ω₀) × k_ext ≈ 0.015 m²/kg
    k_liq_ext = 150.0     # total extinction (m²/kg)
    ω₀_liq    = 0.9999    # single-scatter albedo

    # Ice cloud: r_eff ≈ 30 μm, slightly absorbing in near-IR.
    # k_ext ≈ 50 m²/kg, ω₀ ≈ 0.99
    k_ice_ext = 50.0
    ω₀_ice    = 0.990

    # ── Rayleigh scattering cross section ─────────────────────────────────
    # σ_Ray = (8π³/3) (n²-1)² / (N² λ⁴) where n is air refractive index,
    # N is number density. At 550 nm: σ_Ray ≈ 5.1×10⁻³¹ m² per molecule.
    # In mass units (dividing by m_air = 28.97×1.66e-27 kg):
    #   k_Ray = 5.1e-31 / (28.97×1.66e-27) ≈ 1.06e-5 m²/kg
    # Rayleigh scattering is conservative (ω₀ = 1 exactly).
    k_ray = 1.06e-5   # m²/kg at 550 nm; scales as λ⁻⁴

    # ── Reference pressure ────────────────────────────────────────────────
    P_ref = 101325.0   # Pa

    for ip = 1:npoin
        T    = atmos_data.t_lay[ip]
        P    = atmos_data.p_lay[ip]     # Pa
        vmr  = atmos_data.vmr_h2o[ip]  # mol/mol
        vmro = atmos_data.vmr_o3[ip]   # mol/mol  (add this field if not present)
        qliq = atmos_data.q_liq[ip]    # kg/kg
        qice = atmos_data.q_ice[ip]    # kg/kg
        ρ    = atmos_data.rho[ip]      # kg/m³

        # VMR → mass mixing ratio
        qv   = vmr  * r_h2o
        qo3  = vmro * r_o3

        # Number density (molecules/m³) for Rayleigh — needed for pressure scaling
        # n = P/(kB T); Rayleigh ∝ N ∝ P/T
        p_scale = P / P_ref

        # ── Gas contributions ─────────────────────────────────────────────
        # Rayleigh: purely scattering, scales with air density directly
        σ_ray = ρ * k_ray   # m⁻¹ — conservative scattering

        # H₂O absorption: weak in SW, slight pressure dependence
        κ_vap = ρ * k_vap * qv * p_scale

        # O₃ absorption: independent of pressure (upper-atmosphere UV sink)
        κ_o3  = ρ * k_o3 * qo3

        # ── Cloud contributions ───────────────────────────────────────────
        κ_liq_abs = ρ * k_liq_ext * qliq * (1.0 - ω₀_liq)
        σ_liq_sca = ρ * k_liq_ext * qliq * ω₀_liq

        κ_ice_abs = ρ * k_ice_ext * qice * (1.0 - ω₀_ice)
        σ_ice_sca = ρ * k_ice_ext * qice * ω₀_ice

        # ── Total coefficients ────────────────────────────────────────────
        κ[ip] = κ_vap + κ_o3 + κ_liq_abs + κ_ice_abs
        σ[ip] = σ_ray + σ_liq_sca + σ_ice_sca

        κ_abs_min = 1e-10
        κ[ip]  = max(ρ * κ[ip], κ_abs_min)
        g_eff[ip] = if σ[ip] > 1e-30
            (0.0 * σ_ray + 0.85 * σ_liq_sca + 0.80 * σ_ice_sca) / σ[ip]
        else
            0.0
        end

    end

    κ_ext = κ .+ σ
    ω₀    = σ ./ max.(κ_ext, 1e-30)

    @info "SW extinction   extrema: $(extrema(κ_ext))"
    @info "SW absorption   extrema: $(extrema(κ))"
    @info "SW scattering   extrema: $(extrema(σ))"
    @info "SW SSA          extrema: $(extrema(ω₀))"
    @info "SW effective g  extrema: $(extrema(g_eff))"

    return κ, σ, g_eff
end

"""
    verify_optical_depth(mesh, κ, σ, ω, Je, ngl;
                         sample_xy = [(0.5, 0.5)])

Verify column optical depth computed from extinction coefficients κ (absorption,
m⁻¹) and σ (scattering, m⁻¹) by integrating vertically through the LGL mesh
at a set of requested (x,y) sample locations.

For each sample point the nearest mesh column is found and the total optical
depth τ = ∫ (κ+σ) dz is computed two ways:

  - **GL quadrature**: exact for the polynomial representation; uses the LGL
    weights ω and element Jacobians Je to integrate within each element, then
    sums over elements in the column.

  - **Trapezoidal rule**: integrates over the sorted unique z-nodes in the
    column; serves as an independent check and degrades gracefully if κ_ext
    has inter-element jumps.

A summary table is printed on rank 0 showing both values and their relative
difference. Columns where the two methods disagree by more than 5% are flagged.

# Arguments
- `mesh`  : mesh struct with fields `x, y, z` (Vector, length npoin),
            `connijk[iel,i,j,k]`, `nelem`, `npoin`,
            `xmin, xmax, ymin, ymax, zmin, zmax`
- `κ`     : absorption coefficient (m⁻¹), length npoin
- `σ`     : scattering coefficient (m⁻¹), length npoin
- `ω`     : 1D LGL quadrature weights, length ngl
- `Je`    : element Jacobian `[iel, i, j, k]`
- `ngl`   : number of GL points per spatial direction

# Keyword Arguments
- `sample_xy`: vector of `(x, y)` tuples; the nearest column centroid is used
               for each entry

# Returns
Named tuple with fields:
- `τ_GL`      : GL optical depth per sample column
- `τ_trap`    : trapezoidal optical depth per sample column
- `column_xy` : (x,y) centroid of the matched column for each sample point

# Physical reference ranges
- LW clear-sky:  τ ≈ 2 – 5
- LW cloudy:     τ ≈ 5 – 20
- SW clear-sky:  τ ≈ 0.1 – 0.3  (Rayleigh dominated)
- SW cloudy:     τ ≈ 5 – 50
"""
function verify_optical_depth(mesh, κ, σ, ω, Je, ngl;
                               sample_xy = [(0.5, 0.5)])

    rank  = MPI.Comm_rank(MPI.COMM_WORLD)
    κ_ext = κ .+ σ

    # ── Scale-aware column tolerance ─────────────────────────────────────────
    # Use a fraction of the minimum element width rather than an absolute value.
    # This correctly groups LGL nodes whose x,y coords differ only by floating
    # point rounding from the reference-to-physical mapping.
    x_range = mesh.xmax - mesh.xmin
    y_range = mesh.ymax - mesh.ymin
    tol_xy  = min(x_range, y_range) * 1e-6   # 1 ppm of domain width

    # ── Build column map ──────────────────────────────────────────────────────
    # Round to nearest tol_xy to group nodes at the same (x,y) position
    column_nodes = Dict{Tuple{Float64,Float64}, Vector{Int}}()
    for ip = 1:mesh.npoin
        key = (round(mesh.x[ip] / tol_xy) * tol_xy,
               round(mesh.y[ip] / tol_xy) * tol_xy)
        push!(get!(column_nodes, key, Int[]), ip)
    end

    # ── Sanity check: all columns should have the same number of nodes ────────
    col_sizes  = [length(v) for v in values(column_nodes)]
    expected_n = maximum(col_sizes)   # most common column size
    n_thin     = count(s -> s < expected_n * 0.9, col_sizes)
    if n_thin > 0 && rank == 0
        @warn "$(n_thin) columns have fewer nodes than expected ($expected_n). " *
              "tol_xy = $(tol_xy) may be too tight — LGL coordinate rounding " *
              "is splitting columns. Try increasing tol_xy by 10×."
    end

    col_keys = collect(keys(column_nodes))

    # ── Match sample points ───────────────────────────────────────────────────
    n_samples    = length(sample_xy)
    matched_cols = Vector{Tuple{Float64,Float64}}(undef, n_samples)

    for (s, (x_req, y_req)) in enumerate(sample_xy)
        best_dist = Inf
        best_key  = col_keys[1]
        for key in col_keys
            d = (key[1] - x_req)^2 + (key[2] - y_req)^2
            if d < best_dist
                best_dist = d
                best_key  = key
            end
        end
        matched_cols[s] = best_key
        if rank == 0 && sqrt(best_dist) > 0.1*(mesh.xmax - mesh.xmin)
            @warn "Sample $s: requested ($x_req,$y_req) matched to " *
                  "($(round(best_key[1],digits=4)), $(round(best_key[2],digits=4))); " *
                  "distance $(round(sqrt(best_dist),sigdigits=3)) m"
        end
    end

    # ── Build element-column lookup ───────────────────────────────────────────
    # For each element find which column it belongs to by checking all
    # horizontal nodes (i,j at fixed k=1), not just i=j=1, and taking
    # a majority vote. This is more robust than checking a single node.
    elem_col = Dict{Int, Tuple{Float64,Float64}}()
    for iel = 1:mesh.nelem
        # Collect (x,y) keys for all horizontal nodes in this element
        key_votes = Dict{Tuple{Float64,Float64}, Int}()
        for i = 1:ngl, j = 1:ngl
            ip_ref = mesh.connijk[iel, i, j, 1]
            key    = (round(mesh.x[ip_ref] / tol_xy) * tol_xy,
                      round(mesh.y[ip_ref] / tol_xy) * tol_xy)
            key_votes[key] = get(key_votes, key, 0) + 1
        end
        # Pick the key with the most votes
        best_key = argmax(key_votes)
        elem_col[iel] = best_key
    end

    # ── Verify element-column assignment ─────────────────────────────────────
    # Each column should have nelem / n_columns elements.
    # If any column has zero elements the tol_xy is still too tight.
    n_cols        = length(col_keys)
    elems_per_col = mesh.nelem / n_cols
    col_elem_count = Dict{Tuple{Float64,Float64}, Int}()
    for iel = 1:mesh.nelem
        k = elem_col[iel]
        col_elem_count[k] = get(col_elem_count, k, 0) + 1
    end
    n_empty = count(k -> !haskey(col_elem_count, k), col_keys)
    if n_empty > 0 && rank == 0
        @warn "$n_empty columns have no elements assigned. " *
              "Increase tol_xy by 10× (current: $tol_xy)."
    end

    # ── Integrate each matched column ─────────────────────────────────────────
    τ_GL_all   = zeros(Float64, n_samples)
    τ_trap_all = zeros(Float64, n_samples)

    for (s, col_key) in enumerate(matched_cols)

        # ── GL quadrature ─────────────────────────────────────────────────
        τ_GL  = 0.0
        i_fix = ceil(Int, ngl/2)   # use middle i,j to avoid boundary nodes
        j_fix = ceil(Int, ngl/2)   # which may have degenerate metrics

        n_elem_in_col = 0
        for iel = 1:mesh.nelem
            elem_col[iel] == col_key || continue
            n_elem_in_col += 1
            for k = 1:ngl
                ip        = mesh.connijk[iel, i_fix, j_fix, k]
                dz_weight = ω[k] * Je[iel, i_fix, j_fix, k] /
                            (ω[i_fix] * ω[j_fix])
                τ_GL     += κ_ext[ip] * dz_weight
            end
        end

        if n_elem_in_col == 0 && rank == 0
            @warn "Column $s at $(col_key): no elements found. " *
                  "tol_xy mismatch between column_nodes and elem_col."
        end

        # ── Trapezoidal rule ──────────────────────────────────────────────
        nodes_in_col = column_nodes[col_key]

        if isempty(nodes_in_col)
            τ_trap_all[s] = 0.0
            continue
        end

        z_vals       = mesh.z[nodes_in_col]
        sort_idx     = sortperm(z_vals)
        nodes_sorted = nodes_in_col[sort_idx]
        z_sorted     = z_vals[sort_idx]

        # Deduplicate shared interface nodes
        tol_z       = (mesh.zmax - mesh.zmin) * 1e-8
        z_unique    = Float64[]
        κext_unique = Float64[]
        for (ip, z_val) in zip(nodes_sorted, z_sorted)
            if !isempty(z_unique) && abs(z_val - z_unique[end]) < tol_z
                κext_unique[end] = 0.5*(κext_unique[end] + κ_ext[ip])
            else
                push!(z_unique,    z_val)
                push!(κext_unique, κ_ext[ip])
            end
        end

        τ_trap = 0.0
        for i = 1:length(z_unique)-1
            τ_trap += 0.5*(κext_unique[i] + κext_unique[i+1]) * (z_unique[i+1] - z_unique[i])
        end

        τ_GL_all[s]   = τ_GL
        τ_trap_all[s] = τ_trap
    end

    # ── Print summary ─────────────────────────────────────────────────────────
    if rank == 0
        w = 74
        println("=" ^ w)
        println("  Optical Depth Verification")
        println("=" ^ w)
        @printf("  %-5s  %-12s  %-12s  %-14s  %-14s  %-9s\n",
                "Col", "x (m)", "y (m)", "τ_GL", "τ_trap", "Δτ/τ (%)")
        println("-" ^ w)
        for s = 1:n_samples
            x_c, y_c = matched_cols[s]
            τg        = τ_GL_all[s]
            τt        = τ_trap_all[s]
            rel_diff  = abs(τg - τt) / max(τg, 1e-30) * 100.0
            flag      = rel_diff > 5.0 ? " !" : ""
            @printf("  %-5d  %-12.4f  %-12.4f  %-14.6f  %-14.6f  %-9.2f%s\n",
                    s, x_c, y_c, τg, τt, rel_diff, flag)
        end
        println("=" ^ w)
        println()
        println("  Physical reference ranges:")
        println("  LW clear-sky:  τ ≈ 2 – 5")
        println("  LW cloudy:     τ ≈ 5 – 20")
        println("  SW clear-sky:  τ ≈ 0.1 – 0.3")
        println("  SW cloudy:     τ ≈ 5 – 50")
        println()
        for s = 1:n_samples
            τg       = τ_GL_all[s]
            τt       = τ_trap_all[s]
            rel_diff = abs(τg - τt) / max(τg, 1e-30) * 100.0
            if rel_diff > 5.0
                @warn "Column $s: GL and trapezoidal differ by " *
                      "$(round(rel_diff,digits=1))%."
            end
            if τg < 1e-6
                @warn "Column $s: optical depth ≈ 0. Check tol_xy and " *
                      "that atmos_to_rad_* was called."
            end
            if τg > 100.0
                @warn "Column $s: very large optical depth τ=$(round(τg,digits=1)). " *
                      "Check units of κ and σ (should be m⁻¹)."
            end
        end
    end

    return (τ_GL      = τ_GL_all,
            τ_trap    = τ_trap_all,
            column_xy = matched_cols)
end

# ─────────────────────────────────────────────────────────────────────────────
# Shared geometry helpers
# ─────────────────────────────────────────────────────────────────────────────

"""
    make_boundary_predicates(mesh)

Returns a named tuple of boundary predicate functions closed over `mesh` extents.
Call once at setup and pass the result to the BC and RHS functions.
"""
function make_boundary_predicates(mesh)
    zmin = mesh.zmin; zmax = mesh.zmax
    return (
        is_top    = z -> abs(z - zmax) < 1e-5,
        is_bottom = z -> abs(z - zmin) < 1e-5,
    )
end

function check_beam_flux(extra_mesh, sw, ngl)

    connijk_ang = extra_mesh.extra_connijk
    nop_ang     = extra_mesh.extra_nop
    nelem_ang   = extra_mesh.extra_nelem
    coords_ang  = extra_mesh.extra_coords
    Je_ang      = extra_mesh.extra_metrics.Je
    ωθ_weights  = extra_mesh.ωθ
    ωϕ_weights  = extra_mesh.ωθ

    sw.μ₀ > 0.0 || error("μ₀ must be positive")
    θ_sun = π - acos(clamp(sw.μ₀, 0.0, 1.0))

    flux_integrated  = 0.0
    solid_angle_full = 0.0
    solid_angle_in   = 0.0

    for e = 1:nelem_ang
        nop = nop_ang[e]
        for iθ = 1:nop+1, iϕ = 1:nop+1
            ip_ang = connijk_ang[e, iθ, iϕ]
            θ      = coords_ang[1, ip_ang]
            ϕ      = coords_ang[2, ip_ang]
            Je     = Je_ang[e, iθ, iϕ]
            ω_θ    = ωθ_weights[iθ]
            ω_ϕ    = ωϕ_weights[iϕ]

            dΩ = Je * ω_θ * ω_ϕ
            solid_angle_full += dΩ

            cos(θ) < 0.0 || continue   # inflow at TOA

            solid_angle_in += dΩ

            Δθ = θ - θ_sun
            Δϕ = ϕ - sw.φ₀
            Δϕ = Δϕ - 2π * round(Δϕ / (2π))
            ang_dist² = Δθ^2 + (sin(θ_sun) * Δϕ)^2

            # Radiance: S₀/(π δ²), no μ₀ in numerator
            I_beam = (sw.S₀_flux / (π * sw.δ_beam^2)) *
                     exp(-ang_dist² / sw.δ_beam^2)

            # Flux = ∫ I |cosθ| dΩ; the μ₀ projection comes from cosθ naturally
            flux_integrated += I_beam * abs(cos(θ)) * dΩ
        end
    end

    expected  = sw.S₀_flux * sw.μ₀
    rel_err   = abs(flux_integrated - expected) / expected * 100.0

    @info "Beam flux check:"
    @info "  δ_beam              : $(sw.δ_beam) rad  ($(round(rad2deg(sw.δ_beam),digits=1))°)"
    @info "  θ_sun               : $(round(θ_sun,digits=4)) rad  ($(round(rad2deg(θ_sun),digits=1))°)"
    @info "  ∫ dΩ full sphere    : $(round(solid_angle_full,digits=6))  (exact 4π = $(round(4π,digits=6)))"
    @info "  ∫ dΩ inflow hemi    : $(round(solid_angle_in,  digits=6))  (exact 2π = $(round(2π,digits=6)))"
    @info "  Integrated flux     : $(round(flux_integrated, digits=4)) W/m²"
    @info "  Expected  (S₀ μ₀)  : $(round(expected,        digits=4)) W/m²"
    @info "  Relative error      : $(round(rel_err,         digits=2))%"

    # Residual error after removing the μ₀² double-projection
    flux_old_convention = flux_integrated / sw.μ₀   # what old check would give
    @info "  [Diagnostic] flux/μ₀ = $(round(flux_old_convention,digits=4)) W/m²  " *
          "(should equal S₀ = $(sw.S₀_flux) if old double-μ₀ issue was present)"

    if rel_err > 5.0
        @warn "Flux error > 5%. The beam is under-resolved on this angular mesh. " *
              "Recommended minimum δ_beam ≈ $(round(rad2deg(π/(4*ngl)), digits=1))° " *
              "for a $(round(Int, π/maximum(nop_ang .+ 1)))×$(round(Int,2π/maximum(nop_ang .+ 1))) " *
              "element mesh of order $(maximum(nop_ang))."
    end

    return (flux = flux_integrated, expected = expected, rel_err = rel_err,
            solid_angle_full = solid_angle_full, solid_angle_in = solid_angle_in)
end

"""
    build_sw_lateral_bc_profile(mesh, κ_ext, ω, Je, ngl)

Precompute the cumulative vertical optical depth profile τ(z) from TOA
downward, averaged horizontally across the domain. Used to set the
attenuated direct beam lateral BC for the shortwave.

Returns a linear interpolant: τ_from_TOA(z) giving the vertical optical
depth between z and zmax. The slant-path optical depth at solar zenith
angle θ_sun is τ(z)/μ₀.
"""
function build_sw_lateral_bc_profile(mesh, κ_ext, ngl)

    # ── Find a vertical stack of elements through the domain centre ───────────
    # In a structured hex mesh connijk[iel,i,j,k] gives the global node index.
    # Nodes at the same (i,j) across elements stacked in z form an exact column.
    # We find one such stack by:
    # 1. Picking a reference node near the domain centre at z=zmin
    # 2. Finding all nodes with the same (x,y) by exact connijk lookup

    # Fix i,j at the middle of the reference element
    i_fix = ceil(Int, ngl/2)
    j_fix = ceil(Int, ngl/2)

    # Find the element whose centre is closest to (x_mid, y_mid, zmin)
    x_mid = 0.5*(mesh.xmin + mesh.xmax)
    y_mid = 0.5*(mesh.ymin + mesh.ymax)

    best_dist = Inf
    ref_iel   = 1
    for iel = 1:mesh.nelem
        ip = mesh.connijk[iel, i_fix, j_fix, 1]
        dx = mesh.x[ip] - x_mid
        dy = mesh.y[ip] - y_mid
        dz = mesh.z[ip] - mesh.zmin
        d  = dx^2 + dy^2 + 100*dz^2   # weight z heavily to prefer bottom elements
        if d < best_dist
            best_dist = d
            ref_iel   = iel
        end
    end

    # ── Collect all nodes in this column via connectivity ─────────────────────
    # The reference x,y coordinates come from connijk — exact by construction
    ref_ip = mesh.connijk[ref_iel, i_fix, j_fix, 1]
    x_col  = mesh.x[ref_ip]
    y_col  = mesh.y[ref_ip]

    # Collect every node ip where x[ip] ≈ x_col AND y[ip] ≈ y_col
    # Use a tolerance relative to the element size, not the domain size
    # Element size estimate: domain / number of elements per direction
    n_elem_per_dir = round(Int, mesh.nelem^(1/3))
    elem_size      = min(mesh.xmax-mesh.xmin,
                         mesh.ymax-mesh.ymin) / max(n_elem_per_dir, 1)
    tol_xy         = elem_size * 1e-4   # 0.01% of element size — tight but safe

    nodes_in_col = Int[]
    for ip = 1:mesh.npoin
        if abs(mesh.x[ip] - x_col) < tol_xy &&
           abs(mesh.y[ip] - y_col) < tol_xy
            push!(nodes_in_col, ip)
        end
    end

    # ── Validate ──────────────────────────────────────────────────────────────
    if isempty(nodes_in_col)
        @error "Column finder: no nodes found near " *
               "($(round(x_col,digits=4)), $(round(y_col,digits=4))). " *
               "This should not happen since ref_ip was found via connijk."
        return [mesh.zmax, mesh.zmin], [0.0, 0.0]
    end

    z_span = maximum(mesh.z[nodes_in_col]) - minimum(mesh.z[nodes_in_col])
    if z_span < 0.5*(mesh.zmax - mesh.zmin)
        @warn "Column spans only $(round(z_span,digits=2)) m of " *
              "$(round(mesh.zmax-mesh.zmin,digits=2)) m total. " *
              "Column may be incomplete."
    end

    @info "SW lateral BC profile:"
    @info "  Reference element: $ref_iel"
    @info "  Column (x,y): ($(round(x_col,digits=4)), $(round(y_col,digits=4)))"
    @info "  Nodes found: $(length(nodes_in_col))"
    @info "  z range: [$(round(minimum(mesh.z[nodes_in_col]),digits=2)), " *
                      "$(round(maximum(mesh.z[nodes_in_col]),digits=2))]"

    # ── Sort by z descending (TOA first), deduplicate interfaces ─────────────
    z_vals   = mesh.z[nodes_in_col]
    sort_idx = sortperm(z_vals, rev=true)
    nodes_s  = nodes_in_col[sort_idx]
    z_s      = z_vals[sort_idx]

    tol_z       = (mesh.zmax - mesh.zmin) * 1e-8
    z_unique    = Float64[]
    κext_unique = Float64[]
    for (ip, zv) in zip(nodes_s, z_s)
        if !isempty(z_unique) && abs(zv - z_unique[end]) < tol_z
            κext_unique[end] = 0.5*(κext_unique[end] + κ_ext[ip])
        else
            push!(z_unique,    zv)
            push!(κext_unique, κ_ext[ip])
        end
    end

    # ── Cumulative τ from TOA downward (trapezoidal) ──────────────────────────
    n          = length(z_unique)
    τ_from_TOA = zeros(Float64, n)
    for i = 2:n
        dz             = z_unique[i-1] - z_unique[i]   # positive, descending z
        τ_from_TOA[i]  = τ_from_TOA[i-1] +
                         0.5*(κext_unique[i-1] + κext_unique[i]) * dz
    end

    @info "  Total optical depth τ(surface): $(round(τ_from_TOA[end], digits=4))"

    return z_unique, τ_from_TOA
end


function interp_optical_depth(z, z_prof, τ_from_TOA)
    # z_prof is sorted from zmax (TOA) down to zmin (surface)
    # τ_from_TOA increases from 0 at TOA to τ_total at surface
    z < z_prof[end] && return τ_from_TOA[end]   # below surface
    z > z_prof[1]   && return 0.0               # above TOA

    # z_prof is descending so use reverse search
    idx = searchsortedfirst(z_prof, z, rev=true) - 1
    idx = clamp(idx, 1, length(z_prof)-1)
    frac = (z_prof[idx] - z) / (z_prof[idx] - z_prof[idx+1] + 1e-30)
    return τ_from_TOA[idx] + frac * (τ_from_TOA[idx+1] - τ_from_TOA[idx])
end