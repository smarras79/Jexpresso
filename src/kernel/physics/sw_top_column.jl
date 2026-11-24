# ── SW top boundary from the atmosphere above the lid (:RT_sw_top_bc => :column, the default) ───────────────────────
#
# :RT_sw_top_bc => :vacuum is the old behaviour: the direct beam enters the lid at S₀ and no diffuse radiance enters
# through the top. With :column (default since 2026-11-24) the air above the lid is the column of lw_top_column.jl
# (same profile settings, :RT_lw_top_profile etc.), with the gas SW optics of atmos_to_rad_shortwave (Rayleigh
# k_ray, conservative, g = 0, plus water vapour k_vap·qv·p/p_ref; the model carries no ozone):
#   - the direct beam enters the lid as S₀ e^{-τ_above/μ₀} (τ_above = extinction optical depth above the lid);
#   - the diffuse radiance entering through the top is the plane-parallel radiance at the lid of the column made of
#     the air above the lid on top of the domain's SW reference column (build_sw_column_optics), from
#     build_sw_lateral_plane_parallel (exact single scattering, which carries the azimuth relative to the sun, plus
#     the azimuth-averaged multiple scattering of the 1-D solver, with the surface albedo);
#   - real lateral walls (:RT_sw_lateral_bc :plane_parallel / :single_scatter) read the same extended-column object.
# At μ₀ = 0.5 over a 3 km lid the beam loses ~10 % (mostly Rayleigh, about half of which returns as diffuse light);
# the vapour absorption above the lid is a few percent.

const SW_TOP = Ref{Any}(nothing)        # nothing (:vacuum) or (lat = callable (z, θ, ϕ), z_lid, τ_above)

"""
    sw_top_setup!(mesh, atmos_data, bdy, inputs, comm, κ_ext_sw, σ_sw, g_eff, S0, μ0, φ0, albedo) -> attenuation

Collective. Sets SW_TOP for this RT call and returns the factor e^{-τ_above/μ₀} by which the direct beam is reduced at
the lid (1 for :vacuum).
"""
function sw_top_setup!(mesh, atmos_data, bdy, inputs, comm, κ_ext_sw, σ_sw, g_eff, S0, μ0, φ0, albedo)
    mode = Symbol(get(inputs, :RT_sw_top_bc, :column))
    mode in (:vacuum, :column) || error(":RT_sw_top_bc must be :vacuum or :column")
    SW_TOP[] = nothing
    (mode == :vacuum || μ0 <= 0) && return 1.0
    zl, T_lid, p_lid, qv_lid = rt_lid_state(mesh, atmos_data, bdy, comm)
    z, T, p, qv = rt_top_profile(inputs, zl, T_lid, p_lid, qv_lid)          # ascending from the lid
    k_ray = 1.06e-5; k_vap = 0.006; P_ref = 101325.0                        # atmos_to_rad_shortwave
    na = length(z)
    ea = zeros(na); ωa = zeros(na)
    for k = 1:na
        ρ = p[k] / (287.04 * T[k] * (1 + 0.61 * qv[k]))
        σr = ρ * k_ray; κv = ρ * k_vap * qv[k] * p[k] / P_ref
        ea[k] = σr + κv; ωa[k] = σr / ea[k]
    end
    dom = build_sw_column_optics(mesh, κ_ext_sw, σ_sw, g_eff, mesh.ngl)     # collective; z descending from the lid
    # extended column, TOA first: the air above the lid (without the lid level itself), then the domain column
    ze = vcat(reverse(z[2:end]), dom.z); ee = vcat(reverse(ea[2:end]), dom.e)
    ωe = vcat(reverse(ωa[2:end]), dom.ω); ge = vcat(zeros(na - 1), dom.g)
    τe = zeros(length(ze))
    for i = 2:length(ze); τe[i] = τe[i-1] + 0.5 * (ee[i-1] + ee[i]) * (ze[i-1] - ze[i]); end
    τ_above = τe[na]                                                        # at the lid
    lat = build_sw_lateral_plane_parallel((z = ze, τ = τe, ω = ωe, g = ge, e = ee), S0, μ0, φ0, albedo)
    SW_TOP[] = (lat = lat, z_lid = Float64(dom.z[1]), τ_above = τ_above)
    att = exp(-τ_above / μ0)
    if MPI.Comm_rank(comm) == 0
        x, w = gausslegendre_01(16)
        Fd = 0.0                                                            # diffuse flux entering through the lid
        for (m, wm) in zip(x, w), ϕ in range(0, 2π; length=73)[1:72]
            Fd += wm * m * lat(dom.z[1], π - acos(m), ϕ) * (2π / 72)
        end
        println(@sprintf(" # SW top boundary (:column): extinction optical depth above the lid %.4f; direct beam at the lid %.1f W/m² (×%.4f), diffuse entering through the top %.2f W/m²",
                         τ_above, S0 * μ0 * att, att, Fd))
    end
    return att
end
