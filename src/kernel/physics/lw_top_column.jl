# ── LW top boundary from the atmosphere above the lid (:RT_lw_top_bc => :column, the default) ──────────────────────
#
# :RT_lw_top_bc => :vacuum gives the old vacuum-like isotropic σT_space⁴/π (T_space = :RT_T_space, default 0 K). That
# is right only at the real top of the atmosphere: with the gray LW optics of atmos_to_rad_longwave the air above a low
# lid is opaque (BOMEX, 3 km lid: optical depth ≈ 5 above it, downward flux ≈ 340 W/m², σT_lid⁴ ≈ 359 W/m²), so the
# vacuum lid cools the top of the domain strongly (≈ −90 K/day at the BOMEX lid).
#
# With :RT_lw_top_bc => :column (default since 2026-11-24), every RT call builds a 1-D column above the lid and uses its downward radiance
# I↓(μ) at the lid as the top inflow, for each direction cosine μ = |cos θ|:
#   - lid state: domain mean of T, p and water vapour over the top-boundary nodes (collective);
#   - profile above the lid (:RT_lw_top_profile):
#       :standard (default)  AFGL tropical standard atmosphere (Anderson et al. 1986) up to 40 km, matched to the lid:
#                            T + ΔT e^{-(z - z_lid)/H} and qv · r^{e^{-(z - z_lid)/H}}, with ΔT, r the lid mismatches and
#                            H = :RT_lw_top_blend_m (default 3000 m); the usual practice of LES/CRM radiation codes
#                            (a background profile above the model top)
#       :lapse               T falls by :RT_lw_top_lapse (K/m, default 6.5e-3) to :RT_lw_top_tmin (default 200 K),
#                            qv decays with scale height :RT_lw_top_qv_scale_m (default 2000 m), up to 40 km
#     pressure is hydrostatic from the lid; levels every :RT_lw_top_dz_m (default 25 m);
#   - optics: the gray LW absorption of atmos_to_rad_longwave (dry air + vapour, pressure scaled), clear sky;
#   - radiance: the same exact layer recursion as the lateral columns (lw_column_radiance), from I = 0 at 40 km.
# The same I↓(μ) also feeds the lateral and initial-guess columns (lw_lateral.jl), so all LW inflows agree.
# SW is unchanged (the beam enters the lid at S₀; the clear-sky SW absorption above a low lid is a few percent).

const LW_GRAY_K_DRY = 0.002        # m²/kg, dry air (× p/p_ref), as atmos_to_rad_longwave
const LW_GRAY_K_VAP = 0.05         # m²/kg per kg/kg of vapour (× p/p_ref)
const LW_GRAY_P_REF = 101325.0

# top-inflow radiance as a function of μ = |cos θ| of a downward direction; nothing = the vacuum default (T_space)
const LW_TOP_FN = Ref{Any}(nothing)

# AFGL tropical standard atmosphere (Anderson et al. 1986, AFGL-TR-86-0110): z [km], T [K], H2O [ppmv]
const AFGL_TROP_Z  = [0.0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25,
                      27.5, 30, 32.5, 35, 37.5, 40]
const AFGL_TROP_T  = [299.7, 293.7, 287.7, 283.7, 277.0, 270.3, 263.6, 257.0, 250.3, 243.6, 237.0, 230.1, 223.6, 217.0,
                      210.3, 203.7, 197.0, 194.8, 198.8, 202.7, 206.7, 210.7, 214.6, 217.0, 219.2, 221.4, 227.0, 232.3,
                      237.7, 243.1, 248.5, 254.0]
const AFGL_TROP_H2O = [2.593e4, 1.949e4, 1.534e4, 8.600e3, 4.441e3, 3.346e3, 2.101e3, 1.289e3, 7.637e2, 4.098e2,
                       1.912e2, 7.306e1, 2.905e1, 9.900, 6.220, 4.000, 3.000, 2.900, 2.750, 2.600, 2.600, 2.650, 2.800,
                       2.900, 3.200, 3.250, 3.600, 4.000, 4.300, 4.600, 4.900, 5.200]

function lw_interp1(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    x <= xs[1] && return ys[1]
    x >= xs[end] && return ys[end]
    k = searchsortedlast(xs, x)
    t = (x - xs[k]) / (xs[k+1] - xs[k])
    return ys[k] + t * (ys[k+1] - ys[k])
end

"""
    lw_top_profile(z_lid, T_lid, p_lid, qv_lid; mode, blend, lapse, tmin, qv_scale, dz, ztop) -> (z, T, p, qv)

Column above the lid (ascending from z_lid to ztop), see the header of lw_top_column.jl.
"""
function lw_top_profile(z_lid::Float64, T_lid::Float64, p_lid::Float64, qv_lid::Float64; mode::Symbol=:standard,
                        blend::Float64=3000.0, lapse::Float64=6.5e-3, tmin::Float64=200.0, qv_scale::Float64=2000.0,
                        dz::Float64=25.0, ztop::Float64=40000.0)
    ztop > z_lid || error("lw_top_profile: the lid ($(z_lid) m) is above the profile top ($(ztop) m)")
    n = max(2, ceil(Int, (ztop - z_lid) / dz) + 1)
    z = collect(range(z_lid, ztop; length=n))
    T = similar(z); qv = similar(z)
    zk = AFGL_TROP_Z .* 1000.0
    qstd(zz) = lw_interp1(zk, AFGL_TROP_H2O, zz) * 1e-6 * 0.622
    if mode == :standard
        ΔT = T_lid - lw_interp1(zk, AFGL_TROP_T, z_lid)
        r = max(qv_lid, 1e-12) / qstd(z_lid)
        for k in eachindex(z)
            f = exp(-(z[k] - z_lid) / blend)
            T[k] = lw_interp1(zk, AFGL_TROP_T, z[k]) + ΔT * f
            qv[k] = qstd(z[k]) * r^f
        end
    elseif mode == :lapse
        for k in eachindex(z)
            T[k] = max(T_lid - lapse * (z[k] - z_lid), tmin)
            qv[k] = qv_lid * exp(-(z[k] - z_lid) / qv_scale)
        end
    else
        error(":RT_lw_top_profile must be :standard or :lapse")
    end
    p = similar(z); p[1] = p_lid
    Rd = 287.04; g = 9.80665
    for k = 2:n
        Tm = 0.5 * (T[k] * (1 + 0.61 * qv[k]) + T[k-1] * (1 + 0.61 * qv[k-1]))
        p[k] = p[k-1] * exp(-g * (z[k] - z[k-1]) / (Rd * Tm))
    end
    return z, T, p, qv
end

"""
    lw_top_column(z, T, p, qv) -> LWColumn

Gray extinction (atmos_to_rad_longwave's dry + vapour terms) and Planck radiance at the levels of the column above
the lid; zero downward radiance at its top.
"""
function lw_top_column(z, T, p, qv)
    κ = similar(z); B = similar(z)
    for k in eachindex(z)
        ρ = p[k] / (287.04 * T[k] * (1 + 0.61 * qv[k]))
        s = p[k] / LW_GRAY_P_REF
        κ[k] = ρ * (LW_GRAY_K_DRY * s + LW_GRAY_K_VAP * qv[k] * s)
        B[k] = LW_σ_SB * T[k]^4 / π
    end
    return LWColumn(z, κ, B, 0.0, 0.0)
end

"""
    lw_top_setup!(mesh, atmos_data, bdy, inputs, comm)

Collective. Sets LW_TOP_FN for this RT call: the downward radiance at the lid of the column above it (memoised per μ;
:RT_lw_top_bc => :column, the default), or nothing for :vacuum (σT_space⁴/π).
"""
function lw_top_setup!(mesh, atmos_data, bdy, inputs, comm)
    mode = Symbol(get(inputs, :RT_lw_top_bc, :column))
    mode in (:vacuum, :column) || error(":RT_lw_top_bc must be :vacuum or :column")
    if mode == :vacuum
        LW_TOP_FN[] = nothing
        return nothing
    end
    zl, T_lid, p_lid, qv_lid = rt_lid_state(mesh, atmos_data, bdy, comm)
    z, T, p, qv = rt_top_profile(inputs, zl, T_lid, p_lid, qv_lid)
    col = lw_top_column(z, T, p, qv)
    memo = Dict{Float64,Float64}()
    f = μ -> get!(memo, μ) do
        μ < 1e-6 ? col.B[1] : lw_column_radiance(col, -μ; top_fn=false)[1]
    end
    LW_TOP_FN[] = f
    if MPI.Comm_rank(comm) == 0
        τ = sum(0.5 * (col.κ[k] + col.κ[k+1]) * (z[k+1] - z[k]) for k = 1:length(z)-1)
        x, w = gausslegendre_01(16)
        Fd = 2π * sum(w[i] * x[i] * f(x[i]) for i in eachindex(x))
        println(@sprintf(" # LW top boundary (:column, %s): lid z %.0f m, T %.2f K, p %.0f Pa, qv %.2e; optical depth above %.2f; downward flux %.1f W/m² (σT_lid⁴ %.1f)",
                         get(inputs, :RT_lw_top_profile, :standard), zl, T_lid, p_lid, qv_lid, τ, Fd, LW_σ_SB * T_lid^4))
    end
    return nothing
end

# Domain-mean state at the lid (collective): height, T, p, water-vapour mixing ratio over the top-boundary nodes
function rt_lid_state(mesh, atmos_data, bdy, comm)
    r_mol = 18.015e-3 / 28.964e-3
    acc = zeros(4); ztop = -Inf
    for ip = 1:mesh.npoin
        bdy.is_top(mesh.z[ip]) || continue
        acc[1] += atmos_data.t_current[ip]; acc[2] += atmos_data.p_lay[ip]
        acc[3] += atmos_data.vmr_h2o[ip] * r_mol; acc[4] += 1; ztop = max(ztop, mesh.z[ip])
    end
    s = MPI.Allreduce(acc, MPI.SUM, comm); zl = MPI.Allreduce(ztop, MPI.MAX, comm)
    s[4] > 0 || error("rt_lid_state: no top-boundary nodes found")
    return zl, s[1] / s[4], s[2] / s[4], s[3] / s[4]
end

# The column above the lid shared by the LW and SW top conditions (settings :RT_lw_top_profile, :RT_lw_top_*)
rt_top_profile(inputs, zl, T_lid, p_lid, qv_lid) =
    lw_top_profile(zl, T_lid, p_lid, qv_lid;
                   mode = Symbol(get(inputs, :RT_lw_top_profile, :standard)),
                   blend = Float64(get(inputs, :RT_lw_top_blend_m, 3000.0)),
                   lapse = Float64(get(inputs, :RT_lw_top_lapse, 6.5e-3)),
                   tmin = Float64(get(inputs, :RT_lw_top_tmin, 200.0)),
                   qv_scale = Float64(get(inputs, :RT_lw_top_qv_scale_m, 2000.0)),
                   dz = Float64(get(inputs, :RT_lw_top_dz_m, 25.0)))

# Gauss-Legendre nodes and weights on [0, 1]
function gausslegendre_01(n::Int)
    β = [k / sqrt(4k^2 - 1) for k = 1:n-1]
    E = eigen(SymTridiagonal(zeros(n), β))
    return 0.5 .* (E.values .+ 1), E.vectors[1, :] .^ 2
end
