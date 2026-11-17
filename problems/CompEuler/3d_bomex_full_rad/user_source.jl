#----------------------------------------------------------------------------------------------
# BOMEX forcing (Siebesma et al. 2003; values cross-checked against the MicroHH BOMEX case).
# Each term is switched on from user_inputs (:bomex_<term> => true); initialize() copies the
# switches into BOMEX_FORCING because user_source! has no access to `inputs`.
#----------------------------------------------------------------------------------------------
Base.@kwdef mutable struct St_BomexForcing
    wind::Bool = false            # initial wind + Coriolis + geostrophic forcing
    subsidence::Bool = false      # large-scale subsidence on hl and qt
    radiation::Bool  = false      # prescribed radiative cooling (stand-in while RT is off)
    qt_ls::Bool      = false      # large-scale moisture tendency (drying)
    sponge::Bool     = false      # top buffer layer (MicroHH BOMEX: z >= 2250 m)
end
const BOMEX_FORCING = St_BomexForcing()

function bomex_set_forcing!(inputs)
    BOMEX_FORCING.wind = get(inputs, :bomex_wind, false)
    BOMEX_FORCING.subsidence = get(inputs, :bomex_subsidence, false)
    BOMEX_FORCING.radiation  = get(inputs, :bomex_radiation, false)
    BOMEX_FORCING.qt_ls      = get(inputs, :bomex_qt_ls, false)
    BOMEX_FORCING.sponge     = get(inputs, :bomex_sponge, false)
    return BOMEX_FORCING
end

const BOMEX_F0 = 0.376e-4                       # Coriolis parameter [1/s]

# top buffer layer (MicroHH BOMEX: zstart = 2250 m, sigma = 0.00223 1/s, beta = 2)
const BOMEX_SPONGE_ZSTART = 2250.0
const BOMEX_SPONGE_SIGMA  = 0.00223
const BOMEX_SPONGE_ZTOP   = 3000.0

# Profiles below are piecewise in z (and the reference-table slopes are piecewise constant), so
# they are evaluated at z rounded to 0.1 mm: mesh generators leave ~1e-8 m noise in node heights,
# and a node level that sits exactly on a breakpoint (e.g. z = 1500 m) would otherwise fall on
# different branches in different columns, making the forcing horizontally non-uniform.
@inline bomex_zlevel(z) = round(z; digits = 4)

# geostrophic wind [m/s]
@inline bomex_geostrophic(z) = (-10.0 + 1.8e-3 * z, 0.0)

# radiative cooling dθl/dt [K/s]: -2 K/day below 1500 m, linear to 0 at 3000 m
@inline bomex_radcool(z) = z <= 1500.0 ? -2.0 / 86400.0 :
                           z <= 3000.0 ? (-2.0 + 2.0 * (z - 1500.0) / 1500.0) / 86400.0 : 0.0

# large-scale moisture tendency dqt/dt [kg/kg/s]: -1.2e-8 below 300 m, linear to 0 at 500 m
@inline bomex_qtls(z) = z <= 300.0 ? -1.2e-8 :
                        z <= 500.0 ? -1.2e-8 + 1.2e-8 * (z - 300.0) / 200.0 : 0.0

# large-scale subsidence w_s [m/s] and dw_s/dz [1/s]: 0 -> -0.65 cm/s at 1500 m -> 0 at 2100 m
@inline function bomex_subsidence(z)
    if z <= 1500.0
        return -0.0065 * z / 1500.0, -0.0065 / 1500.0
    elseif z <= 2100.0
        return -0.0065 + 0.0065 * (z - 1500.0) / 600.0, 0.0065 / 600.0
    else
        return 0.0, 0.0
    end
end

function user_source!(S,
                      q,
                      qe,
                      npoin,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)

    PhysConst = PhysicalConst{Float64}()

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]

    # Coriolis (f-plane, BOMEX latitude ≈ 15°N)
    u_geo = 0.0
    f0 = 0.376e-4  # s⁻¹
    buc = -f0 * ρv
    bvc =  f0 * (ρu - u_geo*ρ)

    # Top Rayleigh sponge (2400–3000 m)
    z_sponge = 2400.0
    z_top    = 3000.0
    α_top    = 0.75
    γ        = 2.0
    β_top    = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / 2)^γ
    end


    β = β_top

    S[1] = 0.0
    S[2] = buc - β * ρu
    S[3] = bvc - β * ρv
    S[4] = -ρ * PhysConst.g - β * ρw
    S[5] = 0.0
    S[6] = 0.0
    S[7] = 0.0
end

function user_source!(S,
                      q,
                      qe,
                      npoin,
                      ::CL, ::PERT;
                      neqs=1,
                      x=0.0, y=0.0, z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)

    PhysConst = PhysicalConst{Float64}()
    z = bomex_zlevel(z)            # see bomex_zlevel

    ρ  = q[1]                      # perturbation density (gravity term below)
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]

    # Coriolis + geostrophic forcing on the TOTAL momentum (f-plane):
    #   du/dt = f (v - v_g),   dv/dt = -f (u - u_g)
    buc = 0.0
    bvc = 0.0
    if BOMEX_FORCING.wind
        ρt  = q[1] + qe[1]
        ug, vg = bomex_geostrophic(z)
        buc =  BOMEX_F0 * ((q[3] + qe[3]) - ρt * vg)
        bvc = -BOMEX_F0 * ((q[2] + qe[2]) - ρt * ug)
    end

    # Top buffer layer (MicroHH BOMEX): above BOMEX_SPONGE_ZSTART every prognostic field is
    # relaxed to its initial state with rate β = σ ((z - zstart)/(ztop - zstart))^2.
    # In PERT form the initial state is ρu = ρ u_init(z), ρv = ρw = 0, ρhl' = ρqt' = 0
    # (the noise only acts below 300 m), so those are the targets; ρ' is not damped.
    β  = 0.0
    u0 = 0.0
    if BOMEX_FORCING.sponge && z >= BOMEX_SPONGE_ZSTART
        r  = (z - BOMEX_SPONGE_ZSTART) / (BOMEX_SPONGE_ZTOP - BOMEX_SPONGE_ZSTART)
        β  = BOMEX_SPONGE_SIGMA * r^2
        u0 = BOMEX_FORCING.wind ? bomex_wind(z)[1] : 0.0
    end
    ρtot = q[1] + qe[1]

    S[1] = 0.0
    S[2] = buc - β * (ρu - ρtot * u0)
    S[3] = bvc - β * ρv
    S[4] = -ρ * PhysConst.g - β * ρw
    S[5] = -β * q[5]
    S[6] = -β * q[6]
    S[7] = 0.0

    # prescribed radiative cooling, applied to the static energy as ρ cp dθl/dt
    if BOMEX_FORCING.radiation
        S[5] += (q[1] + qe[1]) * PhysConst.cp * bomex_radcool(z)
    end
    # large-scale moisture tendency
    if BOMEX_FORCING.qt_ls
        S[6] += (q[1] + qe[1]) * bomex_qtls(z)
    end

    # Subsidence  -ρ w_s ∂φ/∂z  for φ = hl, qt, split as φ = φ_ref(z) + φ' so no large terms cancel:
    #   reference part: -ρ w_s dφ_ref/dz                  (pointwise, from the reference table)
    #   deviation part: -∂(w_s ρ_ref φ')/∂z               (flux, in user_flux!)
    #                   + φ' (ρ_ref dw_s/dz + w_s dρ_ref/dz)  (source correction)
    if BOMEX_FORCING.subsidence
        ws, dws = bomex_subsidence(z)
        if ws != 0.0 || dws != 0.0
            R   = BOMEX_REF[]
            ρt  = q[1] + qe[1]
            ρr  = qe[1]
            c   = ρr * dws + ws * ref_drho_dz(R, z)
            dhl, dqt = ref_dphi_dz(R, z)
            hlp = (q[5] + qe[5]) / ρt - qe[5] / ρr
            qtp = (q[6] + qe[6]) / ρt - qe[6] / ρr
            S[5] += hlp * c - ρt * ws * dhl
            S[6] += qtp * c - ρt * ws * dqt
        end
    end
end

function user_source_gpu(q, qe, x, y, z, PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    T = eltype(q)
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]

    f0 = T(0.376e-4)
    buc = T(-f0 * ρv)
    bvc = T( f0 * ρu)

    # Top sponge
    z_sponge = T(2400.0)
    z_top    = T(3000.0)
    α_top    = T(0.75)
    β_top    = T(0.0)
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / T(2))^T(2)
    end

    # Lateral sponge
    α_lat    = T(0.25)
    lat_frac = T(0.15)
    β_x      = T(0.0)
    β_y      = T(0.0)
    dx = lat_frac * (xmax - xmin)
    dy = lat_frac * (ymax - ymin)
    if x < xmin + dx
        r = (xmin + dx - x) / dx
        β_x = α_lat * sinpi(r / T(2))^T(2)
    elseif x > xmax - dx
        r = (x - (xmax - dx)) / dx
        β_x = α_lat * sinpi(r / T(2))^T(2)
    end
    if y < ymin + dy
        r = (ymin + dy - y) / dy
        β_y = α_lat * sinpi(r / T(2))^T(2)
    elseif y > ymax - dy
        r = (y - (ymax - dy)) / dy
        β_y = α_lat * sinpi(r / T(2))^T(2)
    end
    β_lat = max(β_x, β_y)

    β = β_top + β_lat

    return T(0.0), T(buc - β*ρu), T(bvc - β*ρv),
           T(-ρ*PhysConst.g - β*ρw), T(0.0), T(0.0), T(0.0)
end

function user_scattering_functions(θ, θ1, ϕ, ϕ1, g)
    cos_Θ = sin(θ)*sin(θ1)*cos(ϕ - ϕ1) + cos(θ)*cos(θ1)
    cos_Θ = clamp(cos_Θ, -1.0, 1.0)
    return (1 - g^2) / ((4*π) * (1 + g^2 - 2*g*cos_Θ))^(3/2)
end
