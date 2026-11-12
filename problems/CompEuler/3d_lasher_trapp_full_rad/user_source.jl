# Top buffer layer (no lateral sponges: the domain is periodic and the case does not need them).
#   From z = 5.5 km (above the ~5 km maximum cloud top) to the lid at 8.1 km, rate ramping as
#   sin^2 to LT_SPONGE_TOP at the lid (timescale 1/LT_SPONGE_TOP = 20 s).
#   It relaxes u, v, w, hl, qt and ρqp toward the reference state. hl and qt are relaxed as SPECIFIC
#   quantities: the target for ρhl is hl_ref·ρ, not the reference ρhl. Relaxing the conserved ρhl
#   itself pins ρhl (hence p ≈ (R/cp)·ρhl) regardless of ρ, so the buffer removes the pressure
#   response to density changes; the air above 7 km then evacuated (ρ 0.55 → 0.13, T → 1200 K) and
#   the run blew up at t ≈ 1260 s. Damping momentum only (0.01 1/s) let moisture lifted to
#   7–8 km accumulate as ice/snow/graupel at the lid until the run blew up at t ≈ 4080 s.
#   At rest the perturbation is zero, so the buffer does nothing.
const LT_SPONGE_TOP = 0.05    # 1/s
const LT_Z_SPONGE   = 5500.0  # m

@inline function lt_sponge_rate(z, z_top)
    z >= LT_Z_SPONGE || return zero(z)
    r = (z - LT_Z_SPONGE) / (z_top - LT_Z_SPONGE)
    return LT_SPONGE_TOP * sinpi(r / 2)^2
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
    β = lt_sponge_rate(z, LT_DOMAIN[][5])   # lid height set in initialize()

    # q is the total state here. Reference momentum and ρqp are zero; hl, qt relaxed per unit mass.
    S[1] = 0.0
    S[2] = -β * q[2]
    S[3] = -β * q[3]
    S[4] = -q[1] * PhysConst.g - β * q[4]
    S[5] = -β * (q[5] - q[1] * qe[5] / qe[1])
    S[6] = -β * (q[6] - q[1] * qe[6] / qe[1])
    S[7] = -β * q[7]
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
    β = lt_sponge_rate(z, LT_DOMAIN[][5])   # lid height set in initialize()

    # q is the perturbation: ρ(hl - hl_ref) = ρhl' - hl_ref ρ', likewise for qt; ρqp_ref = 0.
    S[1] = 0.0
    S[2] = -β * q[2]
    S[3] = -β * q[3]
    S[4] = -q[1] * PhysConst.g - β * q[4]
    S[5] = -β * (q[5] - q[1] * qe[5] / qe[1])
    S[6] = -β * (q[6] - q[1] * qe[6] / qe[1])
    S[7] = -β * q[7]
end

function user_source_gpu(q, qe, x, y, z, PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    T = eltype(q)
    z_sponge = T(LT_Z_SPONGE)
    z_top    = T(8100.0)
    β        = T(0.0)
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β = T(LT_SPONGE_TOP) * sinpi(r / T(2))^T(2)
    end

    # ρ(hl - hl_ref) and ρ(qt - qt_ref): same expression in both forms; ρqp_ref = 0.
    d5 = q[5] - q[1] * qe[5] / qe[1]
    d6 = q[6] - q[1] * qe[6] / qe[1]
    d7 = q[7]

    return T(0.0), T(-β*q[2]), T(-β*q[3]),
           T(-q[1]*PhysConst.g - β*q[4]), T(-β*d5), T(-β*d6), T(-β*d7)
end

function user_scattering_functions(θ, θ1, ϕ, ϕ1, g)
    cos_Θ = sin(θ)*sin(θ1)*cos(ϕ - ϕ1) + cos(θ)*cos(θ1)
    cos_Θ = clamp(cos_Θ, -1.0, 1.0)
    return (1 - g^2) / ((4*π) * (1 + g^2 - 2*g*cos_Θ))^(3/2)
end
