# Sponge layers (momentum perturbations only).
#   Top: from z = 5.5 km (above the ~5 km maximum cloud top) to the lid at 8.1 km, rate ramping to
#        LT_SPONGE_TOP (timescale 1/LT_SPONGE_TOP = 100 s at the lid).
# No lateral sponges: the domain is periodic and the case does not need them.
# Previously 0.75 1/s at the top (with 0.25 1/s lateral sponges).
const LT_SPONGE_TOP = 0.01    # 1/s

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

    # Top Rayleigh sponge (5.5–8.1 km)
    z_sponge = 5500.0
    z_top    = LT_DOMAIN[][5]       # set in initialize()
    α_top    = LT_SPONGE_TOP
    γ        = 2.0
    β_top    = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / 2)^γ
    end

    β = β_top

    S[1] = 0.0
    S[2] = -β * ρu
    S[3] = -β * ρv
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

    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]

    # Top Rayleigh sponge (5.5–8.1 km)
    z_sponge = 5500.0
    z_top    = LT_DOMAIN[][5]       # set in initialize()
    α_top    = LT_SPONGE_TOP
    γ        = 2.0
    β_top    = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / 2)^γ
    end

    β = β_top

    S[1] = 0.0
    S[2] = -β * ρu
    S[3] = -β * ρv
    S[4] = -ρ * PhysConst.g - β * ρw
    S[5] = 0.0
    S[6] = 0.0
    S[7] = 0.0
end

function user_source_gpu(q, qe, x, y, z, PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)
    T = eltype(q)
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρw = q[4]

    # Top sponge
    z_sponge = T(5500.0)
    z_top    = T(8100.0)
    α_top    = T(LT_SPONGE_TOP)
    β_top    = T(0.0)
    if z >= z_sponge
        r = (z - z_sponge) / (z_top - z_sponge)
        β_top = α_top * sinpi(r / T(2))^T(2)
    end

    β = β_top

    return T(0.0), T(-β*ρu), T(-β*ρv),
           T(-ρ*PhysConst.g - β*ρw), T(0.0), T(0.0), T(0.0)
end

function user_scattering_functions(θ, θ1, ϕ, ϕ1, g)
    cos_Θ = sin(θ)*sin(θ1)*cos(ϕ - ϕ1) + cos(θ)*cos(θ1)
    cos_Θ = clamp(cos_Θ, -1.0, 1.0)
    return (1 - g^2) / ((4*π) * (1 + g^2 - 2*g*cos_Θ))^(3/2)
end
