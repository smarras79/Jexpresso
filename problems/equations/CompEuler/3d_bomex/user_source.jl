import Thermodynamics as TD
function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1)
    
    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0 #Y is the vertical direction in 3D
    S[4] = -ρ*PhysConst.g
    S[5] = 0.0
   
end

function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1)

    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    # Coriolis
    u_geostrophic::T = -10.0
    u_slope::T = 1.8e-3
    u_geo = u_geostrophic + u_slope * z

    f0 = T(0.376e-4)
    ρu = q[2] - u_geo * q[1]
    ρv = q[3]
    ρw = q[4]
    buc::T= -f0*ρv
    bvc::T= f0*ρu
    bwc::T= 0.0

    # sponge layer
    z_sponge::T = 2400
    z_max::T = 3000
    α_max::T = 0.75
    γ::T = 2
    β_sponge::T = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α_max * sinpi(r / 2)^γ
    end
    ρu_sponge::T = - β_sponge * (ρu)

    S[1] = 0.0
    S[2] = buc + ρu_sponge
    S[3] = bvc #Y is the vertical direction in 3D
    S[4] = -ρ*PhysConst.g
    S[5] = 0.0

end

function user_source_gpu(q,qe,x,y,z,PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)

    T = eltype(q)
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    # Coriolis
    # u_geostrophic::T = -10.0
    # u_slope::T = 1.8e-3
    # u_geo = u_geostrophic + u_slope * z
    u_geo::T = 0.0
    f0 = T(0.376e-4)
    ρu = q[2] - u_geo * q[1]
    ρv = q[3]
    ρw = q[4]
    buc::T= -f0*ρv
    bvc::T= f0*ρu
    bwc::T= 0.0

    # sponge layer
    z_sponge::T = 2400
    z_max::T = 3000
    α_max::T = 0.75
    γ::T = 2
    β_sponge::T = 0.0
    if z >= z_sponge
        r = (z - z_sponge) / (z_max - z_sponge)
        β_sponge = α_max * sinpi(r / 2)^γ
    end
    ρu_sponge::T = - β_sponge * (ρu)
    ρv_sponge::T = - β_sponge * (ρv)
    ρw_sponge::T = - β_sponge * (ρw)


    return T(0.0), T(-buc + ρu_sponge), T(-bvc + ρv_sponge), T(-ρ*PhysConst.g + ρw_sponge), T(0.0), T(0.0), T(0.0)
end


function user_saturation_adjustment(q, qe, z, param_set, lpert)

    T = eltype(q)
    @inbounds begin
        ρ = q[1]
        ρu = q[2]
        ρv = q[3]
        ρw = q[4]
        ρθ = q[5]

        u = ρu/ρ
        v = ρv/ρ
        w = ρw/ρ
        θ = ρθ/ρ
        qt = q[6]/ρ
        ql = q[7]/ρ
        e_kin = T(1 // 2) * (u^2 + v^2 + w^2)
        _grav = T(TP.grav(param_set))
        e_pot = _grav * z
        e_int = θ - e_kin - e_pot
        TS = TD.PhaseEquil_ρeq(param_set, ρ, e_int, qt)

        ρ = TD.air_density(param_set, TS)
        q_pt = TD.PhasePartition(param_set, TS)
        e_tot = TD.total_energy(param_set, TS, e_kin, e_pot)
    end
    return ρ, ρ*u, ρ*v, ρ*w, ρ*e_tot, ρ*qt, ρ*q_pt.liq
    # return T(q[1]), T(q[2]), T(q[3]), T(q[4]), T(q[5]), T(q[6]), T(q[7])
    

end
