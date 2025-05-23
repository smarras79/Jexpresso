function user_source!(S,
                      q, 
                      qe,
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

function user_source!(S,
                      q, 
                      qe,
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


    return T(0.0), T(-buc + ρu_sponge), T(-bvc + ρv_sponge), T(-ρ*PhysConst.g + ρw_sponge), T(0.0)
end
