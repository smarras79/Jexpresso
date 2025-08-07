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

    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0 #Y is the vertical direction in 3D
    S[4] = -ρ*PhysConst.g
    S[5] = 0.0

end

function user_source_gpu(q,qe,x,y,z,PhysConst, xmax, xmin, ymax, ymin, zmax, zmin, lpert)

    T = eltype(q)
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]

    return T(0.0), T(0.0), T(0.0), T(-ρ*PhysConst.g), T(0.0)
end

function user_extinction!(κ, x, y, connijk, iel, extra_mesh, ngl)

    for j=1:ngl
        for i=1:ngl
            ip = connijk[iel,i,j]
            r = 2/6 - sqrt((x[ip]-5*3/8)^2+ (y[ip]-3*2/8)^2)
            κ[i,j] = 5.5/(1+exp(-15*r))
        end
    end

end

function user_extinction(x,y)

    r = 2/6 - sqrt((x-5*3/8)^2+ (y-3*2/8)^2)
    return 5.5/(1+exp(-15*r))
end

function user_scattering_coef(x,y)
  
    r = 2/6 - sqrt((x-5*3/8)^2+ (y-3*2/8)^2)
    return 0.7*5.5/(1+exp(-15*r))

end

function user_scattering_functions(θ,θ1, HG)

    g = 0.7
    return (1/HG)*(1-g^2)/((1+g^2-2*g*cos(θ-θ1))^(3/2))

end

function user_scattering!(σ, Φ, x, y, connijk, iel, nelem_ang, nop_ang, connijk_ang, coords_ang, ngl)
    for j=1:ngl
        for i=1:ngl
            ip = connijk[iel,i,j]
            r = 2/6 - sqrt((x[ip]-5*3/8)^2+ (y[ip]-3*2/8)^2)
            σ[i,j] = 0.7*5.5/(1+exp(-15*r))
        end
    end
    for e=1:nelem_ang
        for i=1:nop_ang[e]
            ip = connijk_ang[e,i]
            for e1=1:nelem_ang
                for i1=1:nop_ang[e1]
                    ip1 = connijk_ang[e1,i1]
                    θ = coords_ang[1,ip]
                    θ1 = coords_ang[1,ip1]
                    g = 0.7
                    Φ[ip,ip1] = (1/4*π)*(1-g^2)/((1+g^2-2*g*cos(θ-θ1))^(3/2)) 
                end
            end
        end
    end
end

function user_rhs(x,y,θ)

    return 0.0
end

function user_rad_bc(x,y,θ)
    if (x == 0.0 || y == 2.0)
        return exp(-((192/(2*π))*(θ-7*π/4))^2)
    else
        return 0.0
    end
end

