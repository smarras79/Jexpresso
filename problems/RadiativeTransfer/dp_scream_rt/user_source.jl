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
            κ[i,j] = 10*exp(-((x[ip]-3/2)/3)^2)*exp(-y[ip]/2)
        end
    end

end

function user_extinction(x, y, z)

            return 10*exp(-((x-3/2)/3)^2)*exp(-y/2)

end


function user_scattering!(σ, Φ, x, y, connijk, iel, nelem_ang, nop_ang, connijk_ang, coords_ang, ngl)
    for j=1:ngl
        for i=1:ngl
            ip = connijk[iel,i,j]
            σ[i,j] = 0.1*10*exp(-((x[ip]-3/2)/3)^2)*exp(-y[ip]/2)
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
                    Φ[ip,ip1] = 1.0#(1/(3*π))*(1 + (cos(θ - θ1))^2)   
                end
            end
        end
    end
end

function user_scattering_coef(x,y,z)

    return 0.1*10*exp(-((x-3/2)/3)^2)*exp(-y/2)

end

function user_scattering_functions(θ,θ1,ϕ,ϕ1,HG)

    return (1/(3*π))*(1 + (cos(θ - θ1))^2)*(1+(cos(ϕ-ϕ1))^2)

end

function user_rhs(x,y,z,θ,ϕ)
    return 0.0
end

function user_rhs_sphere(x,y,z,θ,ϕ)

    return 0.0
end

function user_rad_bc(x,y,z,θ,ϕ)
    if ((abs(z-20000))<1e-5) #&& abs(x-2500)>1 && abs(x-47500)>1 && abs(y-2500)>1 && abs(y-47500)>1)
        return exp(-((92/(2*π))*(θ-π))^2)
    else
        return 0.0
    end
end
