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

function user_extinction(x, y)

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
                    Φ[ip,ip1] = (1/(3*π))*(1 + (cos(θ - θ1))^2)   
                end
            end
        end
    end
end

function user_scattering_coef(x,y)

    return 0.1*10*exp(-((x-3/2)/3)^2)*exp(-y/2)

end

function user_scattering_functions(θ,θ1,HG)

    return (1/(3*π))*(1 + (cos(θ - θ1))^2)

end

function user_rhs(x,y,θ)
    κip = 10*exp(-((x-3/2)/3)^2)*exp(-y/2)
    σip = 0.1*κip
    gip = exp(-((1. / 3) * (x - (3 / 3.)))^2)#exp(-((x-3/3)/3)^2)
    dgip = -(2. / 3^2) * (x - (3 / 3.)) * gip
    hip = exp(-4. * (2 - y) / 2)#exp(-4*(2-y)/2)
    dhip = (4. / 2) * hip
    sip = exp(-((96 / (2. * π)) * (θ - (7. * π / 5.)))^2)#exp(-((96/(2*π))*(θ-7*π/5))^2)
    uip = gip*hip*sip
    propip = (cos(θ)*dgip*hip+sin(θ)*gip*dhip)*sip
    scatter, error = quadgk(θ1 -> sip*(1/(3*π))*(1 + (cos(θ - θ1))^2), 0, 2*π, rtol=1e-13, atol = 1e-13) 
    return (-gip*hip*scatter*σip + κip*uip +  propip)
end

function user_rad_bc(x,y,θ)
    gip = exp(-((1. / 3) * (x - (3 / 3.)))^2)#exp(-((x-3/3)/3)^2)
    hip = exp(-4. * (2 - y) / 2)#exp(-4*(2-y)/2)
    sip = exp(-((96 / (2. * π)) * (θ - (7. * π / 5.)))^2)
    return gip*hip*sip
end
