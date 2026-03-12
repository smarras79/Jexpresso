function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,
                      x=0.0,
                      y=0.0,
                      z=0.0,
                      xmin=0.0,xmax=0.0,
                      ymin=0.0,ymax=0.0,
                      zmin=0.0,zmax=0.0)
    
    PhysConst = PhysicalConst{Float64}()
        
    #
    # S(q(x)) = -ρg
    #
    ρ  = q[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = 0.0
    S[4] = -ρ*PhysConst.g
    S[5] = 0.0
   
end

function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1,
                      x=0.0,
                      y=0.0,
                      z=0.0,
                      xmin=0.0, xmax=0.0,
                      ymin=0.0, ymax=0.0,
                      zmin=0.0, zmax=0.0)

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

function user_scattering_functions(θ,θ1,ϕ,ϕ1,g)
    # cos of scattering angle between directions (θ,ϕ) and (θ1,ϕ1)
    cos_Θ = sin(θ)*sin(θ1)*cos(ϕ - ϕ1) + cos(θ)*cos(θ1)
    cos_Θ = clamp(cos_Θ, -1.0, 1.0)   # guard against floating point outside [-1,1]

    # Henyey-Greenstein
    return (1 - g^2) / ((4*π)*(1 + g^2 - 2*g*cos_Θ))^(3/2)
end
