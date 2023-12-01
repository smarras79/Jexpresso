function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, x=1)

    physConst = PhysicalConst{Float64}()
    γ  = physConst.γ
    cv = physConst.cv
    R  = physConst.Rair
        
    ρ  = q[1]
    ρu = q[2]
    u  = ρu/ρ
    p = q[end]

    T = p/(R*ρ)
    e = cv*T
    
    A = 1.0 + 2.2*(x - 1.5)^2   # A/Athroat
        
    #=F[1] = q[2]
    F[2] = ((q[2]^2)/q[1]) + ((γ - 1.0)/γ)*(q[3] - (γ/2)*((q[2]^2)/q[1]))
    F[3] = ((γ*q[2]*q[3])/q[1]) - (γ*(γ - 1.0)*((q[2]^3)/(2*q[1]^2)))
    =#

    U1 = q[2]
    U2 = q[2]^2/q[1]
    U3 = q[3]
    
    F[1] = U2
    F[2] = U2^2/U1 + (γ - 1.0)*(U3 - 0.5*U2^2/U1)
    F[3] = U2*(U3 + (γ - 1.0)*(U3 - 0.5*U2^2/U1))/U1
    
end
