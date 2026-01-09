function user_source!(S, q, qe, npoin, ::CL, ::TOTAL; 
                      neqs=6, x=0.0, y=0.0, ymin=0.0, ymax=24000.0, 
                      ngl=5, nely=10, xmin=-50000.0, xmax=50000.0)
   
    PhysConst = PhysicalConst{Float64}()
    
    # Gravity source term
    ρ = q[1] - qe[1]
    
    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    S[5] = 0.0
    S[6] = 0.0

    # SPONGE LAYER (top only)
    zs = 16000.0  # Start sponge at 16 km
    
    if y > zs
        betay_coe = sinpi(0.5*(y - zs)/(ymax - zs))
    else
        betay_coe = 0.0
    end
    
    ctop = 0.1 * betay_coe
    
    # Apply damping
    S[1] -= ctop * (q[1] - qe[1])
    S[2] -= ctop * (q[2] - qe[2])
    S[3] -= ctop * (q[3] - qe[3])
    S[4] -= ctop * (q[4] - qe[4])
    S[5] -= ctop * (q[5] - qe[5])
    S[6] -= ctop * (q[6] - qe[6])
    
    return S
end 

function user_source!(S, q, qe, npoin, ::CL, ::PERT; 
                      neqs=6, x=0.0, y=0.0, ymin=0.0, ymax=24000.0, 
                      ngl=5, nely=10, xmin=-50000.0, xmax=50000.0)

    PhysConst = PhysicalConst{Float64}()

    # Gravity source term
    ρ = q[1]

    S[1] = 0.0
    S[2] = 0.0
    S[3] = -ρ*PhysConst.g
    S[4] = 0.0
    S[5] = 0.0
    S[6] = 0.0
    
    # SPONGE LAYER (top only)
    zs = 16000.0  # Start sponge at 16 km
    
    if y >= zs
        betay_coe = sinpi(0.5*(y - zs)/(ymax - zs))^2
    else
        betay_coe = 0.0
    end
    
    ctop = 0.1 * min(betay_coe, 1.0)
    
    # Apply damping
    S[1] -= ctop * q[1]
    S[2] -= ctop * q[2]
    S[3] -= ctop * q[3]
    S[4] -= ctop * q[4]
    S[5] -= ctop * q[5]
    S[6] -= ctop * q[6]
    
    return S
end

function user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)

    T = eltype(x)
    
    # SPONGE LAYER (top only)
    zs = T(16000.0)
    
    if y >= zs
        betay_coe = T(sinpi(T(0.5)*(y - zs)/(ymax - zs))^2)
    else
        betay_coe = T(0.0)
    end
    
    ctop = T(0.1) * min(betay_coe, T(1.0))
    
    # Gravity + damping
    ρ = q[1]
    
    return T(-ctop*q[1]), 
           T(-ctop*q[2]), 
           T(-ctop*q[3] - ρ*PhysConst.g), 
           T(-ctop*q[4]), 
           T(-ctop*q[5]), 
           T(-ctop*q[6])
end