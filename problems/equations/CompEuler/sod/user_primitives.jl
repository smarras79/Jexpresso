function user_primitives!(q::SubArray{TFloat},qe::SubArray{TFloat},uprimitive::SubArray{TFloat},::TOTAL)

    ρ = q[1]
    u = q[2]/ρ
    T = q[3]/ρ - 0.5*u^2
    E = q[3]/ρ

    Pr = 0.1
    γ  = 1.4
    μ  = 0.005
    κ  = Pr/(γ - 1.0)
    
    uprimitive[1] = ρ
    uprimitive[2] = u
    uprimitive[3] = T
    uprimitive[4] = E
    
end
