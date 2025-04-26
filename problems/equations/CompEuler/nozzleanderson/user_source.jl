function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1, x=1, y=1, xmax=1, xmin=1, ymax=1)

    γ = 1.4
    γm1 = 0.4

    A  = 1.0 + 2.2*(x - 1.5)^2
    dAdx = 2.2*(2*x - 3)
    
    ρ  = q[1]/A
    u  = q[2]/q[1]
    T  = γm1*(q[3]/q[1] - 0.5*γ*u*u)
    p  = ρ*T

    S[1] = 0.0
    S[2] = p*dAdx/γ
    S[3] = 0.0
    
end
