function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1, x=0)

    
    U1 = q[2]
    U2 = q[2]^2/q[1]
    U3 = q[3]
    
    S[1] = U2
    S[2] = U2^2/U1
    S[3] = U2*(U3 + (γ - 1.0)*(U3 - 0.5*U2^2/U1))/U1

end
