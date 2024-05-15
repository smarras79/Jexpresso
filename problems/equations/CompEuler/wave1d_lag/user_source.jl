function user_source!(S::SubArray{Float64},
                      q::SubArray{Float64}, 
                      qe::SubArray{Float64},
                      npoin::Int64,
                      ::CL, ::TOTAL;
                      neqs=1,x=0.0, y=0.0, ymin=0.0, ymax=30000.0, xmin = -120000, xmax =120000)
    
    PhysConst = PhysicalConst{Float64}()

    #
    # S(q(x)) = -ρg
    #

    if (x >= 2.49)#nsponge_points * dsy) #&& dbl >= 0.0)
        sponge_coe = 2.0/(1+exp((0.3*(xmax-2.5)-x+2.5)/((xmax)/18)))
    elseif (x <=-2.49)
        sponge_coe = 2.0/(1+exp(-(0.3*(xmin+2.5)-x-2.5)/((xmax)/18)))
    else
        sponge_coe = 0.0
    end
    cs = min(sponge_coe,1)
    S[1] = -(cs)*(q[1])
    S[2] = -(cs)*(q[2])
    return  S
end

function user_source_gpu(q,qe,x,PhysConst, xmax, xmin,lpert)

    T = eltype(q)
    if (x >= 2.49)#nsponge_points * dsy) #&& dbl >= 0.0)
        sponge_coe = 2.0/(1+exp((0.3*(xmax-2.5)-x+2.5)/((xmax)/18)))
    elseif (x <=-2.49)
        sponge_coe = 2.0/(1+exp(-(0.3*(xmin+2.5)-x-2.5)/((xmax)/18)))
    else
        sponge_coe = 0.0
    end
    cs = min(sponge_coe,1)
    return T(-cs*q[1]), T(-cs*q[2])

end
