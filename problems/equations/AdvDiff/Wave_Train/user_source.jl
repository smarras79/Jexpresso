function user_source!(S,
                      q, 
                      qe,
                      npoin::Int64,
                      ::CL, ::PERT;
                      neqs=1, x=0.0, y=0.0, xmin=0.0, xmax=1.0)
    
    PhysConst = PhysicalConst{Float64}()
    xsponge = 5000.0
    if (x >= 5000.0)#nsponge_points * dsy) #&& dbl >= 0.0)
        sponge_coe = 2.0/(1+exp((0.3*(xmax-2.5)-x+2.5)/((xmax)/18)))
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
