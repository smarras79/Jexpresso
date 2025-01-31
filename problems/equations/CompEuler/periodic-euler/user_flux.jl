function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_2D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ip=1)
    
    u = 1.0
    v = 0.0
    
    qu  = u*q[1]
    qv  = v*q[1]
        
    F[1] = qu
    G[1] = qv
    
end


function user_flux_gpu(q,qe,PhysConst,lpert)
    
    T = eltype(q) 
    return T(q[1]), T(0.0)

end
