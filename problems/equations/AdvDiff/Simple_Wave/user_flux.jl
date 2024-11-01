function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    F[1] = 1 * q[1]
    
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4, ip=1)
    F[1] = 1 * q[1]
end

function user_flux_gpu(q,qe,PhysConst,lpert)
    T = eltype(q)
    return T(q[1])
end
