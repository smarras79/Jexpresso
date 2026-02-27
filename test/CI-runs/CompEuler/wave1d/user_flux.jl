function user_flux!(F::SubArray{TFloat}, G::SubArray{TFloat}, SD::NSD_1D,
                    q::SubArray{TFloat},
                    qe::SubArray{TFloat},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4, ip=1)

    F[1] = q[2]
    F[2] = q[1]
    
end
