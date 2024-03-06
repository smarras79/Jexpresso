function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=4)

    U = qe[2]
    H = qe[1]
    u = q[2]
    h = q[1]-qe[1]
    g = 9.81
    F[1] = U*h + H*u
    F[2] = g*h + U*u
    
end

function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::PERT; neqs=4)

    U = qe[2]
    H = qe[1]
    u = q[2]
    h = q[1]
    g = 9.81
    F[1] = U*h + H*u
    F[2] = g*h + U*u
end
