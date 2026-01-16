# spatial
function user_flux!(F::SubArray{Float64}, G::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh::St_mesh,
                    ::CL, ::TOTAL; neqs=1, ind, ip, ip_extra, mesh_extra::St_mesh)

    v = mesh_extra.x
    F[1] = v[ip_extra]*q[1]
    
end

# velocity
function user_flux_extra!(F_extra::SubArray{Float64}, G_extra::SubArray{Float64}, SD::NSD_1D,
                    q::SubArray{Float64},
                    qe::SubArray{Float64},
                    mesh_extra::St_mesh,
                    ::CL, ::TOTAL; neqs, ind, ip, ip_extra, ∇V)

    velocity = mesh_extra.x
    v = velocity[ip_extra]
    v_avg = 2.4
    mu = ( exp(-(v-v_avg)^2/2) + exp(-(v+v_avg)^2/2) )/(2*sqrt(2*pi))

    F_extra[1] = -∇V[1]*(q[1] + mu)
    
end

