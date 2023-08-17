function user_flux!(F::SubArray{Float64}, SD::NSD_1D, q::SubArray{Float64}, mesh::St_mesh; neqs=1)
    F[1] = 1.0*q[1]
end
