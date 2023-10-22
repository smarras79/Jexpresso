<<<<<<< HEAD
function user_flux!(F::SubArray{Float64}, T, SD::NSD_1D, q::SubArray{Float64}, mesh::St_mesh; neqs=1)
    
    F .= 1.0*q
=======
function user_flux!(F::SubArray{Float64}, SD::NSD_1D, q::SubArray{Float64}, mesh::St_mesh; neqs=1)
    F[1] = 1.0*q[1]
>>>>>>> 64430ce5c650911b66b7ed8723529e7a6c591c7f
end
