#=function user_flux(T, SD::NSD_1D, q::Vector, mesh::St_mesh; neqs=1)
    #
    # F(q(x)) = 0.8*q
    #
    F = @MVector zeros(T, neqs)
    
    F[1] = +1.0*q[1]
    
    return F
    
end
=#
function user_flux(T, SD::NSD_1D, q::Array, mesh::St_mesh; neqs=1)
    #
    # F(q(x)) = 0.8*q
    #
    F = zeros(T, mesh.npoin, 1:neqs)
    
    F[:, 1] .= +1.0*q[:, 1]
    
    return F
    
end
