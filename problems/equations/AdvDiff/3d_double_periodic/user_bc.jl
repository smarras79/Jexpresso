function user_bc_dirichlet!(q,
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::TOTAL)

end

function user_bc_dirichlet!(q,
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe, ::PERT)

end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end

function user_bc_dirichlet_gpu(q,qe,x,y,z,t,nx,ny,nz,qbdy,lpert)
    return T(qbdy[1])
end
