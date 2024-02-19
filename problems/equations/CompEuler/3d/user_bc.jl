function user_bc_dirichlet!(q::SubArray{Float64},
                            x::AbstractFloat, y::AbstractFloat, z::AbstractFloat,
                            t::AbstractFloat, tag,
                            qbdy::AbstractArray,
                            nx, ny, nz,
                            xmin, xmax,
                            ymin, ymax,
                            zmin, zmax,
                            qe::SubArray{Float64}, ::TOTAL)

    e = 1.0e-4
    
    xmine = xmin + e
    xmaxe = xmax - e
    ymine = ymin + e
    ymaxe = ymax - e
    zmine = zmin + e
    zmaxe = zmax - e

    if x < xmine || x > xmaxe 
        qbdy[2] = 0.0
    end
    
    if y < ymine || y > ymaxe 
        qbdy[3] = 0.0
    end

    if z < zmine || z > zmaxe 
        qbdy[4] = 0.0
    end
    
end

function user_bc_neumann(q::AbstractArray, gradq::AbstractArray, x::AbstractFloat, y::AbstractFloat, z::AbstractFloat, t::AbstractFloat, tag::String, inputs::Dict)
    flux = zeros(size(q,2),1)
    return flux
end
