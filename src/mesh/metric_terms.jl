include("../mesh/mesh.jl")
include("../basis/basis_structs.jl")


mutable struct St_metrics{TFloat}
    dxdξ::Union{Array{TFloat}, Missing}
    dydξ::Union{Array{TFloat}, Missing}
    dzdξ::Union{Array{TFloat}, Missing}
    dξdx::Union{Array{TFloat}, Missing}
    dηdx::Union{Array{TFloat}, Missing}
    dζdx::Union{Array{TFloat}, Missing}
    J   ::Array{TFloat}
end


function build_metric_terms(SD::NSD_2D, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ)
    
    metrics = St_metrics{T}(dxdξ = zeros(2, Q+1, Q+1, mesh.nelem), #∂x/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                            J    = zeros(   Q+1, Q+1, mesh.nelem)) #    J[1:Nq, 1:Nq, 1:nelem]
    
    
    for iel = 1:mesh.nelem
        for l = 1:Q+1
            for k = 1:Q+1
                for j = 1:N+1
                    for i = 1:N+1
                        nothing
                    end
                end
            end
        end
    end
        
    
end



function build_metric_terms(SD::NSD_3D, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ)
    
    metrics = St_metrics{T}(zeros(3, Q+1, Q+1, Q+1, mesh.nelem), #∂x/∂ξ[3, 1:Nq, 1:Nq, 1:Nq, 1:nelem]
                            zeros(   Q+1, Q+1, Q+1, mesh.nelem)) #    J[   1:Nq, 1:Nq, 1:Nq, 1:nelem]
    
    
end
