include("../mesh/mod_mesh.jl")
include("../basis/basis_structs.jl")


mutable struct St_metrics{TFloat}
    dxdξ::Array{TFloat}
    J   ::Array{TFloat}
end


function build_metric_terms(SD::NSD_2D, mesh::St_mesh, basis::St_Lagrange, ξ)
    
    N, Q = size(ξ, 1)-1, size(ξ, 2)-1

    metrics = St_metrics{T}(zeros(Q+1,Q+1, mesh.nelem), zeros(Q+1,Q+1, mesh.nelem))
    
    
end
