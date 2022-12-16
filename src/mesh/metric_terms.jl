include("../mesh/mesh.jl")
include("../basis/basis_structs.jl")


abstract type AbstractMetricForm end
struct COVAR <: AbstractMetricForm end
struct CNVAR <: AbstractMetricForm end

Base.@kwdef mutable struct St_metrics{TFloat}
    dxdξ::Union{Array{TFloat}, Missing} = zeros(1)
    dxdη::Union{Array{TFloat}, Missing} = zeros(1)
    dxdζ::Union{Array{TFloat}, Missing} = zeros(1)
    
    dydξ::Union{Array{TFloat}, Missing} = zeros(1)
    dydη::Union{Array{TFloat}, Missing} = zeros(1)
    dydζ::Union{Array{TFloat}, Missing} = zeros(1)

    dzdξ::Union{Array{TFloat}, Missing} = zeros(1)
    dzdη::Union{Array{TFloat}, Missing} = zeros(1)
    dzdζ::Union{Array{TFloat}, Missing} = zeros(1)
    
    dξdx::Union{Array{TFloat}, Missing} = zeros(1)
    dξdy::Union{Array{TFloat}, Missing} = zeros(1)
    dξdz::Union{Array{TFloat}, Missing} = zeros(1)
    
    dηdx::Union{Array{TFloat}, Missing} = zeros(1)
    dηdy::Union{Array{TFloat}, Missing} = zeros(1)
    dηdz::Union{Array{TFloat}, Missing} = zeros(1)

    dζdx::Union{Array{TFloat}, Missing} = zeros(1)
    dζdy::Union{Array{TFloat}, Missing} = zeros(1)
    dζdz::Union{Array{TFloat}, Missing} = zeros(1)
    
    Je  ::Array{TFloat} = zeros(1)
end


function build_metric_terms(SD::NSD_2D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(Q+1, Q+1, mesh.nelem), #∂x/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                            dxdη = zeros(Q+1, Q+1, mesh.nelem), #∂x/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                            dydξ = zeros(Q+1, Q+1, mesh.nelem), #∂y/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                            dydη = zeros(Q+1, Q+1, mesh.nelem), #∂y/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                            Je   = zeros(Q+1, Q+1, mesh.nelem)) #   Je[1:Nq, 1:Nq, 1:nelem]

    ψ  = basis.ψ
    dψ = basis.dψ

    @info " metric terms WIP"
    
    for iel = 1:mesh.nelem
        for l = 1:Q+1
            for k = 1:Q+1
                for j = 1:N+1
                    for i = 1:N+1
                        ip = mesh.connijk[i,j,iel] 
                        xij = mesh.x[ip]
                        metrics.dxdξ[k, l, iel] = metrics.dxdξ[1, k, l, iel] + dψ[i,k]*ψ[j,l]*xij
                        metrics.dxdη[k, l, iel] = metrics.dxdη[1, k, l, iel] + ψ[i,k]*dψ[j,l]*xij
                        
                        metrics.dydξ[k, l, iel] = metrics.dydξ[1, k, l, iel] + dψ[i,k]*ψ[j,l]*xij
                        metrics.dydη[k, l, iel] = metrics.dydη[1, k, l, iel] + ψ[i,k]*dψ[j,l]*xij
                    end
                end
            end
            Je[k, l, iel] = metrics.dxdξ[k, l, iel]*metrics.dydη[k, l, iel] - metrics.dydξ[k, l, iel]*metrics.dxdη[k, l, iel]
        end
    end

    
    show(stdout, "text/plain", Je)
    
    return metrics
end

function build_metric_terms(SD::NSD_3D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ)

     metrics = St_metrics{T}(dxdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dydζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             Je   = zeros(Q+1, Q+1, Q+1, mesh.nelem))  #   Je[1:Nq, 1:Nq, 1:nelem]
    
    return metrics
end
