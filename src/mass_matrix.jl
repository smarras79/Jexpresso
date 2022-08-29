using Test
using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using Plots; inspectdr()
using SparseArrays

include("./mesh/mod_mesh.jl")
include("./basis/basis_structs.jl")

abstract type AbstractIntegrationType end
struct Exact <: AbstractIntegrationType end
struct Inexact <: AbstractIntegrationType end


abstract type AbstractMassType end
mutable struct St_Mass_Exact{TFloat} <: AbstractMassType
    M::Array{TFloat, 2} #Dagonal matrix for inexact integration
end

mutable struct St_Mass_Inexact{TFloat} <: AbstractMassType
    M::Array{TFloat, 1} #Dagonal matrix for inexact integration
end
    
#Exact mass matrix

function build_element_matrices!(TP::Exact, ψ, ω, nelem::Int64, N::Int64, Q::Int64, T)

    el_matrices = St_Mass_Exact{T}(zeros(N+1, N+1))

    plt = plot() #Clear plot
    M::Array{TFloat, 2} = zeros(N+1, Q+1)
    #for iel=1:nelem
        for k=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    M[i, j] = M[i, j] + (1.0/8.0)*ω[k]*ψ[i,k]*ψ[j,k]
                    el_matrices.M[i, j] = el_matrices.M[i, j] + 0.5*ω[k]*ψ[i,k]*ψ[j,k]
                end
            end
        end
    #end

    @show el_matrices.M
    
    #display(plot!(spy(M)))
    return el_matrices
end

function build_element_matrices!(TP::Inexact, ψ, ω, nelem::Int64, N::Int64, Q::Int64, T)

    el_matrices = St_Mass_Inexact{T}(zeros(N+1))

    for iel=1:nelem
        for k=1:Q+1
            for i=1:N+1
                for j=1:N+1
                    if (i == j)
                        el_matrices.M[i] = el_matrices.M[i] + 0.5*ω[k]*ψ[i,k]*ψ[j,k]
                    end
                end
            end
        end
    end
    @show el_matrices.M
    
    #display(plot!(spy(M)))
    return el_matrices
    
end

