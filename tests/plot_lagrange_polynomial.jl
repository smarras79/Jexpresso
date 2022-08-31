using Plots; gr()
using LaTeXStrings
using ColorSchemes

include("../src/basis/basis_structs.jl")

"""
    Test driver to show how to build and plot a
    Lagrange polynomial of order N inside reference element [-1,1]
"""


"""
    basis has size (N+1, Q+1)
    ξ  is the array of N+1 Interpolation points (e.g. LGL)
    ξq is the array of Q+1 quadrature points
"""
function plot_curve(ξ, ψ, title)
    
    plt = plot() #Clear plot
    display(plot!(ξ, ψ, title = title, legend=false, lw = 3,
                  xtickfontsize=16, ytickfontsize=16, reuse=false,
                  xlabel="ξ",ylabel="ψ(ξ)"))
    
end

function plot_basis(ψ, ξ, ξq)

    N = size(ψ, 1) - 1
    
    plt = plot() #Clear plot
    for i=1:N+1
        display(plot!(ξq[:], ψ[i,:], title = "Bases", legend=false, lw = 3,
                      xtickfontsize=16, ytickfontsize=16, reuse=false,
                      xlabel="ξ",ylabel="ψᵢ(ξ)"))
        
        scatter!(ξ, zeros(N+1), markersizes=4)
    end    
end


function plot_basis22(ξ, ω_barycentric, Q, N, TFloat)

    ξr  = range(-1, 1, length=Q+1)
    l   = zeros(TFloat, Q+1)
    lgl = St_lgl{TFloat}(zeros(TFloat, N+1),
                         zeros(TFloat, N+1))
    
    plt = plot() #Clear plot
    for k=1:N+1
        
        for i=1:Q+1
            ψ = LagrangeInterpolatingPolynomials(ξr[i], ξ, ω_barycentric)
            l[i] = ψ[k]
        end        
        display(plot!(ξr, l, title = "Bases", legend=false, lw = 3,
                      xtickfontsize=16, ytickfontsize=16, reuse=false,
                      xlabel="ξ",ylabel="ψᵢ(ξ)"))

    end
    
end


function plot_basis3(ξ, ω_barycentric, Q, N, TFloat)

    ξr  = range(-1, 1, length=Q+1)
    l   = zeros(TFloat, Q+1)
    lgl = St_lgl{TFloat}(zeros(TFloat, N+1),
                         zeros(TFloat, N+1))
    
    plt = plot() #Clear plot
    for k=1:N+1
        
        for i=1:Q+1
            ψ = LagrangeInterpolatingPolynomials(ξr[i], ξ, ω_barycentric)
            l[i] = ψ[k]
        end        
        display(plot!(ξr, l, title = "Bases", legend=false, lw = 3,
                      xtickfontsize=16, ytickfontsize=16, reuse=false,
                      xlabel="ξ",ylabel="ψᵢ(ξ)"))

    end
    
end


function plot_basis4(ξ, ξq, N, Q, TFloat)
    
    ψ = zeros(TFloat, N+1, Q+1)
    dψdx = ψ
    
    plt = plot() #Clear plot
    for k=1:Q+1        
        for i=1:N+1
            (ψ, dψdx) = LagrangeInterpolatingPolynomials_classic(ξ, ξq, N, Q, TFloat)
        end      
    end

    for i=1:N+1
        display(plot!(ξq[:], ψ[i,:], title = "Bases", legend=false, lw = 3,
                      xtickfontsize=16, ytickfontsize=16, reuse=false,
                      xlabel="ξ",ylabel="ψᵢ(ξ)"))
    end
    
end

# 
# Run it:
#
#N = 6
#plot_lagrange_polynomial(N, Float64)
