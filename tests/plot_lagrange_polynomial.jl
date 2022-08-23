include("../src/basis/basis_structs.jl")
include("../src/IO/plotting/jeplots.jl")


"""
    Test driver to show how to build and plot a
    Lagrange polynomial of order N inside reference element [-1,1]
"""



function plot_lagrange_polynomial(N, TFloat)

    nx  = 100
    ξr  = range(-1, 1, length=nx)
    l   = zeros(TFloat, nx)
    lgl = St_lgl{TFloat}(zeros(TFloat, N+1),
                         zeros(TFloat, N+1))
    
    build_Integration_points!(lgl, N)
    ξ = lgl.ξ
    ω = BarycentricWeights(ξ)

    for i=1:nx
        ψ = LagrangeInterpolatingPolynomials(ξr[i], ξ, ω)
        l[i] = ψ[1]
    end
    
    plot_curve(ξr, l, "Basis", "ξ", "ψ", [""], :none)
    
end

#Run it:
N = 8

plot_lagrange_polynomial(N, Float64)
