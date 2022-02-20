"""
Jexpresso, grid-based approximation of PDEs in the Julia programming language
This module provides rich set of tools for the numerical solution of PDE, mainly based
on finite element methods.
The module is structured in the following sub-modules:
- [`Jexpresso.Helpers`](@ref)
- [`Jexpresso.Polynomials`](@ref)
The exported names are:
"""

module Jexpresso

include("./basis/build_lgl.jl")

end # module
