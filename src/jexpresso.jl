"""
jexpresso, grid-based approximation of PDEs in the Julia programming language

This module provides rich set of tools for the numerical solution of PDE, mainly based
on finite element methods.

The module is structured in the following sub-modules:

- [`jexpresso.Io`](@ref)

The exported names are:
$(EXPORTS)
"""
module jexpresso

using DocStringExtensions

include("run.jl")

#include("Exports.jl")

end # module
