"""
Jexpresso is a toy software for the numerical solution of 1D, 2D, and 3D PDEs 
using spectral methods and spectral element methods on CPUs and GPUs. 
DISCLAIMER: this is WIP and we have barely started implementing this code.
The module is structured in the following sub-modules:

- [`Jexpresso.arrays`](@ref)
- [`Jexpresso.auxiliary`](@ref)
- [`Jexpresso.io`](@ref)
- [`Jexpresso.kernel`](@ref)
- [`Jexpresso.macros`](@ref)
- [`Jexpresso.problems`](@ref)

The exported names are:
$(EXPORTS)
"""
module Jexpresso

using DocStringExtensions

include("./run.jl")

end
