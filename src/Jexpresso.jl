"""
A research software for the numerical solution of a system of an arbitrary number of conservation laws using continuous spectral elements. DISCLAIMER: this is WIP and only 2D is being maintained until parallelization is complete.

If you are interested in contributing, please get in touch.
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)
"""
module Jexpresso

using Revise
using BenchmarkTools
using Dates
using DelimitedFiles
using DataStructures
using LoopVectorization
using ElasticArrays
using InternedStrings
using LinearAlgebra
using StaticArrays
using StaticArrays: SVector, MVector, MArray, SMatrix, @SMatrix
using DiffEqBase
using DiffEqDevTools
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
using SnoopCompile
using SciMLBase: CallbackSet, DiscreteCallback,
                 ODEProblem, ODESolution, ODEFunction,
                 SplitODEProblem
import SciMLBase: get_du, get_tmp_cache, u_modified!,
                  AbstractODEIntegrator, init, step!, check_error,
                  get_proposed_dt, set_proposed_dt!,
                  terminate!, remake

using UnicodePlots
using Printf

const TInt   = Int64
const TFloat = Float64

using DocStringExtensions

include(joinpath( "..", "problems", "equations", "AbstractEquations.jl"))

include(joinpath( "macros", "je_macros.jl"))

include(joinpath( "kernel", "abstractTypes.jl"))

include(joinpath( "kernel", "globalStructs.jl"))

include(joinpath( "kernel", "physics", "globalConstantsPhysics.jl"))

include(joinpath( "kernel", "physics", "constitutiveLaw.jl"))

include(joinpath( "kernel", "mesh", "mesh.jl"))

include(joinpath( "kernel", "bases", "basis_structs.jl"))

include(joinpath( "kernel", "mesh", "metric_terms.jl"))

include(joinpath( "kernel", "infrastructure", "element_matrices.jl"))

include(joinpath( "kernel", "infrastructure", "sem_setup.jl"))

include(joinpath( "kernel", "infrastructure", "Kopriva_functions.jl"))

include(joinpath( "kernel", "boundaryconditions", "BCs.jl"))

include(joinpath( "kernel", "operators", "rhs.jl"))

include(joinpath( "kernel", "solvers", "TimeIntegrators.jl"))

include(joinpath("kernel", "operators", "rhs_laguerre.jl"))

include(joinpath("kernel", "operators", "filter.jl"))

include(joinpath( "kernel", "solvers", "Axb.jl"))

include(joinpath( "io", "mod_inputs.jl"))

include(joinpath( "io", "write_output.jl"))

include(joinpath( "io", "diagnostics.jl"))

include(joinpath( "auxiliary", "checks.jl"))

include("./run.jl")

end
