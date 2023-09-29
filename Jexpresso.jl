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

#using DocStringExtensions

include(joinpath("problems", "equations", "AbstractEquations.jl"))

include(joinpath("src", "kernel", "abstractTypes.jl"))

include(joinpath("src", "kernel", "globalStructs.jl"))

include(joinpath("src", "kernel", "physics", "globalConstantsPhysics.jl"))

include(joinpath("src", "kernel", "physics", "constitutiveLaw.jl"))

include(joinpath("src", "kernel", "mesh", "mesh.jl"))

include(joinpath("src", "kernel", "bases", "basis_structs.jl"))

include(joinpath("src", "kernel", "mesh", "metric_terms.jl"))

include(joinpath("src", "kernel", "infrastructure", "element_matrices.jl"))

include(joinpath("src", "kernel", "infrastructure", "sem_setup.jl"))

include(joinpath("src", "kernel", "infrastructure", "Kopriva_functions.jl"))

include(joinpath("src", "kernel", "boundaryconditions", "BCs.jl"))

include(joinpath("src", "kernel", "operators", "rhs.jl"))

include(joinpath("src", "kernel", "solvers", "TimeIntegrators.jl"))

include(joinpath("src", "kernel", "solvers", "Axb.jl"))

include(joinpath("src", "io", "mod_inputs.jl"))

include(joinpath("src", "io", "write_output.jl"))

include(joinpath("src", "io", "diagnostics.jl"))

include("./src/run.jl")

end
