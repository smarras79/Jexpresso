"""
A research software for the numerical solution of a system of an arbitrary number of conservation laws using continuous spectral elements. DISCLAIMER: this is WIP and only 2D is being maintained until parallelization is complete.

If you are interested in contributing, please get in touch.
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)
"""
module Jexpresso

using MPI
using KernelAbstractions
using Revise
using BenchmarkTools
using Dates
using DelimitedFiles
using CSV, DataFrames
using DataStructures
using LoopVectorization
using ElasticArrays
using Geodesy
using InternedStrings
using LinearAlgebra
using SpecialFunctions
using StaticArrays
using StaticArrays: SVector, MVector
using DiffEqBase
using DiffEqDevTools
using OrdinaryDiffEq
using OrdinaryDiffEq: solve
using SnoopCompile
using LinearSolve
using LinearSolve: solve
using SciMLBase: CallbackSet, DiscreteCallback,
                 ODEProblem, ODESolution, ODEFunction,
                 SplitODEProblem
using HDF5
import SciMLBase: get_du, get_tmp_cache, u_modified!,
                  AbstractODEIntegrator, init, step!, check_error,
                  get_proposed_dt, set_proposed_dt!,
                  terminate!, remake
import ClimaParams as CP
import Thermodynamics as TD
import Thermodynamics.Parameters as TP


import ClimaComms
@static pkgversion(ClimaComms) >= v"0.6" && ClimaComms.@import_required_backends

using RRTMGP
using RRTMGP.Vmrs
using RRTMGP.LookUpTables
using RRTMGP.AtmosphericStates
using RRTMGP.Optics
using RRTMGP.Sources
using RRTMGP.BCs
using RRTMGP.Fluxes
using RRTMGP.AngularDiscretizations
using RRTMGP.RTE
using RRTMGP.RTESolver
import RRTMGP.Parameters.RRTMGPParameters
using RRTMGP.ArtifactPaths

using UnicodePlots
using Printf
using NCDatasets

TInt   = Int64
TFloat = Float64
cpu    = true

using DocStringExtensions

include(joinpath( "..", "problems", "equations", "AbstractEquations.jl"))

include(joinpath( "macros", "je_macros.jl"))

include(joinpath( "kernel", "abstractTypes.jl"))

include(joinpath( "kernel", "elementLearningStructs.jl"))

include(joinpath( "kernel", "globalStructs.jl"))

include(joinpath( "kernel", "ArtificialViscosity","viscousStructs.jl"))

include(joinpath( "kernel", "ArtificialViscosity","Wall_model.jl"))

include(joinpath( "kernel", "physics", "microphysicsStructs.jl"))

include(joinpath( "kernel", "physics", "microphysics.jl"))

include(joinpath( "kernel", "physics", "saturation.jl"))

include(joinpath( "kernel", "physics", "soundSpeed.jl"))

include(joinpath( "kernel", "physics", "globalConstantsPhysics.jl"))

include(joinpath( "kernel", "physics", "constitutiveLaw.jl"))

include(joinpath( "kernel", "physics", "large_scale.jl"))

include(joinpath( "kernel", "physics", "largescaleStructs.jl"))

include(joinpath( "kernel", "physics", "turbul.jl"))

include(joinpath( "kernel", "physics", "CM_MOST.jl"))

include(joinpath( "kernel", "mesh", "Geom.jl"))

include(joinpath( "kernel", "mesh", "mesh.jl"))

include(joinpath( "kernel", "bases", "basis_structs.jl"))

include(joinpath( "kernel", "mesh", "metric_terms.jl"))

include(joinpath( "kernel", "infrastructure", "element_matrices.jl"))

include(joinpath( "kernel", "mesh", "phys_grid.jl"))

include(joinpath( "kernel", "infrastructure", "params_setup.jl"))

include(joinpath( "kernel", "infrastructure", "sem_setup.jl"))

include(joinpath( "kernel", "infrastructure", "Kopriva_functions.jl"))

include(joinpath( "kernel", "infrastructure", "convert_to_gpu.jl"))

include(joinpath( "kernel", "boundaryconditions", "surface_integral.jl"))

include(joinpath( "kernel", "boundaryconditions", "BCs.jl"))

include(joinpath( "kernel", "operators", "operators.jl"))

include(joinpath( "kernel", "operators", "rhs.jl"))

include(joinpath( "kernel", "operators", "rhs_2point.jl"))

include(joinpath( "kernel", "operators", "rhs_gpu.jl"))

include(joinpath( "kernel", "operators", "rhs_laguerre_gpu.jl"))

include(joinpath( "kernel", "operators", "imex2d.jl"))

include(joinpath( "kernel", "operators", "imex.jl"))

include(joinpath( "kernel", "operators", "rhs_laguerre.jl"))

include(joinpath( "kernel", "operators", "filter.jl"))

include(joinpath( "kernel", "solvers", "TimeIntegrators.jl"))

include(joinpath( "kernel", "solvers", "Axb.jl"))

include(joinpath( "kernel", "Adaptivity", "Projection.jl"))

include(joinpath( "io", "mod_inputs.jl"))

include(joinpath( "io", "mod_print_io.jl"))

include(joinpath( "io", "write_output.jl"))

include(joinpath( "io", "diagnostics.jl"))

include(joinpath( "io", "Extract_topo.jl"))

include(joinpath( "io", "print_matrix.jl"))

include(joinpath( "io", "soundings.jl"))

include(joinpath( "auxiliary", "auxiliary_functions.jl"))

include(joinpath( "auxiliary", "checks.jl"))

include("./run.jl")

# Run the test
# test_create_2d_projection_matrices_numa2d()
end
