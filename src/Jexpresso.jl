"""
A research software for the numerical solution of a system of an arbitrary number of conservation laws using continuous spectral elements. DISCLAIMER: this is WIP and only 2D is being maintained until parallelization is complete.

If you are interested in contributing, please get in touch.
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)
"""
module Jexpresso

# Disable Julia's incremental precompilation for this module.
#
# Jexpresso loads several packages (GridapGmsh, GridapP4est, ClimaComms backends)
# that initialise native C/C++ shared libraries.  On macOS, the dynamic linker
# (dyld) crashes with a heap-corruption segfault (signal 11 in lsl::Allocator::free)
# when these libraries are loaded inside Julia's precompilation subprocess.
#
# __precompile__(false) stops Julia from spawning that subprocess.  The module
# is loaded directly instead, which is slightly slower on a cold start but has
# no effect when running with the --sysimage flag (the normal production path).
#
# Exception: when make_sysimage.jl is building the sysimage it sets
# JEXPRESSO_BUILDING_SYSIMAGE=1 and that env var is inherited by every child
# process PackageCompiler spawns (including the --output-o compilation
# subprocess).  We skip __precompile__(false) so PackageCompiler can include
# Jexpresso in the sysimage.  The resulting sysimage is then used with
# --sysimage, so the precompile path is never taken at runtime.
if !haskey(ENV, "JEXPRESSO_BUILDING_SYSIMAGE")
    __precompile__(false)
end

using MPI
using KernelAbstractions
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

# Core SciML Types and Functions
using OrdinaryDiffEq: OrdinaryDiffEq, solve, Tsit5, CarpenterKennedy2N54
using LinearSolve: LinearSolve, solve
using SciMLBase
using SciMLBase: CallbackSet, DiscreteCallback,
    ODEProblem, ODESolution, ODEFunction,
    SplitODEProblem, get_du, get_tmp_cache, 
    u_modified!, AbstractODEIntegrator, init, 
    step!, check_error, get_proposed_dt, 
    set_proposed_dt!, terminate!, remake

# Solvers and Linear Algebra
using OrdinaryDiffEq: OrdinaryDiffEq, solve, Tsit5 # Specific solvers help
using LinearSolve: LinearSolve, solve

# Utilities
using HDF5
using JLD2
using SnoopCompile # Only keep if you are actively profiling latency
using JLD2

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

using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Adaptivity
using Gridap.Geometry: GridMock
using GridapDistributed
using PartitionedArrays
using GridapGmsh
using GridapP4est

using DocStringExtensions

TInt   = Int64
TFloat = Float64
cpu    = true

include(joinpath( "..", "problems", "AbstractEquations.jl"))

include(joinpath( "macros", "je_macros.jl"))

include(joinpath( "auxiliary", "timing.jl"))

include(joinpath( "kernel", "globalStructs.jl"))

include(joinpath( "kernel", "abstractTypes.jl"))

include(joinpath( "kernel", "elementLearningStructs.jl"))

include(joinpath( "kernel", "ArtificialViscosity","viscousStructs.jl"))

include(joinpath( "kernel", "ArtificialViscosity","Wall_model.jl"))

include(joinpath( "kernel", "coupling", "couplingStructs.jl"))

include(joinpath( "kernel", "physics", "microphysicsStructs.jl"))

include(joinpath( "kernel", "physics", "microphysics.jl"))

include(joinpath( "kernel", "physics", "saturation.jl"))

include(joinpath( "kernel", "physics", "soundSpeed.jl"))

include(joinpath( "kernel", "physics", "globalConstantsPhysics.jl"))

include(joinpath( "kernel", "physics", "constitutiveLaw.jl"))

include(joinpath( "kernel", "physics", "large_scale.jl"))

include(joinpath( "kernel", "physics", "largescaleStructs.jl"))

include(joinpath( "kernel", "physics", "turbul.jl"))

include(joinpath( "kernel", "physics", "SGS.jl"))

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

include(joinpath( "kernel", "operators", "rhs_laguerre.jl"))

include(joinpath( "kernel", "operators", "filter.jl"))

include(joinpath( "kernel", "solvers", "TimeIntegrators.jl"))

include(joinpath( "kernel", "solvers", "Axb.jl"))

include(joinpath( "kernel", "Adaptivity", "Projection.jl"))

include(joinpath( "kernel", "mpi", "mpi_communications.jl"))

include(joinpath( "io", "mod_io.jl"))

include(joinpath( "io", "mod_print_io.jl"))

include(joinpath( "io", "write_output.jl"))

include(joinpath( "io", "diagnostics.jl"))

include(joinpath( "io", "Extract_topo.jl"))

include(joinpath( "io", "print_matrix.jl"))

include(joinpath( "io", "soundings.jl"))

include(joinpath( "auxiliary", "auxiliary_functions.jl"))

include(joinpath( "auxiliary", "checks.jl"))

include("./run.jl")

export @timers

end
