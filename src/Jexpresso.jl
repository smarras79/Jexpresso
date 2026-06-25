"""
A research software for the numerical solution of a system of an arbitrary number of conservation laws using continuous spectral elements. DISCLAIMER: this is WIP and only 2D is being maintained until parallelization is complete.

If you are interested in contributing, please get in touch.
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)
"""
module Jexpresso

using QuadGK
using MPI
using KernelAbstractions
# PERF: the following five `using`s eagerly load big package binaries
# on every MPI rank (Revise ~150 MB, BenchmarkTools, SnoopCompile,
# UnicodePlots, Geodesy) and were not referenced anywhere in the source
# tree. Removed to cut the per-rank baseline. Re-add at the REPL if
# you need them interactively.
using Revise
# using BenchmarkTools
using Dates
using CSV, DataFrames
# PERF: removed `using DataStructures`, `using DelimitedFiles`,
# `using InternedStrings`, `using ElasticArrays` — none of these
# packages were actually referenced anywhere in src/. They each cost
# load time on every Julia start; re-add only when an actual user
# turns up.
using LoopVectorization
# using Geodesy
using LinearAlgebra
# PERF: `using LinearOperators` moved into _ensure_rt_loaded!() —
# the only call site is src/kernel/operators/Axb_rad_mpi.jl, which is
# part of the radiative-transfer code path.
using SpecialFunctions
using SparseArrays
using StaticArrays
using StaticArrays: SVector, MVector
using OrdinaryDiffEq
using OrdinaryDiffEq: solve
# using SnoopCompile
using LinearSolve
using LinearSolve: solve
using SciMLBase: CallbackSet, DiscreteCallback,
                 ODEProblem, ODESolution, ODEFunction,
                 SplitODEProblem, FullSpecialize
# PERF: `using HDF5` moved into Jexpresso._ensure_hdf5_loaded!() —
# only `write_hdf5`/`read_hdf5` in src/io/write_output.jl reference
# h5* functions, and both call the loader before touching them. The
# local sentinel `HDF5()` (abstractTypes.jl) is just a dispatch tag
# and does not need the package.
import SciMLBase: get_du, get_tmp_cache, u_modified!,
                  AbstractODEIntegrator, init, step!, check_error,
                  get_proposed_dt, set_proposed_dt!,
                  terminate!, remake
import ClimaParams as CP
import Thermodynamics as TD
import Thermodynamics.Parameters as TP

# PERF: the radiative-transfer dependency tree (RRTMGP + 11
# submodules, ClimaComms with its MPI-backend macro, NCDatasets) used
# to be `using`-d at the top of this file. They are
# only referenced from:
#
#   * src/kernel/mesh/phys_grid.jl :: compute_radiative_fluxes!(...)
#   * src/io/read_dp_scream.jl     :: read_atmospheric_data + helpers
#   * src/kernel/operators/Axb_rad_mpi.jl :: LinearOperators
#
# All three are gated by `inputs[:lRT_problem]` (set in user_inputs.jl).
# city2d and every other non-RT case pays the load cost for nothing.
#
# Defer to `_ensure_rt_loaded!()` (see bottom of this file) which is
# called from drivers.jl when an RT run is detected. Function bodies
# parse fine without the names in scope — Julia resolves them at first
# call, by which time the lazy loader will have run.
using Serialization

# PERF: UnicodePlots was only referenced from a commented-out
# `UnicodePlots.heatmap` debug call in Projection.jl. Removed from the
# eager load path — the import alone adds ~50 MB per rank.
# using UnicodePlots
using Printf

using Gridap
using Gridap.Arrays
using Gridap.Arrays: Table
using Gridap.Geometry
using Gridap.Fields
using Gridap.ReferenceFEs
using Gridap.CellData
using Gridap.Adaptivity
# PERF: `using Gridap.Geometry: GridMock` was unused (only re-imported
# in src/kernel/mesh/Geom.jl, never referenced). Removed.
using GridapDistributed
using PartitionedArrays
using GridapGmsh
# PERF: `using GridapP4est` (and its companion P4est_wrapper) moved
# into Jexpresso._ensure_amr_loaded!(). They are only touched by AMR
# / initial-refinement / restart-AMR code paths, all gated by
# `inputs[:ladapt]`, `inputs[:linitial_refine]`, `inputs[:lamr]` or
# `inputs[:lrestart_amr]`. city2d and other plain runs skip the load.

using PrecompileTools


TInt   = Int64
TFloat = Float64
cpu    = true

using DocStringExtensions

include(joinpath( "..", "problems", "AbstractEquations.jl"))

include(joinpath( "macros", "je_macros.jl"))

include(joinpath( "auxiliary", "timing.jl"))

include(joinpath( "kernel", "abstractTypes.jl"))

include(joinpath( "kernel", "mesh", "meshStructs.jl"))

include(joinpath( "kernel", "elementLearningStructs.jl"))

include(joinpath( "kernel", "globalStructs.jl"))

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

include(joinpath( "kernel", "physics", "atmos_to_rad.jl"))

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

include(joinpath( "kernel", "solvers", "IMEXTimeIntegrators.jl"))

include(joinpath("kernel", "operators", "Axb_rad_mpi.jl"))

include(joinpath( "kernel", "solvers", "Axb.jl"))

include(joinpath("kernel", "operators", "build_rad_2d.jl"))

include(joinpath("kernel", "operators", "build_rad_3d.jl"))

include(joinpath( "kernel", "operators", "angular_comms.jl"))

include(joinpath( "kernel", "operators", "extra_amr_matrices.jl"))

include(joinpath( "kernel", "operators", "debug_amr_parallel.jl"))

include(joinpath( "kernel", "operators", "mass_assembly_jacc.jl"))

include(joinpath( "kernel", "Adaptivity", "Projection.jl"))

include(joinpath( "kernel", "mpi", "mpi_communications.jl"))

include(joinpath( "io", "mod_inputs.jl"))

include(joinpath( "io", "les_statistics.jl"))

include(joinpath( "io", "mod_print_io.jl"))

include(joinpath( "io", "write_output.jl"))

include(joinpath( "io", "diagnostics.jl"))

include(joinpath( "io", "Extract_topo.jl"))

include(joinpath( "io", "print_matrix.jl"))

include(joinpath( "io", "soundings.jl"))

include(joinpath( "io", "read_dp_scream.jl"))

include(joinpath( "auxiliary", "auxiliary_functions.jl"))

include(joinpath( "auxiliary", "checks.jl"))

# PERF: only evaluate run.jl at module-body time when this file was
# invoked as a Julia SCRIPT (`julia src/Jexpresso.jl CompEuler 3d`).
# When the package is being precompiled the @compile_workload block at
# the bottom of this file already includes run.jl with sod1d args, and
# including it from here too would define `parse_commandline` twice
# during precompile and violate Julia ≥ 1.10's no-method-overwrite
# rule:
#
#   WARNING: Method definition parse_commandline() in module Jexpresso
#   ... overwritten on the same line (check for duplicate calls to
#   `include`).
#   ERROR: Method overwriting is not permitted during Module
#   precompilation.
#
# When the package is loaded via `using Jexpresso` (precompile cache
# hit) the module body doesn't re-evaluate at all, so this branch is
# moot — REPL users invoke `Jexpresso.run_case(eqs, eqs_case)` to
# start a case (defined further down).
#
# The script-form `julia src/Jexpresso.jl CompEuler 3d` still works:
# `abspath(PROGRAM_FILE)` matches this file's path AND
# `jl_generating_output` returns 0 (we're not generating precompile
# output), so the include fires and the historic auto-run path is
# preserved.
if abspath(PROGRAM_FILE) == abspath(@__FILE__) &&
   ccall(:jl_generating_output, Cint, ()) == 0
    include("./run.jl")
end

export @timers

# ──────────────────────────────────────────────────────────────────────
# Lazy dependency loaders.
#
# Heavy packages whose symbols are only touched by a specific feature
# (radiative transfer, HDF5 I/O, NetCDF I/O, AMR) are kept out of the
# top-level `using` block and brought into Jexpresso's namespace on
# demand by these helpers. The contract is that each `_ensure_*!()`
# call MUST happen before any code that references the corresponding
# package symbols — i.e. at the top of the entry-point function for
# that feature. Function bodies parse fine without the names in scope:
# Julia resolves them at the first call, by which time the loader has
# run.
#
# Why this design (as opposed to Julia 1.9 package extensions): the
# feature code already lives inside Jexpresso (phys_grid.jl,
# Axb_rad_mpi.jl, write_output.jl, mesh.jl). Splitting them into ext/
# would require moving those files and re-wiring the include chain.
# Deferred `@eval using` keeps the source layout intact and saves the
# same load time for the common non-feature run (city2d et al.).
#
# Macros are an exception: macroexpansion happens at parse time, so any
# `Package.@macro` reference in a function body would fail to expand
# even if the call site is dead code. None of the deferred packages
# are referenced this way at present; if a future caller needs e.g.
# `Infiltrator.@exfiltrate` they must either eager-load it or wrap the
# call in `@eval`.
# ──────────────────────────────────────────────────────────────────────

# ─── Radiative-transfer dependency tree ───────────────────────────────
# Touched only by:
#   * src/kernel/mesh/phys_grid.jl  :: compute_radiative_fluxes!
#   * src/kernel/operators/Axb_rad_mpi.jl :: LinearOperators wrappers
#   * src/kernel/operators/build_rad_*.jl
# All gated by `inputs[:lRT_problem]` (drivers.jl).
const _RT_LOADED = Ref(false)
"""
    Jexpresso._ensure_rt_loaded!()

Lazily bring RRTMGP + its submodules, ClimaComms and LinearOperators
into Jexpresso's namespace. No-op after the first call.
"""
function _ensure_rt_loaded!()
    _RT_LOADED[] && return nothing
    @eval Jexpresso begin
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
        using RRTMGP.ArtifactPaths
        import RRTMGP.Parameters.RRTMGPParameters
        import RRTMGP: get_artifact_path
        using LinearOperators
    end
    # NCDatasets is also touched from this code path (Dataset(...) in
    # compute_radiative_fluxes!), so pull it in via its dedicated loader.
    _ensure_netcdf_loaded!()
    _RT_LOADED[] = true
    return nothing
end

# ─── NetCDF I/O (via NCDatasets) ──────────────────────────────────────
# Touched only by:
#   * src/io/write_output.jl     :: write_NetCDF (outformat=NETCDF())
#   * src/io/read_dp_scream.jl   :: read_atmospheric_data (RT data)
#   * src/kernel/mesh/phys_grid.jl :: compute_radiative_fluxes! (RT)
const _NETCDF_LOADED = Ref(false)
"""
    Jexpresso._ensure_netcdf_loaded!()

Lazily `using NCDatasets`. Safe to call from any NetCDF entry point
(write_NetCDF, read_atmospheric_data, compute_radiative_fluxes!).
"""
function _ensure_netcdf_loaded!()
    _NETCDF_LOADED[] && return nothing
    @eval Jexpresso using NCDatasets
    _NETCDF_LOADED[] = true
    return nothing
end

# ─── HDF5 I/O ─────────────────────────────────────────────────────────
# Touched only by:
#   * src/io/write_output.jl :: write_hdf5, read_hdf5
# Reached when `outformat::HDF5` dispatch fires (diagnostic write or
# restart). The local sentinel `HDF5()` (abstractTypes.jl) is unrelated
# to the package import — it's just a dispatch tag — so it stays even
# without `using HDF5`.
const _HDF5_LOADED = Ref(false)
"""
    Jexpresso._ensure_hdf5_loaded!()

Lazily `using HDF5`. Called from write_hdf5/read_hdf5 before the first
`h5open`/`h5read`.
"""
function _ensure_hdf5_loaded!()
    _HDF5_LOADED[] && return nothing
    @eval Jexpresso using HDF5
    _HDF5_LOADED[] = true
    return nothing
end

# ─── AMR / p4est forest ──────────────────────────────────────────────
# GridapP4est + P4est_wrapper are touched only by the AMR-style paths:
#   * mesh.jl  :: UniformlyRefinedForestOfOctreesDiscreteModel,
#                 OctreeDistributedDiscreteModel, GridapP4est.pXest_copy
#   * mesh.jl  :: load_p4est_checkpoint_model (uses P4est_wrapper.p8est_*)
#   * sem_setup.jl :: AMR restart branch
#   * TimeIntegrators.jl :: write_p4est_checkpoint
# All gated by `inputs[:ladapt]`, `inputs[:linitial_refine]`,
# `inputs[:lamr]`, or `inputs[:lrestart_amr]`.
const _AMR_LOADED = Ref(false)
"""
    Jexpresso._ensure_amr_loaded!()

Lazily bring GridapP4est and P4est_wrapper into Jexpresso's namespace.
Called from every entry point that touches the p4est forest. No-op
after the first call.
"""
function _ensure_amr_loaded!()
    _AMR_LOADED[] && return nothing
    @eval Jexpresso begin
        using GridapP4est
        using P4est_wrapper
    end
    _AMR_LOADED[] = true
    return nothing
end

# Run the test
# test_create_2d_projection_matrices_numa2d()


# ──────────────────────────────────────────────────────────────────────
# REPL entry point.
#
# Lets the user run a case after a single `using Jexpresso`:
#
#     julia> using Jexpresso                       # ← once, hits precompile cache
#     julia> Jexpresso.run_case("CompEuler", "3d") # ← every subsequent run, fast
#
# This replaces the older `include("./src/Jexpresso.jl")` pattern which
# evaluated the whole `module Jexpresso ... end` body again and
# triggered `WARNING: replacing module Jexpresso` on every iteration —
# along with throwing away the precompile cache for the orchestration
# layer.
#
# What's actually cheap on the second call: every dependency is already
# loaded, every Jexpresso source file is already compiled, and the
# typed-barrier RHS specialisations from the previous call stay hot in
# the method cache. The only re-evaluated code is the case-specific
# `user_inputs.jl` / `user_flux.jl` / … plus the per-call orchestration
# in `run.jl` (no module redefinition, so no warning).
#
# CI_MODE=true points the case loader at test/CI-runs/<eqs>/<eqs_case>
# instead of problems/<eqs>/<eqs_case>; matches the third positional
# arg of the historical command-line form.
# ──────────────────────────────────────────────────────────────────────
# Record of which case's user_*.jl files are currently loaded in this
# session (case dir + per-file mtimes). run.jl consults this to skip
# re-`include`ing them when the SAME case is re-run unchanged. Re-including
# redefines user_flux!/user_source!/… which invalidates the compiled rhs!,
# forcing the RHS + SciML integrator to recompile inside the warm-up on
# every run_case call (the ~30 s "Precompile warm-up" freeze on a re-run).
# Skipping the redundant include keeps an unchanged re-run launch-cost-only;
# switching cases or editing any user_*.jl bumps the check so changes still
# take effect.
const _LOADED_CASE_DIR  = Ref{String}("")
const _CASE_FILE_MTIMES = Dict{String,Float64}()

"""
    Jexpresso.run_case(eqs, eqs_case; CI_MODE=false)

Run a single Jexpresso case from the REPL. `eqs` and `eqs_case` are
the directory names under `problems/<eqs>/<eqs_case>/` (e.g.
`run_case("CompEuler", "3d")`). Returns `nothing`.

Prefer this to `include("./src/Jexpresso.jl")` — the include path
re-defines the module on every call (the `WARNING: replacing module
Jexpresso` you see) and discards the precompile cache for the
orchestration layer. With `run_case`, `using Jexpresso` happens
exactly once per session and the second case onward is launch-cost-only.
"""
function run_case(eqs::AbstractString, eqs_case::AbstractString;
                  CI_MODE::Bool=false)
    push!(empty!(ARGS), String(eqs), String(eqs_case), string(CI_MODE))
    # run.jl is a script — no `module …` declaration — so re-`include`ing
    # it from inside the already-loaded Jexpresso module just re-runs
    # the orchestration top-to-bottom in this module's scope. The
    # function/case-file redefinitions it performs are silent (only
    # module redefinitions print the WARNING).
    Base.include(@__MODULE__, joinpath(@__DIR__, "run.jl"))
    return nothing
end


# PERF: precompile workload — runs `test/CI-runs/CompEuler/sod1d`
# one driver pass through during `Pkg.precompile()` so PrecompileTools
# records every method that fires inside the integrator + RHS chain
# and bakes it into the package cache. Without this the very first
# `Jexpresso.run_case(...)` in a fresh REPL pays the full JIT cost
# for the integrator step + callback-specialised RHS — tens of
# seconds. With it, the first call is launch-cost-only.
#
# Why sod1d: it's the smallest case in the tree that exercises the
# full pipeline (1D mesh built by Jexpresso so no GMSH dependency,
# `nelx=100`, `nop=4`, `lvisc=true`, `tend=0.2` with `Δt=1e-4` →
# 2000 timesteps). Trades a few seconds of precompile time for
# baked-in code for every downstream case that uses the same RHS
# kernels.
#
# Files exercised: test/CI-runs/CompEuler/sod1d/user_inputs.jl plus
# the five sibling user_*.jl + initialize.jl files (all six are
# included unconditionally by src/run.jl, so all six must exist for
# the workload to succeed — this is what tripped us up before sod1d
# was added back).
@setup_workload begin
    # The workload below calls MPI.Init() (via run.jl). If precompilation
    # is triggered from a process that was itself launched by mpiexec —
    # e.g. a cold-cache `mpiexec -n 3 julia -e 'using Jexpresso; ...'` —
    # the precompile worker inherits the launcher's PMI/hydra environment
    # variables (PMI_FD, PMI_RANK, ...). MPICH then attempts a full PMI
    # handshake against a socket it doesn't own and aborts with
    # `PMI_Get_appnum returned -1`, failing precompilation on every rank.
    # Scrub those variables so MPI.Init falls back to singleton init
    # inside the precompile sandbox, exactly as it does when
    # `Pkg.precompile()` runs from a plain serial REPL. This only runs
    # during precompilation (@setup_workload is a no-op otherwise), so
    # the actual mpiexec-launched compute processes are unaffected.
    for k in collect(keys(ENV))
        if startswith(k, "PMI_")    || startswith(k, "PMIX_")  ||
           startswith(k, "HYDRA_")  || startswith(k, "HYDI_")  ||
           startswith(k, "MPIEXEC_")|| startswith(k, "OMPI_")  ||
           startswith(k, "MPI_LOCALRANKID") || startswith(k, "MPI_LOCALNRANKS")
            delete!(ENV, k)
        end
    end

    # ── MPI-driven driver workload ──────────────────────────────────────
    # Running the sod1d driver below bakes the hot-path JIT (RHS, the SciML
    # integrator, the callback-specialized warm-up) into the precompile
    # cache, so the first `run_case` in a fresh process is launch-cost-only
    # instead of paying ~tens of seconds of JIT. The workload is kept *lean*
    # — drivers.jl caps the run to 3 timesteps while `jl_generating_output`
    # is set — so precompilation stays fast; a full 2000-step sod1d pass is
    # NOT run here.
    #
    # The driver calls MPI.Init(), which on an InfiniBand cluster brings up
    # the libfabric (OFI) fabric. On a login node that step is hostile:
    #
    #   * verbs/mlx5 (the IB default) tries to allocate an RDMA queue pair
    #     and aborts — verbs are disabled / locked memory too low:
    #         Failed to modify UD QP to INIT on mlx5_0: Operation not permitted
    #         create_vni_context: Cannot allocate memory
    #   * tcp enumerates every NIC (IPoIB, bonded, …) and does reverse-DNS
    #     during MPI.Init — this can stall for many minutes.
    #
    # So for the single-process precompile worker we steer libfabric onto
    # the shared-memory `shm` provider (unless the user pinned FI_PROVIDER),
    # which needs no network and is the singleton-init happy path, and
    # restore the previous state afterwards so nothing leaks past
    # precompilation.
    #
    # The workload is OFF by default: on an InfiniBand login node MPI.Init
    # cannot reliably bring up a fabric (verbs aborts, tcp/shm can stall for
    # many minutes), and we will not let precompilation hang on that. The
    # package still precompiles fully without it; you only lose the warm-up,
    # so the first solve of a given problem shape pays its JIT at runtime.
    #
    # Opt in — only on a compute node where the fabric is healthy, or any
    # machine where MPI.Init is cheap (e.g. a laptop) — with:
    #
    #     JEXPRESSO_PRECOMPILE_WORKLOAD=1   (also: true / yes / on)
    #
    # When opted in the run is lean (drivers.jl caps it to 3 steps while
    # generating precompile output) and uses FI_PROVIDER=shm for the
    # single-process worker.
    _run_workload = lowercase(get(ENV, "JEXPRESSO_PRECOMPILE_WORKLOAD", "0")) in
                    ("1", "true", "yes", "on")
    _fi_provider_was_set = haskey(ENV, "FI_PROVIDER")
    _fi_provider_prev    = get(ENV, "FI_PROVIDER", "")
    if _run_workload && !_fi_provider_was_set
        ENV["FI_PROVIDER"] = "shm"
    end

    @compile_workload begin
        if _run_workload
            push!(empty!(ARGS), "CompEuler", "sod1d", "true")
            include(joinpath(@__DIR__, "run.jl"))   # lean driver pass (3 steps)
        end
    end

    # Undo the temporary FI_PROVIDER override (no-op if we never set it, or
    # if the user had pinned it — we left theirs untouched above).
    if _run_workload && !_fi_provider_was_set
        delete!(ENV, "FI_PROVIDER")
    elseif _run_workload
        ENV["FI_PROVIDER"] = _fi_provider_prev
    end
end
end
