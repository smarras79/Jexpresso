using ArgParse

#--------------------------------------------------------
# The problem name is a command line argument:
#
# 1. Launch Julia:
# >> julia --project=.
#
# 2. Push equations name to ARGS
#    You need this only when you run a new equations
#
#    julia > push!(empty!(ARGS), BENCHMARK::String, CASE_NAME::String);
#    julia > include(./src/Jexpresso.jl)
#
#    BENCHMARK is the name of the user's directory that contains a user-defined CASE_NAME
#    CASE_NAME is the name of the user's subdirectory $JEXPRESSO/problems/BENCHMARK/CASE_NAME (e.g. theta)
#
# Ex. To run the rising thermal bubble benchmark: $JEXPRESSO/problems/CompEuler/theta
# 
#  julia > push!(empty!(ARGS), "CompEuler", "theta");
#  julia > include(./src/Jexpresso.jl)
#
# To create a new case:
#
# mkdir $JEXPRESSO/problems/USER_DEFINED_DIR/
# mkdir $JEXPRESSO/problems/USER_DEFINED_DIR/USER_DEFINED_CASE_NAME
#
# ex.:
# mkdir $JEXPRESSO/problems/acoustics
# mkdir $JEXPRESSO/problems/acoustics/acoustics2d
#
#  julia > push!(empty!(ARGS), "acoustics", "acoustics2d");
#  julia > include(./src/Jexpresso.jl)
#
#--------------------------------------------------------
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "eqs"
        help = "Directoy that contains some user-defined cases"
        default = "CompEuler"
        required = false
        
        "eqs_case"
        help = "case name in equations directory"
        default = "wave1d"
        required = false
        
        "CI_MODE"
        help = "CI_MODE: true or false"
        default = "false"
        required = false
    end

    return parse_args(s)
end


#--------------------------------------------------------
# MPI initialisation.
#
# By default this is the same one-shot init used historically: a single
# COMM_WORLD shared with no other code. When the env var JEXPRESSO_COUPLED
# is set (1 / true / yes), we instead use je_init_mpi_and_split_comm()
# which splits COMM_WORLD by APPID and lets us detect a coupled run with
# another MPI code such as Alya. In standalone runs the two paths are
# functionally equivalent; the env-gated path is opt-in to guarantee no
# behavioural drift when JEXPRESSO_COUPLED is unset.
#--------------------------------------------------------
const _JEXPRESSO_COUPLED_ENV = lowercase(get(ENV, "JEXPRESSO_COUPLED", ""))
const JEXPRESSO_COUPLING_ENABLED = _JEXPRESSO_COUPLED_ENV in ("1", "true", "yes", "on")

if JEXPRESSO_COUPLING_ENABLED
    _world, _local_comm, _wsize, _wrank, _lsize, _lrank, _is_coupled =
        je_init_mpi_and_split_comm()
    comm   = _local_comm
    rank   = _lrank
    nparts = _lsize
else
    MPI.Init()
    _world      = MPI.COMM_WORLD
    _local_comm = MPI.COMM_WORLD
    _is_coupled = false
    comm   = MPI.COMM_WORLD
    rank   = MPI.Comm_rank(comm)
    nparts = MPI.Comm_size(comm)
end

#--------------------------------------------------------
# Parse command line args:
#--------------------------------------------------------
parsed_args                = parse_commandline()
parsed_equations           = string(parsed_args["eqs"])
parsed_equations_case_name = string(parsed_args["eqs_case"])
parsed_CI_mode             = string(parsed_args["CI_MODE"])
driver_file                = string(dirname(@__DIR__()), "/problems/drivers.jl")

# Check if running under CI environment and set directory accordingly
if parsed_CI_mode == "true"
    case_name_dir = string(dirname(@__DIR__()), "/test/CI-runs", "/", parsed_equations, "/", parsed_equations_case_name)
else
    case_name_dir = string(dirname(@__DIR__()), "/problems", "/", parsed_equations, "/", parsed_equations_case_name)
end

user_input_file      = string(case_name_dir, "/user_inputs.jl")
user_flux_file       = string(case_name_dir, "/user_flux.jl")
user_source_file     = string(case_name_dir, "/user_source.jl")
user_bc_file         = string(case_name_dir, "/user_bc.jl")
user_initialize_file = string(case_name_dir, "/initialize.jl")
user_primitives_file = string(case_name_dir, "/user_primitives.jl")
# OPTIONAL per-case file: a closed-form/reference solution the 1D plotter
# overlays on its output (see problems/CompEuler/sod1d/user_analytic.jl).
# Most cases have none, so this one is included only when present.
user_analytic_file   = string(case_name_dir, "/user_analytic.jl")

# PERF: only (re-)include the driver + the case's user_*.jl files when
# something actually changed since the last run in this session. Re-running
# the SAME case unchanged must NOT re-`include` these, because that
# redefines user_flux!/user_source!/user_bc!/… and redefining a method
# Julia has already specialized into the compiled `rhs!` invalidates it —
# so the next `solve(...)` (the warm-up step in time_loop!) recompiles the
# entire RHS + integrator from scratch, which is the ~30 s freeze at
# "# Precompile warm-up (1 step solve)" on a second identical run_case.
#
# Reload when the case directory changed (different case → its functions
# really must be redefined) or when any of these files was edited (mtime
# bump → pick up the user's change). Otherwise reuse the already-compiled
# code: an unchanged re-run is then launch-cost-only.
_case_load_files = [driver_file, user_input_file, user_flux_file,
                    user_source_file, user_bc_file, user_initialize_file,
                    user_primitives_file]
isfile(user_analytic_file) && push!(_case_load_files, user_analytic_file)
_need_case_reload = (_LOADED_CASE_DIR[] != case_name_dir) ||
    any(f -> get(_CASE_FILE_MTIMES, f, -1.0) != mtime(f), _case_load_files)
if _need_case_reload
    for _f in _case_load_files
        include(_f)
        _CASE_FILE_MTIMES[_f] = mtime(_f)
    end
    _LOADED_CASE_DIR[] = case_name_dir
end
#--------------------------------------------------------
# Read User Inputs:
#--------------------------------------------------------
mod_inputs_print_welcome(rank)
inputs = Dict{}()

inputs = user_inputs()
# Make the case directory available to the mesh/SEM cache helpers so cache
# files live next to user_inputs.jl (per-case), not next to the shared
# *.msh file.  This is what lets two cases that happen to point at the same
# gmsh file keep separate caches and what makes "running a new case"
# automatically miss the cache without any user-visible flag.
inputs[:_case_dir]            = case_name_dir
# Whether THIS case ships a reference solution. Checked by the 1D plotter
# instead of a bare isdefined(): user_analytic_solution is a method on the
# Jexpresso module, so once sod1d has been run in a session the definition
# survives a switch to another 1D case, and a bare isdefined() would then
# happily overlay Sod's exact solution on an unrelated problem.
inputs[:_has_analytic]        = isfile(user_analytic_file)
inputs[:_parsed_equations]    = parsed_equations
inputs[:_parsed_case_name]    = parsed_equations_case_name
inputs[:_user_input_file]     = user_input_file

#--------------------------------------------------------
# CI output conventions.
#
# A CI run is only useful if the comparison can find its result: the CI
# machinery reads the .h5 files in test/CI-runs/<eqs>/<case>/output/ and
# checks them against test/CI-ref/<eqs>/<case>/output/. That requires
# three things of the deck, and a deck copied straight out of problems/
# typically satisfies none of them — it asks for VTK, in an output
# directory of its own, in a timestamped subdirectory:
#
#   :outformat         => "hdf5"   the comparison reads HDF5
#   :output_dir        => "none"   put output next to the case inputs
#   :loverwrite_output => true     "output", not "output-05Aug2025-181233"
#
# So CI mode forces them, and says so, rather than leaving the user to
# discover an empty output directory an hour later. Everything else in
# the deck (tend, diagnostics, mesh, …) is the case author's business.
#
# Set JEXPRESSO_CI_OUTPUT=0 to keep the deck's own settings — e.g. to get
# VTK out of a CI deck for visualisation.
#--------------------------------------------------------
if parsed_CI_mode == "true" &&
   !(lowercase(get(ENV, "JEXPRESSO_CI_OUTPUT", "1")) in ("0", "false", "no", "off"))

    for (ci_key, ci_value) in (:outformat         => "hdf5",
                               :output_dir        => "none",
                               :loverwrite_output => true)
        previous = get(inputs, ci_key, nothing)
        unchanged = previous isa AbstractString && ci_value isa AbstractString ?
                    lowercase(previous) == ci_value : previous == ci_value
        inputs[ci_key] = ci_value
        if !unchanged && rank == 0
            println(" # CI_MODE: forcing :", ci_key, " => ", repr(ci_value),
                    previous === nothing ? " (deck left it unset)" :
                                           string(" (deck asked for ", repr(previous), ")"))
        end
    end

    #----------------------------------------------------------
    # Output cadence: one write, at the end of the run.
    #
    # A reference solution is the final state; intermediate dumps cost
    # wall time and, since write_hdf5 names its files var_<i>_<rank>.h5
    # with no time index, every one of them is overwritten by the next
    # anyway. A deck asking for `:diagnostics_at_times => (0:10:1000)`
    # would therefore write 101 times to produce the one snapshot that
    # survives — and add 101 integrator tstops on the way.
    #
    # Setting :diagnostics_at_times makes mod_inputs_user_inputs! zero
    # :ndiagnostics_outputs (its else branch), which is why the deck's
    # value is dropped here rather than overwritten.
    #----------------------------------------------------------
    if haskey(inputs, :tend) && inputs[:tend] isa Number
        ci_final_time = [Float64(inputs[:tend])]
        if get(inputs, :diagnostics_at_times, nothing) != ci_final_time && rank == 0
            println(" # CI_MODE: forcing a single output at t=", inputs[:tend],
                    " (deck asked for ",
                    repr(get(inputs, :diagnostics_at_times,
                             get(inputs, :ndiagnostics_outputs, "nothing"))), ")")
        end
        inputs[:diagnostics_at_times] = ci_final_time
        delete!(inputs, :ndiagnostics_outputs)
    end
    inputs[:lwrite_initial] = false

    # NOTE: CI mode deliberately does NOT switch on :lstep_heartbeat. With a
    # single write at :tend nothing is printed during the solve, so a slow run
    # looks like a hung one — but per-step output is noise in a healthy run.
    # Set JEXPRESSO_STEP_HEARTBEAT=1 when you need to see progress or measure
    # the step rate (first 5 steps, then every 100th; see TimeIntegrators.jl).
end

if rank == 0
    print(" # mod_inputs_user_inputs! (filling defaults) ......... ")
    flush(stdout)
end
_t_defaults = time_ns()
mod_inputs_user_inputs!(inputs, rank)
if rank == 0
    @printf("DONE (%.2f s)\n", (time_ns() - _t_defaults) / 1e9)
    flush(stdout)
end

#--------------------------------------------------------
# Create output directory if it doesn't exist:
#--------------------------------------------------------
user_defined_output_dir = inputs[:output_dir]

if inputs[:loverwrite_output]
    outstring = string("output")
else        
    outstring = rank == 0 ? string("output-",  Dates.format(now(), "dduyyyy-HHMMSS")) : ""
    outstring = MPI.bcast(outstring, 0, comm)
end
if user_defined_output_dir == "none"
    OUTPUT_DIR = joinpath(case_name_dir, outstring)
    inputs[:output_dir] = OUTPUT_DIR
else
    OUTPUT_DIR = joinpath(user_defined_output_dir, parsed_equations, parsed_equations_case_name, outstring)
    inputs[:output_dir] = OUTPUT_DIR
end
if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

#--------------------------------------------------------
# Create restart output/inupt directory if it doesn't exist:
#--------------------------------------------------------
if (!haskey(inputs, :restart_output_file_path))
    inputs[:restart_output_file_path] = joinpath(OUTPUT_DIR,string("restart"))
end

if (haskey(inputs, :lrestart))
    if(inputs[:lrestart] == true && !haskey(inputs, :restart_input_file_path))
        inputs[:restart_input_file_path] = inputs[:restart_output_file_path]
    end
else
    inputs[:lrestart] = false
end

#--------------------------------------------------------
# Save a copy of user_inputs.jl for the case being run
#--------------------------------------------------------
if rank == 0
    cp(user_input_file, joinpath(OUTPUT_DIR, basename(user_input_file)); force = true)
end

#--------------------------------------------------------
# Typed cache of :energy_equation for hot-path callers (user_flux!,
# user_primitives!, user_uout!). Reading the non-const module-global
# `inputs` from inside per-quadrature-point loops costs ~38 ns per
# access and inhibits inlining; a Ref{Bool} const reads in ~2 ns and
# is type-stable.
#--------------------------------------------------------
if !@isdefined(ENERGY_EQUATION_THETA)
    const ENERGY_EQUATION_THETA = Ref{Bool}(false)
end
ENERGY_EQUATION_THETA[] = (inputs[:energy_equation] == "theta")

#--------------------------------------------------------
# Convert inputs Dict → NamedTuple if the user opted in.
#
# Carrying inputs as a NamedTuple removes the per-call dictionary lookups
# from hot paths (RHS, time loop, BCs) and lets the compiler specialise on
# the concrete field types.  All historical `inputs[:key]` reads continue
# to work; only WRITES (which would mutate the container) are forbidden,
# so any `inputs[:key] = value` left in downstream code must be expressed
# as a rebind `inputs = (; inputs..., key = value)` (see drivers.jl).
#
# All inputs mutations above this point operate on the Dict, so they are
# safe; the conversion happens here, just before the coupling handshake
# and the with_mpi block.
#
# val_lsaturation is a Val-wrapped boolean made available so RHS kernels
# can dispatch on it without paying a runtime branch + dictionary lookup
# on every step (matches the ab/hacky pattern).
#--------------------------------------------------------
if get(inputs, :use_named_tuples, false) == true
    inputs = NamedTuple(inputs)
end

val_lsaturation = Val(get(inputs, :lsaturation, false))
inputs = inputs isa NamedTuple ?
    (; inputs..., comm = MPI.COMM_WORLD, val_lsaturation = val_lsaturation) :
    inputs

#--------------------------------------------------------
# Coupling handshake (must happen OUTSIDE the with_mpi block so the
# coupling MPI calls execute before with_mpi's closure JIT begins; see
# couplingStructs.jl for the rationale). When JEXPRESSO_COUPLED is unset
# this branch is skipped entirely and the standalone code path is
# byte-for-byte the historical one.
#--------------------------------------------------------
if JEXPRESSO_COUPLING_ENABLED
    _is_coupled = je_perform_coupling_handshake(_world, nparts)
    if _is_coupled
        je_receive_alya_data(_world, nparts)
        je_prefetch_caches!(inputs, nparts, _local_comm, _world)
        je_early_coupling_sync!(_local_comm, _world)
    end
end

#--------------------------------------------------------
# use Metal (for apple) or CUDA (non apple) if we are on GPU
#--------------------------------------------------------
# UX: this is the boundary where PartitionedArrays' with_mpi block
# JIT-compiles its closure on first call (cold mpirun). Tens of
# seconds is expected here on cold start; subsequent runs in the
# SAME Julia session are instant. Pre-baking this JIT into the
# package precompile cache is what PrecompileTools.@compile_workload
# would do — see the note at the end of this file's output.
if rank == 0
    print(" # Entering with_mpi block (PartitionedArrays JIT on first call) ......... ")
    flush(stdout)
end
_t_withmpi = time_ns()
if JEXPRESSO_COUPLING_ENABLED
    with_mpi(; comm = _local_comm) do distribute
        if rank == 0
            @printf("DONE (%.2f s) → calling driver()\n", (time_ns() - _t_withmpi) / 1e9)
            flush(stdout)
        end
        driver(nparts,
               distribute,
               inputs,
               OUTPUT_DIR,
               TFloat;
               world      = _world,
               is_coupled = _is_coupled)

        # In coupled (MPMD) mode Fortran waits on a final
        # MPI_Barrier(MPI_COMM_WORLD) for clean shutdown. The matching
        # Julia barrier must be issued from inside the with_mpi block
        # while `_world` is still in scope.
        if _is_coupled
            MPI.Barrier(_world)
        end
    end
else
    with_mpi() do distribute
        if rank == 0
            @printf("DONE (%.2f s) → calling driver()\n", (time_ns() - _t_withmpi) / 1e9)
            flush(stdout)
        end
        driver(nparts,
               distribute,
               inputs,
               OUTPUT_DIR,
               TFloat)
    end
end
