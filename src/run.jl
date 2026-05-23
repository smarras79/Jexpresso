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

include(driver_file)

include(user_input_file)
include(user_flux_file)
include(user_source_file)
include(user_bc_file)
include(user_initialize_file)
include(user_primitives_file)
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
inputs[:_parsed_equations]    = parsed_equations
inputs[:_parsed_case_name]    = parsed_equations_case_name
inputs[:_user_input_file]     = user_input_file
mod_inputs_user_inputs!(inputs, rank)
println("[init][rank=$rank] mod_inputs_user_inputs! returned"); flush(stdout)

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
println("[init][rank=$rank] OUTPUT_DIR ready: $OUTPUT_DIR"); flush(stdout)
MPI.Barrier(comm)
println("[init][rank=$rank] passed barrier after OUTPUT_DIR"); flush(stdout)

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
println("[init][rank=$rank] inputs container finalised (Dict/NamedTuple)"); flush(stdout)
MPI.Barrier(comm)
println("[init][rank=$rank] about to enter with_mpi block ..."); flush(stdout)

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
if JEXPRESSO_COUPLING_ENABLED
    with_mpi(; comm = _local_comm) do distribute
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
    # Probe MPI state before calling PartitionedArrays.with_mpi(), to isolate
    # whether the hang originates in our own MPI handling or inside
    # PartitionedArrays' closure setup (which typically does Comm_dup +
    # threadlevel checks).
    println("[init][rank=$rank] MPI.Initialized() = $(MPI.Initialized())"); flush(stdout)
    if MPI.Initialized()
        try
            _thread_level = MPI.Query_thread()
            println("[init][rank=$rank] MPI.Query_thread() = $_thread_level"); flush(stdout)
        catch e
            println("[init][rank=$rank] MPI.Query_thread() threw: $e"); flush(stdout)
        end
    end
    println("[init][rank=$rank] probing MPI.Comm_dup(comm) directly ..."); flush(stdout)
    _probe_dup = MPI.Comm_dup(comm)
    println("[init][rank=$rank] MPI.Comm_dup(comm) returned ok"); flush(stdout)
    # Skip explicit Comm_free — MPI.jl frees duplicated comms at finalize anyway,
    # and the symbol name differs across MPI.jl versions.

    println("[init][rank=$rank] >>> calling with_mpi(; comm=MPI.COMM_WORLD) now"); flush(stdout)
    with_mpi(; comm = MPI.COMM_WORLD) do distribute
        local _r = MPI.Comm_rank(distribute.comm)
        local _s = MPI.Comm_size(distribute.comm)
        println("[init][rank=$_r] inside with_mpi: distribute.comm acquired (size=$_s)"); flush(stdout)
        MPI.Barrier(distribute.comm)
        println("[init][rank=$_r] inside with_mpi: passed barrier, calling driver()"); flush(stdout)
        driver(nparts,
               distribute,
               inputs,
               OUTPUT_DIR,
               TFloat)
    end
end
