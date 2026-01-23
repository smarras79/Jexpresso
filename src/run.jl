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

        "--coupling"
        help = "Enable MPI coupling mode (intercommunicator approach)"
        action = :store_true

        "--code-id"
        help = "Code ID for coupling (1-indexed integer)"
        arg_type = Int
        default = 1

        "--n-codes"
        help = "Total number of coupled codes"
        arg_type = Int
        default = 1

        "--code-name"
        help = "Descriptive name for this code"
        arg_type = String
        default = "Jexpresso"

        "--gather-coupling"
        help = "Enable gather-based MPI coupling mode (shared COMM_WORLD)"
        action = :store_true
    end

    return parse_args(s)
end


MPI.Init()
comm::MPI.Comm = MPI.COMM_WORLD
rank::Int = MPI.Comm_rank(comm)
nparts::Int = MPI.Comm_size(comm)

#--------------------------------------------------------
# Parse command line args:
#--------------------------------------------------------
parsed_args                = parse_commandline()
parsed_equations           = string(parsed_args["eqs"])
parsed_equations_case_name = string(parsed_args["eqs_case"])
parsed_CI_mode             = string(parsed_args["CI_MODE"])
parsed_coupling            = parsed_args["coupling"]
parsed_code_id             = parsed_args["code-id"]
parsed_n_codes             = parsed_args["n-codes"]
parsed_code_name           = parsed_args["code-name"]
parsed_gather_coupling     = parsed_args["gather-coupling"]
driver_file                = string(dirname(@__DIR__()), "/problems/drivers.jl")

#--------------------------------------------------------
# Initialize coupling if enabled:
#--------------------------------------------------------
coupling_ctx::Union{JexpressoCoupling.CouplingContext, Nothing} = nothing

# Check for mutually exclusive coupling modes
if parsed_coupling && parsed_gather_coupling
    if rank == 0
        @error "Cannot enable both --coupling and --gather-coupling simultaneously"
    end
    MPI.Abort(comm, 1)
end

if parsed_coupling
    if rank == 0
        @info "Initializing MPI coupling mode (intercommunicator): code_id=$parsed_code_id, n_codes=$parsed_n_codes, name=$parsed_code_name"
    end

    coupling_ctx = JexpressoCoupling.initialize_coupling(
        comm,
        parsed_code_id,
        parsed_n_codes;
        code_name=parsed_code_name
    )

    # Use local communicator for Jexpresso's internal operations
    comm = coupling_ctx.comm_local
    rank = coupling_ctx.local_rank
    nparts = coupling_ctx.local_size

    if rank == 0
        @info "Jexpresso running with $(nparts) local ranks (world_rank=$(coupling_ctx.world_rank))"
    end
elseif parsed_gather_coupling
    if rank == 0
        @info "Initializing gather-based MPI coupling mode (shared COMM_WORLD)"
    end

    # Participate in MPI_Comm_split with MPI_UNDEFINED
    # This creates a NULL communicator but allows participation in the collective operation
    MPI_UNDEFINED = Int32(-32766)
    par_comm_final = MPI.Comm_split(comm, MPI_UNDEFINED, rank)

    if rank == 0
        @info "Participated in MPI_Comm_split with MPI_UNDEFINED"
    end

    # Prepare application name for MPI_Gather (128 bytes, MPI_CHARACTER compatible)
    app_name_str = parsed_code_name
    send_buf = fill(Int8(' '), 128)
    nb = min(length(codeunits(app_name_str)), 128)
    @inbounds for i in 1:nb
        send_buf[i] = Int8(codeunits(app_name_str)[i])
    end

    # Participate in MPI_Gather
    # In MPMD mode, Jexpresso ranks are typically not global rank 0 (external code gets rank 0)
    # Jexpresso sends its application name to the global root
    MPI.Gather!(send_buf, nothing, 0, comm)

    if rank == 0
        @info "Gather-based coupling initialization complete (sent '$app_name_str' to global root). Continuing with normal Jexpresso execution on $(nparts) ranks."
    end

    # Continue with original communicator (no change to comm/rank/nparts)
    # Jexpresso operates on full MPI_COMM_WORLD
end

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
mod_inputs_user_inputs!(inputs, rank)

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
# use Metal (for apple) or CUDA (non apple) if we are on GPU
#--------------------------------------------------------
# Pass the correct communicator (comm_local when coupling, COMM_WORLD otherwise)
with_mpi(comm=comm) do distribute

    driver(nparts,
           distribute,
           inputs, # input parameters from src/user_input.jl
           OUTPUT_DIR,
           TFloat)

end

#--------------------------------------------------------
# Cleanup coupling if it was initialized:
#--------------------------------------------------------
if !isnothing(coupling_ctx)
    JexpressoCoupling.finalize_coupling(coupling_ctx)
end
