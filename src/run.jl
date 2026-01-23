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
# Main execution function
# For coupling mode: call after setting communicator with set_mpi_comm()
# For standalone mode: called automatically at module load
#
# NOTE: This must be a function (not top-level code) because:
# 1. Setup code needs to use the custom communicator (for MPI.bcast, etc)
# 2. The communicator can only be set AFTER the module is loaded
# 3. Therefore setup must be delayed until after module load
# 4. This requires wrapping setup in a function
# 5. Many functions expect certain variables as module globals (legacy design)
# 6. Hence the 'global' declarations - these make local variables into module globals
#--------------------------------------------------------
function jexpresso_main()

    # Initialize MPI if not already done
    if !MPI.Initialized()
        MPI.Init()
    end

    # Get communicator - custom one if set (coupling mode), otherwise COMM_WORLD
    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    nparts = MPI.Comm_size(comm)

    #--------------------------------------------------------
    # Parse command line args:
    #--------------------------------------------------------
    parsed_args = parse_commandline()

    # These must be global because other functions access them by name (legacy design)
    global parsed_equations           = string(parsed_args["eqs"])
    global parsed_equations_case_name = string(parsed_args["eqs_case"])
    global parsed_CI_mode             = string(parsed_args["CI_MODE"])
    global driver_file                = string(dirname(@__DIR__()), "/problems/drivers.jl")

    # Check if running under CI environment and set directory accordingly
    if parsed_CI_mode == "true"
        case_name_dir = string(dirname(@__DIR__()), "/test/CI-runs", "/", parsed_equations, "/", parsed_equations_case_name)
    else
        case_name_dir = string(dirname(@__DIR__()), "/problems", "/", parsed_equations, "/", parsed_equations_case_name)
    end

    # These must be global because other functions access them by name (legacy design)
    global user_input_file      = string(case_name_dir, "/user_inputs.jl")
    global user_flux_file       = string(case_name_dir, "/user_flux.jl")
    global user_source_file     = string(case_name_dir, "/user_source.jl")
    global user_bc_file         = string(case_name_dir, "/user_bc.jl")
    global user_initialize_file = string(case_name_dir, "/initialize.jl")
    global user_primitives_file = string(case_name_dir, "/user_primitives.jl")

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
    # Use Base.invokelatest to handle world age issue when dynamically loading functions
    Base.invokelatest(mod_inputs_print_welcome, rank)

    # inputs must be global because many functions access it by name (legacy design)
    global inputs = Dict{}()
    inputs = Base.invokelatest(user_inputs)
    Base.invokelatest(mod_inputs_user_inputs!, inputs, rank)

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
    # Use Base.invokelatest for dynamically loaded driver function
    # IMPORTANT: Pass custom communicator to with_mpi for coupling mode
    with_mpi(; comm=comm) do distribute

        Base.invokelatest(driver,
                         nparts,
                         distribute,
                         inputs, # input parameters from src/user_input.jl
                         OUTPUT_DIR,
                         TFloat)

    end
end

#--------------------------------------------------------
# Auto-execute if this file is run directly (not in coupling mode)
# Skip auto-execution if JEXPRESSO_COUPLING_MODE env var is set
#--------------------------------------------------------
if !haskey(ENV, "JEXPRESSO_COUPLING_MODE")
    jexpresso_main()
end
