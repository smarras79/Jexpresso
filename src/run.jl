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
    
    #-----------------------------------------------------------------
    # Initialize MPI and build communicators
    #-----------------------------------------------------------------
    comm, rank, nparts = je_mpi_init()

    
    #-----------------------------------------------------------------
    # Parse command line args:
    #-----------------------------------------------------------------
    je_parse_args()

    
    #-----------------------------------------------------------------
    # Include user_*_file.jl
    #-----------------------------------------------------------------
    include(driver_file)
    include(user_input_file)
    include(user_flux_file)
    include(user_source_file)
    include(user_bc_file)
    include(user_initialize_file)
    include(user_primitives_file)

    
    #-----------------------------------------------------------------
    # Read User Inputs:
    #-----------------------------------------------------------------
    # Use Base.invokelatest to handle world age issue when dynamically loading functions
    Base.invokelatest(mod_inputs_print_welcome, rank)

    # inputs must be global because many functions access it by name (legacy design)
    global inputs = Dict{}()
    inputs        = Base.@invokelatest user_inputs()
    Base.invokelatest(mod_inputs_user_inputs!, inputs, rank)

    
    #-----------------------------------------------------------------
    # Create output directory if it doesn't exist:
    #-----------------------------------------------------------------
    OUTPUT_DIR = mod_io_mkoutdir!(inputs)

    
    #-----------------------------------------------------------------
    # IMPORTANT: Pass custom comm to with_mpi when coupling codes
    #-----------------------------------------------------------------
    with_mpi(; comm=comm) do distribute

        Base.@invokelatest driver(nparts,
                                  distribute,
                                  inputs,
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
