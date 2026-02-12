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
# Coupled run with Alya:
#  mpirun --tag-output -np 2 ./alya/Alya_enhanced.x : -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
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
    # Initialize MPI
    #-----------------------------------------------------------------
    if !MPI.Initialized()
        MPI.Init()
    end

    #-----------------------------------------------------------------
    # Coupling detection and communicator setup.
    #
    # If JEXPRESSO_MPI_COMM was already set externally
    # (e.g. by Jexpresso-mini-coupled.jl), skip this block.
    #
    # Otherwise, split COMM_WORLD by APPID.  This is a collective
    # operation: ALL ranks (Alya + Julia) must call Comm_split.
    #   - Standalone mode: all ranks share the same APPID => lsize == wsize
    #   - Coupled mode:    Alya has a different APPID     => lsize < wsize
    #-----------------------------------------------------------------
    if JEXPRESSO_MPI_COMM[] === nothing
        world = MPI.COMM_WORLD
        wrank = MPI.Comm_rank(world)
        wsize = MPI.Comm_size(world)

        appid = try parse(Int, get(ENV, "APPID", "2")) catch; 2 end

        local_comm = MPI.Comm_split(world, appid, wrank)
        lrank = MPI.Comm_rank(local_comm)
        lsize = MPI.Comm_size(local_comm)

        if lsize < wsize
            #-------------------------------------------------------------
            # Coupled mode detected: we share COMM_WORLD with another code
            #-------------------------------------------------------------
            set_mpi_comm(local_comm)
            set_mpi_comm_world(world)

            if lrank == 0
                println("[Jexpresso] Coupled mode: world_size=$wsize, local_size=$lsize, appid=$appid")
                flush(stdout)
            end

            nranks_other = wsize - lsize

            # Exchange app identity (Gather names to world rank 0)
            local_chars = Vector{UInt8}(rpad("JEXPRESSO", 128, ' '))
            MPI.Gather!(local_chars, nothing, 0, world)

            #-------------------------------------------------------------
            # Receive remote grid metadata from Alya via Bcast on COMM_WORLD
            # Alya broadcasts from world rank 0.
            # Use Int32/Float32 to match Fortran MPI types.
            #-------------------------------------------------------------
            ndime_buf = Vector{Int32}(undef, 1)
            MPI.Bcast!(ndime_buf, 0, world)
            ndime = ndime_buf[1]

            rem_min = Vector{Float32}(undef, 3)
            rem_max = Vector{Float32}(undef, 3)
            rem_nx  = Vector{Int32}(undef, 3)
            for idime in 1:3
                MPI.Bcast!(@view(rem_min[idime:idime]), 0, world)
                MPI.Bcast!(@view(rem_max[idime:idime]), 0, world)
                MPI.Bcast!(@view(rem_nx[idime:idime]),  0, world)
            end

            alya2world_l = zeros(Int32, nranks_other)
            alya2world   = MPI.Allreduce(alya2world_l, MPI.SUM, world)

            a_l = zeros(Int32, wsize, wsize)
            a   = MPI.Allreduce(a_l, MPI.SUM, world)

            # Store coupling data for je_couplingSetup() to consume later
            set_coupling_data(Dict{Symbol,Any}(
                :ndime         => ndime,
                :rem_min       => rem_min,
                :rem_max       => rem_max,
                :rem_nx        => rem_nx,
                :alya2world    => alya2world,
                :couple_matrix => a,
            ))

            if lrank == 0
                println("[Jexpresso] Received from Alya: ndime=$ndime, " *
                        "rem_min=$rem_min, rem_max=$rem_max, rem_nx=$rem_nx")
                flush(stdout)
            end
        end
        # If lsize == wsize we are in standalone mode: nothing extra to do.
        # je_mpi_init() will use COMM_WORLD via get_mpi_comm().
    end

    #-----------------------------------------------------------------
    # Build communicators (uses local_comm if coupled, COMM_WORLD otherwise)
    #-----------------------------------------------------------------
    comm, rank, nparts = je_mpi_init()

    #-----------------------------------------------------------------
    # Parse command line args:
    #-----------------------------------------------------------------
    mod_io_parse_args()

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
        Base.@invokelatest driver(nparts, distribute, inputs, OUTPUT_DIR, TFloat)
    end

end

#--------------------------------------------------------
# Auto-execute if this file is run directly (not in coupling mode)
# Skip auto-execution if JEXPRESSO_COUPLING_MODE env var is set
#--------------------------------------------------------
if !haskey(ENV, "JEXPRESSO_COUPLING_MODE")
    jexpresso_main()
end
