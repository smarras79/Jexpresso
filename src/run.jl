#--------------------------------------------------------
# run.jl (Modified for Late Handshake)
#--------------------------------------------------------
function jexpresso_main()

    # Initialize MPI and split communicator:
    world, local_comm, wsize, wrank, lsize, lrank, is_coupled = je_init_mpi_and_split_comm()
    
    # Parse args and include files
    # (Note: These are fast, but Julia will JIT compile them)
    mod_io_parse_args()

    include(driver_file)
    include(user_input_file)
    include(user_flux_file)
    include(user_source_file)
    include(user_bc_file)
    include(user_initialize_file)
    include(user_primitives_file)

    # Load User Inputs
    Base.invokelatest(mod_inputs_print_welcome, lrank)

    global inputs = Dict{}()
    inputs = Base.@invokelatest user_inputs()
    Base.invokelatest(mod_inputs_user_inputs!, inputs, lrank)

    OUTPUT_DIR = mod_io_mkoutdir!(inputs)

    # Launch Driver with both Local and World communicators
    # We use Base.invokelatest to avoid world-age issues with included files
    with_mpi(; comm=local_comm) do distribute
        Base.@invokelatest driver(lsize, distribute, inputs, OUTPUT_DIR, TFloat, world)

        # In coupled (MPMD) mode Fortran calls MPI_Barrier(MPI_COMM_WORLD) after
        # its time loop to coordinate clean shutdown.  Julia must call the matching
        # barrier from inside the with_mpi block, where all Julia ranks are still
        # executing in parallel and `world` is in scope.  Without this, Fortran
        # blocks forever at its barrier while Julia races ahead to MPI.Finalize.
        if is_coupled
            println("[Jexpresso rank $lrank] time loop done — entering world barrier"); flush(stdout)
            MPI.Barrier(world)
            println("[Jexpresso rank $lrank] world barrier passed"); flush(stdout)
        end
    end
end

if !haskey(ENV, "JEXPRESSO_COUPLING_MODE")
    jexpresso_main()
else
    jexpresso_main()
    MPI.Finalize()
end
