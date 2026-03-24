#--------------------------------------------------------
# run.jl
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

    # When building a PackageCompiler sysimage we only need 1 timestep to
    # exercise every code path and capture all JIT specialisations.
    # The full-length run is the job of the actual production execution.
    if haskey(ENV, "JEXPRESSO_WARMUP")
        inputs[:tend]                = inputs[:tinit] + inputs[:Δt]
        inputs[:ndiagnostics_outputs] = 0
        inputs[:loutput_pert]        = false
        inputs[:lwrite_initial]      = false
    end

    # COUPLING HANDSHAKE — must happen OUTSIDE the with_mpi block.
    #
    # Julia JIT-compiles a closure when it is first *called*, not when it is
    # defined.  If the handshake and je_receive_alya_data were placed inside
    # the with_mpi do-block, Julia would JIT-compile the entire closure body
    # (including je_receive_alya_data and all its transitive callees) before
    # executing a single line of it.  Alya blocks at its MPI_Bcast(ndime)
    # waiting for Julia to call MPI_Bcast — if that call is inside a closure
    # that is still being JIT-compiled, Alya waits for the full JIT (~15 s).
    #
    # Keeping the handshake and je_receive_alya_data here — as plain top-level
    # calls — means Julia JIT-compiles and *executes* them immediately,
    # completing Alya's Bcasts before the with_mpi closure JIT begins.
    #
    # je_receive_alya_data is idempotent: the second call from inside
    # setup_coupling_and_mesh is a no-op.
    is_coupled = je_perform_coupling_handshake(world, lsize)
    if is_coupled
        je_receive_alya_data(world, lsize)
    end

    # Launch Driver with both Local and World communicators.
    # We use Base.@invokelatest to avoid world-age issues with included files.
    with_mpi(; comm=local_comm) do distribute
        Base.@invokelatest driver(lsize, distribute, inputs, OUTPUT_DIR, TFloat, world, is_coupled)

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

if iszero(ccall(:jl_generating_output, Cint, ())) && !haskey(ENV, "JEXPRESSO_SYSIMAGE_BUILD")
    if !haskey(ENV, "JEXPRESSO_COUPLING_MODE")
        jexpresso_main()
    else
        jexpresso_main()
        MPI.Finalize()
    end
end
