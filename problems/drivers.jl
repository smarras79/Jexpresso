# PERF: `using Distributions` removed — the file never references
# Distributions.* (the EL sample loop uses Base.rand). Distributions
# pulls in StatsBase / StatsFuns / SpecialFunctions for nothing.
#
# `using ONNXRunTime` is kept inside `elementLearningStructs.jl`,
# which is the only place that actually calls into it. drivers.jl
# itself never references the symbol, so we don't need to pay for
# the ~150 MB ONNX runtime binary on every rank just to import this
# file for non-EL runs (city2d, etc).

function driver(nparts,
                distribute,
                inputs,
                OUTPUT_DIR::String,
                TFloat;
                world      = nothing,
                is_coupled::Bool = false)

    comm = distribute.comm
    rank = MPI.Comm_rank(comm)

    # PERF/UX: narrow progress prints so the user can see where the
    # silent-wall first-call JIT time is going (sem_setup, mesh read,
    # metric build, …). This prints every phase boundary on rank 0; if
    # one of them takes obviously longer than the others, that's the
    # JIT hotspot. Cheap: just five flushed prints + a counter.
    if rank == 0
        print(" # driver() entered, calling sem_setup ......... ")
        flush(stdout)
    end
    _t_sem = time_ns()

    #---------------------------------------------------------
    # Time span (shared by both standalone and coupled paths).
    #---------------------------------------------------------
    if inputs[:lamr] == true
        amr_freq = inputs[:amr_freq]
        Δt_amr   = amr_freq * inputs[:Δt]
        tspan    = [TFloat(inputs[:tinit]), TFloat(inputs[:tinit] + Δt_amr)]
    else
        tspan    = [TFloat(inputs[:tinit]), TFloat(inputs[:tend])]
    end

    # PERF: during package precompilation (PrecompileTools @compile_workload)
    # this driver runs only to bake the hot-path JIT — RHS, the SciML
    # integrator, the callback-specialized warm-up — into the precompile
    # cache. The actual time integration is throwaway there, so cap the run
    # to a handful of steps. A full sod1d pass is 2000 steps (tend=0.2,
    # Δt=1e-4); running all of them during precompile is what made it take
    # minutes. The warm-up pre-pass plus the first few real steps still
    # trigger every specialization we want cached, so cold runs stay fast
    # without precompilation paying for a full simulation. `jl_generating_
    # output` is 1 only while generating precompile output, so real runs are
    # never shortened.
    if ccall(:jl_generating_output, Cint, ()) == 1 && haskey(inputs, :Δt)
        _dt_pc = TFloat(inputs[:Δt])
        tspan  = [TFloat(inputs[:tinit]), TFloat(inputs[:tinit]) + 3 * _dt_pc]
    end

    #---------------------------------------------------------
    # Mesh + initial state (+ coupling object in MPMD mode).
    #
    #   - Coupled path:  setup_coupling_and_mesh does sem_setup +
    #                    convert_mesh_arrays! + initialize internally and
    #                    returns the populated CouplingData.
    #   - Standalone:    historical sem_setup → (lRT_problem early exit) →
    #                    initialize. lwarmup pre-pass kept for memory-bound
    #                    runs that need a coarse-mesh JIT warmup.
    #---------------------------------------------------------
    coupling = nothing
    if is_coupled
        @assert world !== nothing "world communicator must be supplied when is_coupled=true"
        coupling, sem, partitioned_model, qp = setup_coupling_and_mesh(world, nparts, inputs, nparts, distribute, rank, OUTPUT_DIR, TFloat)
    else
        if inputs[:lwarmup] == true
            if rank == 0 println(BLUE_FG(string(" # JIT pre-compilation of large problem ..."))) end

            input_mesh = inputs[:gmsh_filename]
            # inputs may be a Dict (mutable) or a NamedTuple (immutable, when
            # :use_named_tuples => true).  Swap the mesh field via the
            # container's idiomatic update path.
            if inputs isa NamedTuple
                inputs_warm = (; inputs..., gmsh_filename = inputs[:gmsh_filename_c])
                sem_dummy   = sem_setup(inputs_warm, nparts, distribute)
            else
                inputs[:gmsh_filename] = inputs[:gmsh_filename_c]
                sem_dummy = sem_setup(inputs, nparts, distribute)
                inputs[:gmsh_filename] = input_mesh
            end
            sem_dummy = nothing

            if rank == 0 println(BLUE_FG(string(" # JIT pre-compilation of large problem ... DONE"))) end
        end

        sem, partitioned_model = sem_setup(inputs, nparts, distribute)

        if rank == 0
            @printf("DONE (%.2f s)\n", (time_ns() - _t_sem) / 1e9)
            print(" # initialize() ......... ")
            flush(stdout)
        end
        _t_init = time_ns()

        if inputs[:backend] != CPU()
            convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs)
        end

        # lRT_problem builds its own problem and exits early — it does not
        # use params_setup or time_loop!. Standalone-only branch.
        if inputs[:lRT_problem]
            # PERF: bring RRTMGP + ClimaComms + LinearOperators +
            # NCDatasets into scope. These were eagerly
            # loaded at the top of src/Jexpresso.jl, costing every
            # non-RT run (city2d, sod1d, …) ~tens of seconds of load
            # time for code they never call. _ensure_rt_loaded!() is a
            # cheap no-op once the first RT run has triggered it.
            Jexpresso._ensure_rt_loaded!()
            build_rad(sem, inputs)
            return
        end

        qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)
        if rank == 0
            @printf("DONE (%.2f s)\n", (time_ns() - _t_init) / 1e9)
            flush(stdout)
        end
    end

    #---------------------------------------------------------
    # Parameters setup (shared).
    #---------------------------------------------------------
    if rank == 0 println(" # Params_setup ..................................") end

    params, u = params_setup(sem, qp, inputs, OUTPUT_DIR, TFloat, tspan; coupling = coupling)

    if rank == 0 println(" # Params_setup .................................. DONE") end

    # test of projection matrix for solutions from old to new, i.e., coarse to fine, fine to coarse
    # test_projection_solutions(sem.mesh, qp, sem.partitioned_model, inputs, nparts, sem.distribute)

    #---------------------------------------------------------
    # Time integration / linear solve.
    #---------------------------------------------------------
    if !inputs[:llinsolve]
        #-----------------------------------------------------------------------------------
        # Hyperbolic/parabolic problems that lead to Mdq/dt = RHS.
        # is_coupled / coupling are no-ops for standalone runs.
        #-----------------------------------------------------------------------------------

        # Precompile warm-up: run a 1-step solve OUTSIDE the @time block
        # so its JIT-compile allocations are not counted in the displayed
        # total. The real `time_loop!` call below then starts with all
        # the hot-path methods already specialized. Disable via
        # JEXPRESSO_PRECOMPILE_WARMUP=0 (env), --no-precompile-warmup
        # (CLI), or :lprecompile_warmup => false (user_inputs.jl).
        precompile_warmup_run!(inputs, params, u, partitioned_model,
                               is_coupled, coupling)

        @time solution = time_loop!(inputs, params, u, partitioned_model,
                                    is_coupled, coupling)

        # PLOT NOTICE: Plotting is called from inside time_loop using callbacks.

    else
        #-----------------------------------------------------------------------------------
        # Problems that lead to Lx = RHS (standalone only).
        #-----------------------------------------------------------------------------------
        if (inputs[:backend] == CPU())
            if inputs[:lelementLearning]
                element_learning_linsolve!(sem, params, qp, inputs, OUTPUT_DIR, TFloat, rank)
            else
                standard_linsolve!(sem, params, qp, inputs, OUTPUT_DIR)
            end
        else
            println( " ")
            println( " WARNING!!! drivers.jl:L114")
            println( " WARNING: CHECK IF THIS GPU IMPLEMENTATION OF Ax=b still works")
            println( " ")
            nothing
        end
    end
end
