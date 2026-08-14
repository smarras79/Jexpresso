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
    # Spherical shell (2D manifold embedded in 3D).
    #
    # The GRID is built by the ordinary gmsh path, exactly like every other
    # case: mod_mesh_mesh_driver → mod_mesh_read_gmsh!. That reader recognises
    # a 2D manifold embedded in 3D from the model itself (see `lmanifold` in
    # src/kernel/mesh/mesh.jl), keeps z through the high-order node placement,
    # and snaps the LGL points onto the shell. Gridap's topology already gives
    # one entry per UNIQUE edge, so the panel seams of a watertight cubed
    # sphere stitch together with no special handling.
    #
    # What stays specific to the manifold is downstream of the grid:
    #
    #   metrics (sphere_metrics.jl)  →  initial condition (the case's
    #   initialize.jl)  →  the :ode_solver time loop (sphere_time_loop.jl), whose RHS
    #   is the SEM surface divergence in sphere_rhs.jl fed by the case's
    #   user_flux! / user_source!.
    #
    # Two early exits, for working on the grid or the initial state alone:
    #
    #   :lgrid_only => true   build the grid, write it, return.
    #   :linit_only => true   build the grid AND the initial condition, write
    #                         both (into ONE file: the fields ride on the grid),
    #                         return without integrating.
    #
    # :lgrid_only wins if both are set.
    #---------------------------------------------------------
    if get(inputs, :lspherical_shell, false) == true

        #
        # A live Julia session keeps the ALREADY-COMPILED Jexpresso module.
        # run.jl re-includes THIS file on every run_case, but the src/ files are
        # only evaluated when the module itself is (re)loaded. Pulling new code
        # into a running session therefore leaves a fresh drivers.jl calling
        # functions the stale module has never seen — which surfaces as a bare
        # `UndefVarError: <name> not defined in Jexpresso` from the line below.
        # Say what it actually means.
        #
        for _w in (:write_vtk_sphere_grid,
                   :project_momentum_to_sphere!, :sphere_normal_momentum,
                   :build_sphere_metrics, :build_sphere_params, :sphere_time_loop!)
            isdefined(@__MODULE__, _w) && continue
            error(" # ERROR drivers.jl: `", _w, "` is not defined in the loaded Jexpresso module.\n",
                  " #   The module in this Julia session is older than the source tree on disk.\n",
                  " #   Restart Julia and run the case in a fresh session:\n",
                  " #     julia --project=.\n",
                  " #     julia> using Jexpresso\n",
                  " #     julia> Jexpresso.run_case(\"ShallowWater\", \"SWsphere\")")
        end

        if MPI.Comm_size(comm) > 1
            error(" # ERROR drivers.jl: the spherical-shell path is serial for now. Run it on a single MPI rank.")
        end

        if rank == 0
            println()
            println(" # :lspherical_shell => true — reading the shell grid through the ordinary gmsh path")
            flush(stdout)
        end

        smesh, _ = mod_mesh_mesh_driver(inputs, nparts, distribute)

        smesh.lmanifold || error(" # ERROR drivers.jl: :lspherical_shell => true but ",
                                 inputs[:gmsh_filename], " is not a curved 2D surface embedded in 3D.")

        # ONE file per run, as (ngl-1)² sub-elements per spectral element — the
        # same convention as write_vtk_grid_only for the flat cases. When the
        # initial condition has been built its fields ride on that same file.
        if get(inputs, :lgrid_only, false) == true
            rank == 0 && write_vtk_sphere_grid(smesh, "sphere_grid_ho", OUTPUT_DIR)
            if rank == 0
                println(" # :lgrid_only => true — grid built and written to ", OUTPUT_DIR, ". Stopping here.")
            end
            return smesh
        end

        #
        # Metric terms of the manifold: the contravariant basis that turns the
        # (F,G,H) of user_flux.jl into a SURFACE divergence, plus the diagonal
        # mass matrix. This is what the flat 2D metric machinery cannot supply.
        #
        smetrics = build_sphere_metrics(smesh, inputs; verbose = (rank == 0))
        if get(inputs, :lcheck_grid, true) == true
            ok = check_sphere_metrics(smesh, smetrics; verbose = (rank == 0))
            if !ok && get(inputs, :lstop_on_bad_grid, true) == true
                error(" # ERROR drivers.jl: the spherical shell metrics failed their consistency checks.")
            end
        end

        qsphere = initialize(smesh.SD, 0, smesh, inputs, OUTPUT_DIR, TFloat)

        #
        # Lagrange-multiplier projection, Marras/Kopera/Giraldo Eq. (9)-(11).
        #
        # In a real run this belongs at the end of every time step (or RK
        # stage): it removes the momentum component NORMAL to the shell that
        # the discrete operators accumulate, which is what keeps the fluid on
        # the sphere. There is no time loop yet, so applying it here does two
        # useful things instead — it reports how far off the manifold the
        # initial state is (it should be at round-off, since the Galewsky jet
        # is built from tangential unit vectors), and it guarantees that
        # whatever is written to VTK is exactly tangential.
        #
        if get(inputs, :llagrange_projection, true) == true

            drift = project_momentum_to_sphere!(qsphere.qn, smesh; ivar = 2)

            # qout (the primitives written to VTK) was derived from qn BEFORE
            # the projection, so refresh it through the case's own user_uout!
            # rather than shipping output that disagrees with the state.
            for ip = 1:smesh.npoin
                user_uout!(ip, inputs[:SOL_VARS_TYPE],
                           @view(qsphere.qout[ip, :]),
                           @view(qsphere.qn[ip, :]),
                           @view(qsphere.qe[ip, :]))
            end

            if rank == 0
                @printf(" # Lagrange projection P = I - xxᵀ/r²: removed max|(φu)·x̂| = %.3e\n", drift)
            end
        end

        if get(inputs, :linit_only, false) == true
            rank == 0 && write_vtk_sphere_grid(smesh, "sphere_grid_ho", OUTPUT_DIR; q = qsphere)
            if rank == 0
                println(" # :linit_only => true — grid + initial condition written to ", OUTPUT_DIR, ". Stopping here.")
            end
            return smesh, qsphere
        end

        #
        # Time integration. sphere_time_loop! writes the initial condition and
        # then one VTK file per output time, and applies the Lagrange projection
        # after every RK stage.
        #
        sparams = build_sphere_params(smesh, smetrics, inputs; neqs = qsphere.neqs)

        @time tfinal = sphere_time_loop!(smesh, smetrics, sparams, qsphere,
                                         inputs, OUTPUT_DIR; verbose = (rank == 0))

        return smesh, qsphere, tfinal
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
            flush(stdout)
        end

        # Grid-only run: dump the high-order grid and stop before the initial
        # condition. Same switch as the spherical-shell branch above, for the
        # cases that DO go through sem_setup.
        if get(inputs, :lgrid_only, false) == true
            write_vtk_grid_only(sem.mesh.SD, sem.mesh, "grid_ho", OUTPUT_DIR,
                                distribute(LinearIndices((nparts,))), nparts)
            if rank == 0
                println(" # :lgrid_only => true — grid built and written to ", OUTPUT_DIR, ". Stopping here.")
            end
            return sem
        end

        if rank == 0
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
            build_rad!(sem, partitioned_model, inputs, nparts, distribute)

            return
        end

        qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)
        if rank == 0
            @printf("DONE (%.2f s)\n", (time_ns() - _t_init) / 1e9)
            flush(stdout)
        end
    end
    #---------------------------------------------------------
    # Time span (shared by both standalone and coupled paths).
    #---------------------------------------------------------
    # Build tspan after initialize() so VTK/HDF5 restarts that update inputs[:tinit] are reflected.
    tspan = build_tspan(inputs, TFloat)

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
