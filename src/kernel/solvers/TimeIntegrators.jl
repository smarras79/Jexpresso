# Enhanced TimeIntegrators.jl with IMEX support for Jexpresso

# Import IMEX module functions
include("../operators/imex.jl")

"""
    amr_strategy_imex!(inputs, params, u, t, partitioned_model)
    
AMR strategy adapted for IMEX problems
"""
function amr_strategy_imex!(inputs, params, u, t, partitioned_model)
    # Call existing AMR but create SplitODEProblem instead of ODEProblem
    # This is a placeholder - implement based on your existing AMR code
    
    new_params, new_partitioned_model = amr_strategy!(inputs, params, u, t, partitioned_model)
    
    # Create new IMEX functions
    f_ex!, f_im!, jac_prototype = create_imex_functions(new_params, 
                                                       get(inputs, :imex_acoustic_implicit, true))
    
    # Create new SplitODEProblem
    new_prob = SplitODEProblem(f_ex!, f_im!, u, (t, inputs[:tend]), new_params)
    
    return new_prob, new_partitioned_model
end

function time_loop!(inputs, params, u, args...)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    partitioned_model = args[1]
    println_rank(" # Solving ODE  ................................ "; msg_rank = rank)
    
    # Check if IMEX solver is requested
    use_imex = get(inputs, :use_imex_solver, false)
    #imex_acoustic_implicit = get(inputs, :imex_acoustic_implicit, true)
    imex_acoustic_implicit = false
    if use_imex
        println_rank(" # Using IMEX time integration ................"; msg_rank = rank)
        
        # Create IMEX problem
        rhs_ex!, rhs_im!, jac_prototype = create_imex_functions(params, imex_acoustic_implicit)
        
        # Pass jac_prototype to SplitODEProblem, not to algorithm
        prob = SplitODEProblem(rhs_ex!, rhs_im!, u, params.tspan, params;) #jac_prototype=(jac_prototype, jac_prototype))
        
        # Select IMEX algorithm
        imex_algorithm = get_imex_algorithm(inputs, jac_prototype)
        
    else
        println_rank(" # Using explicit time integration ............"; msg_rank = rank)
        
        # Standard explicit problem
        prob = ODEProblem(rhs!, u, params.tspan, params)
        
        # Use specified explicit algorithm
        imex_algorithm = inputs[:ode_solver]
    end

    #------------------------------------------------------------------------
    # Runtime callbacks (unchanged from original)
    #------------------------------------------------------------------------
    dosetimes = inputs[:diagnostics_at_times]
    idx_ref   = Ref{Int}(0)
    c         = Float64(0.0)
    rad_time  = inputs[:radiation_time_step]
    lnew_mesh = true   
    
    function two_stream_condition(u, t, integrator)
        if (rem(t,rad_time) < 1e-3)
            return true
        else
            return false
        end
    end

    function do_radiation!(integrator)
        println(" doing two stream radiation heat flux calculations at t=", integrator.t)
        @info "doing rad test"
        compute_radiative_fluxes!(lnew_mesh, params.mesh, params.uaux, params.qp.qe, 
                                 params.mp, params.phys_grid, params.inputs[:backend], 
                                 params.SOL_VARS_TYPE)
    end
    
    ret_dosetime_ref  = Ref{Bool}(false)

    function condition(u, t, integrator)
        idx  = findfirst(x -> x == t, dosetimes)
        if idx !== nothing
            idx_ref[] = idx
            ret_dosetime_ref[] = true
        else
            ret_dosetime_ref[] = false
        end

        tol = 1e-6
        return ret_dosetime_ref[]
    end
    
    function affect!(integrator)
        idx          = idx_ref[]
        ret_dosetime = ret_dosetime_ref[]
        if ret_dosetime == true
            println_rank(" #  t=", integrator.t; msg_rank = rank)

            #CFL - modify for IMEX if needed
            if (inputs[:ladapt] == false && params.SOL_VARS_TYPE == TOTAL())
                if use_imex
                    # IMEX allows larger time steps due to implicit acoustic treatment
                    computeCFL_IMEX(params.mesh.npoin, integrator.p.qp.neqs,
                                   inputs[:Δt],
                                   params.mesh.Δeffective_s,
                                   integrator,
                                    params.SD; visc=inputs[:μ])
                else
                    computeCFL(params.mesh.npoin, integrator.p.qp.neqs,
                             inputs[:Δt],
                             params.mesh.Δeffective_s,
                             integrator,
                             params.SD; visc=inputs[:μ])
                end
            end
            
            write_output(integrator.p.SD, integrator.u, integrator.p.uaux, integrator.t, idx,
                         integrator.p.mesh, integrator.p.mp,
                         integrator.p.connijk_original, integrator.p.poin_in_bdy_face_original,
                         integrator.p.x_original, integrator.p.y_original, integrator.p.z_original,
                         inputs[:output_dir], inputs,
                         integrator.p.qp.qvars,
                         integrator.p.qp.qoutvars,
                         inputs[:outformat];
                         nvar=integrator.p.qp.neqs, qexact=integrator.p.qp.qe)
        end
    end
    
    cb_rad = DiscreteCallback(two_stream_condition, do_radiation!)
    cb     = DiscreteCallback(condition, affect!)    
    cb_amr = DiscreteCallback(condition, affect!)
    CallbackSet(cb)#,cb_rad)
    #------------------------------------------------------------------------
    # END runtime callbacks
    #------------------------------------------------------------------------

    #
    # Write initial conditions:
    #
    if rank == 0 
        solver_type = use_imex ? "IMEX" : "explicit"
        println(" # Write initial condition to ",  typeof(inputs[:outformat]), 
                " (", solver_type, " solver) .........")
    end
    
    write_output(params.SD, u, params.uaux, inputs[:tinit], 0,
                 params.mesh, params.mp,
                 params.connijk_original, params.poin_in_bdy_face_original,
                 params.x_original, params.y_original, params.z_original,
                 inputs[:output_dir], inputs,
                 params.qp.qvars, params.qp.qoutvars,
                 inputs[:outformat];
                 nvar=params.qp.neqs, qexact=params.qp.qe)
    
    if rank == 0  
        println(" # Write initial condition to ",  typeof(inputs[:outformat]), " ......... END") 
    end
    
    #
    # Simulation - Choose solver path
    #
    if use_imex
        # IMEX solver with adaptive time stepping based on advective CFL
        #solution = solve_imex_problem(prob, imex_algorithm, inputs, cb, dosetimes, rank)
        solution = solve_imex_problem(prob,  inputs[:ode_solver], inputs, cb, dosetimes, rank)
    else
        # Standard explicit solver (original path)
        solution = solve(prob,
                         inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                         callback = CallbackSet(cb), tstops = dosetimes,
                         save_everystep = false,
                         adaptive=false,
                         #adaptive=inputs[:ode_adaptive_solver],
                         saveat = range(inputs[:tinit], inputs[:tend], 
                                        length=inputs[:ndiagnostics_outputs]))
    end
    
    # AMR support (works with both explicit and IMEX)
    if inputs[:lamr] == true
        while solution.t[end] < inputs[:tend]
            if use_imex
                @time prob, partitioned_model = amr_strategy_imex!(inputs, prob.p, 
                                                                  solution.u[end][:], 
                                                                  solution.t[end], 
                                                                  partitioned_model)
                @time solution = solve_imex_problem(prob, imex_algorithm, inputs, 
                                                   cb_amr, dosetimes, rank)
            else
                @time prob, partitioned_model = amr_strategy!(inputs, prob.p, 
                                                             solution.u[end][:], 
                                                             solution.t[end], 
                                                             partitioned_model)
                @time solution = solve(prob,
                                     inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                                     callback = cb_amr, tstops = dosetimes,
                                     save_everystep = false,
                                     adaptive=inputs[:ode_adaptive_solver],
                                     saveat = [])
            end
        end
    end
    
    solver_type = use_imex ? "IMEX" : "explicit"
    println_rank(" # Solving ODE (", solver_type, ") .............. DONE"; msg_rank = rank)
    
    return solution
end
