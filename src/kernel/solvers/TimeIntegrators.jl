using BenchmarkTools

function time_loop!(inputs, params, u)

    println(" # Solving ODE  ................................ ")
    
    prob = ODEProblem(rhs!,
                      u,
                      params.tspan,
                      params);
    
#=
    #analysis_interval = 100
    #analysis_callback = AnalysisCallback(semi, interval = analysis_interval,
    #                                     extra_analysis_integrals = (entropy,))

    #alive_callback = AliveCallback(analysis_interval = analysis_interval)

    save_solution = SaveSolutionCallback(interval = 100,
                                         save_initial_solution = true,
                                         save_final_solution = true,
                                         solution_variables = cons2prim)

    amr_controller = ControllerThreeLevel(semi, IndicatorMax(semi, variable = first),
                                          base_level = 4,
                                          med_level = 5, med_threshold = 0.1,
                                          max_level = 6, max_threshold = 0.6)
    amr_callback = AMRCallback(semi, amr_controller,
                               interval = 5,
                               adapt_initial_condition = true,
                               adapt_initial_condition_only_refine = true)

    stepsize_callback = StepsizeCallback(cfl = 1.6)

    callbacks = CallbackSet(summary_callback,
                            analysis_callback, alive_callback,
                            save_solution,
                            amr_callback, stepsize_callback);=#
                            
    # Initialize the solution variable
    solution = nothing

    try
        @time solution = solve(prob,
                            inputs[:ode_solver], dt=Float32(inputs[:Δt]),
                            save_everystep = false,
                            adaptive=inputs[:ode_adaptive_solver],
                            saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    catch e
        # println("Instability detected: ", e)
        if isa(e, SciMLBase.InstabilityException)
            solution = e.sol
        else
            println("Error detected: ", e)
        end
    end
    
    if solution === nothing
        error("Solution could not be obtained due to an error.")
    end

    
  #=  summary_callback = SummaryCallback()

    =#


    
    println(" # Solving ODE  ................................ DONE")
    
    return solution
end
