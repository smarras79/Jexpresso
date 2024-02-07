using BenchmarkTools


function time_loop!(inputs, params, u)
    
    prob = ODEProblem(rhs!,
                      u,
                      params.tspan,
                      params);

    output_range = floor((inputs[:tend] - inputs[:tinit])/inputs[:ndiagnostics_outputs])

    if typeof(params.mesh.SD) == NSD_2D
        if(inputs[:lexact_integration])
            N = params.mesh.ngl
            Q = N + 1
            mass_ini   = compute_mass!(params.uaux, u, params.qp.qe, params.mesh, params.metrics, params.ω, params.qp.neqs, params.QT, Q, params.basis.ψ)
        else
            mass_ini   = compute_mass!(params.uaux, u, params.qp.qe, params.mesh, params.metrics, params.ω, params.qp.neqs, params.QT, params.SOL_VARS_TYPE)
        end
        println(" # Initial Mass  :   ", mass_ini)
        energy_ini = 0.0
    end
    
    @time solution = solve(prob,
                           inputs[:ode_solver], dt=inputs[:Δt],
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    println(" # Solving ODE  ................................ DONE")

    
    return solution
end
