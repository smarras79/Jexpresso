include("../abstractTypes.jl")

function time_loop!(QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M,
                    De, Le,
                    Δt,
                    inputs::Dict,
                    OUTPUT_DIR::String,
                    T)
    
    #
    # ODE: solvers come from DifferentialEquations.j;
    #
    # Initialize
    println(" # Solving ODE ................................")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
    u = zeros(T, mesh.npoin*qp.neqs);
    global q1 = zeros(T, mesh.npoin, qp.neqs);
    global q2 = zeros(T, mesh.npoin, qp.neqs);
    global q3 = zeros(T, mesh.npoin, qp.neqs);
    global zb = zeros(T, mesh.npoin);
    
    for i=1:qp.neqs
        idx = (i-1)*mesh.npoin
        u[idx+1:i*mesh.npoin] .= qp.qn[:,i]
        global q1[:,i] .= qp.qn[:,i]
        global q2[:,i] .= qp.qn[:,i]
        global q3[:,i] .= qp.qn[:,i]
    end
    if (typeof(PT) == ShallowWater)
        for i=1:mesh.npoin
            global zb[i] = bathymetry(mesh.x[i])
        end
    end
    deps = zeros(1,1)
    tspan  = (inputs[:tinit], inputs[:tend])
    params = (; T, SD=mesh.SD, QT, PT, neqs=qp.neqs, basis, ω, mesh, metrics, inputs, M, De, Le, Δt, deps)
    prob   = ODEProblem(rhs!,
                        u,
                        tspan,
                        params);
    
    @time    solution = solve(prob,
                              inputs[:ode_solver],
                              dt = inputs[:Δt],
                              save_everystep=false,
                              saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]),
                              progress=true,
                              progress_message = (dt, u, p, t) -> t)
    
    println(" # Solving ODE  ................................ DONE")
    
    return solution
    
end
