using LinearAlgebra
using DiffEqBase
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
using SnoopCompile
import SciMLBase
using WriteVTK

include("../abstractTypes.jl")

function time_loop!(QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M,
                    De, Le,
                    Nt, Δt,
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
    for i=1:qp.neqs
        idx = (i-1)*mesh.npoin
        u[idx+1:i*mesh.npoin] .= qp.qn[:,i]
        global q1[:,i] .= qp.qn[:,i]
        global q2[:,i] .= qp.qn[:,i]
    end
    #@info qp.neqs
    tspan  = (inputs[:tinit], inputs[:tend])
    params = (; T, SD=mesh.SD, QT, PT, neqs=qp.neqs, basis, ω, mesh, metrics, inputs, M, De, Le, Δt)
    prob   = ODEProblem(rhs!,
                        u,
                        tspan,
                        params);
    
    @time    solution = solve(prob,
                              inputs[:ode_solver],
                              dt = Δt,
                              save_everystep=false,
                              saveat = range(T(0.), Nt*T(Δt), length=inputs[:ndiagnostics_outputs]),
                              progress = true,
                              progress_message = (dt, u, p, t) -> t)
    println(" # Solving ODE  ................................ DONE")
    
    return solution
    
end
