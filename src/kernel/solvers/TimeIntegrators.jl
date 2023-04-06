using LinearAlgebra
using DiffEqBase
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
using SnoopCompile
import SciMLBase
using WriteVTK

include("../abstractTypes.jl")
include("../../io/write_output.jl")

function time_loop!(SD,
                    QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M,
                    De, Le,
                    Nt, Δt,
                    neqs, 
                    inputs::Dict,
                    OUTPUT_DIR::String,
                    T)
    
    #
    # ODE: solvers come from DifferentialEquations.j;
    #
    # Initialize
    println(" # Solving ODE ................................")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
    u = zeros(T, mesh.npoin*neqs);
    for i=1:neqs
        idx = (i-1)*mesh.npoin
        u[idx+1:i*mesh.npoin] .= qp.qn[:,i]
    end
    #@info neqs
    tspan  = (inputs[:tinit], inputs[:tend])    
    params = (; T, SD, QT, PT, neqs, basis, ω, mesh, metrics, inputs, M, De, Le, Δt)
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
