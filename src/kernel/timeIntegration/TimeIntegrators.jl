using LinearAlgebra
using DiffEqBase
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
using LinearSolve
using SnoopCompile
import SciMLBase

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
                    neqns, 
                    inputs::Dict,
                    BCT,
                    OUTPUT_DIR::String,
                    T)
    
    #
    # ODE: solvers come from DifferentialEquations.j;
    #
    #Initialize
    u      = zeros(T, mesh.npoin);
    u     .= qp.qn[:,1];
    tspan  = (inputs[:tinit], inputs[:tend])    
    params = (; T, SD, QT, PT, BCT, neqns, basis, ω, mesh, metrics, inputs, M, De, Le)
    prob   = ODEProblem(rhs!,
                        u,
                        tspan,
                        params);
    
    println(" # Solving ODE with ................................")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
    
    @time    sol = solve(prob,
                         inputs[:ode_solver],
                         dt = Δt,
                         save_everystep=false,
                         saveat = range(T(0.), Nt*T(Δt), length=inputs[:ndiagnostics_outputs]),
                         progress = true,
                         progress_message = (dt, u, p, t) -> t)
    println(" # Solving ODE with  ................................ DONE")
    
    write_output(sol, SD, mesh, OUTPUT_DIR, inputs, inputs[:outformat])
    
end

function solveAx!(SD,
                  QT,
                  PT,
                  mesh::St_mesh,
                  metrics::St_metrics,
                  basis, ω,
                  qp::St_SolutionVars,
                  M, L,
                  neqns,
                  inputs::Dict,
                  BCT,
                  OUTPUT_DIR::String,
                  T)
    
    #
    # Ax=b solve come from LinearSolve.jl;
    #
    println(" # Solve Ax=b ................................")
    println(" #  ", inputs[:ode_solver])
    # Naive B.C.
    ϵ = eps(Float32)
    for ip=1:mesh.npoin
        x, y = mesh.x[ip], mesh.y[ip]
        if( (x > 1.0 - ϵ) || (x < -1.0 + ϵ))
            qp.qn[ip,1] = sinpi(2*y)
            for jp=1:mesh.npoin
                L[ip,jp] = 0.0
            end
            L[ip,ip] = 1.0
        end
        if( (y > 1.0 - ϵ) || (y < -1.0 + ϵ))
            qp.qn[ip,1] = 0.0
            for jp=1:mesh.npoin
                L[ip,jp] = 0.0
            end
            L[ip,ip] = 1.0
        end        
    end
    #END Naive B.C..
    
    RHS = build_rhs_source(SD,
                          QT,
                          PT,
                          qp.qn,
                          mesh,
                          M, #M is sparse for exact integration
                          T)

    
    prob = LinearProblem(L, RHS);
    
    @time sol = solve(prob, inputs[:ode_solver])

    println(" # Solving Ax=b ................................ DONE")
    
    write_output(sol, SD, mesh, OUTPUT_DIR, inputs, inputs[:outformat])
    
    
end
