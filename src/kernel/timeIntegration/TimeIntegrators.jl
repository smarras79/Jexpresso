using LinearAlgebra
using DiffEqBase
using OrdinaryDiffEq
using OrdinaryDiffEq: SplitODEProblem, solve, IMEXEuler
import SciMLBase


include("../abstractTypes.jl")
include("../infrastructure/element_matrices.jl")
include("../../io/plotting/jeplots.jl")

function rhs_old!(du, u, params, t)

    T       = Float64
    SD      = params.SD
    QT      = params.QT
    PT      = params.PT
    BCT     = params.BCT
    neqns   = params.neqns
    basis   = params.basis
    mesh    = params.mesh
    metrics = params.metrics
    inputs  = params.inputs
    ω       = params.ω
    M       = params.M
    
    u[mesh.npoin_linear] = 0.0
    #
    # rhs[ngl,ngl,nelem]
    #
    rhs_el      = build_rhs(SD, QT, PT, neqns, u, basis.ψ, basis.dψ, ω, mesh, metrics, T)
    rhs_diff_el = build_rhs_diff(SD, QT, PT, neqns, u, basis.ψ, basis.dψ, ω, inputs[:νx], inputs[:νy], mesh, metrics, T)
    
    apply_boundary_conditions!(SD, rhs_el, u, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω,t, BCT, neqns)
    
    #du = DSSijk_rhs(SD, rhs_el + inputs[:δvisc]*rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, mesh.nop, T)
    du =     DSS(SD, QT, rhs_el + inputs[:δvisc]*rhs_diff_el, mesh.connijk, mesh.nelem, mesh.npoin, mesh.nop, T)

    divive_by_mass_matrix!(du, M, QT)
    
    apply_periodicity!(SD, rhs_el, u, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, t, BCT, neqns)
    
    @info "" maximum(u) maximum(du)
    #error("sasa")
    return du
end


function rhs!(du, u, params, time)

    #SD::NSD_1D, QT::Inexact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, De, u)
    T       = Float64
    SD      = params.SD
    QT      = params.QT
    PT      = params.PT
    BCT     = params.BCT
    neqns   = params.neqns
    basis   = params.basis
    mesh    = params.mesh
    metrics = params.metrics
    inputs  = params.inputs
    ω       = params.ω
    M       = params.M
    De      = params.De

    RHS = build_rhs(SD, QT, PT, BCT, neqns, u, basis.ψ, basis.dψ, ω, mesh, metrics, M, De, time, inputs, T)    
    du .= RHS
    
    return du #This is already DSSed
end



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
    # ODE
    #
    u = zeros(T, mesh.npoin);
    u .= qp.qn[:,1];
    params = (; T, SD, QT, PT, BCT, neqns, basis, ω, mesh, metrics, inputs, M, De)
    tspan = (inputs[:tinit], inputs[:tend])
    
    prob = ODEProblem(rhs!,
                      u,
                      tspan,
                      params);
    
    
    println(" # Solving ODE with ................................\n")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
    
    @time    sol = solve(prob,
                inputs[:ode_solver],
                dt = Δt,
                saveat = range(T(0.), Nt*T(Δt), length=inputs[:ndiagnostics_outputs]),
                progress = true,
                progress_message = (dt, u, p, t) -> t)
    println(" # Solving ODE with  ................................ DONE\n")

    
    
    for iout = 1: inputs[:ndiagnostics_outputs]
        title = @sprintf "Tracer: final solution at t=%6.4f" sol.t[iout]
        plot_results(SD, mesh.x, mesh.y, sol.u[iout], title, string(OUTPUT_DIR, "/it.", iout, ".png"))
    end
    
end
