include("../abstractTypes.jl")
using BenchmarkTools

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

    #-----------------------------------------------------------------
    # Initialize:
    # u     -> solution array
    # uaux  -> solution auxiliary array (Eventually to be removed)
    # F,G   -> physical flux arrays
    # S     -> source array
    # rhs* -> inviscid and viscous ELEMENT rhs
    # RHS* -> inviscid and viscous GLOBAL  rhs
    #-----------------------------------------------------------------
    u           = zeros(T, mesh.npoin*qp.neqs);
    uaux        = zeros(T, mesh.npoin, qp.neqs)
    F           = zeros(T, mesh.ngl, mesh.ngl, qp.neqs)
    G           = zeros(T, mesh.ngl, mesh.ngl, qp.neqs)
    S           = zeros(T, mesh.ngl, mesh.ngl, qp.neqs)
    rhs_el      = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem, qp.neqs)
    rhs_diff_el = zeros(T, mesh.ngl, mesh.ngl, mesh.nelem, qp.neqs)
    RHS         = zeros(T, mesh.npoin, qp.neqs)
    RHS_visc    = zeros(T, mesh.npoin, qp.neqs)

    #The following are currently used by B.C.
    gradu       = zeros(2, 1, 1) #zeros(2,mesh.npoin,nvars)
    ubdy        = zeros(qp.neqs)
    bdy_flux    = zeros(qp.neqs,1)
    #-----------------------------------------------------------------
    
    for i=1:qp.neqs
        idx = (i-1)*mesh.npoin
        u[idx+1:i*mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1 .= qp.qn[:,i]
        qp.qnm2 .= qp.qn[:,i]
        
    end
    
    deps = zeros(1,1)
    tspan  = (inputs[:tinit], inputs[:tend])
    
    params = (T, F, G, S,
              uaux, gradu, ubdy, bdy_flux,
              rhs_el, rhs_diff_el,
              RHS, RHS_visc, 
              SD=mesh.SD, QT, PT,
              neqs=qp.neqs,
              basis, ω, mesh, metrics, inputs,
              M, De, Le,
              Δt, deps,
              qp.qnm1, qp.qnm2, qp.μ)
    
    prob = ODEProblem(rhs!,
                      u,
                      tspan,
                      params);
    
    @time solution = solve(prob,
                           inputs[:ode_solver], dt=inputs[:Δt],
                           save_everystep = false,
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    
    println(" # Solving ODE  ................................ DONE")
    
    return solution
    
end
