using BenchmarkTools

function time_loop!(QT,
                    PT,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M,
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
    u            = zeros(T, mesh.npoin*qp.neqs)
    uaux         = zeros(T, mesh.npoin, qp.neqs)
    uaux_el      = zeros(T, mesh.nelem, mesh.ngl, mesh.ngl, qp.neqs)
    rhs_el       = zeros(T, mesh.nelem, mesh.ngl, mesh.ngl, qp.neqs)
    rhs_diff_el  = zeros(T, mesh.nelem, mesh.ngl, mesh.ngl, qp.neqs)
    rhs_diffξ_el = zeros(T, mesh.nelem, mesh.ngl, mesh.ngl, qp.neqs)
    rhs_diffη_el = zeros(T, mesh.nelem, mesh.ngl, mesh.ngl, qp.neqs)
    F            = zeros(T, mesh.ngl, mesh.ngl, qp.neqs)
    G            = zeros(T, mesh.ngl, mesh.ngl, qp.neqs)
    S            = zeros(T, mesh.ngl, mesh.ngl, qp.neqs)
    RHS          = zeros(T, mesh.npoin, qp.neqs)
    RHS_visc     = zeros(T, mesh.npoin, qp.neqs)

    #The following are currently used by B.C.
    gradu       = zeros(2, 1, 1) #zeros(2,mesh.npoin,nvars)
    ubdy        = zeros(qp.neqs)
    bdy_flux    = zeros(qp.neqs,1)
    
    uprimitive = zeros(T, mesh.ngl, mesh.ngl, qp.neqs)
        
    #-----------------------------------------------------------------
    
    for i=1:qp.neqs
        idx = (i-1)*mesh.npoin
        u[idx+1:i*mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1[:,i] = @view(qp.qn[:,i])
        qp.qnm2[:,i] = @view(qp.qn[:,i])
        
    end
    
    deps = zeros(1,1)
    tspan  = (inputs[:tinit], inputs[:tend])    
    visc_coeff = (inputs[:νρ], inputs[:νx], inputs[:νy], inputs[:κ])
    
    params = (T, F, G, S,
              uaux, uaux_el,
              ubdy, gradu, bdy_flux, #for B.C.
              rhs_el, rhs_diff_el,
              rhs_diffξ_el, rhs_diffη_el,
              uprimitive,
              RHS, RHS_visc, 
              SD=mesh.SD, QT, PT,
              neqs=qp.neqs,
              basis, ω, mesh, metrics,
              inputs, visc_coeff,              
              M, Δt, deps,
              qp.qnm1, qp.qnm2, qp.μ)
    
    prob = ODEProblem(rhs!,
                      u,
                      tspan,
                      params);

    output_range = floor((inputs[:tend] - inputs[:tinit])/inputs[:ndiagnostics_outputs])
    @time solution = solve(prob,
                           inputs[:ode_solver], dt=inputs[:Δt],
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
                           #save_start=false, save_end=true);
    
    println(" # Solving ODE  ................................ DONE")

    return solution
end
