using BenchmarkTools

function time_loop!(QT,
                    PT,
                    CL,
                    mesh::St_mesh,
                    metrics::St_metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M, Minv, 
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
    lexact_integration = inputs[:lexact_integration]
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
    gradu      = zeros(T, 2, 1, 1) #zeros(2,mesh.npoin,nvars)
    ubdy       = zeros(qp.neqs)
    bdy_flux   = zeros(qp.neqs,1)    
    uprimitive = zeros(T, mesh.ngl, mesh.ngl, qp.neqs+1)
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
              SD=mesh.SD, QT, PT, CL,
              neqs=qp.neqs,
              basis, ω, mesh, metrics,
              inputs, visc_coeff,              
              M, Minv,
              Δt, deps,
              qp.qe, qp.qnm1, qp.qnm2, qp.μ)
    
    prob = ODEProblem(rhs!,
                      u,
                      tspan,
                      params);

    output_range = floor((inputs[:tend] - inputs[:tinit])/inputs[:ndiagnostics_outputs])

    if(lexact_integration)
       N = mesh.ngl
       Q = N + 1
       mass_ini   = compute_mass!(uaux, u, mesh, metrics, ω,qp.neqs,QT,Q,basis.ψ)
    else
       mass_ini   = compute_mass!(uaux, u, mesh, metrics, ω,qp.neqs,QT)
    end

    energy_ini = compute_energy!(uaux, u, mesh, metrics, ω,qp.neqs)
    println(" # Initial Mass  :   ", mass_ini)
    println(" # Initial Energy: ", energy_ini)
    
    @time solution = solve(prob,
                           inputs[:ode_solver], dt=inputs[:Δt],
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    println(" # Solving ODE  ................................ DONE")

    println(" # Diagnostics  ................................ ")
    print_diagnostics(mass_ini, energy_ini, uaux, solution, mesh, metrics, ω, qp.neqs,QT,basis.ψ)
    println(" # Diagnostics  ................................ DONE")
    
    return solution
end
