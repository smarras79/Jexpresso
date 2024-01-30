using BenchmarkTools

function time_loop!(QT,            #Quadrature type: Inexact() vs Exaxt()
                    PT,            #Problem type: this is the name of the directory such as CompEuler
                    SOL_VARS_TYPE, #TOTAL() vs PERT() for total vs perturbation solution variables
                    CL,            #law type: CL() vs NCL() for conservative vs not-conservative forms
                    AD,            #AbstractDiscretization: ContGal(), DiscGal(), FD()
                    mesh::St_mesh,
                    metrics,
                    basis, ω,
                    qp::St_SolutionVars,
                    M, Minv, 
                    Δt,
                    inputs::Dict,
                    OUTPUT_DIR::String,
                    T;fx=zeros(Float64,1,1), fy = zeros(Float64,1,1),fy_lag = zeros(Float64,1,1))
    
    #
    # ODE: solvers come from DifferentialEquations.j;
    #
    # Initialize
    println(" # Solving ODE ................................")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
    lexact_integration = inputs[:lexact_integration]
    # extract backend
    backend = inputs[:backend]
    #-----------------------------------------------------------------
    # Initialize:
    # u     -> solution array
    # uaux  -> solution auxiliary array (Eventually to be removed)
    # F,G   -> physical flux arrays
    # S     -> source array
    # rhs* -> inviscid and viscous ELEMENT rhs
    # RHS* -> inviscid and viscous GLOBAL  rhs
    #-----------------------------------------------------------------
    u            = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin)*qp.neqs)
    uaux         = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs) 
    uaux_el      = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    rhs_el       = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    rhs_diff_el  = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    rhs_diffξ_el = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    rhs_diffη_el = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    F            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    G            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    S            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    RHS          = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
    RHS_visc     = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
    vaux         = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin)) #generic auxiliary array for general use
    
    #The following are currently used by B.C.
    gradu      = KernelAbstractions.zeros(backend, T, 2, 1, 1) #zeros(2,Int64(mesh.npoin),nvars)
    ubdy       = KernelAbstractions.zeros(backend, T, qp.neqs)
    bdy_flux   = KernelAbstractions.zeros(backend, T, qp.neqs,1)    
    uprimitive = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs+1)
    #filter arrays
    q_t = KernelAbstractions.zeros(backend, T,qp.neqs,Int64(mesh.ngl),Int64(mesh.ngl))
    q_ti = KernelAbstractions.zeros(backend, T,Int64(mesh.ngl),Int64(mesh.ngl))
    #fy_t = KernelAbstractions.allocate(backend, T, Int64(mesh.ngl),Int64(mesh.ngl))
    #KernelAbstractions.copyto!(backend, fy_t, transpose(fy))
    #fy_t_lag = KernelAbstractions.allocate(backend, T, Int64(mesh.ngr), Int64(mesh.ngl))
    #KernelAbstractions.copyto!(backend, fy_t_lag, transpose(fy_lag))
    fy_t = transpose(fy)
    fy_t_lag = transpose(fy_lag)
    fqf = KernelAbstractions.zeros(backend, T,qp.neqs,Int64(mesh.ngl),Int64(mesh.ngl))
    b = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngl), qp.neqs)
    B = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs) 
    #store grid limits to save time
    xmax = maximum(mesh.x)
    xmin = minimum(mesh.x)
    ymax = maximum(mesh.y)
    ymin = minimum(mesh.y)
    
    #The following are only built and active if Laguerre boundaries are to be used
    if ( "Laguerre" in mesh.bdy_edge_type)
        uaux_el_lag      = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        rhs_el_lag       = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        rhs_diff_el_lag  = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        rhs_diffξ_el_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        rhs_diffη_el_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        F_lag            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        G_lag            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        S_lag            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        RHS_lag          = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
        RHS_visc_lag     = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
        uprimitive_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs+1)
        q_t_lag = KernelAbstractions.zeros(backend, T,qp.neqs,Int64(mesh.ngl),Int64(mesh.ngr))
        q_ti_lag = KernelAbstractions.zeros(backend, T,Int64(mesh.ngl),Int64(mesh.ngr))
        fqf_lag = KernelAbstractions.zeros(backend, T,qp.neqs,Int64(mesh.ngl),Int64(mesh.ngr))
        b_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngl), Int64(mesh.ngr), qp.neqs)
        B_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
    end
    if (inputs[:llaguerre_1d])
        uaux_el_lag      = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        rhs_el_lag       = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        rhs_diff_el_lag  = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        rhs_diffξ_el_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        rhs_diffη_el_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        F_lag            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        G_lag            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        S_lag            = KernelAbstractions.zeros(backend, T, Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        RHS_lag          = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
        RHS_visc_lag     = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
        uprimitive_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs+1)
        q_t_lag = KernelAbstractions.zeros(backend, T, qp.neqs,Int64(mesh.ngr),Int64(mesh.ngl))
        q_ti_lag = KernelAbstractions.zeros(backend, T,Int64(mesh.ngr),Int64(mesh.ngl))
        fqf_lag = KernelAbstractions.zeros(backend, T,qp.neqs,Int64(mesh.ngr),Int64(mesh.ngl))
        b_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.nelem), Int64(mesh.ngr), Int64(mesh.ngl), qp.neqs)
        B_lag = KernelAbstractions.zeros(backend, T, Int64(mesh.npoin), qp.neqs)
    end

    #-----------------------------------------------------------------
    for i=1:qp.neqs
        idx = (i-1)*mesh.npoin
        u[idx+1:i*mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1[:,i] = @view(qp.qn[:,i])
        qp.qnm2[:,i] = @view(qp.qn[:,i])
        
    end
    
    deps  = zeros(1,1)
    tspan = (inputs[:tinit], inputs[:tend])    
    visc_coeff = inputs[:μ]#(inputs[:νρ], inputs[:νx], inputs[:νy], inputs[:κ], inputs[:κ], inputs[:κ], inputs[:κ])
    ivisc_equations = inputs[:ivisc_equations]   
 
    if ("Laguerre" in mesh.bdy_edge_type || inputs[:llaguerre_1d])
 
        params = (T, F, G, S,
                  uaux, uaux_el, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el,
                  uprimitive,
                  q_t, q_ti, fqf, b, B,
                  q_t_lag, q_ti_lag, fqf_lag, b_lag, B_lag,
                  RHS, RHS_visc,
                  F_lag, G_lag, S_lag, 
                  uaux_el_lag,
                  rhs_el_lag,
                  rhs_diff_el_lag,
                  rhs_diffξ_el_lag, rhs_diffη_el_lag,
                  RHS_lag, RHS_visc_lag, uprimitive_lag, 
                  SD=mesh.SD, QT, CL, PT, AD,
                  SOL_VARS_TYPE,
                  neqs=qp.neqs,
		  basis=basis[1], basis_lag = basis[2], ω = ω[1], ω_lag = ω[2], mesh, metrics = metrics[1], metrics_lag = metrics[2], 
                  inputs, visc_coeff, ivisc_equations,
                  M, Minv,
                  Δt, deps, xmax, xmin, ymax, ymin,
                  qp.qe, qp.qnm1, qp.qnm2, qp.μ,fx,fy, fy_t, fy_lag, fy_t_lag, laguerre=true)
    else
        params = (T, F, G, S,
                  uaux, uaux_el, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el,
                  uprimitive,
                  q_t, q_ti, fqf, b, B,
                  RHS, RHS_visc,
                  SD=mesh.SD, QT, CL, PT, AD, 
                  SOL_VARS_TYPE, 
                  neqs=qp.neqs,
                  basis, ω, mesh, metrics,
                  inputs, visc_coeff, ivisc_equations,
                  M, Minv,
                  Δt, deps, xmax, xmin, ymax, ymin,
                  qp.qe, qp.qnm1, qp.qnm2, qp.μ,fx,fy,fy_t,laguerre=false)
    end 
    
    prob = ODEProblem(rhs!,
                      u,
                      tspan,
                      params);

    output_range = floor((inputs[:tend] - inputs[:tinit])/inputs[:ndiagnostics_outputs])

    if typeof(mesh.SD) == NSD_2D
        if(lexact_integration)
            N = mesh.ngl
            Q = N + 1
            mass_ini   = compute_mass!(uaux, u, params.qe, mesh, metrics, ω, qp.neqs, QT, Q, basis.ψ)
        else
            mass_ini   = compute_mass!(uaux, u, params.qe, mesh, metrics, ω, qp.neqs, QT, SOL_VARS_TYPE)
        end
        println(" # Initial Mass  :   ", mass_ini)
        energy_ini = 0.0
        #energy_ini = compute_energy!(uaux, u, params.qe, mesh, metrics, ω,qp.neqs)
        #println(" # Initial Energy: ", energy_ini)
    end

   
    @time solution = solve(prob,
                           inputs[:ode_solver], dt=inputs[:Δt],
                           save_everystep = false,
                           adaptive=inputs[:ode_adaptive_solver],
                           saveat = range(inputs[:tinit], inputs[:tend], length=inputs[:ndiagnostics_outputs]));
    println(" # Solving ODE  ................................ DONE")


    if typeof(mesh.SD) == NSD_2D
        println(" # Diagnostics  ................................ ")
        if ("Laguerre" in mesh.bdy_edge_type)
            #    print_diagnostics(mass_ini, energy_ini, uaux, solution, mesh, metrics, ω, qp.neqs,QT,basis[1].ψ)
        else
            #   print_diagnostics(mass_ini, energy_ini, uaux, solution, mesh, metrics, ω, qp.neqs,QT,basis.ψ)
        end
        println(" # Diagnostics  ................................ DONE")
    end

    return solution
end
