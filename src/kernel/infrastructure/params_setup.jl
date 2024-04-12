function params_setup(sem,
                      qp::St_SolutionVars,
                      inputs::Dict,
                      OUTPUT_DIR::String,
                      T)

    #
    # ODE: solvers come from DifferentialEquations.j;
    #
    # Initialize
    println(" # Build arrays and params ................................ ")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
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
    u            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin)*Int64(qp.neqs))
    uaux         = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))

    if sem.mesh.nsd == 1
        uaux_el      = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_el       = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diff_el  = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffξ_el = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffη_el = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffζ_el = KernelAbstractions.zeros(backend, T, 0)
        F            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(qp.neqs))
        G            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(qp.neqs))
        H            = KernelAbstractions.zeros(backend, T, 0)
        S            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(qp.neqs))
        uprimitive   = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(qp.neqs)+1)
        flux_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem),Int64(sem.mesh.ngl), 2*qp.neqs)
        source_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), qp.neqs)
        qbdy_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nedges_bdy), Int64(sem.mesh.ngl), qp.neqs)
    elseif  sem.mesh.nsd == 2
        uaux_el      = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_el       = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diff_el  = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffξ_el = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffη_el = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffζ_el = KernelAbstractions.zeros(backend, T, 0)
        F            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        G            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        H            = KernelAbstractions.zeros(backend, T, 0)
        S            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        uprimitive   = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs)+1)
        flux_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem),Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), 2*qp.neqs)
        source_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem),Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), qp.neqs)
        qbdy_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nedges_bdy), Int64(sem.mesh.ngl), qp.neqs)
    elseif  sem.mesh.nsd == 3
        uaux_el      = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_el       = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diff_el  = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffξ_el = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffη_el = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffζ_el = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        F            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        G            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        H            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        S            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
        uprimitive   = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs)+1)
        flux_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), 2*qp.neqs)
        source_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), qp.neqs)
        qbdy_gpu = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nedges_bdy), Int64(sem.mesh.ngl), qp.neqs)
    end
    
    RHS      = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
    RHS_visc = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
    vaux     = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin)) #generic auxiliary array for general use
    
    #The following are currently used by B.C.
    gradu    = KernelAbstractions.zeros(backend, T, 2, 1, 1) #KernelAbstractions.zeros(2,Int64(sem.mesh.npoin),nvars)
    ubdy     = KernelAbstractions.zeros(backend, T, Int64(qp.neqs))
    bdy_flux = KernelAbstractions.zeros(backend, T, Int64(qp.neqs),1)    
   
    #filter arrays
    q_t  = KernelAbstractions.zeros(backend, T,Int64(qp.neqs),Int64(sem.mesh.ngl),Int64(sem.mesh.ngl))
    q_ti = KernelAbstractions.zeros(backend, T,Int64(sem.mesh.ngl),Int64(sem.mesh.ngl))
    fy_t = transpose(sem.fy)
    fy_t_lag = transpose(sem.fy_lag)
    fqf = KernelAbstractions.zeros(backend, T,Int64(qp.neqs),Int64(sem.mesh.ngl),Int64(sem.mesh.ngl))
    b = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngl), Int64(qp.neqs))
    B = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs)) 
    #store grid limits to save time
    xmax = maximum(sem.mesh.x)
    xmin = minimum(sem.mesh.x)
    ymax = maximum(sem.mesh.y)
    ymin = minimum(sem.mesh.y)

    
    #The following are only built and active if Laguerre boundaries are to be used
    if ( "Laguerre" in sem.mesh.bdy_edge_type)
        uaux_el_lag      = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        rhs_el_lag       = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        rhs_diff_el_lag  = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        rhs_diffξ_el_lag = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        rhs_diffη_el_lag = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        F_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        G_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        H_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        S_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        RHS_lag          = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
        RHS_visc_lag     = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
        uprimitive_lag = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs)+1)
        q_t_lag = KernelAbstractions.zeros(backend, T,Int64(qp.neqs),Int64(sem.mesh.ngl),Int64(sem.mesh.ngr))
        q_ti_lag = KernelAbstractions.zeros(backend, T,Int64(sem.mesh.ngl),Int64(sem.mesh.ngr))
        fqf_lag = KernelAbstractions.zeros(backend, T,Int64(qp.neqs),Int64(sem.mesh.ngl),Int64(sem.mesh.ngr))
        b_lag = KernelAbstractions.zeros(Int64(sem.mesh.nelem), Int64(sem.mesh.ngl), Int64(sem.mesh.ngr), Int64(qp.neqs))
        B_lag = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
    end
    if (inputs[:llaguerre_1d_right]||inputs[:llaguerre_1d_left])
        uaux_el_lag      = KernelAbstractions.zeros(backend, T, sem.mesh.nelem_semi_inf, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_el_lag       = KernelAbstractions.zeros(backend, T, sem.mesh.nelem_semi_inf, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diff_el_lag  = KernelAbstractions.zeros(backend, T, sem.mesh.nelem_semi_inf, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffξ_el_lag = KernelAbstractions.zeros(backend, T, sem.mesh.nelem_semi_inf, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        rhs_diffη_el_lag = KernelAbstractions.zeros(backend, T, sem.mesh.nelem_semi_inf, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        F_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        G_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        H_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        S_lag            = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        RHS_lag          = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
        RHS_visc_lag     = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
        uprimitive_lag = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs)+1)
        q_t_lag = KernelAbstractions.zeros(backend, T,Int64(qp.neqs),Int64(sem.mesh.ngr),Int64(sem.mesh.ngl))
        q_ti_lag = KernelAbstractions.zeros(backend, T,Int64(sem.mesh.ngr),Int64(sem.mesh.ngl))
        fqf_lag = KernelAbstractions.zeros(backend, T,Int64(qp.neqs),Int64(sem.mesh.ngr),Int64(sem.mesh.ngl))
        b_lag = KernelAbstractions.zeros(Int64(sem.mesh.nelem), Int64(sem.mesh.ngr), Int64(sem.mesh.ngl), Int64(qp.neqs))
        B_lag = KernelAbstractions.zeros(backend, T, Int64(sem.mesh.npoin), Int64(qp.neqs))
    end

    #-----------------------------------------------------------------
    for i=1:qp.neqs
        idx = (i-1)*sem.mesh.npoin
        u[idx+1:i*sem.mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1[:,i] = @view(qp.qn[:,i])
        qp.qnm2[:,i] = @view(qp.qn[:,i])
        
    end
    
    deps  = KernelAbstractions.zeros(backend, T, 1,1)
    Δt = inputs[:Δt]
    tspan = (T(inputs[:tinit]), T(inputs[:tend]))    
    if (backend == CPU())
        visc_coeff = inputs[:μ]#(inputs[:νρ], inputs[:νx], inputs[:νy], inputs[:κ], inputs[:κ], inputs[:κ], inputs[:κ])
    else
        coeffs = zeros(TFloat,qp.neqs)
        coeffs .= inputs[:μ]
        visc_coeff = KernelAbstractions.allocate(backend,TFloat,qp.neqs)
        KernelAbstractions.copyto!(backend,visc_coeff,coeffs)
    end
    ivisc_equations = inputs[:ivisc_equations]   
    
    if ("Laguerre" in sem.mesh.bdy_edge_type || inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])
        
        params = (backend, T, F, G, H, S,
                  uaux, uaux_el, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el,rhs_diffζ_el,
                  uprimitive,
                  flux_gpu, source_gpu, qbdy_gpu,
                  q_t, q_ti, fqf, b, B,
                  q_t_lag, q_ti_lag, fqf_lag, b_lag, B_lag,
                  RHS, RHS_visc,
                  F_lag, G_lag, S_lag, 
                  uaux_el_lag,
                  rhs_el_lag,
                  rhs_diff_el_lag,
                  rhs_diffξ_el_lag, rhs_diffη_el_lag,
                  RHS_lag, RHS_visc_lag, uprimitive_lag, 
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD,
                  sem.SOL_VARS_TYPE,
                  neqs=qp.neqs,
		  basis=sem.basis[1], basis_lag = sem.basis[2], ω = sem.ω[1], ω_lag = sem.ω[2], sem.mesh, metrics = sem.metrics[1], metrics_lag = sem.metrics[2], 
                  inputs, visc_coeff, ivisc_equations,
                  sem.matrix.M, sem.matrix.Minv,tspan,
                  Δt, deps, xmax, xmin, ymax, ymin,
                  qp, sem.fx, sem.fy, fy_t, sem.fy_lag, fy_t_lag, laguerre=true)
    else
        params = (backend, T, F, G, H, S,
                  uaux, uaux_el, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el,rhs_diffζ_el,
                  uprimitive,
                  flux_gpu, source_gpu, qbdy_gpu,
                  q_t, q_ti, fqf, b, B,
                  RHS, RHS_visc,
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD, 
                  sem.SOL_VARS_TYPE, 
                  neqs=qp.neqs,
                  sem.basis, sem.ω, sem.mesh, sem.metrics,
                  inputs, visc_coeff, ivisc_equations,
                  sem.matrix.M, sem.matrix.Minv,tspan,
                  Δt, deps, xmax, xmin, ymax, ymin,
                  qp, sem.fx, sem.fy, fy_t,laguerre=false)
    end 


    println(" # Build arrays and params ................................ DONE")

    return params, u
    
end

