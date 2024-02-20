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
    
    #-----------------------------------------------------------------
    # Initialize:
    # u     -> solution array
    # uaux  -> solution auxiliary array (Eventually to be removed)
    # F,G   -> physical flux arrays
    # S     -> source array
    # rhs* -> inviscid and viscous ELEMENT rhs
    # RHS* -> inviscid and viscous GLOBAL  rhs
    #-----------------------------------------------------------------
    u            = zeros(T, sem.mesh.npoin*qp.neqs)
    uaux         = zeros(T, sem.mesh.npoin, qp.neqs)

    if sem.mesh.nsd == 1
        uaux_el      = zeros(T, sem.mesh.nelem, sem.mesh.ngl, qp.neqs)
        rhs_el       = zeros(T, sem.mesh.nelem, sem.mesh.ngl, qp.neqs)
        rhs_diff_el  = zeros(T, sem.mesh.nelem, sem.mesh.ngl, qp.neqs)
        rhs_diffξ_el = zeros(T, sem.mesh.nelem, sem.mesh.ngl, qp.neqs)
        rhs_diffη_el = zeros(T, sem.mesh.nelem, sem.mesh.ngl, qp.neqs)
        rhs_diffζ_el = zeros(T, 0)
        F            = zeros(T, sem.mesh.ngl, qp.neqs)
        G            = zeros(T, sem.mesh.ngl, qp.neqs)
        H            = zeros(T, 0)
        S            = zeros(T, sem.mesh.ngl, qp.neqs)
        uprimitive   = zeros(T, sem.mesh.ngl, qp.neqs+1)
    elseif  sem.mesh.nsd == 2
        uaux_el      = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_el       = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diff_el  = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diffξ_el = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diffη_el = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diffζ_el = zeros(T, 0)
        F            = zeros(T, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        G            = zeros(T, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        H            = zeros(T, 0)
        S            = zeros(T, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        uprimitive   = zeros(T, sem.mesh.ngl, sem.mesh.ngl, qp.neqs+1)
    elseif  sem.mesh.nsd == 3
        uaux_el      = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_el       = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diff_el  = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diffξ_el = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diffη_el = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        rhs_diffζ_el = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        F            = zeros(T, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        G            = zeros(T, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        H            = zeros(T, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        S            = zeros(T, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
        uprimitive   = zeros(T, sem.mesh.ngl, sem.mesh.ngl, sem.mesh.ngl, qp.neqs+1)
    end
    
    RHS      = zeros(T, sem.mesh.npoin, qp.neqs)
    RHS_visc = zeros(T, sem.mesh.npoin, qp.neqs)
    vaux     = zeros(T, sem.mesh.npoin) #generic auxiliary array for general use
    
    #The following are currently used by B.C.
    gradu    = zeros(T, 2, 1, 1) #zeros(2,sem.mesh.npoin,nvars)
    ubdy     = zeros(qp.neqs)
    bdy_flux = zeros(qp.neqs,1)    
   
    #filter arrays
    q_t  = zeros(Float64,qp.neqs,sem.mesh.ngl,sem.mesh.ngl)
    q_ti = zeros(Float64,sem.mesh.ngl,sem.mesh.ngl)
    fy_t = transpose(sem.fy)
    fy_t_lag = transpose(sem.fy_lag)
    fqf = zeros(Float64,qp.neqs,sem.mesh.ngl,sem.mesh.ngl)
    b = zeros(sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngl, qp.neqs)
    B = zeros(Float64, sem.mesh.npoin, qp.neqs) 
    #store grid limits to save time
    xmax = maximum(sem.mesh.x)
    xmin = minimum(sem.mesh.x)
    ymax = maximum(sem.mesh.y)
    ymin = minimum(sem.mesh.y)

    
    #The following are only built and active if Laguerre boundaries are to be used
    if ( "Laguerre" in sem.mesh.bdy_edge_type)
        uaux_el_lag      = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        rhs_el_lag       = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        rhs_diff_el_lag  = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        rhs_diffξ_el_lag = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        rhs_diffη_el_lag = zeros(T, sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        F_lag            = zeros(T, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        G_lag            = zeros(T, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        H_lag            = zeros(T, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        S_lag            = zeros(T, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        RHS_lag          = zeros(T, sem.mesh.npoin, qp.neqs)
        RHS_visc_lag     = zeros(T, sem.mesh.npoin, qp.neqs)
        uprimitive_lag = zeros(T, sem.mesh.ngl, sem.mesh.ngr, qp.neqs+1)
        q_t_lag = zeros(Float64,qp.neqs,sem.mesh.ngl,sem.mesh.ngr)
        q_ti_lag = zeros(Float64,sem.mesh.ngl,sem.mesh.ngr)
        fqf_lag = zeros(Float64,qp.neqs,sem.mesh.ngl,sem.mesh.ngr)
        b_lag = zeros(sem.mesh.nelem, sem.mesh.ngl, sem.mesh.ngr, qp.neqs)
        B_lag = zeros(Float64, sem.mesh.npoin, qp.neqs)
    end
    if (inputs[:llaguerre_1d_right]||inputs[:llaguerre_1d_left])
        uaux_el_lag      = zeros(T, sem.mesh.nelem_semi_inf, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        rhs_el_lag       = zeros(T, sem.mesh.nelem_semi_inf, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        rhs_diff_el_lag  = zeros(T, sem.mesh.nelem_semi_inf, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        rhs_diffξ_el_lag = zeros(T, sem.mesh.nelem_semi_inf, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        rhs_diffη_el_lag = zeros(T, sem.mesh.nelem_semi_inf, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        F_lag            = zeros(T, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        G_lag            = zeros(T, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        H_lag            = zeros(T, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        S_lag            = zeros(T, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        RHS_lag          = zeros(T, sem.mesh.npoin, qp.neqs)
        RHS_visc_lag     = zeros(T, sem.mesh.npoin, qp.neqs)
        uprimitive_lag = zeros(T, sem.mesh.ngr, sem.mesh.ngl, qp.neqs+1)
        q_t_lag = zeros(Float64,qp.neqs,sem.mesh.ngr,sem.mesh.ngl)
        q_ti_lag = zeros(Float64,sem.mesh.ngr,sem.mesh.ngl)
        fqf_lag = zeros(Float64,qp.neqs,sem.mesh.ngr,sem.mesh.ngl)
        b_lag = zeros(sem.mesh.nelem, sem.mesh.ngr, sem.mesh.ngl, qp.neqs)
        B_lag = zeros(Float64, sem.mesh.npoin, qp.neqs)
    end

    #-----------------------------------------------------------------
    for i=1:qp.neqs
        idx = (i-1)*sem.mesh.npoin
        u[idx+1:i*sem.mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1[:,i] = @view(qp.qn[:,i])
        qp.qnm2[:,i] = @view(qp.qn[:,i])
        
    end
    
    deps  = zeros(1,1)
    Δt = inputs[:Δt]
    tspan = (inputs[:tinit], inputs[:tend])    
    visc_coeff = inputs[:μ]#(inputs[:νρ], inputs[:νx], inputs[:νy], inputs[:κ], inputs[:κ], inputs[:κ], inputs[:κ])
    ivisc_equations = inputs[:ivisc_equations]   
    
    if ("Laguerre" in sem.mesh.bdy_edge_type || inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])
        
        params = (T, F, G, H, S,
                  uaux, uaux_el, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el,rhs_diffζ_el,
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
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD,
                  sem.SOL_VARS_TYPE,
                  neqs=qp.neqs,
		  basis=sem.basis[1], basis_lag = sem.basis[2], ω = sem.ω[1], ω_lag = sem.ω[2], sem.mesh, metrics = sem.metrics[1], metrics_lag = sem.metrics[2], 
                  inputs, visc_coeff, ivisc_equations,
                  sem.matrix.M, sem.matrix.Minv,tspan,
                  Δt, deps, xmax, xmin, ymax, ymin,
                  qp, sem.fx, sem.fy, fy_t, sem.fy_lag, fy_t_lag, laguerre=true)
    else
        params = (T, F, G, H, S,
                  uaux, uaux_el, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el,rhs_diffζ_el,
                  uprimitive,
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

