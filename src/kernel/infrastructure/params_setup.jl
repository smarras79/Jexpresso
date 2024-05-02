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
    u        = zeros(T, sem.mesh.npoin*qp.neqs)
    uaux     = zeros(T, sem.mesh.npoin, qp.neqs)
    rhs      = allocate_rhs(sem.mesh.SD, sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T; neqs=qp.neqs)
    fluxes   = allocate_fluxes(sem.mesh.SD, sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T; neqs=qp.neqs)
    vaux     = zeros(T, sem.mesh.npoin) #generic auxiliary array for general use
            
    #The following are currently used by B.C.
    gradu    = zeros(T, 2, 1, 1)
    ubdy     = zeros(qp.neqs)
    bdy_flux = zeros(qp.neqs,1)    

    #
    # filter arrays
    #
    filter =  allocate_filter(sem.mesh.SD,
                              sem.mesh.nelem,
                              sem.mesh.npoin,
                              sem.mesh.ngl,
                              TFloat;
                              neqs=qp.neqs,
                              lfilter=inputs[:lfilter])
    
    fy_t     = transpose(sem.fy)
     
    #q_t = filter.q_t
    #fqf = filter.fqf
    #q_ti = filter.q_ti
    #b = filter.b
    #B = filter.B
    
    #store grid limits to save time
    xmax = maximum(sem.mesh.x)
    xmin = minimum(sem.mesh.x)
    ymax = maximum(sem.mesh.y)
    ymin = minimum(sem.mesh.y)

    
    #The following are only built and active if Laguerre boundaries are to be used
    if ( "Laguerre" in sem.mesh.bdy_edge_type ||
        inputs[:llaguerre_1d_right] == true   ||
        inputs[:llaguerre_1d_left]  == true )
        
        rhs_lag = allocate_rhs_lag(sem.mesh.SD,
                                   sem.mesh.nelem_semi_inf,
                                   sem.mesh.npoin,
                                   sem.mesh.ngl,
                                   sem.mesh.ngr,
                                   TFloat; neqs=qp.neqs)

        fluxes_lag = allocate_fluxes_lag(sem.mesh.SD,
                                         sem.mesh.ngl,
                                         sem.mesh.ngr,
                                         TFloat; neqs=qp.neqs)

        filter_lag =  allocate_filter_lag(sem.mesh.SD,
                                          sem.mesh.nelem_semi_inf,
                                          sem.mesh.npoin,
                                          sem.mesh.ngl,
                                          sem.mesh.ngr,
                                          TFloat;
                                          neqs=qp.neqs,
                                          lfilter=inputs[:lfilter])
        fy_t_lag = transpose(sem.fy_lag)
        
        #q_t_lag  =  filter_lag.q_t_lag
        #q_ti_lag =  filter_lag.q_ti_lag
        #fqf_lag  =  filter_lag.fqf_lag
        #b_lag    =  filter_lag.b_lag
        #B_lag    =  filter_lag.B_lag

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
        
        params = (T, rhs, fluxes,
                  rhs_lag, fluxes_lag,
                  filter, filter_lag,
                  uaux, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                 
                  #q_t, q_ti, fqf, b, B,
                  #q_t_lag, q_ti_lag, fqf_lag, b_lag, B_lag,
                  
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD,
                  sem.SOL_VARS_TYPE,
                  neqs=qp.neqs,
		  basis=sem.basis[1], basis_lag = sem.basis[2], ω = sem.ω[1], ω_lag = sem.ω[2], sem.mesh, metrics = sem.metrics[1], metrics_lag = sem.metrics[2], 
                  inputs, visc_coeff, ivisc_equations,
                  sem.matrix.M, sem.matrix.Minv,tspan,
                  Δt, deps, xmax, xmin, ymax, ymin,
                  qp, sem.fx, sem.fy, fy_t, sem.fy_lag, fy_t_lag, laguerre=true)
    else
        params = (rhs, fluxes,
                  filter,
                  T,
                  uaux, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  #q_t, q_ti, fqf, b, B,
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

