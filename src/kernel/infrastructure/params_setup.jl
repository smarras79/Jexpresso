function params_setup(sem,
                      qp::St_SolutionVars,
                      inputs::Dict,
                      OUTPUT_DIR::String,
                      T)

    
    println(" # Build arrays and params ................................ ")
    @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]

    backend = inputs[:backend]
    
    uODE = allocate_uODE(sem.mesh.SD,
                         sem.mesh.npoin,
                         T, backend;
                         neqs=qp.neqs)
    
    rhs    = allocate_rhs(sem.mesh.SD,
                          sem.mesh.nelem,
                          sem.mesh.npoin,
                          sem.mesh.ngl,
                          T, backend;
                          neqs=qp.neqs)
    
    fluxes = allocate_fluxes(sem.mesh.SD,
                             sem.mesh.npoin,
                             sem.mesh.ngl,
                             T, backend;
                             neqs=qp.neqs)
    
    gpuAux = allocate_gpuAux(sem.mesh.SD,
                             sem.mesh.nelem,
                             sem.mesh.nedges_bdy,
                             sem.mesh.nfaces_bdy,
                             sem.mesh.ngl,
                             T, backend;
                             neqs=qp.neqs)

    u            = uODE.u
    uaux         = uODE.uaux
    vaux         = uODE.vaux
    F            = fluxes.F
    G            = fluxes.G
    H            = fluxes.H
    S            = fluxes.S
    RHS          = rhs.RHS
    RHS_visc     = rhs.RHS_visc
    rhs_el       = rhs.rhs_el
    rhs_diff_el  = rhs.rhs_diff_el
    rhs_diffξ_el = rhs.rhs_diffξ_el
    rhs_diffη_el = rhs.rhs_diffη_el
    rhs_diffζ_el = rhs.rhs_diffζ_el

    #
    # GPU arrays
    #
    flux_gpu       = gpuAux.flux_gpu
    source_gpu     = gpuAux.source_gpu
    qbdy_gpu       = gpuAux.qbdy_gpu
    uprimitive     = fluxes.uprimitive
    
    #
    # filter arrays
    #
    filter = allocate_filter(sem.mesh.SD, sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T, backend; neqs=qp.neqs, lfilter=inputs[:lfilter])
    fy_t   = transpose(sem.fy)
    q_t    = filter.q_t
    q_ti   = filter.q_ti
    fqf    = filter.fqf
    b      = filter.b
    B      = filter.B

    #    
    # B.C. arrays
    #
    gradu    = KernelAbstractions.zeros(backend, T, 2, 1, 1) #KernelAbstractions.zeros(2,Int64(sem.mesh.npoin),nvars)
    ubdy     = KernelAbstractions.zeros(backend, T, Int64(qp.neqs))
    bdy_flux = KernelAbstractions.zeros(backend, T, Int64(qp.neqs),1)    
    
    xmax = maximum(sem.mesh.x); xmin = minimum(sem.mesh.x)
    ymax = maximum(sem.mesh.y); ymin = minimum(sem.mesh.y)
    zmax = maximum(sem.mesh.z); zmin = minimum(sem.mesh.z)
    
   
    #
    # Laguerre arrays
    #
    if ( "Laguerre" in sem.mesh.bdy_edge_type ||
        inputs[:llaguerre_1d_right] == true   ||
        inputs[:llaguerre_1d_left]  == true )
        
        rhs_lag = allocate_rhs_lag(sem.mesh.SD,
                                   sem.mesh.nelem_semi_inf,
                                   sem.mesh.npoin,
                                   sem.mesh.ngl,
                                   sem.mesh.ngr,
                                   T,
                                   backend;
                                   neqs = qp.neqs)


        fluxes_lag = allocate_fluxes_lag(sem.mesh.SD,
                                         sem.mesh.ngl,
                                         sem.mesh.ngr,
                                         T,
                                         backend;
                                         neqs = qp.neqs)

        filter_lag =  allocate_filter_lag(sem.mesh.SD,
                                          sem.mesh.nelem_semi_inf,
                                          sem.mesh.npoin,
                                          sem.mesh.ngl,
                                          sem.mesh.ngr,
                                          T,
                                          backend;
                                          neqs = qp.neqs,
                                          lfilter = inputs[:lfilter])

        gpuAux_lag = allocate_gpuAux_lag(sem.mesh.SD,
                                         sem.mesh.nelem_semi_inf,
                                         sem.mesh.nedges_bdy,
                                         sem.mesh.nfaces_bdy,
                                         sem.mesh.ngl,
                                         sem.mesh.ngr,
                                         T, backend;
                                         neqs=qp.neqs)

        

        RHS_lag          = rhs_lag.RHS_lag
        RHS_visc_lag     = rhs_lag.RHS_visc_lag
        rhs_el_lag       = rhs_lag.rhs_el_lag
        rhs_diff_el_lag  = rhs_lag.rhs_diff_el_lag
        rhs_diffξ_el_lag = rhs_lag.rhs_diffξ_el_lag
        rhs_diffη_el_lag = rhs_lag.rhs_diffη_el_lag
        rhs_diffζ_el_lag = rhs_lag.rhs_diffζ_el_lag
        
        F_lag            = fluxes_lag.F_lag
        G_lag            = fluxes_lag.G_lag
        H_lag            = fluxes_lag.H_lag
        S_lag            = fluxes_lag.S_lag
        uprimitive_lag   = fluxes_lag.uprimitive_lag
        
        flux_lag_gpu   = gpuAux_lag.flux_lag_gpu
        source_lag_gpu = gpuAux_lag.source_lag_gpu
        qbdy_lag_gpu   = gpuAux_lag.qbdy_lag_gpu
        
        fy_t_lag = transpose(sem.fy_lag)
        q_t_lag  = filter_lag.q_t_lag
        q_ti_lag = filter_lag.q_ti_lag
        fqf_lag  = filter_lag.fqf_lag
        b_lag    = filter_lag.b_lag
        B_lag    = filter_lag.B_lag
        
    end
    
    
    for i=1:qp.neqs
        idx = (i-1)*sem.mesh.npoin
        u[idx+1:i*sem.mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1[:,i] = @view(qp.qn[:,i])
        qp.qnm2[:,i] = @view(qp.qn[:,i])
        
    end
    
    deps  = KernelAbstractions.zeros(backend, T, 1,1)
    Δt    = inputs[:Δt]
    tspan = [T(inputs[:tinit]), T(inputs[:tend])]
    if (backend == CPU())
        visc_coeff = inputs[:μ]
    else
        coeffs     = zeros(TFloat,qp.neqs)
        coeffs    .= inputs[:μ]
        visc_coeff = KernelAbstractions.allocate(backend,TFloat,qp.neqs)
        KernelAbstractions.copyto!(backend,visc_coeff,coeffs)
    end
    ivisc_equations = inputs[:ivisc_equations]   
    
    params_base = (backend,
                   T, inputs,
                   uaux, vaux,
                   ubdy, gradu, bdy_flux,                   
                   RHS, RHS_visc,
                   rhs_el, rhs_diff_el,
                   rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                   uprimitive,
                   F, G, H, S,
                   flux_gpu, source_gpu, qbdy_gpu,
                   q_t, q_ti, fqf, b, B,
                   SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD, 
                   sem.SOL_VARS_TYPE, 
                   neqs=qp.neqs,
                   sem.basis, sem.ω, sem.mesh, sem.metrics,
                   visc_coeff, ivisc_equations,
                   sem.matrix.M, sem.matrix.Minv,tspan,
                   Δt, deps, xmax, xmin, ymax, ymin,
                   qp, sem.fx, sem.fy, fy_t, laguerre=false)

    if ("Laguerre" in sem.mesh.bdy_edge_type ||
        inputs[:llaguerre_1d_right] ||
        inputs[:llaguerre_1d_left])
        
        params_lag = (flux_lag_gpu, source_lag_gpu,
                      qbdy_lag_gpu,
                      uprimitive_lag,
                      F_lag, G_lag, S_lag, 
                      rhs_el_lag,
                      rhs_diff_el_lag,
                      rhs_diffξ_el_lag, rhs_diffη_el_lag,
                      RHS_lag, RHS_visc_lag,
                      q_t_lag, q_ti_lag, fqf_lag,
                      b_lag, B_lag, 
                      sem.fy_lag, fy_t_lag, laguerre=true)
        
        params = (params_base..., params_lag)
        params_base = nothing
        params_lag  = nothing
        dump(params)
        @info params.tspan

        @mystop
    else
        params = params_base
        params_base = nothing
    end

    println(" # Build arrays and params ................................ DONE")

    return params, u
    
end

