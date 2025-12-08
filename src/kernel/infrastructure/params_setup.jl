function params_setup(sem,
                      qp::St_SolutionVars,
                      inputs::Dict,
                      OUTPUT_DIR::String,
                      T,
                      tspan = [T(inputs[:tinit]), T(inputs[:tend])])

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Build arrays and params ................................ "; msg_rank = rank, suppress = sem.mesh.msg_suppress)
    if rank == 0 && tspan[1] == T(inputs[:tinit])
        @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
    elseif rank == 0 && tspan[1] != T(inputs[:tinit])
        @info " " tspan[1] tspan[end]
    end

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

    fijk   = allocate_fijk(sem.mesh.SD,
                           sem.mesh.ngl,
                           T, backend;
                           neqs=qp.neqs)

    ∇f     = allocate_∇f(sem.mesh.SD,
                         sem.mesh.nelem,
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

    gpuMoist = allocate_gpuMoist(sem.mesh.SD,
                             sem.mesh.npoin,
                             sem.mesh.nelem,
                             sem.mesh.ngl,
                             T, backend, inputs[:lmoist];
                             neqs=qp.neqs)
    
    ncf_arrays = allocate_ncfArrays(sem.mesh.SD,
                                    sem.mesh.num_ncf_pg,
                                    sem.mesh.num_ncf_cg,
                                    sem.mesh.ngl,
                                    T, backend;
                                    neqs=qp.neqs)

    if inputs[:lvisc] == true && inputs[:visc_model] != AV()
        viscsgs = allocate_visc(sem.mesh.SD, sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T, backend; neqs=qp.neqs)
    else
        viscsgs = allocate_visc(sem.mesh.SD, 1, 1, 1, T, backend; neqs=1)
    end
    
    u            = uODE.u
    uaux         = uODE.uaux
    vaux         = uODE.vaux
    utmp         = uODE.utmp
    F            = fluxes.F
    G            = fluxes.G
    H            = fluxes.H
    S            = fluxes.S
    fijk         = fijk.fijk
    ∇f_el        = ∇f.∇f_el
    RHS          = rhs.RHS
    RHS_visc     = rhs.RHS_visc
    rhs_el       = rhs.rhs_el
    rhs_diff_el  = rhs.rhs_diff_el
    rhs_diffξ_el = rhs.rhs_diffξ_el
    rhs_diffη_el = rhs.rhs_diffη_el
    rhs_diffζ_el = rhs.rhs_diffζ_el
    μsgs         = viscsgs.μ
    
    rhs_el_tmp   = rhs.rhs_el_tmp
    
    #------------------------------------------------------------------------------------
    # non conforming faces arrays and mpi cache
    #------------------------------------------------------------------------------------
    q_el      = ncf_arrays.q_el
    q_el_pro  = ncf_arrays.q_el_pro 
    q_ghost_p = ncf_arrays.q_ghost_p 
    q_ghost_c = ncf_arrays.q_ghost_c 
    # mpi cache for ncf
    cache_ghost_p = SendReceiveCache(comm, @view(q_ghost_p[:, 1]), sem.mesh.pgip_owner)
    cache_ghost_c = SendReceiveCache(comm, q_ghost_c, sem.mesh.cgip_owner)
    #------------------------------------------------------------------------------------
    # boundary flux arrays
    #------------------------------------------------------------------------------------
    bdy_fluxes = allocate_bdy_fluxes(sem.mesh.SD,
                          sem.mesh.nfaces_bdy,
                          sem.mesh.nedges_bdy,
                          sem.mesh.npoin,
                          sem.mesh.ngl,
                          T, backend;
                          neqs=qp.neqs)
    
    F_surf = bdy_fluxes.F_surf
    S_face = bdy_fluxes.S_face
    S_flux = bdy_fluxes.S_flux
    #------------------------------------------------------------------------------------
    # GPU arrays
    #------------------------------------------------------------------------------------
    flux_gpu       = gpuAux.flux_gpu
    source_gpu     = gpuAux.source_gpu
    qbdy_gpu       = gpuAux.qbdy_gpu
    uprimitive     = fluxes.uprimitive
    flux_micro     = gpuMoist.flux_micro
    source_micro   = gpuMoist.source_micro
    adjusted       = gpuMoist.adjusted
    Pm             = gpuMoist.Pm
    #------------------------------------------------------------------------------------
    # filter arrays
    #------------------------------------------------------------------------------------
    filter = allocate_filter(sem.mesh.SD, sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T, backend; neqs=qp.neqs, lfilter=inputs[:lfilter])
    fy_t   = transpose(sem.fy)
    fz_t   = transpose(sem.fz)
    q_t    = filter.q_t
    q_ti   = filter.q_ti
    q_tij  = filter.q_tij
    fqf    = filter.fqf
    b      = filter.b
    B      = filter.B

    #------------------------------------------------------------------------------------  
    # B.C. arrays
    #------------------------------------------------------------------------------------
    gradu    = KernelAbstractions.zeros(backend, T, 2, 1, 1)
    ubdy     = KernelAbstractions.zeros(backend, T, Int64(qp.neqs))
    bdy_flux = KernelAbstractions.zeros(backend, T, Int64(qp.neqs),1)    

    #------------------------------------------------------------------------------------
    # Some domain parameters
    #------------------------------------------------------------------------------------
    xmax = sem.mesh.xmax; xmin = sem.mesh.xmin
    ymax = sem.mesh.ymax; ymin = sem.mesh.ymin
    zmax = sem.mesh.zmax; zmin = sem.mesh.zmin
        
    #------------------------------------------------------------------------------------
    # Laguerre arrays
    #------------------------------------------------------------------------------------
    if ( sem.mesh.lLaguerre ||
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
        
        flux_lag_gpu     = gpuAux_lag.flux_lag_gpu
        source_lag_gpu   = gpuAux_lag.source_lag_gpu
        qbdy_lag_gpu     = gpuAux_lag.qbdy_lag_gpu
        
        fy_t_lag         = transpose(sem.fy_lag)
        q_t_lag          = filter_lag.q_t_lag
        q_ti_lag         = filter_lag.q_ti_lag
        fqf_lag          = filter_lag.fqf_lag
        b_lag            = filter_lag.b_lag
        B_lag            = filter_lag.B_lag
    end
    
    #------------------------------------------------------------------------------------
    # Allocate micophysics arrays
    #------------------------------------------------------------------------------------
    mp = allocate_SamMicrophysics(sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T, backend, sem.mesh.SD; lmoist=inputs[:lmoist])
    #------------------------------------------------------------------------------------
    # Allocate large scale tendencies arrays
    #------------------------------------------------------------------------------------
    LST = allocate_LargeScaleTendencies(sem.mesh.npoin, sem.mesh, inputs, T, backend; lLST=inputs[:LST])
    #------------------------------------------------------------------------------------
    # Allocate wall model arrays
    #------------------------------------------------------------------------------------
    WM = allocate_Wall_model(sem.mesh.nfaces_bdy, sem.mesh.ngl, T, backend; lwall_model=inputs[:lwall_model])
    #------------------------------------------------------------------------------------
    # Allocate Thermodynamic params for bomex case
    #------------------------------------------------------------------------------------
    PhysConst = PhysicalConst{TFloat}()
    thermo_params = create_updated_TD_Parameters(PhysConst.potential_temperature_reference_pressure)
    
    #------------------------------------------------------------------------------------
    # Populate solution arrays
    #------------------------------------------------------------------------------------
    if (sem.mesh.SD != NSD_1D()) && !(sem.mesh.lLaguerre)
        if rank == 0
            @info "start conformity4ncf_q!"
        end
        g_dss_cache_qp = setup_assembler(sem.mesh.SD, qp.qn, sem.mesh.ip2gip, sem.mesh.gip2owner)
        conformity4ncf_q!(qp.qn, rhs_el_tmp, @view(utmp[:,:]), vaux, 
                          g_dss_cache_qp,
                          sem.mesh.SD, 
                          sem.QT, sem.mesh.connijk,
                          sem.mesh, sem.matrix.Minv, 
                          sem.metrics.Je, sem.ω, sem.AD, 
                          qp.neqs,
                          q_el, q_el_pro,
                          cache_ghost_p, q_ghost_p,
                          cache_ghost_c, q_ghost_c,
                          sem.interp; ladapt = inputs[:ladapt])
        conformity4ncf_q!(qp.qe, rhs_el_tmp, @view(utmp[:,:]), vaux, 
                          g_dss_cache_qp,
                          sem.mesh.SD, 
                          sem.QT, sem.mesh.connijk, 
                          sem.mesh, sem.matrix.Minv, 
                          sem.metrics.Je, sem.ω, sem.AD, 
                          qp.neqs,
                          q_el, q_el_pro,
                          cache_ghost_p, q_ghost_p,
                          cache_ghost_c, q_ghost_c,
                          sem.interp; ladapt = inputs[:ladapt])
        MPI.Barrier(comm)
        if rank == 0
            @info "end conformity4ncf_q!"
        end
    end
    for i=1:qp.neqs
        idx = (i-1)*sem.mesh.npoin
        u[idx+1:i*sem.mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1[:,i] = @view(qp.qn[:,i])
        qp.qnm2[:,i] = @view(qp.qn[:,i])
    end
    
    deps  = KernelAbstractions.zeros(backend, T, 1,1)
    Δt    = inputs[:Δt]
    #if (backend == CPU())
    #    visc_coeff = zeros(TFloat, qp.neqs)
    #    if inputs[:lvisc]
    #        visc_coeff .= inputs[:μ]
    #    end
    #else
   
    if inputs[:lvisc]
        coeffs = zeros(TFloat, qp.neqs)
        if size(inputs[:μ]) > size(coeffs)
            coeffs .= inputs[:μ][1:qp.neqs]
        elseif size(inputs[:μ]) <= size(coeffs)
            coeffs .= inputs[:μ]
        end
        visc_coeff = KernelAbstractions.allocate(backend,TFloat, qp.neqs)
        KernelAbstractions.copyto!(backend,visc_coeff,coeffs)
    else
        visc_coeff = KernelAbstractions.allocate(backend, TFloat, 1)
        visc_coeff = [0.0]
    end

    # setup timer
    timers = Dict{String, MPIFunctionTimer}()
    #------------------------------------------------------------------------------------
    # Populate params tuple to carry global arrays and constants around
    #------------------------------------------------------------------------------------
    if (sem.mesh.lLaguerre ||
        inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])
        g_dss_cache = setup_assembler(sem.mesh.SD, RHS, sem.mesh.ip2gip, sem.mesh.gip2owner)
        params = (backend, T, F, G, H, S,
                  uaux, vaux, utmp,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el, rhs_el_tmp,
                  rhs_diffξ_el, rhs_diffη_el,rhs_diffζ_el,
                  uprimitive,
                  flux_gpu, source_gpu, qbdy_gpu,
                  q_t, q_ti, q_tij, fqf, b, B,
                  q_t_lag, q_ti_lag, fqf_lag, b_lag, B_lag,
                  flux_lag_gpu, source_lag_gpu,
                  qbdy_lag_gpu,
                  RHS, RHS_visc,
                  F_lag, G_lag, S_lag, 
                  F_surf, S_face, S_flux, M_surf_inv = sem.matrix.M_surf_inv, M_edge_inv = sem.matrix.M_edge_inv,
                  rhs_el_lag,
                  rhs_diff_el_lag,
                  rhs_diffξ_el_lag, rhs_diffη_el_lag,
                  RHS_lag, RHS_visc_lag, uprimitive_lag, 
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD,
                  sem.SOL_VARS_TYPE,
                  neqs=qp.neqs,
                  sem.mesh,
                  sem.connijk_original, sem.poin_in_bdy_face_original, sem.x_original, sem.y_original, sem.z_original,
		  basis=sem.basis[1], basis_lag = sem.basis[2],
                  ω = sem.ω[1], ω_lag = sem.ω[2],
                  metrics = sem.metrics[1], metrics_lag = sem.metrics[2], 
                  inputs, VT = inputs[:visc_model], visc_coeff,
                  WM,
                  sem.matrix.M, sem.matrix.Minv, g_dss_cache=g_dss_cache, tspan,
                  Δt, deps, xmax, xmin, ymax, ymin, zmin, zmax,
                  qp, mp, sem.fx, sem.fy, fy_t, sem.fy_lag, fy_t_lag, sem.fz, fz_t, laguerre=true,
                  timers)
        
    else
        g_dss_cache = setup_assembler(sem.mesh.SD, RHS, sem.mesh.ip2gip, sem.mesh.gip2owner)
        params = (backend,
                  T, inputs,
                  uaux, vaux, utmp,
                  ubdy, gradu, bdy_flux,                   
                  RHS, RHS_visc,
                  fijk, ∇f_el,
                  rhs_el, rhs_diff_el, rhs_el_tmp,
                  rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                  uprimitive,
                  F, G, H, S,
                  F_surf, S_face, S_flux, M_surf_inv = sem.matrix.M_surf_inv, M_edge_inv = sem.matrix.M_edge_inv,
                  flux_gpu, source_gpu, qbdy_gpu,
                  flux_micro, source_micro, adjusted, Pm,
                  q_el, q_el_pro, q_ghost_p, q_ghost_c,
                  cache_ghost_p, cache_ghost_c,
                  q_t, q_ti, q_tij, fqf, b, B,
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD, 
                  sem.SOL_VARS_TYPE, 
                  neqs=qp.neqs,
                  sem.connijk_original, sem.poin_in_bdy_face_original, sem.x_original, sem.y_original, sem.z_original,
                  sem.basis, sem.ω, sem.mesh, sem.metrics,
                  thermo_params, VT = inputs[:visc_model], visc_coeff,
                  sem.matrix.M, sem.matrix.Minv, g_dss_cache=g_dss_cache,
                  tspan, Δt, xmax, xmin, ymax, ymin, zmin, zmax,
                  WM,
                  phys_grid = sem.phys_grid,
                  qp, mp, LST, sem.fx, sem.fy, fy_t, sem.fz, fz_t, laguerre=false,
                  OUTPUT_DIR,
                  timers,
                  sem.interp, sem.project, sem.nparts, sem.distribute)
    end

    println_rank(" # Build arrays and params ................................ DONE"; msg_rank = rank, suppress = sem.mesh.msg_suppress)

    return params, u
    
end
