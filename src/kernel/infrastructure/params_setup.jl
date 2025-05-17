function params_setup(sem,
                      qp::St_SolutionVars,
                      inputs::Dict,
                      OUTPUT_DIR::String,
                      T,
                      tspan = [T(inputs[:tinit]), T(inputs[:tend])])

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    println_rank(" # Build arrays and params ................................ "; msg_rank = rank, suppress = sem.mesh.msg_suppress)
    if rank == 0
        @info " " inputs[:ode_solver] inputs[:tinit] inputs[:tend] inputs[:Δt]
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

    μdsgs  = allocate_musgs(sem.mesh.SD,
                              sem.mesh.nelem,
                              sem.mesh.npoin,
                              sem.mesh.ngl,
                              inputs[:visc_model],
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

    u            = uODE.u
    uaux         = uODE.uaux
    vaux         = uODE.vaux
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
    
    #------------------------------------------------------------------------------------
    # boundary flux arrays
    #------------------------------------------------------------------------------------
    bdy_fluxes = allocate_bdy_fluxes(sem.mesh.SD,
                          sem.mesh.nfaces_bdy,
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
    xmin = MPI.Allreduce(minimum(sem.mesh.x), MPI.MIN, comm)
    xmax = MPI.Allreduce(maximum(sem.mesh.x), MPI.MAX, comm)
    ymin = MPI.Allreduce(minimum(sem.mesh.y), MPI.MIN, comm)
    ymax = MPI.Allreduce(maximum(sem.mesh.y), MPI.MAX, comm)
    zmin = MPI.Allreduce(minimum(sem.mesh.z), MPI.MIN, comm)
    zmax = MPI.Allreduce(maximum(sem.mesh.z), MPI.MAX, comm)

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
    mp = allocate_SamMicrophysics(sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T, backend; lmoist=inputs[:lmoist])
    #------------------------------------------------------------------------------------
    # Allocate large scale tendencies arrays
    #------------------------------------------------------------------------------------
    LST = allocate_LargeScaleTendencies(sem.mesh.npoin, sem.mesh, inputs, T, backend; lLST=inputs[:LST]) 
    #------------------------------------------------------------------------------------
    # Allocate Thermodynamic params for bomex case
    #------------------------------------------------------------------------------------
    thermo_params = create_updated_TD_Parameters(TFloat(101325.0))
        

    #------------------------------------------------------------------------------------
    # Populate solution arrays
    #------------------------------------------------------------------------------------
    for i=1:qp.neqs
        idx = (i-1)*sem.mesh.npoin
        u[idx+1:i*sem.mesh.npoin] = @view qp.qn[:,i]
        qp.qnm1[:,i] = @view(qp.qn[:,i])
        qp.qnm2[:,i] = @view(qp.qn[:,i])
        
    end
    
    deps  = KernelAbstractions.zeros(backend, T, 1,1)
    Δt    = inputs[:Δt]
    # tspan = [T(inputs[:tinit]), T(inputs[:tend])]
    if (backend == CPU())
        visc_coeff = inputs[:μ]
    else
        coeffs     = zeros(TFloat,qp.neqs)
        coeffs    .= inputs[:μ]
        visc_coeff = KernelAbstractions.allocate(backend,TFloat,qp.neqs)
        KernelAbstractions.copyto!(backend,visc_coeff,coeffs)
    end

    #------------------------------------------------------------------------------------
    # Populate params tuple to carry global arrays and constants around
    #------------------------------------------------------------------------------------
    if (sem.mesh.lLaguerre ||
        inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])
        pM = setup_assembler(sem.mesh.SD, RHS, sem.mesh.ip2gip, sem.mesh.gip2owner)
        params = (backend, T, F, G, H, S,
                  uaux, vaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el,rhs_diffζ_el,
                  uprimitive,
                  flux_gpu, source_gpu, qbdy_gpu,
                  q_t, q_ti, q_tij, fqf, b, B,
                  q_t_lag, q_ti_lag, fqf_lag, b_lag, B_lag,
                  flux_lag_gpu, source_lag_gpu,
                  qbdy_lag_gpu,
                  RHS, RHS_visc,
                  F_lag, G_lag, S_lag,
                  F_surf, S_face, S_flux, M_surf_inv = sem.matrix.M_surf_inv,
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
                  μdsgs,
                  sem.matrix.M, sem.matrix.Minv, pM=pM, tspan,
                  Δt, deps, xmax, xmin, ymax, ymin, zmin, zmax,
                  qp, mp, sem.fx, sem.fy, fy_t, sem.fy_lag, fy_t_lag, sem.fz, fz_t, laguerre=true)
        
    else
        pM = setup_assembler(sem.mesh.SD, RHS, sem.mesh.ip2gip, sem.mesh.gip2owner)
        params = (backend,
                  T, inputs,
                  uaux, vaux,
                  ubdy, gradu, bdy_flux,                   
                  RHS, RHS_visc,
                  fijk, ∇f_el,
                  rhs_el, rhs_diff_el,
                  rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el,
                  uprimitive,
                  F, G, H, S,
                  F_surf, S_face, S_flux, M_surf_inv = sem.matrix.M_surf_inv,
                  flux_gpu, source_gpu, qbdy_gpu,
                  flux_micro, source_micro, adjusted, Pm,
                  q_t, q_ti, q_tij, fqf, b, B,
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.PT, sem.AD, 
                  sem.SOL_VARS_TYPE, 
                  neqs=qp.neqs,
                  sem.connijk_original, sem.poin_in_bdy_face_original, sem.x_original, sem.y_original, sem.z_original,
                  sem.basis, sem.ω, sem.mesh, sem.metrics,
                  thermo_params, VT = inputs[:visc_model], visc_coeff,
                  μdsgs,
                  sem.matrix.M, sem.matrix.Minv, pM=pM,
                  tspan, Δt, xmax, xmin, ymax, ymin, zmin, zmax,
                  phys_grid = sem.phys_grid,
                  qp, mp, LST, sem.fx, sem.fy, fy_t, sem.fz, fz_t, laguerre=false,
                  OUTPUT_DIR,
                  sem.interp, sem.project, sem.partitioned_model, sem.nparts, sem.distribute)
    end

    println_rank(" # Build arrays and params ................................ DONE"; msg_rank = rank, suppress = sem.mesh.msg_suppress)

    return params, u
    
end
