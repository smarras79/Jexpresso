function params_setup(sem,
                      qp::St_SolutionVars,
                      inputs::Dict,
                      OUTPUT_DIR::String,
                      T,
                      tspan = [T(inputs[:tinit]), T(inputs[:tend])];
                      coupling = nothing)

    comm = get_mpi_comm()
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
    fluxaux      = uODE.fluxaux
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
    μ_max        = viscsgs.μ_max
    
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
    filter = allocate_filter(sem.mesh.SD, sem.mesh.nelem, sem.mesh.npoin, sem.mesh.ngl, T, backend; neqs=qp.neqs, lfilter=inputs[:lfilter], ladapt=inputs[:ladapt])
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
    WM = allocate_Wall_model(sem.mesh.nfaces_bdy, sem.mesh.ngl, T, backend; lwall_model=inputs[:lwall_model], lmoist=inputs[:lmoist])
    #------------------------------------------------------------------------------------
    # Allocate Thermodynamic params for bomex case
    #------------------------------------------------------------------------------------
    PhysConst = PhysicalConst{TFloat}()
    thermo_params = create_updated_TD_Parameters(PhysConst.potential_temperature_reference_pressure)
    # nelem/neqs/SD are read only by the DSGS allocator: DynSGS carries one
    # coefficient per ELEMENT plus neqs-long reduction buffers, and it returns
    # a closure struct only in 3D (in 1D/2D viscous_rhs_el! has its own DSGS
    # branch that never looks at params.sgs -- see sgsStructs.jl).
    sgs        = allocate_SGS(sem.mesh.npoin, TFloat, backend, PhysConst, inputs[:visc_model];
                              C_s   = inputs[:C_s],
                              nelem = sem.mesh.nelem,
                              neqs  = qp.neqs,
                              SD    = sem.mesh.SD,
                              C1    = get(inputs, :dsgs_C1, 1.0),
                              C2    = get(inputs, :dsgs_C2, 0.5))
    if sgs isa AbstractSGSModel
        # mod_inputs.jl always populates :lrichardson, so read it rather than
        # supplying a second, unreachable default here — the two disagreed, and
        # the mod_inputs one (false) silently won for every deck that did not
        # set the key, leaving the buoyancy correction off in runs that assumed
        # it was on.
        sgs.lrichardson   = inputs[:lrichardson]
        sgs.ltheta_eqn    = !(haskey(inputs, :energy_equation) && inputs[:energy_equation] == "energy")
        # Both closures now honour it. Vreman's limiter is a no-op when this is
        # false, so the key means the same thing for either model.
        sgs.lwall_damping = inputs[:lwall_damping] == true
        # DynSGS only. Adding the Smagorinsky viscosity to the residual one is
        # opt-in and off by default; see the field's comment in sgsStructs.jl
        # for why a wall-modelled PBL is the case where it is worth asking for.
        if sgs isa SGS_DSGS
            sgs.ladd_smagorinsky = get(inputs, :dsgs_add_smagorinsky, false) == true
        end

        if sgs.lwall_damping
            # Distance to the lower wall at every node. The wall-normal
            # coordinate is the last one: z in 3D, y in 2D.
            coords_h = Array(sem.mesh.coords)
            idir     = sem.mesh.SD == NSD_3D() ? 3 : 2
            # Reduce over ranks rather than reading mesh.zmin/ymin: those carry a
            # -1.0 sentinel until the mesh driver fills them, and a silently
            # shifted wall distance is worse than the cost of one Allreduce.
            wall0    = MPI.Allreduce(minimum(@view coords_h[idir, :]), MPI.MIN, comm)
            zw       = max.(@view(coords_h[idir, :]) .- wall0, 0.0)
            KernelAbstractions.copyto!(backend, sgs.zwall, TFloat.(zw))
            if inputs[:lwarp] && rank == 0
                @warn(":lwall_damping uses height above the domain floor, but :lwarp is on. " *
                      "Over terrain that is not the distance to the wall, so the near-wall " *
                      "limit will under-damp above the hill.")
            end
        end

        # :μ means two different things depending on :visc_model. Under AV() it
        # is the constant kinematic viscosity in m²/s; under a dynamic model it
        # multiplies the eddy viscosity that the closure has already computed,
        # so anything other than 0 or 1 silently rescales C_s by sqrt(:μ). Decks
        # carried AV-tuned values (5, 10, 15) into Smagorinsky runs this way,
        # which is a 2.2-3.9x inflation of the Smagorinsky constant.
        if rank == 0 && inputs[:lvisc]
            bad = [(i, m) for (i, m) in enumerate(inputs[:μ]) if m != 0.0 && m != 1.0]
            if !isempty(bad)
                @warn("visc_model $(typeof(inputs[:visc_model])) computes its own eddy viscosity, " *
                      "but :μ is not a 0/1 mask: entries $(bad) multiply it. " *
                      "Effective C_s = C_s*sqrt(:μ) = " *
                      "$(round.(sgs.C_s .* sqrt.(Float64[m for (_, m) in bad]), digits=3)). " *
                      "Set :μ => [0.0, 1, 1, ...] and tune :C_s instead.")
            end
        end
    end
    sgs_stress = zeros(TFloat, Int64(sem.mesh.npoin), 12)
    
    #------------------------------------------------------------------------------------
    # Populate solution arrays
    #------------------------------------------------------------------------------------
    if (sem.mesh.SD != NSD_1D()) && !(sem.mesh.lLaguerre)
        if rank == 0
            println(" # start conformity4ncf_q!")
        end
        g_dss_cache_qp = setup_assembler(sem.mesh.SD, qp.qn, sem.mesh.ip2gip, sem.mesh.gip2owner)
        conformity4ncf_q!(qp.qn, rhs_el_tmp, @view(utmp[:,:]), vaux, 
                          g_dss_cache_qp,
                          sem.mesh.SD, 
                          sem.QT, sem.mesh.connijk,
                          sem.mesh, sem.matrix.Minv, 
                          sem.metrics.Je, sem.ω, sem.AD, 
                          q_el, q_el_pro,
                          cache_ghost_p, q_ghost_p,
                          cache_ghost_c, q_ghost_c,
                          sem.interp; ladapt=inputs[:ladapt], neqs=qp.neqs)
        conformity4ncf_q!(qp.qe, rhs_el_tmp, @view(utmp[:,:]), vaux, 
                          g_dss_cache_qp,
                          sem.mesh.SD, 
                          sem.QT, sem.mesh.connijk, 
                          sem.mesh, sem.matrix.Minv, 
                          sem.metrics.Je, sem.ω, sem.AD, 
                          q_el, q_el_pro,
                          cache_ghost_p, q_ghost_p,
                          cache_ghost_c, q_ghost_c,
                          sem.interp; ladapt=inputs[:ladapt], neqs=qp.neqs)
        MPI.Barrier(comm)
        if rank == 0
            println(" # end conformity4ncf_q!")
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

    # Per-element, per-equation DSGS viscosity buffer.
    #   μ_dsgs[iel, ieq] is filled by compute_dsgs_viscosity! every RHS
    #   call (allocation-free). For 1D E-form Marras gives a single μ
    #   shared by all equations; the same value is written to every
    #   column so the VTU / PNG output sees per-equation slots. For 2D
    #   θ-form the columns carry distinct per-equation indicators:
    #     [:,1] = ν_ρ (diagnostic, not applied)
    #     [:,2] = μ_ρu          [:,3] = μ_ρv
    #     [:,4] = κ_θ  (already scaled by Pr/(γ-1))
    #
    # In 3D the matrix is DIAGNOSTIC ONLY. What the RHS applies there comes
    # from sgs.μ_el through SGS_diffusion (the same path Smagorinsky takes),
    # and these columns just mirror that split so write_output.jl has a
    # mu_dsgs_<var> field per equation:
    #     [:,1]   = :μ[1]·μ/Sc_t          (0 with the usual :μ[1] = 0)
    #     [:,2:4] = :μ[ieq]·μ             momentum
    #     [:,5]   = :μ[5]·μ/Pr_t          θ
    ldsgs     = inputs[:lvisc] == true && inputs[:visc_model] == DSGS()
    ldsgs_mhd = inputs[:lvisc] == true && inputs[:visc_model] == DSGS_MHD()

    # DYNSGS NEEDS TOTAL VARIABLES, and failing that it does not misbehave
    # mildly -- it blows up on the first step, which looks like a time-step or
    # a physics problem and is neither.
    #
    # The residual coefficient is normalised by the domain spread of the
    # solution (SGS.jl:729-731):
    #
    #     mu_res|e = C1 * D^2 * max_i ||R_i||inf,e / ||q_i - <q_i>||inf,Omega
    #
    # and that denominator is taken from `q` AS STORED. Under PERT() `q` is
    # already the deviation from the reference state, so for the usual quiescent
    # start the momentum slots are IDENTICALLY ZERO at t = 0: the denominator
    # collapses to `eps`, mu_res saturates its own limiter, and the run gets
    # mu_max = C2*D*(|u|+c) applied everywhere from step one. On CompEuler/theta
    # that is ~17,000 m^2/s, a diffusion number of 29 against an explicit limit
    # near 0.5 -- instantaneous instability.
    #
    # Marras et al. (JCP 2015) define the method on the total variables, and
    # every DSGS deck in problems/ uses TOTAL(). Refusing here rather than
    # warning: there is no partially-correct outcome to fall back to.
    if (ldsgs || ldsgs_mhd) && inputs[:SOL_VARS_TYPE] == PERT()
        error("DynSGS (:visc_model => ", inputs[:visc_model], ") requires ",
              ":SOL_VARS_TYPE => TOTAL(), but this case sets PERT(). The ",
              "residual coefficient is normalised by the domain spread of the ",
              "solution variables, and under PERT() those are perturbations -- ",
              "identically zero in the momentum slots at t = 0 for a quiescent ",
              "start -- so the normalisation collapses and the run is handed ",
              "its maximum viscosity on the first step. Set TOTAL(), or use a ",
              "non-dynamic :visc_model (SMAG(), VREM(), AV()), which are all ",
              "fine with PERT(). CompEuler/theta_dsgs is CompEuler/theta set up ",
              "for this model.")
    end
    if ldsgs || ldsgs_mhd
        μ_dsgs       = KernelAbstractions.zeros(backend, TFloat,
                                                Int64(sem.mesh.nelem), Int64(qp.neqs))
        μ_dsgs_pnode = KernelAbstractions.zeros(backend, TFloat,
                                                Int64(sem.mesh.npoin), Int64(qp.neqs))
    else
        μ_dsgs       = KernelAbstractions.zeros(backend, TFloat, 1, 1)
        μ_dsgs_pnode = KernelAbstractions.zeros(backend, TFloat, 1, 1)
    end

    # DynSGS-MHD extras.
    #
    # dsgs_qnm1/qnm2 are the BDF2 history the residual is built on. They
    # exist separately from qp.qnm1/qnm2 because those are advanced on
    # every RK *stage* — fine as generic scratch, useless as a time
    # derivative. rhs! advances this pair exactly once per time step, and
    # dsgs_thist records when it last did. Both start at the initial state
    # so the very first residual is identically zero rather than 3q/(2Δt).
    #
    # dsgs_avg / dsgs_denom are the per-equation domain-reduction scratch,
    # preallocated so compute_dsgs_viscosity! stays allocation-free.
    if ldsgs || ldsgs_mhd
        # Shaped like qp.qn / uaux, NOT (npoin, neqs): uaux carries one
        # extra trailing column (pressure) beyond the neqs solution slots,
        # which is why qp.qnm1/qnm2 are allocated from dims1 too. Sizing
        # these to neqs makes `dsgs_qnm2 .= uaux` a DimensionMismatch.
        dsgs_qnm1 = KernelAbstractions.zeros(backend, TFloat,
                                             Int64(size(qp.qn,1)), Int64(size(qp.qn,2)))
        dsgs_qnm2 = KernelAbstractions.zeros(backend, TFloat,
                                             Int64(size(qp.qn,1)), Int64(size(qp.qn,2)))
        for i = 1:size(qp.qn,2)
            dsgs_qnm1[:,i] = @view(qp.qn[:,i])
            dsgs_qnm2[:,i] = @view(qp.qn[:,i])
        end
        dsgs_avg   = KernelAbstractions.zeros(backend, TFloat, Int64(qp.neqs))
        dsgs_denom = KernelAbstractions.zeros(backend, TFloat, Int64(qp.neqs))
    else
        dsgs_qnm1  = KernelAbstractions.zeros(backend, TFloat, 1, 1)
        dsgs_qnm2  = KernelAbstractions.zeros(backend, TFloat, 1, 1)
        dsgs_avg   = KernelAbstractions.zeros(backend, TFloat, 1)
        dsgs_denom = KernelAbstractions.zeros(backend, TFloat, 1)
    end
    dsgs_thist = Ref{Float64}(-1.0e30)
    # Two flags set by _build_rhs! on every RHS call and read by
    # viscous_rhs_el!. Neither is ever snapshotted: both are recomputed from
    # `time` on every call, so the integrator warm-up, the precompile pass and
    # the restart path need no restore site for them.
    #
    #   dsgs_fresh  this is the FIRST RK stage of a step, so `uaux` is qⁿ and
    #               the two buffers hold qⁿ⁻¹ and qⁿ⁻²: the three levels the
    #               BDF2 stencil needs are consistent and the coefficient may
    #               be rebuilt. False on every other stage, where `uaux` is a
    #               stage state and rebuilding would measure the tendency
    #               rather than the residual.
    #   dsgs_hist   the buffers hold two GENUINELY DISTINCT past steps. False
    #               for the first two steps, when they still carry the initial
    #               state and the stencil cannot be a time derivative at all.
    #
    # See the comment at `ldsgs_step` in rhs.jl for what each of them is
    # protecting against.
    dsgs_fresh = Ref{Bool}(false)
    dsgs_hist  = Ref{Bool}(false)

    #---------------------------------------------------------------------
    # dsgs_wres -- WHICH NODES THE RESIDUAL MAY BE TAKEN AT. 1 interior,
    # 0 on a domain boundary.
    #
    # DynSGS measures how badly the SEMI-DISCRETE PDE is satisfied. At a node
    # whose value is set by a strongly-imposed boundary condition, it is not
    # that equation that advances the solution, and the residual there is not
    # a discretisation error -- it is the boundary condition.
    #
    # Concretely: apply_boundary_conditions_dirichlet! (rhs.jl, at the TOP of
    # every RHS call) projects the wall-normal momentum out of uaux at every
    # free-slip node. The inviscid RHS then puts it straight back, and the next
    # call projects it out again. So at a wall node
    #
    #     BDF2(qⁿ, qⁿ⁻¹, qⁿ⁻²)   sees the PROJECTED states
    #     M⁻¹·RHS                does NOT contain the projection
    #
    # and they differ by the whole projected flux, every step, for ever. It is
    # not small: measured on CompEuler/rtb2d_dsgs the worst residual in the
    # domain was on a boundary node at EVERY step of the run, always on the
    # wall-normal momentum equation, and it drove the element coefficient to
    # ~700 m²/s against ~14 in the interior. Since the coefficient is constant
    # per element, every element touching a wall lit up -- a picture of red
    # squares along the walls where there is no gradient at all, and the
    # bubble the model is supposed to be tracking two colour-bar decades down.
    #
    # Masking those nodes out of the L∞ costs nothing: an element keeps its
    # interior nodes (16 of 25 in 2D at nop = 4, 75 of 125 in 3D), and it is
    # only the max over the element that is being taken.
    #
    # The DENOMINATORS are not masked. They are domain norms of the SOLUTION,
    # which is perfectly well defined on a boundary.
    #---------------------------------------------------------------------
    dsgs_wres = ones(TFloat, Int64(sem.mesh.npoin))
    if ldsgs || ldsgs_mhd
        for arr in (sem.mesh.poin_in_bdy_face, sem.mesh.poin_in_bdy_edge)
            length(arr) == 0 && continue
            for k in eachindex(arr)
                ip = Int(arr[k])
                1 <= ip <= length(dsgs_wres) && (dsgs_wres[ip] = zero(TFloat))
            end
        end
        DSGS_NOMASK[] && fill!(dsgs_wres, one(TFloat))
        nint = count(!iszero, dsgs_wres)
        println_rank(" #   DynSGS residual taken on ", nint, " of ",
                     length(dsgs_wres), " nodes (",
                     length(dsgs_wres) - nint, " on the boundary, masked out)";
                     msg_rank = rank, suppress = sem.mesh.msg_suppress)
    end
    # JEXPRESSO_DSGS_MONITOR=1 -> one line per step with what the model
    # produced (kernel/physics/SGS.jl, _dsgs_monitor). Read once here rather
    # than per RHS call: `haskey(ENV, ...)` on the hot path is a dictionary
    # lookup per step for a value that cannot change.
    DSGS_MONITOR[] = get(ENV, "JEXPRESSO_DSGS_MONITOR", "0") ∉ ("0", "false", "")
    # :dsgs_residual => :tendency (default) | :strict. See the long comment on
    # DSGS_STRICT in kernel/physics/SGS.jl for the measurement that decided the
    # default -- the two choices measure different things and :strict, which is
    # the one that is literally BDF2, reads the model's own viscous term.
    _dsgs_res = Symbol(get(inputs, :dsgs_residual, :tendency))
    _dsgs_res in (:tendency, :strict) ||
        error(" # ERROR params_setup.jl: :dsgs_residual must be :tendency ",
              "(the default) or :strict; got $(_dsgs_res).")
    DSGS_STRICT[] = _dsgs_res === :strict
    DSGS_NOMASK[]      = get(ENV, "JEXPRESSO_DSGS_NOMASK",     "0") ∉ ("0","false","")

    # Per-equation scratch the 2D DSGS path uses to pack the
    # per-element coefficient before calling _expansion_visc!:
    #   visc_coeff_dsgs[1] = 0                          (mass)
    #   visc_coeff_dsgs[2..end-1] = μ_dsgs[iel]         (momentum)
    #   visc_coeff_dsgs[end]      = Pr/(γ-1)·μ_dsgs[iel] (energy/θ)
    visc_coeff_dsgs = KernelAbstractions.zeros(backend, TFloat, Int64(qp.neqs))

    # setup timer
    timers = Dict{String, MPIFunctionTimer}()

    # LES statistics z-level cache (computed once, shared across timesteps)
    nprofiles        = length(inputs[:lesprofile_vars])
    nstress          = length(inputs[:lesstress_vars])
    les_stat_cache    = build_les_stat_cache(sem.mesh, nprofiles, nstress, TFloat, backend)
    les_cross_section = build_les_cross_section(sem.mesh, nprofiles, nstress, TFloat)
    les_bottom_cache  = build_les_bottom_cache(sem.mesh, sem.metrics, inputs)

    #------------------------------------------------------------------------------------
    # Populate params tuple to carry global arrays and constants around
    #------------------------------------------------------------------------------------
    if (sem.mesh.lLaguerre ||
        inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])
        g_dss_cache = setup_assembler(sem.mesh.SD, RHS, sem.mesh.ip2gip, sem.mesh.gip2owner)
        params = (backend, T, F, G, H, S,
                  uaux, vaux, utmp, fluxaux,
                  ubdy, gradu, bdy_flux, #for B.C.
                  rhs_el, rhs_diff_el, rhs_el_tmp,
                  rhs_diffξ_el, rhs_diffη_el,rhs_diffζ_el, μ_max,
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
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.AD,
                  sem.SOL_VARS_TYPE, sem.volume_flux,
                  neqs=qp.neqs,
                  sem.mesh,
                  sem.connijk_original, sem.poin_in_bdy_face_original, sem.x_original, sem.y_original, sem.z_original,
		  basis=sem.basis[1], basis_lag = sem.basis[2],
                  ω = sem.ω[1], ω_lag = sem.ω[2],
                  metrics = sem.metrics[1], metrics_lag = sem.metrics[2], 
                  inputs, VT = inputs[:visc_model], visc_coeff, μ_dsgs, μ_dsgs_pnode, visc_coeff_dsgs,
                  dsgs_qnm1, dsgs_qnm2, dsgs_avg, dsgs_denom, dsgs_thist, dsgs_fresh, dsgs_hist, dsgs_wres,
                  WM,
                  sem.matrix.M, sem.matrix.Minv, g_dss_cache=g_dss_cache, tspan,
                  Δt, deps, xmax, xmin, ymax, ymin, zmin, zmax,
                  qp, mp, sem.fx, sem.fy, fy_t, sem.fy_lag, fy_t_lag, sem.fz, fz_t, laguerre=true,
                  les_stat_cache,
                  les_cross_section,
                  les_bottom_cache,
                  sgs,
                  sgs_stress,
                  timers,
                  coupling = coupling)

    else
        g_dss_cache = setup_assembler(sem.mesh.SD, RHS, sem.mesh.ip2gip, sem.mesh.gip2owner)
        
        params = (backend,
                  T, inputs,
                  uaux, vaux, utmp, fluxaux,
                  ubdy, gradu, bdy_flux,
                  RHS, RHS_visc,
                  fijk, ∇f_el,
                  rhs_el, rhs_diff_el, rhs_el_tmp,
                  rhs_diffξ_el, rhs_diffη_el, rhs_diffζ_el, μ_max,
                  uprimitive,
                  F, G, H, S,
                  F_surf, S_face, S_flux, M_surf_inv = sem.matrix.M_surf_inv, M_edge_inv = sem.matrix.M_edge_inv,
                  flux_gpu, source_gpu, qbdy_gpu,
                  flux_micro, source_micro, adjusted, Pm,
                  q_el, q_el_pro, q_ghost_p, q_ghost_c,
                  cache_ghost_p, cache_ghost_c,
                  q_t, q_ti, q_tij, fqf, b, B,
                  SD=sem.mesh.SD, sem.QT, sem.CL, sem.AD,
                  sem.SOL_VARS_TYPE, sem.volume_flux,
                  neqs=qp.neqs,
                  sem.connijk_original, sem.poin_in_bdy_face_original, sem.x_original, sem.y_original, sem.z_original,
                  sem.basis, sem.ω, sem.mesh, sem.metrics,
                  thermo_params, VT = inputs[:visc_model], visc_coeff, μ_dsgs, μ_dsgs_pnode, visc_coeff_dsgs,
                  dsgs_qnm1, dsgs_qnm2, dsgs_avg, dsgs_denom, dsgs_thist, dsgs_fresh, dsgs_hist, dsgs_wres,
                  sem.matrix.M, sem.matrix.Minv, g_dss_cache=g_dss_cache,
                  tspan, Δt, xmax, xmin, ymax, ymin, zmin, zmax,
                  WM,
                  phys_grid = sem.phys_grid,
                  atmos_data = sem.atmos_data,
                  qp, mp, LST, sem.fx, sem.fy, fy_t, sem.fz, fz_t, laguerre=false,
                  les_stat_cache,
                  les_cross_section,
                  les_bottom_cache,
                  sgs,
                  sgs_stress,
                  OUTPUT_DIR,
                  timers,
                  sem.interp, sem.project, sem.nparts, sem.distribute,
                  coupling = coupling)
    end

    #------------------------------------------------------------------------------------
    # HEVI: the vertical implicit operator reads the mesh, the metrics, the
    # basis and the reference state off `params`, so it can only be built once
    # `params` exists -- and it has to live ON `params`, because the integrator
    # reaches it as `p.hevi`. NamedTuple merge makes a new object; that happens
    # exactly once, here.
    #
    # Built eagerly rather than on first use: setup does global reductions and
    # an all-to-all, and putting a collective behind a lazily-taken branch is
    # how a run deadlocks on the rank that took the other one.
    #------------------------------------------------------------------------------------
    if hevi_enabled(inputs)
        if sem.mesh.SD != NSD_3D()
            error("HEVI is 3D only; this case is $(typeof(sem.mesh.SD)). ",
                  "Use an explicit integrator, or drop :lhevi / :ode_solver => HEVI_ARK.")
        elseif backend != CPU()
            error("HEVI has no GPU path: the column solve is a LAPACK banded LU on the ",
                  "host, and the gather/scatter moves host arrays. Run with CPU(), or ",
                  "use an explicit integrator on the GPU.")
        end
        params = merge(params, (hevi = build_hevi(params, inputs),))
    end

    #------------------------------------------------------------------------------------
    # IMEX3D: same story as HEVI, and for the same reasons -- the 3D acoustic
    # operator, its column preconditioner and the Krylov workspace all read the
    # mesh, the metrics, the basis and the reference state off `params`, and the
    # integrator reaches them as `p.imex`.
    #
    # Eagerly, not lazily: setup does global reductions, an all-to-all and (with
    # :imex_verify on) a full Krylov solve. Putting collectives behind a lazily
    # taken branch is how a run deadlocks on the rank that took the other one.
    #------------------------------------------------------------------------------------
    if imex3d_enabled(inputs)
        if sem.mesh.SD != NSD_3D()
            error("IMEX3D is 3D only; this case is $(typeof(sem.mesh.SD)). ",
                  "Use an explicit integrator, or drop :limex / :ode_solver => IMEX_ARK.")
        elseif backend != CPU()
            error("IMEX3D has no GPU path: the column preconditioner is a LAPACK banded ",
                  "LU on the host and the gather/scatter moves host arrays. Run with ",
                  "CPU(), or use an explicit integrator on the GPU.")
        end
        params = merge(params, (imex = build_imex3d(params, inputs),))
    end

    # SPLIT_EXPLICIT: same story as HEVI -- the fast operator reads the mesh,
    # the metrics, the basis and the reference state off `params`, and the
    # integrator reaches it as `p.substep`.
    if inputs[:ode_solver] isa SPLIT_EXPLICIT
        if sem.mesh.SD != NSD_3D()
            error("SPLIT_EXPLICIT is 3D only; this case is $(typeof(sem.mesh.SD)).")
        elseif backend != CPU()
            error("SPLIT_EXPLICIT has no GPU path.")
        end
        params = merge(params, (substep = build_substep(params, inputs, inputs[:Δt]),))
    end

    println_rank(" # Build arrays and params ................................ DONE"; msg_rank = rank, suppress = sem.mesh.msg_suppress)

    return params, u
    
end
