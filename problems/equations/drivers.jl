function driver(nparts,
                distribute,
                inputs::Dict,
                    OUTPUT_DIR::String,
                TFloat) 
    comm  = distribute.comm
    rank = MPI.Comm_rank(comm)
    sem = sem_setup(inputs, nparts, distribute)
    
    if (inputs[:backend] != CPU())
        convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs)
    end
    
    qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)

    # test of projection matrix for solutions from old to new, i.e., coarse to fine, fine to coarse
    # test_projection_solutions(sem.mesh, qp, sem.partitioned_model, inputs, nparts, sem.distribute)
    if inputs[:ladapt] == true
        if rank == 0
            @info "start conformity4ncf_q!"
        end
        pM = setup_assembler(sem.mesh.SD, qp.qn, sem.mesh.ip2gip, sem.mesh.gip2owner)
        @time conformity4ncf_q!(qp.qn, pM, sem.mesh.SD, sem.QT, sem.mesh.connijk, sem.mesh, sem.matrix.Minv, sem.metrics.Je, sem.ω, sem.AD, qp.neqs+1, sem.interp)
        @time conformity4ncf_q!(qp.qe, pM, sem.mesh.SD, sem.QT, sem.mesh.connijk, sem.mesh, sem.matrix.Minv, sem.metrics.Je, sem.ω, sem.AD, qp.neqs+1, sem.interp)
        
        MPI.Barrier(comm)
        if rank == 0
            @info "end conformity4ncf_q!"
        end
    end

    if (inputs[:amr] == true)
        amr_freq = inputs[:amr_freq]
        Δt_amr   = amr_freq * inputs[:Δt]
        tspan    = [TFloat(inputs[:tinit]), TFloat(inputs[:tinit] + Δt_amr)]
    else
        tspan = [TFloat(inputs[:tinit]), TFloat(inputs[:tend])]
    end
    params, u =  params_setup(sem,
                              qp,
                              inputs,
                              OUTPUT_DIR,
                              TFloat,
                              tspan)
    
    if !inputs[:llinsolve]
        #
        # Hyperbolic/parabolic problems that lead to Mdq/dt = RHS
        #
        @time solution = time_loop!(inputs, params, u)
        # PLOT NOTICE: Plotting is called from inside time_loop using callbacks.
        
    else
        #
        # Problems that lead to Ax = b
        #
        RHS = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin))

        if (inputs[:backend] == CPU())
          
            Minv = diagm(sem.matrix.Minv)
            
            L_temp = Minv * sem.matrix.L
            sem.matrix.L .= L_temp
            
            for ip =1:sem.mesh.npoin
                b = user_source(RHS[ip],
                                params.qp.qn[ip],
                                params.qp.qe[ip],
                                sem.mesh.npoin, inputs[:CL], inputs[:SOL_VARS_TYPE];
                                neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip])
                RHS[ip] = b
            end

            for ip = 1:sem.mesh.npoin
                sem.matrix.L[ip,ip] += inputs[:rconst][1]
            end

            #-----------------------------------------------------
            # Element-learning infrastructure
            #-----------------------------------------------------
            if inputs[:lelementLearning]
                elementLearning_Axb(sem.mesh, sem.matrix.L, RHS)
            end
            #-----------------------------------------------------
            # END Element-learning infrastructure
            #-----------------------------------------------------
            
            apply_boundary_conditions_lin_solve!(sem.matrix.L, 0.0, params.qp.qe,
                                                 params.mesh.x, params.mesh.y, params.mesh.z,
                                                 params.metrics.nx,
                                                 params.metrics.ny,
                                                 params.metrics.nz,
                                                 sem.mesh.npoin, params.mesh.npoin_linear, 
                                                 params.mesh.poin_in_bdy_edge,
                                                 params.mesh.poin_in_bdy_face,
                                                 params.mesh.nedges_bdy,
                                                 params.mesh.nfaces_bdy,
                                                 params.mesh.ngl, params.mesh.ngr,
                                                 params.mesh.nelem_semi_inf,
                                                 params.basis.ψ, params.basis.dψ,
                                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                 RHS, 0.0, params.ubdy,
                                                 params.mesh.connijk_lag, params.mesh.bdy_edge_in_elem,
                                                 params.mesh.bdy_edge_type,
                                                 params.ω, qp.neqs, params.inputs, params.AD, sem.mesh.SD)
      
         
        else
            k = lin_solve_rhs_gpu_2d!(inputs[:backend])
            k(RHS, qp.qn, qp.qe, sem.mesh.x, sem.mesh.y, qp.neqs; ndrange = sem.mesh.npoin)
            KernelAbstractions.synchronize(inputs[:backend])
            
            Minv = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin), Int64(sem.mesh.npoin))
            k =  diagm_gpu!(inputs[:backend])
            k(Minv , sem.matrix.Minv; ndrange = sem.mesh.npoin)
            L_temp = Minv * sem.matrix.L
            sem.matrix.L .= L_temp

            k = apply_boundary_conditions_gpu_lin_solve!(inputs[:backend])
            k(RHS, sem.matrix.L, sem.mesh.poin_in_bdy_edge, sem.mesh.npoin; ndrange = (sem.mesh.nedges_bdy*sem.mesh.ngl), workgroupsize = (sem.mesh.ngl))
            KernelAbstractions.synchronize(inputs[:backend])
            if ("Laguerre" in sem.mesh.bdy_edge_type)
                k = apply_boundary_conditions_lag_gpu_lin_solve!(inputs[:backend])
                k(RHS, sem.matrix.L, sem.mesh.connijk_lag, sem.mesh.ngl, sem.mesh.ngr, sem.mesh.npoin, sem.mesh.nelem_semi_inf, inputs[:lperiodic_laguerre];
                  ndrange = (sem.mesh.nelem_semi_inf*sem.mesh.ngl,sem.mesh.ngr), workgroupsize = (sem.mesh.ngl,sem.mesh.ngr))
                KernelAbstractions.synchronize(inputs[:backend])
            end
            k = add_to_diag!(inputs[:backend])
            k(sem.matrix.L, TFloat(10.0); ndrange = sem.mesh.npoin)
            KernelAbstractions.synchronize(inputs[:backend])
        end
        
        @time solution = solveAx(sem.matrix.L, RHS, inputs[:ode_solver])
        
        write_output(params.SD, solution.u, params.uaux, 0.0, 1,
                     sem.mesh, nothing,
                     nothing, nothing,
                     0.0, 0.0, 0.0,
                     OUTPUT_DIR, inputs,
                     params.qp.qvars,
                     params.qp.qoutvars,
                     inputs[:outformat];
                     nvar=params.qp.neqs, qexact=params.qp.qe)
        
    end
end

function elementLearning_Axb(mesh::St_mesh, A, RHS)

    @info "∂Oxdd"
    println(mesh.∂O)
    @info "∂τddddd"
    println(mesh.∂τ)

    mesh.lengthO =  mesh.length∂O +  mesh.lengthτO
    
    @info mesh.lengthΓ, mesh.lengthO, mesh.length∂τ, mesh.lengthτO, mesh.length∂O
    
    EL = allocate_elemLearning(mesh.nelem, mesh.ngl,
                               mesh.length∂O,
                               mesh.length∂τ,
                               TFloat, inputs[:backend])

    nelintpoints = (mesh.ngl-2)*(mesh.ngl-2)
    nelpoints = size(mesh.conn)[2]
    intaux = nelpoints - nelintpoints
    
    for iel=1:mesh.nelem

        #
        # A∂oᵥₒₜ
        #
        ii = 1
        for i = intaux+1:nelpoints
            ipo = mesh.conn[iel, i]

            for i1=1:length(mesh.∂O)
                
                iO = mesh.∂O[i1]
                
                EL.A∂Ovo[i1, ii, iel] = A[iO, ipo]
            end

            #
            # Aᵥₒ∂τₜ
            #
            for i1=1:length(mesh.∂τ)
                
                iτ = mesh.∂τ[i1]
                
                EL.Avo∂τ[ii, i1, iel] = A[ipo, iτ]
            end

            #
            # Aᵥₒᵥₒ
            #
            jj = 1
            for j = intaux+1:nelpoints          
                jpo = mesh.conn[iel, j]
                
                EL.Avovo[ii,jj,iel] = A[ipo, jpo]
                println(EL.Avovo[ii, jj, iel])

                jj += 1
                
            end
            ii += 1
        end

        #
        # Hᵥₒᵥₒ[iel] = A⁻¹ᵥₒᵥₒ[iel]
        # 
        EL.Hvovo[:,:,iel] = inv(EL.Avovo[:,:,iel])
        
    end
    
    #
    # A∂O∂τ ⊂ A∂τ∂τ
    #
    for j1=1:length(mesh.∂τ)
        jτ1 = mesh.∂τ[j1]
        
        for i1=1:length(mesh.∂O)
            
            iO1 = mesh.∂O[i1]
            
            EL.A∂O∂τ[i1, j1] = A[iO1, jτ1]
        end
        
        for j2=1:length(mesh.∂τ)
            jτ2 = mesh.∂τ[j2]
            
            EL.A∂τ∂τ[j1, j2] = A[jτ1, jτ2]
        end
            
    end
    
    #
    # B∂O∂τ[:,:] = A∂O∂τ - Sum_{iel} A∂Oᵥₒ[:,:,iel]*A⁻¹ᵥₒᵥ[:,:,iel]*Aᵥₒ∂τ[:,:,iel]
    #
    intermediate_product = zeros(mesh.length∂O, mesh.length∂τ, mesh.nelem)
    for i in 1:size(EL.Avo∂τ, 3)
        intermediate_product[:,:,i] = EL.A∂Ovo[:,:,i]*EL.Hvovo[:,:,i]*EL.Avo∂τ[:,:,i]
    end
    EL.B∂O∂τ = EL.A∂O∂τ - sum(intermediate_product, dims=3)
    
    @info size(EL.B∂O∂τ)
    @mystop
    
end
