function driver(nranks,
                distribute,
                inputs::Dict,
                OUTPUT_DIR::String,
                TFloat,
                world,
                is_coupled::Bool = je_perform_coupling_handshake(world, nranks))

    if rank == 0 @info " Params_setup .................................." end

    tspan = [TFloat(inputs[:tinit]), TFloat(inputs[:tend])]
    if (inputs[:lamr] == true)
        amr_freq = inputs[:amr_freq]
        Δt_amr   = amr_freq * inputs[:Δt]
        tspan    = [TFloat(inputs[:tinit]), TFloat(inputs[:tinit] + Δt_amr)]
    end

    coupling = nothing
    lsize = nranks
    if rank == 0 println(" is_coupled: ", is_coupled) end    
    if is_coupled

        # 2. Complete coupling setup
        coupling, sem, partitioned_model, qp = @time setup_coupling_and_mesh(
            world, lsize, inputs, nranks, distribute, rank, OUTPUT_DIR, TFloat
        )
        
        # Now call params_setup with correct order: sem first, then coupling
        params, u = @time params_setup(
            sem, 
            coupling,  # Pass coupling as second argument
            qp, 
            inputs, 
            OUTPUT_DIR, 
            TFloat, 
            tspan
        )
        
    else
        #---------------------------------------------------------
        # SEM setup
        #---------------------------------------------------------
        sem, partitioned_model = sem_setup(inputs, lsize, distribute, rank)   
        
        #---------------------------------------------------------
        # Initialize.jl is contained in the user's problem case directory
        #---------------------------------------------------------
        qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)

        params, u = params_setup(
            sem,
            nothing,  # No coupling in standalone mode
            qp,
            inputs,
            OUTPUT_DIR,
            TFloat,
            tspan
        )
        
    end
    
    #---------------------------------------------------------
    # Parameters setup
    #---------------------------------------------------------   
  
    if rank == 0 @info " Params_setup .................................. END" end
    
    if !inputs[:llinsolve]
        #---------------------------------------------------------
        # Evolutionary problems that lead to Mdq/dt = RHS
        #---------------------------------------------------------
        
        @time solution = time_loop!(
            inputs,
            params,
            u,
            partitioned_model,
            is_coupled,
            is_coupled ? coupling : nothing  # Pass coupling object
        )
        
    else
        #---------------------------------------------------------
        # Problems that lead to Lx = RHS
        #---------------------------------------------------------
        RHS   = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin))
        Mdiag = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(sem.mesh.npoin))

        EL = allocate_elemLearning(sem.mesh.nelem, sem.mesh.ngl,
                                   sem.mesh.length∂O,
                                   sem.mesh.length∂τ,
                                   sem.mesh.lengthΓ,
                                   TFloat, inputs[:backend];
                                   Nsamp=inputs[:Nsamp])
        
        if (inputs[:backend] == CPU())
            
            bufferin  = Vector{Vector{Float64}}()
            bufferout = Vector{Vector{Float64}}()
            total_cols_writtenin  = 0  # Track how many columns we've written
            total_cols_writtenout = 0  # Track how many columns we've written

            # Clear/initialize file at start
            if isfile("input_tensor.csv")
                rm("input_tensor.csv")
            end
            if isfile("output_tensor.csv")
                rm("output_tensor.csv")
            end
            for isamp=1:inputs[:Nsamp]
                #
                # L*q = M*RHS   See algo 12.18 of Giraldo's book
                #
                #Minv          = diagm(sem.matrix.Minv)
                #sem.matrix.L .= Minv * sem.matrix.L
                
                # 2.a/b
                μ             = 1
                avisc         = zeros(TFloat, sem.mesh.ngl^2)
                ranvisc       = isamp*μ #+ 10*rand()
                avisc[:]     .= ranvisc
                #sem.matrix.L .= ranvisc*sem.matrix.L
                
                for ip =1:sem.mesh.npoin
                    RHS[ip] = user_source!(RHS[ip],
                                           params.qp.qn[ip],
                                           params.qp.qe[ip],
                                           sem.mesh.npoin,
                                           inputs[:CL], inputs[:SOL_VARS_TYPE];
                                           neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip],
                                           xmax=sem.mesh.xmax, xmin=sem.mesh.xmin,
                                           ymax=sem.mesh.ymax, ymin=sem.mesh.ymin)
                end
                RHS = sem.matrix.M.*RHS
                
                if inputs[:lsparse] ==  false
                    for ip = 1:sem.mesh.npoin
                        sem.matrix.L[ip,ip] += inputs[:rconst][1]
                    end
                end
                
                apply_boundary_conditions_lin_solve!(sem.matrix.L,
                                                     0.0, params.qp.qe,
                                                     params.mesh.coords,
                                                     params.metrics.nx,
                                                     params.metrics.ny,
                                                     params.metrics.nz,
                                                     sem.mesh.npoin,
                                                     params.mesh.npoin_linear, 
                                                     params.mesh.poin_in_bdy_edge,
                                                     params.mesh.poin_in_bdy_face,
                                                     params.mesh.nedges_bdy,
                                                     params.mesh.nfaces_bdy,
                                                     params.mesh.ngl, params.mesh.ngr,
                                                     params.mesh.nelem_semi_inf,
                                                     params.basis.ψ, params.basis.dψ,
                                                     0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                     RHS, 0.0, params.ubdy,
                                                     params.mesh.connijk_lag,
                                                     params.mesh.bdy_edge_in_elem,
                                                     params.mesh.bdy_edge_type,
                                                     params.ω, qp.neqs,
                                                     params.inputs, params.AD, sem.mesh.SD)
                
                #-----------------------------------------------------
                # Element-learning infrastructure
                #-----------------------------------------------------
                if inputs[:lelementLearning] == false && inputs[:lsparse] ==  false
                    println(" # Solve x=inv(A)*b: full storage")
                    @time solution = solveAx(sem.matrix.L, RHS, inputs[:ode_solver])
                else
                    if inputs[:lelementLearning]
                        
                        @time elementLearning_Axb!(params.qp.qn, params.uaux, sem.mesh,
                                                   sem.matrix.L, RHS, EL,
                                                   avisc,
                                                   bufferin, bufferout;
                                                   isamp=isamp,
                                                   total_cols_writtenin=total_cols_writtenin,
                                                   total_cols_writtenout=total_cols_writtenout)
                        
                    elseif inputs[:lsparse]
                        println(" # Solve x=inv(A)*b: sparse storage")
                        @time params.qp.qn = sem.matrix.L\RHS
                    end
                end

                
                #usol = inputs[:lelementLearning] ? params.qp.qn : solution.u
                usol = params.qp.qn
                
                args = (params.SD, usol, params.uaux, 1, isamp,
                        sem.mesh, nothing,
                        nothing, nothing,
                        0.0, 0.0, 0.0,
                        OUTPUT_DIR, inputs,
                        params.qp.qvars,
                        params.qp.qoutvars,
                        inputs[:outformat])
                
                write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe)
                
                @info "isamp" isamp
                
            end #isamp loop
            
            #-----------------------------------------------------
            # END Element-learning infrastructure
            #-----------------------------------------------------
            
        else
            println( " ")
            println( " WARNING!!! drivers.jl:L114")
            println( " WARNING: CHECK IF THIS GPU IMPLEMENTATION OF Ax=b still works")
            println( " ")
            nothing
            #=k = lin_solve_rhs_gpu_2d!(inputs[:backend])
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
            KernelAbstractions.synchronize(inputs[:backend])=#
        end
        
        #usol = inputs[:lelementLearning] ? params.qp.qn : solution.u
        if inputs[:lelementLearning] || inputs[:lsparse]
            usol = params.qp.qn
        else
            #usol = params.qp.qn
            usol = solution.u
        end
        args = (params.SD, usol, params.uaux, 0.0, 1,
                sem.mesh, nothing,
                nothing, nothing,
                0.0, 0.0, 0.0,
                OUTPUT_DIR, inputs,
                params.qp.qvars,
                params.qp.qoutvars,
                inputs[:outformat])

        write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe)
        
    end

end
