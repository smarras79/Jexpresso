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
    if (sem.mesh.SD == NSD_2D())
        build_radiative_transfer_problem(sem.mesh, inputs, 1, sem.mesh.ngl, sem.basis.dψ, sem.basis.ψ, sem.ω, sem.metrics.Je, 
                                     sem.metrics.dξdx, sem.metrics.dξdy, sem.metrics.dηdx, sem.metrics.dηdy, 
                                     sem.metrics.nx, sem.metrics.ny, sem.mesh.elem_to_edge, sem.mesh.extra_mesh, sem.QT, NSD_2D(), sem.AD)
    else
        build_radiative_transfer_problem(sem.mesh, inputs, 1, sem.mesh.ngl, sem.basis.dψ, sem.basis.ψ, sem.ω, sem.metrics.Je,
                                     sem.metrics.dξdx, sem.metrics.dξdy, sem.metrics.dξdz, 
                                     sem.metrics.dηdx, sem.metrics.dηdy, sem.metrics.dηdz,
                                     sem.metrics.dζdx, sem.metrics.dζdy, sem.metrics.dζdz,
                                     sem.metrics.nx, sem.metrics.ny, sem.metrics.nz, 
                                     sem.mesh.elem_to_face, sem.mesh.extra_mesh, sem.QT, NSD_3D(), sem.AD)
    
    end
    #=
    qp = initialize(sem.mesh.SD, sem.PT, sem.mesh, inputs, OUTPUT_DIR, TFloat)

    # test of projection matrix for solutions from old to new, i.e., coarse to fine, fine to coarse
    # test_projection_solutions(sem.mesh, qp, sem.partitioned_model, inputs, nparts, sem.distribute)
    if inputs[:ladapt] == true
        if rank == 0
            @info "start conformity4ncf_q!"
        end
        @time conformity4ncf_q!(qp.qn, sem.matrix.pM, sem.mesh.SD, sem.QT, sem.mesh.connijk, sem.mesh, sem.matrix.Minv, sem.metrics.Je, sem.ω, sem.AD, qp.neqs+1, sem.interp)
        @time conformity4ncf_q!(qp.qe, sem.matrix.pM, sem.mesh.SD, sem.QT, sem.mesh.connijk, sem.mesh, sem.matrix.Minv, sem.metrics.Je, sem.ω, sem.AD, qp.neqs+1, sem.interp)
        
        MPI.Barrier(comm)
        if rank == 0
            @info "end conformity4ncf_q!"
        end
    end

    if (inputs[:ladapt] == true) && (inputs[:amr] == true)
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

            
            #-----------------------------------------------------
            # Element-learning infrastructure
            #-----------------------------------------------------
            if inputs[:lelementLearning]
                elementLearning_Axb!(params.qp.qn, params.uaux, sem.mesh, sem.matrix.L, RHS)
            else
                solution = solveAx(sem.matrix.L, RHS, inputs[:ode_solver])
            end
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
        
        usol = inputs[:lelementLearning] ? params.qp.qn : solution.u
        args = (params.SD, usol, params.uaux, 0.0, 1,
                     sem.mesh, nothing,
                     nothing, nothing,
                     0.0, 0.0, 0.0,
                     OUTPUT_DIR, inputs,
                     params.qp.qvars,
                     params.qp.qoutvars,
                     inputs[:outformat])

        write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe)
        
    end=#
end

function elementLearning_Axb!(u, uaux, mesh::St_mesh, A, ubdy)
    
    mesh.lengthO =  mesh.length∂O +  mesh.lengthIo
        
    EL = allocate_elemLearning(mesh.nelem, mesh.ngl,
                               mesh.length∂O,
                               mesh.length∂τ,
                               mesh.lengthΓ,
                               TFloat, inputs[:backend])

    nelintpoints = (mesh.ngl-2)*(mesh.ngl-2)
    nelpoints = size(mesh.conn)[2]
    elnbdypoints = nelpoints - nelintpoints
    for iel=1:mesh.nelem
        #
        # A∂oᵥₒ
        #
        ii = 1
        for i = elnbdypoints+1:nelpoints
            ipo = mesh.conn[iel, i]
            
            for io=1:length(mesh.∂O)
                io1 = mesh.∂O[io]
                EL.A∂Ovo[io, ii, iel] = A[io1, ipo]
            end

            #
            # Aᵥₒ∂τ
            #
            for jτ = 1:mesh.length∂τ
                jτ1 = mesh.∂τ[jτ]
                EL.Avo∂τ[ii, jτ, iel] = A[ipo, jτ1]
            end

            #
            # Aᵥₒᵥₒ
            #
            jj = 1
            for j = elnbdypoints+1:nelpoints          
                jpo = mesh.conn[iel, j]
                
                EL.Avovo[ii, jj, iel] = A[ipo, jpo]
                jj += 1
            end
            
            #
            # Aᵥₒᵥb
            #
            for j = 1:elnbdypoints
                jpb = mesh.conn[iel, j]
                
                EL.Avovb[ii, j, iel] = A[ipo, jpb]
            end
            
            ii += 1
        end
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
    # A∂OIo
    #
    for jo=1:mesh.length∂O
        jo1 = mesh.∂O[jo]
        
        for io=1:mesh.lengthIo            
            io1 = mesh.Io[io]
            
            EL.A∂OIo[jo, io] = A[jo1, io1]
        end
    end
     #
    # AIo∂O
    #
    for jo=1:mesh.length∂O
        jo1 = mesh.∂O[jo]
        
        for io=1:mesh.lengthIo            
            io1 = mesh.Io[io]
            
            EL.AIo∂O[io, jo] = A[io1, jo1]
        end
    end
    #
    # AIoIo
    #
    for io = 1:mesh.lengthIo
        io1 = mesh.Io[io]
        
        for jo = 1:mesh.lengthIo
            jo1 = mesh.Io[jo]

            EL.AIoIo[io, jo] = A[io1, jo1]
        end
    end
    #
    # AIo∂τ
    #
    for jτ = 1:mesh.length∂τ
        jτ1 = mesh.∂τ[jτ]
        
        for io=1:mesh.lengthIo
            io1 = mesh.Io[io]
            
            EL.AIo∂τ[io, jτ] = A[io1, jτ1]
        end
    end
    invAIoIo = similar(EL.AIoIo)
    invAIoIo = inv(EL.AIoIo)

    dims = (mesh.lengthIo, mesh.lengthΓ)
    AIoΓ = similar(EL.AIoIo, dims);
    
    #------------------------------------------------------------------------
    # Eq. (13)
    #------------------------------------------------------------------------    
    #  B∂O∂τ[:,:] = A∂O∂τ - Sum_{iel} A∂Oᵥₒ[:,:,iel]*A⁻¹ᵥₒᵥₒ[:,:,iel]*Aᵥₒ∂τ[:,:,iel] -> A∂O∂τ - Sum_{iel}A⋅B⋅C
    #
    
    #
    # LOCAL VERSION (eq 13)
    #
    ABC = zeros(mesh.length∂O, mesh.length∂τ, mesh.nelem)
    BC  = zeros(size(EL.Avo∂τ)[1], size(EL.Avo∂τ)[2])
    for iel = 1:mesh.nelem
        
        # BC = A⁻¹ᵥₒᵥₒ[:,:,iel]⋅Aᵥₒ∂τ[:,:,iel]
        LinearAlgebra.mul!(BC, inv(EL.Avovo[:,:,iel]), EL.Avo∂τ[:,:,iel])
        
        # ABC = A∂Oᵥₒ[:,:,iel]⋅BC
        LinearAlgebra.mul!(@view(ABC[:,:,iel]), @view(EL.A∂Ovo[:,:,iel]), @view(BC[:,:]))
    end
    ∑el = similar(EL.A∂O∂τ)
    ∑el = sum(ABC, dims=3)
    EL.B∂O∂τ = EL.A∂O∂τ - ∑el # (13)

    #
    # WARNING: for large grids this double loop may be a bottleneck
    #
    for i1=1:length(mesh.∂O)      #row    B[i1][i2]        
        for i2=1:length(mesh.∂O)  #column B[i1][i2]
            
            j2 = findall(x->x==mesh.∂O[i2], mesh.∂τ)[1]
            EL.B∂O∂O[i1, i2] = EL.B∂O∂τ[i1, j2]
        end        
    end
    
    gΓ = zeros(mesh.lengthΓ)
    for iΓ = 1:mesh.lengthΓ
        g1=mesh.Γ[iΓ]
        
        jτ = findall(x->x==mesh.Γ[iΓ], mesh.∂τ)[1]
        EL.B∂O∂Γ[:, iΓ] .= EL.B∂O∂τ[:, jτ]

        gΓ[iΓ] = ubdy[g1, 1]
    end
    
    #------------------------------------------------------------------------
    # Eq. (11)
    #------------------------------------------------------------------------    
    #
    # B∂O∂Γ⋅gΓ
    #
    BOΓg = zeros(mesh.length∂O)
    LinearAlgebra.mul!(BOΓg, EL.B∂O∂Γ, gΓ)

    u∂O = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(mesh.length∂O))
    invB∂O∂O = similar(EL.B∂O∂O)
    invB∂O∂O = inv(EL.B∂O∂O)
    #
    #  u∂O = -B⁻¹∂O∂O⋅(B∂OΓ⋅gΓ)
    #
    LinearAlgebra.mul!(u∂O, -invB∂O∂O, BOΓg)

    #
    #  AIo,Γ
    #
    AIoΓ = similar(A, (mesh.lengthIo, mesh.lengthΓ))
    for iΓ = 1:mesh.lengthΓ
        g1=mesh.Γ[iΓ]        
        for io = 1:mesh.lengthIo
            io1 = mesh.Io[io]
            
            AIoΓ[io, iΓ] = A[io1, g1]
        end
    end

    #
    # Eq (12)
    #
    AIoΓg = similar(AIoΓ, (mesh.lengthIo))
    LinearAlgebra.mul!(AIoΓg, AIoΓ, gΓ)

    AIou∂O = similar(AIoΓg)
    LinearAlgebra.mul!(AIou∂O, EL.AIo∂O, u∂O)

    dims = (mesh.lengthIo)
    uIo = similar(u∂O, dims)
    LinearAlgebra.mul!(uIo, -invAIoIo, (AIou∂O + AIoΓg))

    for io = 1:mesh.lengthIo
        io1 = mesh.Io[io]
        u[io1] = uIo[io]
    end
    for io = 1:mesh.length∂O
        io1 = mesh.∂O[io]
        u[io1] = u∂O[io]
    end
    for io = 1:mesh.lengthΓ
        io1 = mesh.Γ[io]
        u[io1] = gΓ[io]
    end

end

#
    #= GLOBAL VERSION (eq 10)
    # BC = A⁻¹ᵢₒᵢₒ⋅Aᵢₒ∂τ
    dims = (mesh.lengthIo, mesh.length∂τ)
    BC = similar(EL.AIoIo, dims);
    LinearAlgebra.mul!(BC, invAIoIo, EL.AIo∂τ)

    # ABC = A∂Oᵢₒ⋅BC
    dims = (mesh.length∂O, mesh.length∂τ)
    ABC  = similar(EL.AIoIo, dims)
    LinearAlgebra.mul!(ABC, EL.A∂OIo, BC)
    
    EL.B∂O∂τ .= EL.A∂O∂τ .- ABC
    globalB∂O∂τ = copy(EL.B∂O∂τ)
    =#
    
