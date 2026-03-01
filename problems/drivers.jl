using Distributions

function driver(nparts,
                distribute,
                inputs::Dict,
                OUTPUT_DIR::String,
                TFloat)
    
    comm  = distribute.comm
    rank = MPI.Comm_rank(comm)
    
    if inputs[:lwarmup] == true
        if rank == 0
            println(BLUE_FG(string(" # JIT pre-compilation of large problem ...")))
	    end
        input_mesh = inputs[:gmsh_filename]
        inputs[:gmsh_filename] = inputs[:gmsh_filename_c]
        sem_dummy = sem_setup(inputs, nparts, distribute)
        inputs[:gmsh_filename] = input_mesh
        
        #check_memory(" Right after sem_dummy setup.")
        
        # --- MEMORY CLEANUP ---
        # 1. Explicitly drop the reference to the dummy object
        sem_dummy = nothing 

        # 2. Force a full garbage collection run
        #GC.gc()
        
        #check_memory(" At GC() after sem_dummy setup.")
        
        if rank == 0
            println(BLUE_FG(string(" # JIT pre-compilation of large problem ... END")))
        end
    end
    #check_memory(" Before sem_setup.")
                
    sem, partitioned_model = sem_setup(inputs, nparts, distribute)

    #check_memory(" After sem_setup.")
    
    if (inputs[:backend] != CPU())
        convert_mesh_arrays!(sem.mesh.SD, sem.mesh, inputs[:backend], inputs)
    end
    if (inputs[:lRT_problem])
        if (sem.mesh.SD == NSD_2D())
            build_radiative_transfer_problem(sem.mesh, inputs, 1, sem.mesh.ngl, sem.basis.dψ, sem.basis.ψ, sem.ω, sem.metrics.Je, 
                                     sem.metrics.dξdx, sem.metrics.dξdy, sem.metrics.dηdx, sem.metrics.dηdy, 
                                     sem.metrics.nx, sem.metrics.ny, sem.mesh.elem_to_edge, sem.mesh.extra_mesh, sem.QT, NSD_2D(), sem.AD)
        else
            κ = zeros(sem.mesh.npoin)
            σ = zeros(sem.mesh.npoin)
            if (inputs[:lRT_from_data])
                @info "reading atmospheric data to build extinction and scattering coefficients"
                filename = inputs[:RT_data_file]
                data = read_atmospheric_data(filename)
                data_interp = interpolate_atmosphere_to_mesh(data,sem.mesh)
                κ, σ = atmos_to_rad(data_interp,sem.mesh.npoin)
            end

            build_radiative_transfer_problem(sem.mesh, inputs, 1, sem.mesh.ngl, sem.basis.dψ, sem.basis.ψ, sem.ω, sem.metrics.Je,
                                     sem.metrics.dξdx, sem.metrics.dξdy, sem.metrics.dξdz, 
                                     sem.metrics.dηdx, sem.metrics.dηdy, sem.metrics.dηdz,
                                     sem.metrics.dζdx, sem.metrics.dζdy, sem.metrics.dζdz,
                                     sem.metrics.nx, sem.metrics.ny, sem.metrics.nz, 
                                     sem.mesh.elem_to_face, sem.mesh.extra_mesh, κ, σ, sem.QT, NSD_3D(), sem.AD)
        end
    else
        qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)

        #check_memory(" After initialize.")
        # println("Rank $rank: $(Sys.free_memory() / 2^30) GB free")
        # test of projection matrix for solutions from old to new, i.e., coarse to fine, fine to coarse
        # test_projection_solutions(sem.mesh, qp, sem.partitioned_model, inputs, nparts, sem.distribute)

        if (inputs[:lamr] == true)
            amr_freq = inputs[:amr_freq]
            Δt_amr   = amr_freq * inputs[:Δt]
            tspan    = [TFloat(inputs[:tinit]), TFloat(inputs[:tinit] + Δt_amr)]
        else
            tspan = [TFloat(inputs[:tinit]), TFloat(inputs[:tend])]
        end

        if rank == 0
            @info " Params_setup .................................."
        end
        params, u =  params_setup(sem,
                              qp,
                              inputs,
                              OUTPUT_DIR,
                              TFloat,
                              tspan)
    
        if rank == 0
            @info " Params_setup .................................. END"
        end
    
        # test of projection matrix for solutions from old to new, i.e., coarse to fine, fine to coarse
        # test_projection_solutions(sem.mesh, qp, sem.partitioned_model, inputs, nparts, sem.distribute)
    
        if !inputs[:llinsolve]
            #
            # Hyperbolic/parabolic problems that lead to Mdq/dt = RHS
            #
            @time solution = time_loop!(inputs, params, u, partitioned_model)
            # PLOT NOTICE: Plotting is called from inside time_loop using callbacks.
        
        else
            #
            # Problems that lead to Lx = RHS
            #
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
                â             = zeros(TFloat, sem.mesh.ngl, sem.mesh.ngl)
                avisc         = zeros(TFloat, sem.mesh.ngl^2)
                avisc[:]     .= 0.5 + rand() #Uniform distribution between 0.5 and 1.5
                avisc[:]     .= ranvisc
                
                expansion_2d(â, ψ)
                @mystop
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

function elementLearning_Axb!(u, uaux, mesh::St_mesh,
                              A, ubdy, EL,
                              avisc, 
                              bufferin, bufferout;
                              isamp=1,
                              total_cols_writtenin=0,
                              total_cols_writtenout=0)
    
    mesh.lengthO =  mesh.length∂O +  mesh.lengthIo
    
    #=EL = allocate_elemLearning(mesh.nelem, mesh.ngl,
                               mesh.length∂O,
                               mesh.length∂τ,
                               mesh.lengthΓ,
                               TFloat, inputs[:backend];
                               Nsamp=inputs[:Nsamp])
    =#
    
    nelintpoints = (mesh.ngl-2)*(mesh.ngl-2)
    nelpoints    = size(mesh.conn)[2]
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
    # inv(AiIoIo)
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

    u∂O      = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(mesh.length∂O))
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

    for io = 1:mesh.lengthIo # all internal without edges
        io1 = mesh.Io[io]
        u[io1] = uIo[io]
    end
    for io = 1:mesh.length∂O # all skeleton no boundaries
        io1 = mesh.∂O[io]
        u[io1] = u∂O[io]
    end
    for io = 1:mesh.lengthΓ # # all boundaries
        io1 = mesh.Γ[io]
        u[io1] = gΓ[io]
    end
    #@info u
    #@mystop
    # ∂O U Γ
    skeletonAndbdy                                              = zeros(Int64, length(mesh.∂τ))
    skeletonAndbdy[1:mesh.length∂O]                            .= mesh.∂O[:]
    skeletonAndbdy[mesh.length∂O+1:mesh.length∂O+mesh.lengthΓ] .= mesh.Γ[:]

    
#    #for iel=1:mesh.nelem
#        for iedge = 1:4
#            @info iel, mesh.edge2pedge[iedge]
#        end
#    end
#    @mystop
#    bdy_edge_in_elem
#    edge2pedge

    #
    # uvb ⊂ u∂τ  SHUKAI meeting
    #
#=
    @info "mesh.cell_face_ids "
    @info mesh.cell_face_ids, size(mesh.cell_face_ids)
    @info "mesh.facet_cell_ids"
    @info mesh.facet_cell_ids, size(mesh.facet_cell_ids)

    for iedge=1:mesh.nedges
        @info " iedge ", iedge, "belongs to element ",  mesh.facet_cell_ids[iedge]
    end
    for iedge=1:mesh.nedges
        @info " edge ", iedge, " belong to elem ", mesh.edge_in_elem[iedge]
    end

    @info "-----"
    =#
    uvb = zeros(Float64, mesh.nelem, elnbdypoints)
    u∂τ = zeros(Float64, length(mesh.∂τ))
    for iskel = 1:length(mesh.∂τ)
        is = skeletonAndbdy[iskel]
#        @info iskel, is
        u∂τ[iskel] = u[is]
     #x   @info iskel, is, u∂τ[iskel]
    end

    
     for iel=1:mesh.nelem
         #
         # 
         #
         ii = 1
         #
         # Aᵥₒᵥb
         #
         for isk = 1:length(mesh.∂τ)
             ipsk = skeletonAndbdy[isk]
             
             for ibdy = 1:elnbdypoints
                 #jpb = mesh.conn[iel, j]
                 uvb[iel, ibdy] = u∂τ[ipsk]
                # @info iel, ibdy, uvb[iel, ibdy]
             end
             ii += 1
         end
     end

    
    #=for iel = 1:mesh.nelem
        for ibdyel = 1:
        for j1=1:length(mesh.∂τ)
            jτ1 = skeletonAndbdy[j1]
            
            uvb[iel, ibdyel] = u∂τ[jτ1]
        end
    end=#

    
    
    #
    # ML: input/outpute tensors to use in training (?):
    #
    # 1. Set B∂τ∂τ := A∂τ∂τ
    #    
    T2  = zeros(size(EL.Avovo)[1], size(EL.Avovb)[2])
    T1  = zeros(size(EL.Avovb)[2], size(EL.Avovb)[2])
    Bie = similar(T2)
    
    # 2.c
    EL.input_tensor[:, isamp] .= avisc[:]

    # 2.d        
    for iel = 1:1 #mesh.nelem
        
        Avbvo = transpose(EL.Avovb[:,:,iel])
        
        # T2 = -A⁻¹ᵥₒᵥₒ[:,:,iel]⋅Avovb[:,:,iel]
        LinearAlgebra.mul!(T2, -inv(EL.Avovo[:,:,iel]), EL.Avovb[:,:,iel])
        Bie .= -T2
        
        # T1 = Avbvo[:,:,iel]⋅T2 = - Avbvo⋅A⁻¹ᵥₒᵥₒ⋅Avovb
        LinearAlgebra.mul!(@view(T1[:,:]), @view(Avbvo[:,:]), @view(T2[:,:]))
        
        # 2.e
        # Output tensor:
        EL.output_tensor[:, isamp] .= vec(Bie)  # Bie = -T2ie
    end

    #------------------------------------------------------------------------
    # Write input/output_bufferin.csv
    #------------------------------------------------------------------------
    write_MLtensor(@view(EL.input_tensor[:, isamp]), bufferin, total_cols_writtenin, "input_tensor.csv")
    write_MLtensor(@view(EL.output_tensor[:, isamp]), bufferout, total_cols_writtenout, "output_tensor.csv")
    #------------------------------------------------------------------------
    # END write input/output_buffer.csv
    #------------------------------------------------------------------------
    
    #!!!!!  ML code GOES HERE !!!!!
    
end

# Point evaluation: interpolate at a single point (ξ, η)
function expansion_2d(a::Matrix, ψ::Vector)
    
    # Tensor product form: ψᵀ * A * ψ
    return dot(ψ, a * ψ)
    
end
