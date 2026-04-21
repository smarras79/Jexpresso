using Distributions
using ONNXRunTime

function driver(nparts,
                distribute,
                inputs,
                OUTPUT_DIR::String,
                TFloat)
    
    comm  = distribute.comm
    rank = MPI.Comm_rank(comm)
    
    if inputs[:lwarmup] == true

        if rank == 0 println(BLUE_FG(string(" # JIT pre-compilation of large problem ..."))) end
        
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
        
        if rank == 0 println(BLUE_FG(string(" # JIT pre-compilation of large problem ... DONE"))) end
    end
    #check_memory(" Before sem_setup.")
                
    sem, partitioned_model = sem_setup(inputs, nparts, distribute)
    
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
        
        if (inputs[:lamr] == true)
            amr_freq = inputs[:amr_freq]
            Δt_amr   = amr_freq * inputs[:Δt]
            tspan    = [TFloat(inputs[:tinit]), TFloat(inputs[:tinit] + Δt_amr)]
        else
            tspan = [TFloat(inputs[:tinit]), TFloat(inputs[:tend])]
        end

        if rank == 0 println(" # Params_setup ..................................") end
        
        params, u =  params_setup(sem,
                                  qp,
                                  inputs,
                                  OUTPUT_DIR,
                                  TFloat,
                                  tspan)
        
        if rank == 0 println(" # Params_setup .................................. DONE") end
    
        # test of projection matrix for solutions from old to new, i.e., coarse to fine, fine to coarse
        # test_projection_solutions(sem.mesh, qp, sem.partitioned_model, inputs, nparts, sem.distribute)
    
        if !inputs[:llinsolve]
            #-----------------------------------------------------------------------------------
            # Hyperbolic/parabolic problems that lead to Mdq/dt = RHS
            #-----------------------------------------------------------------------------------
            @time solution = time_loop!(inputs, params, u, partitioned_model)
            
            # PLOT NOTICE: Plotting is called from inside time_loop using callbacks.
            
        else
            #-----------------------------------------------------------------------------------
            # Problems that lead to Lx = RHS
            #-----------------------------------------------------------------------------------
            npoin          = sem.mesh.npoin
            nelem          = sem.mesh.nelem
            nelem_semi_inf = params.mesh.nelem_semi_inf
            ngl            = sem.mesh.ngl
            ngr            = sem.mesh.ngr
            
            RHS   = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(npoin))
            Mdiag = KernelAbstractions.zeros(inputs[:backend], TFloat, Int64(npoin))
            
            if (inputs[:backend] == CPU())

                if inputs[:lelementLearning]
                    if rank == 0 println(BLUE_FG(string(" # ALLOCATE FOR ELEMENT LEARNING ......."))) end

                    
                    nelintpoints = (ngl - 2)^2
                    nelpoints    = ngl^2
                    elnbdypoints = nelpoints - nelintpoints
                    
                    EL = @time allocate_elemLearning(nelem, ngl,
                                                     sem.mesh.length∂O,
                                                     sem.mesh.length∂τ,
                                                     sem.mesh.lengthΓ,
                                                     TFloat, inputs[:backend];
                                                     Nsamp=inputs[:Nsamp],
                                                     lEL_Sample=inputs[:lEL_Sample])

                    if rank == 0 println(BLUE_FG(string(" # ALLOCATE FOR ELEMENT LEARNING ....... DONE"))) end
                    
                    BOΓg        = zeros(sem.mesh.length∂O)
                    gΓ          = zeros(sem.mesh.lengthΓ)
                    lvtk_sample = false
                    
                    if EL.lEL_Sample
                        #-----------------------------------------------------
                        # 1. Sampling
                        #-----------------------------------------------------
                        bufferin  = Vector{Vector{Float64}}()
                        bufferout = Vector{Vector{Float64}}()
                        total_cols_writtenin  = 0
                        total_cols_writtenout = 0
                        
                        if isfile("input_tensor.csv");  rm("input_tensor.csv");  end
                        if isfile("output_tensor.csv"); rm("output_tensor.csv"); end

                        # ── Allocate ONCE outside the loop ────────────────────────────────────────
                        A       = sem.matrix.L
                        A_∂τ∂τ  = A[sem.mesh.∂τ, sem.mesh.∂τ]
                        avisc   = zeros(TFloat, 1, ngl^2)          # shape fixed, values change each iter
                        nfeatures = size(avisc, 2)

                        wbuf = EL_WorkBuffers(params.mesh, A, A_∂τ∂τ, nfeatures,
                                              nelintpoints, elnbdypoints,
                                              "./JX_NN_model.onnx")  # load_inference called ONCE here

                        for isamp = 1:inputs[:Nsamp]
                            println(" # --- sample = $isamp")

                            # avisc changes each sample — update values in-place, no reallocation
                            ranvisc      = 0.5 + rand()
                            avisc[1, :] .= ranvisc

                            for ip = 1:npoin
                                user_source!(RHS[ip], params.qp.qn[ip], params.qp.qe[ip],
                                             npoin, inputs[:CL], inputs[:SOL_VARS_TYPE];
                                             neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip],
                                             xmax=sem.mesh.xmax, xmin=sem.mesh.xmin,
                                             ymax=sem.mesh.ymax, ymin=sem.mesh.ymin)
                            end
                            RHS = sem.matrix.M .* RHS

                            apply_boundary_conditions_lin_solve!(sem.matrix.L,
                                                                 0.0, params.qp.qe,
                                                                 params.mesh.coords,
                                                                 params.metrics.nx,
                                                                 params.metrics.ny,
                                                                 params.metrics.nz,
                                                                 npoin,
                                                                 params.mesh.npoin_linear,
                                                                 params.mesh.poin_in_bdy_edge,
                                                                 params.mesh.poin_in_bdy_face,
                                                                 params.mesh.nedges_bdy,
                                                                 params.mesh.nfaces_bdy,
                                                                 ngl, ngr,
                                                                 nelem_semi_inf,
                                                                 params.basis.ψ, params.basis.dψ,
                                                                 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                                                                 RHS, 0.0, params.ubdy,
                                                                 params.mesh.connijk_lag,
                                                                 params.mesh.bdy_edge_in_elem,
                                                                 params.mesh.bdy_edge_type,
                                                                 params.ω, qp.neqs,
                                                                 params.inputs, params.AD, sem.mesh.SD)

                            # wbuf reused — no new allocations, no new ONNX sessions
                            elementLearning_Axb!(params.qp.qn, params.uaux, sem.mesh,
                                                 A, RHS, EL,
                                                 avisc,
                                                 bufferin, bufferout,
                                                 BOΓg, gΓ, wbuf;
                                                 isamp=isamp,
                                                 total_cols_writtenin=total_cols_writtenin,
                                                 total_cols_writtenout=total_cols_writtenout)
                        end # isamp loop

                        total_cols_writtenin  = flush_MLtensor!(bufferin,  total_cols_writtenin,  "input_tensor.csv")
                        total_cols_writtenout = flush_MLtensor!(bufferout, total_cols_writtenout, "output_tensor.csv")

                        if rank == 0 println(BLUE_FG(" # EL SAMPLING .......... DONE")) end
                        
                    else
                        #-----------------------------------------------------
                        # 2. Inference:
                        #-----------------------------------------------------
                        #
                        # L*q = M*RHS   See algo 12.18 of Giraldo's book
                        #
                        # 2.a/b
                        μ        = 1
                        #â        = zeros(TFloat, ngl, ngl)
                        avisc      = zeros(TFloat, 1, ngl^2)
                        avisc[1,:].= 0.5 + rand() #Uniform distribution between 0.5 and 1.5
                        nfeatures  = size(avisc, 2)
                        #ψ        = sem.basis.ψ
                        #expansion_2d!(â, ψ)
                        
                        for ip =1:npoin
                            RHS[ip] = user_source!(RHS[ip],
                                                   params.qp.qn[ip],
                                                   params.qp.qe[ip],
                                                   npoin,
                                                   inputs[:CL], inputs[:SOL_VARS_TYPE];
                                                   neqs=1, x=sem.mesh.x[ip], y=sem.mesh.y[ip],
                                                   xmax=sem.mesh.xmax, xmin=sem.mesh.xmin,
                                                   ymax=sem.mesh.ymax, ymin=sem.mesh.ymin)
                        end
                        RHS = sem.matrix.M.*RHS
                        
                        apply_boundary_conditions_lin_solve!(sem.matrix.L,
                                                             0.0, params.qp.qe,
                                                             params.mesh.coords,
                                                             params.metrics.nx,
                                                             params.metrics.ny,
                                                             params.metrics.nz,
                                                             npoin,
                                                             params.mesh.npoin_linear, 
                                                             params.mesh.poin_in_bdy_edge,
                                                             params.mesh.poin_in_bdy_face,
                                                             params.mesh.nedges_bdy,
                                                             params.mesh.nfaces_bdy,
                                                             ngl, ngr,
                                                             nelem_semi_inf,
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
                        nfeatures    = size(avisc, 2)
                        A            = sem.matrix.L
                        A_∂τ∂τ       = A[sem.mesh.∂τ, sem.mesh.∂τ]   # needed by EL_WorkBuffers constructor

                        wbuf = EL_WorkBuffers(params.mesh, A, A_∂τ∂τ, nfeatures,
                                              nelintpoints, elnbdypoints,
                                              "./JX_NN_model.onnx")
                        
                        total_cols_writtenin  = 0
                        total_cols_writtenout = 0
                        
                        println(GREEN_FG(string(" # INFERENCE: call to elementLearning_Axb! .......... ")))
                        elementLearning_Axb!(params.qp.qn, params.uaux, sem.mesh,
                                             A, RHS, EL,
                                             avisc,
                                             [0.0], [0.0],
                                             BOΓg, gΓ, wbuf;
                                             isamp=1,
                                             total_cols_writtenin=total_cols_writtenin,
                                             total_cols_writtenout=total_cols_writtenout)

                        println(GREEN_FG(string(" # INFERENCE: call to elementLearning_Axb! .......... DONE")))
                        usol = params.qp.qn
                        neqs = params.qp.neqs
                        args = (params.SD, usol, params.uaux, 1, 1,
                                sem.mesh, nothing,
                                nothing, nothing,
                                0.0, 0.0, 0.0,
                                OUTPUT_DIR, inputs,
                                params.qp.qvars,
                                params.qp.qoutvars,
                                inputs[:outformat])
                        
                        write_output(args...; nvar=neqs, qexact=params.qp.qe)
                        #-----------------------------------------------------
                        # END Element-learning infrastructure
                        #-----------------------------------------------------
                    end

                else
                    
                    #-----------------------------------------------------
                    # L*q = M*RHS   See algo 12.18 of Giraldo's book
                    #-----------------------------------------------------
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
                    
                    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage ..............")))
                    sol = @btime solveAx($sem.matrix.L, $RHS, inputs[:ode_solver])
                    println(YELLOW_FG(string(" # Solve x=inv(A)*b: sparse storage .............. DONE")))
                    
                    args = (params.SD, sol.u, params.uaux, 1, 1,
                            sem.mesh, nothing,
                            nothing, nothing,
                            0.0, 0.0, 0.0,
                            OUTPUT_DIR, inputs,
                            params.qp.qvars,
                            params.qp.qoutvars,
                            inputs[:outformat])
                    
                    write_output(args...; nvar=params.qp.neqs, qexact=params.qp.qe)
                end
            else
                println( " ")
                println( " WARNING!!! drivers.jl:L114")
                println( " WARNING: CHECK IF THIS GPU IMPLEMENTATION OF Ax=b still works")
                println( " ")
                nothing
            end
        end
    end
end

# Point evaluation: interpolate at a single point (ξ, η)
function expansion_2d!(a::Matrix, ψ::Matrix)
    
    # Tensor product form: ψᵀ * A * ψ
    return dot(ψ, a * ψ)
    
end
