include("../mesh/restructure_for_periodicity.jl")
include("../mesh/warping.jl")

function sem_setup(inputs::Dict, nparts, distribute, adapt_flags = nothing, partitioned_model_coarse = nothing, omesh = nothing)
    
    comm = distribute.comm
    rank = MPI.Comm_rank(comm)
    
    fx = zeros(Float64,1,1)
    fy = zeros(Float64,1,1)
    fz = zeros(Float64,1,1)
    fy_lag = zeros(Float64,1,1)
    Nξ    = inputs[:nop]
    lexact_integration = inputs[:lexact_integration]    
    PT    = inputs[:equations]
    AD    = inputs[:AD]
    CL    = inputs[:CL]
    phys_grid = zeros(Float64,1,1)
    SOL_VARS_TYPE = inputs[:SOL_VARS_TYPE]
    
    #--------------------------------------------------------
    # Create/read mesh
    # return mesh::St_mesh
    # and Build interpolation nodes
    #             the user decides among LGL, GL, etc. 
    # Return:
    # ξ = ND.ξ.ξ
    # ω = ND.ξ.ω
    #--------------------------------------------------------
    if isnothing(adapt_flags)
        mesh, partitioned_model = mod_mesh_mesh_driver(inputs, nparts, distribute)
    else
        mesh, partitioned_model, n2o_ele_map = mod_mesh_mesh_driver(inputs, nparts, distribute, adapt_flags, partitioned_model_coarse, omesh)
    end
    
    if (inputs[:xscale] != 1.0 && inputs[:xdisp] != 0.0)
        mesh.x .= (mesh.x .+ TFloat(inputs[:xdisp])) .*TFloat(inputs[:xscale]*0.5)
    elseif (inputs[:xscale] != 1.0)
        mesh.x = mesh.x*TFloat(inputs[:xscale]*0.5)
    elseif (inputs[:xdisp] != 0.0)
        mesh.x .= (mesh.x .+ TFloat(inputs[:xdisp]))
    end
    if (inputs[:yscale] != 1.0 && inputs[:ydisp] != 0.0)
        mesh.y .= (mesh.y .+ inputs[:ydisp]) .*inputs[:yscale] * 0.5
    elseif(inputs[:yscale] != 1.0)
        mesh.y .= (mesh.y) .*inputs[:yscale]*0.5
    elseif(inputs[:ydisp] != 0.0)
        mesh.y .= (mesh.y .+ inputs[:ydisp])
    end
    mesh.ymax = maximum(mesh.y)
    
    #--------------------------------------------------------
    # Build interpolation and quadrature points/weights
    #--------------------------------------------------------
    ξω  = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop, inputs[:backend])    
    interp, project = build_projection_1d(ξω.ξ)
    ξ,ω = ξω.ξ, ξω.ω    
    if lexact_integration
        #
        # Exact quadrature:
        # Quadrature order (Q = N+1) ≠ polynomial order (N)
        #
        QT  = Exact() #Quadrature Type
        QT_String = "Exact"
        Qξ  = Nξ + 1
        
        ξωQ   = basis_structs_ξ_ω!(inputs[:quadrature_nodes], Qξ)
        ξq, ω = ξωQ.ξ, ξωQ.ω
    else  
        #
        # Inexact quadrature:
        # Quadrature and interpolation orders coincide (Q = N)
        #
        QT  = Inexact() #Quadrature Type
        QT_String = "Inexact"
        Qξ  = Nξ
        ξωq = ξω
        ξq  = ξ
        ω   = ξω.ω
    end
    SD = mesh.SD
    
    #--------------------------------------------------------
    # Build Lagrange polynomials:
    #
    # Return:
    # ψ     = basis.ψ[N+1, Q+1]
    # dψ/dξ = basis.dψ[N+1, Q+1]
    #--------------------------------------------------------
    if (mesh.nsd > 1)
        #@info "grid size data" "xmax, ymax, xmin, ymin", maximum(mesh.x), maximum(mesh.y), minimum(mesh.x), minimum(mesh.y) 
        if ("Laguerre" in mesh.bdy_edge_type[:])
            basis1 = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            ξω2 = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta], inputs[:backend])
            ξ2,ω2 = ξω2.ξ, ξω2.ω
            basis2 = build_Interpolation_basis!(ScaledLaguerreBasis(), ξ2, ξ2, inputs[:laguerre_beta], TFloat, inputs[:backend])
            basis = (basis1, basis2)
            ω1 = ω
            ω = (ω1,ω2)
            if (inputs[:lfilter])
                ξω3 = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.ngr-1, inputs[:backend])
                ξ3,ω3 = ξω3.ξ, ξω3.ω
                if (inputs[:backend] == CPU())
                    fx = init_filter(mesh.ngl-1,ξ,inputs[:mu_x],mesh,inputs)
                    fy = init_filter(mesh.ngl-1,ξ,inputs[:mu_y],mesh,inputs)
                #fy_lag = init_filter(mesh.ngr-1,ξ3,inputs[:mu_y],mesh,inputs)
                    fy_lag = init_filter(mesh.ngr-1,ξ2,inputs[:mu_y],mesh,inputs)
                else
                    ξ_gl = KernelAbstractions.zeros(CPU(), Float64, Int64(mesh.ngl))
                    KernelAbstractions.copyto!(CPU(),ξ_gl,ξ)
                    ξ_gr = KernelAbstractions.zeros(CPU(), Float64, Int64(mesh.ngr))
                    KernelAbstractions.copyto!(CPU(),ξ_gr,ξ2)
                    fx_1 = init_filter(mesh.ngl-1,ξ_gl,inputs[:mu_x],mesh,inputs)
                    fy_1 = init_filter(mesh.ngl-1,ξ_gl,inputs[:mu_y],mesh,inputs)
                #fy_lag = init_filter(mesh.ngr-1,ξ3,inputs[:mu_y],mesh,inputs)
                    fy_lag_1 = init_filter(mesh.ngr-1,ξ_gr,inputs[:mu_y],mesh,inputs)
                    fx = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    fy = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    fy_lag = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngr), Int64(mesh.ngr))
                    KernelAbstractions.copyto!(inputs[:backend], fx, fx_1)
                    KernelAbstractions.copyto!(inputs[:backend], fy, fy_1)
                    KernelAbstractions.copyto!(inputs[:backend], fy_lag, fy_lag_1)
                end
            end
            
            if (inputs[:lwarp])
                warp_mesh!(mesh,inputs)
            end
            @info " Build metrics ......"
            @time metrics1 = build_metric_terms(SD, COVAR(), mesh, basis1, Nξ, Qξ, ξ, ω1, TFloat; backend = inputs[:backend])
            @time metrics2 = build_metric_terms(SD, COVAR(), mesh, basis1, basis2, Nξ, Qξ, mesh.ngr, mesh.ngr, ξ, ω1, ω2, TFloat; backend = inputs[:backend])
            metrics = (metrics1, metrics2)
            @info " Build metrics ...... DONE"
            @time periodicity_restructure!(mesh,mesh.x,mesh.y,mesh.z,mesh.xmax,
                                           mesh.xmin,mesh.ymax,mesh.ymin,mesh.zmax,mesh.zmin,mesh.poin_in_bdy_face,
                                           mesh.poin_in_bdy_edge,mesh.ngl,mesh.ngr,mesh.nelem,mesh.npoin,mesh.nsd,mesh.bdy_edge_type,
                                           mesh.bdy_face_type,mesh.bdy_face_in_elem,mesh.bdy_edge_in_elem,
                                           mesh.connijk,mesh.connijk_lag,mesh.npoin_linear,mesh.nelem_semi_inf,inputs,inputs[:backend])

            matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat;
                                             ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend])
            
        else
            @info " Build interpolation bases ......"
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            @info " Build interpolation bases ...... END"
            ω1 = ω
            ω = ω1
            if (inputs[:lfilter])
                if (inputs[:backend] == CPU())
                    fx = init_filter(mesh.ngl-1,ξ,inputs[:mu_x],mesh,inputs)
                    fy = init_filter(mesh.ngl-1,ξ,inputs[:mu_y],mesh,inputs)
                    if (mesh.nsd >2)
                        fz = init_filter(mesh.ngl-1,ξ,inputs[:mu_z],mesh,inputs)
                    end
                else
                    ξ_temp = KernelAbstractions.zeros(CPU(), Float64, Int64(mesh.ngl))
                    KernelAbstractions.copyto!(CPU(),ξ_temp,ξ)
                    fx_1 = init_filter(mesh.ngl-1,ξ_temp,inputs[:mu_x],mesh,inputs)
                    fy_1 = init_filter(mesh.ngl-1,ξ_temp,inputs[:mu_y],mesh,inputs)
                    fx = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    fy = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                    KernelAbstractions.copyto!(inputs[:backend], fx, fx_1)
                    KernelAbstractions.copyto!(inputs[:backend], fy, fy_1)
                    if (mesh.nsd > 2)
                        fz_1 = init_filter(mesh.ngl-1,ξ_temp,inputs[:mu_z],mesh,inputs)
                        fz = KernelAbstractions.allocate(inputs[:backend], TFloat, Int64(mesh.ngl), Int64(mesh.ngl))
                        KernelAbstractions.copyto!(inputs[:backend], fz, fz_1)
                    end
                end
            end
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            if (mesh.nsd > 2)
                if (inputs[:lwarp])
                    warp_mesh_3D!(mesh,inputs)
                end
            else
                if (inputs[:lwarp])
                    warp_mesh!(mesh,inputs)
                end
            end
            if rank == 0
                @info " Build metrics ......"
            end
            @time metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, ω, TFloat; backend = inputs[:backend])
            if rank == 0
                @info " Build metrics ...... END"
            end
            if (inputs[:lphysics_grid])
                phys_grid = init_phys_grid(mesh, inputs,inputs[:nlay_pg],inputs[:nx_pg],inputs[:ny_pg],mesh.xmin,mesh.xmax,mesh.ymin,mesh.ymax,mesh.zmin,mesh.zmax,inputs[:backend])
            end 
            if rank == 0
                @info " Build periodicity infrastructure ......"
            end
            @time periodicity_restructure!(mesh,mesh.x,mesh.y,mesh.z,mesh.xmax,
                                           mesh.xmin,mesh.ymax,mesh.ymin,mesh.zmax,mesh.zmin,mesh.poin_in_bdy_face,
                                           mesh.poin_in_bdy_edge,mesh.ngl,mesh.ngr,mesh.nelem,mesh.npoin,mesh.nsd,mesh.bdy_edge_type,
                                           mesh.bdy_face_type,mesh.bdy_face_in_elem,mesh.bdy_edge_in_elem,
                                           mesh.connijk,mesh.connijk_lag,mesh.npoin_linear,mesh.nelem_semi_inf,inputs,inputs[:backend])
            if rank == 0
                @info " Build periodicity infrastructure ...... DONE"
            end

#@mystop(" L 152 sem_setup")
            
            #warp_mesh!(mesh,inputs)
            
            if rank == 0
                @info " Matrix wrapper ......"
            end
            matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend], interp)
            if rank == 0
                @info " Matrix wrapper ...... END"
            end
            
        end
    else
        
        if(inputs[:llaguerre_1d_right] || inputs[:llaguerre_1d_left])

            basis1 = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            ξω2 = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta], inputs[:backend])
            ξ2,ω2 = ξω2.ξ, ξω2.ω
            basis2 = build_Interpolation_basis!(ScaledLaguerreBasis(), ξ2, ξ2, inputs[:laguerre_beta], TFloat, inputs[:backend])
            basis = (basis1, basis2)
            ω1 = ω
            ω = (ω1,ω2)
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
             @info " Build metrics ......"
             @time metrics1 = build_metric_terms(SD, COVAR(), mesh, basis[1], Nξ, Qξ, ξ, ω, TFloat;backend = inputs[:backend])
             @time metrics2 = build_metric_terms_1D_Laguerre(SD, COVAR(), mesh, basis[2], mesh.ngr, mesh.ngr, ξ2, ω2, inputs, TFloat;backend = inputs[:backend])
            metrics = (metrics1, metrics2)
             @info " Build metrics ...... DONE"
            matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend])
        else
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])

            ω1 = ω
            ω = ω1
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
             @time metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, ω, TFloat; backend = inputs[:backend])

            if (inputs[:lperiodic_1d])
                @time periodicity_restructure!(mesh,mesh.x,mesh.y,mesh.z,mesh.xmax,
                                           mesh.xmin,mesh.ymax,mesh.ymin,mesh.zmax,mesh.zmin,mesh.poin_in_bdy_face,
                                           mesh.poin_in_bdy_edge,mesh.ngl,mesh.ngr,mesh.nelem,mesh.npoin,mesh.nsd,mesh.bdy_edge_type,
                                           mesh.bdy_face_type,mesh.bdy_face_in_elem,mesh.bdy_edge_in_elem,
                                           mesh.connijk,mesh.connijk_lag,mesh.npoin_linear,mesh.nelem_semi_inf,inputs,inputs[:backend])
            end
            matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend])
        end
    end


    #--------------------------------------------------------
    # Build matrices
    #--------------------------------------------------------
    if isnothing(adapt_flags)
        return (; QT, PT, CL, AD, SOL_VARS_TYPE, mesh, metrics, basis, ω, matrix, fx, fy, fy_lag, fz, phys_grid, interp, project, partitioned_model, nparts, distribute)
    else
        return (; QT, PT, CL, AD, SOL_VARS_TYPE, mesh, metrics, basis, ω, matrix, fx, fy, fy_lag, fz, phys_grid, interp, project, partitioned_model, nparts, distribute), n2o_ele_map
    end
    
end
