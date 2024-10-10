include("../mesh/restructure_for_periodicity.jl")
include("../mesh/warping.jl")

function sem_setup(inputs::Dict)
    
    fx = zeros(Float64,1,1)
    fy = zeros(Float64,1,1)
    fy_lag = zeros(Float64,1,1)
    Nξ    = inputs[:nop]
    lexact_integration = inputs[:lexact_integration]    
    PT    = inputs[:equations]
    AD    = inputs[:AD]
    CL    = inputs[:CL]
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
    mesh = mod_mesh_mesh_driver(inputs)
    
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
            #@time periodicity_restructure!(mesh,inputs)
            if (inputs[:lwarp])
                warp_mesh!(mesh,inputs)
            end
            metrics1 = build_metric_terms(SD, COVAR(), mesh, basis1, Nξ, Qξ, ξ, ω1, TFloat; backend = inputs[:backend])
            metrics2 = build_metric_terms(SD, COVAR(), mesh, basis1, basis2, Nξ, Qξ, mesh.ngr, mesh.ngr, ξ, ω1, ω2, TFloat; backend = inputs[:backend])
            metrics = (metrics1, metrics2)
            @time periodicity_restructure!(mesh,inputs,inputs[:backend])
            
            matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat;
                                             ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend])
            
        else
            
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])
            ω1 = ω
            ω = ω1
            if (inputs[:lfilter])
                fx = init_filter(mesh.ngl-1,ξ,inputs[:mu_x],mesh,inputs)
                fy = init_filter(mesh.ngl-1,ξ,inputs[:mu_y],mesh,inputs)
            end
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            if (inputs[:lwarp])
                warp_mesh!(mesh,inputs)
            end
            @info " metrics"
            @time metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, ω, TFloat; backend = inputs[:backend])
            @time periodicity_restructure!(mesh,inputs,inputs[:backend]) 
            #warp_mesh!(mesh,inputs)
            matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend])
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
            metrics1 = build_metric_terms(SD, COVAR(), mesh, basis[1], Nξ, Qξ, ξ, ω, TFloat;backend = inputs[:backend])
            metrics2 = build_metric_terms_1D_Laguerre(SD, COVAR(), mesh, basis[2], mesh.ngr, mesh.ngr, ξ2, ω2, inputs, TFloat;backend = inputs[:backend])
            metrics = (metrics1, metrics2) 
            matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend])
        else
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat, inputs[:backend])

            ω1 = ω
            ω = ω1
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, ω, TFloat; backend = inputs[:backend])

            if (inputs[:lperiodic_1d])
                periodicity_restructure!(mesh,inputs,inputs[:backend])
            end
            matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation], backend = inputs[:backend])
        end
    end

    #--------------------------------------------------------
    # Build matrices
    #--------------------------------------------------------
    
    return (; QT, PT, CL, AD, SOL_VARS_TYPE, mesh, metrics, basis, ω, matrix, fx, fy, fy_lag)
end
