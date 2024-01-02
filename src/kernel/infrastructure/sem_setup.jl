include("../mesh/restructure_for_periodicity.jl")
include("../mesh/warping.jl")

function sem_setup(inputs::Dict)
    fx = zeros(Float64,1,1)
    fy = zeros(Float64,1,1)
    fy_lag = zeros(Float64,1,1)
    Nξ = inputs[:nop]
    lexact_integration = inputs[:lexact_integration]    
    PT    = inputs[:equations]
    AD    = inputs[:AD]
    
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
        mesh.x .= (mesh.x .+ inputs[:xdisp]) .*inputs[:xscale]*0.5
    elseif (inputs[:xscale] != 1.0)
        mesh.x = mesh.x*inputs[:xscale]*0.5
    elseif (inputs[:xdisp] != 0.0)
        mesh.x .= (mesh.x .+ inputs[:xdisp])
    end
    if (inputs[:yscale] != 1.0 && inputs[:ydisp] != 0.0)
        mesh.y .= (mesh.y .+ inputs[:ydisp]) .*inputs[:yscale] * 0.5
    elseif(inputs[:yscale] != 1.0)
        mesh.y .= (mesh.y) .*inputs[:yscale]*0.5
    elseif(inputs[:ydisp] != 0.0)
        mesh.y .= (mesh.y .+ inputs[:ydisp])
    end
    mesh.ymax = maximum(mesh.y)
    #@info "xmax, ymax", maximum(mesh.x), maximum(mesh.y)    
    #warp_mesh!(mesh,inputs)    
    #--------------------------------------------------------
    # Build interpolation and quadrature points/weights
    #--------------------------------------------------------
    ξω  = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.nop)    
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
        if ("Laguerre" in mesh.bdy_edge_type[:])
            basis1 = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
            ξω2 = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta])
            ξ2,ω2 = ξω2.ξ, ξω2.ω
            basis2 = build_Interpolation_basis!(ScaledLaguerreBasis(), ξ2, ξ2, inputs[:laguerre_beta], TFloat)
            basis = (basis1, basis2)
            ω1 = ω
            ω = (ω1,ω2)
            if (inputs[:lfilter])
                ξω3 = basis_structs_ξ_ω!(inputs[:interpolation_nodes], mesh.ngr-1)
                ξ3,ω3 = ξω3.ξ, ξω3.ω
                fx = init_filter(mesh.ngl-1,ξ,inputs[:mu_x],mesh,inputs)
                fy = init_filter(mesh.ngl-1,ξ,inputs[:mu_y],mesh,inputs)
                #fy_lag = init_filter(mesh.ngr-1,ξ3,inputs[:mu_y],mesh,inputs)
                fy_lag = init_filter(mesh.ngr-1,ξ2,inputs[:mu_y],mesh,inputs)
            end
            @time periodicity_restructure!(mesh,inputs)
            if (inputs[:lwarp])
                warp_mesh!(mesh,inputs)
            end
            metrics1 = build_metric_terms(SD, COVAR(), mesh, basis1, Nξ, Qξ, ξ, ω1, TFloat)
            metrics2 = build_metric_terms(SD, COVAR(), mesh, basis1, basis2, Nξ, Qξ, mesh.ngr, mesh.ngr, ξ, ω1, ω2, TFloat)
            metrics = (metrics1, metrics2)
            
            matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
        else
            
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
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
            @time metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, ω, TFloat)
            
            @info " periodicity_restructure!"
            @time periodicity_restructure!(mesh,inputs)
            
            #warp_mesh!(mesh,inputs)
            matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
        end
    else 
        if(inputs[:llaguerre_1d])

            basis1 = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
            ξω2 = basis_structs_ξ_ω!(LGR(), mesh.ngr-1,inputs[:laguerre_beta])
            ξ2,ω2 = ξω2.ξ, ξω2.ω
            basis2 = build_Interpolation_basis!(ScaledLaguerreBasis(), ξ2, ξ2, inputs[:laguerre_beta], TFloat)
            basis = (basis1, basis2)
            ω1 = ω
            ω = (ω1,ω2)
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            metrics1 = build_metric_terms(SD, COVAR(), mesh, basis[1], Nξ, Qξ, ξ, ω, TFloat)
            metrics2 = build_metric_terms_1D_Laguerre(SD, COVAR(), mesh, basis[2], mesh.ngr, mesh.ngr, ξ2, ω2, inputs, TFloat)
            metrics = (metrics1, metrics2) 
            matrix = matrix_wrapper_laguerre(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
        else
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
            ω1 = ω
            ω = ω1
            #--------------------------------------------------------
            # Build metric terms
            #--------------------------------------------------------
            metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, ω, TFloat)

            if (inputs[:lperiodic_1d])
                periodicity_restructure!(mesh,inputs)
            end
            matrix = matrix_wrapper(AD, SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
        end
    end
    #--------------------------------------------------------
    # Build matrices
    #--------------------------------------------------------

    return (; QT, PT, mesh, metrics, basis, ω, matrix, fx, fy, fy_lag)
end
