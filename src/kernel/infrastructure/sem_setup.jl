include("../bases/basis_structs.jl")
include("../infrastructure/element_matrices.jl")
include("../mesh/restructure_for_periodicity.jl")
include("../mesh/warping.jl")

function sem_setup(inputs::Dict)
    fx = zeros(Float64,1,1)
    fy = zeros(Float64,1,1) 
    Nξ = inputs[:nop]
    lexact_integration = inputs[:lexact_integration]    
    PT    = inputs[:equations]
    
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
    mesh.x = mesh.x*100000.0
    mesh.y .= (mesh.y .+ 1) .*15000.0
    @info "xmax", maximum(mesh.x),maximum(mesh.y)
    mesh.ymax = maximum(mesh.y)
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
        
        ξωQ   = basis_structs_ξ_ω!(inputs[:quadrature_nodes], mesh.nop)
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
            ξω2 = basis_structs_ξ_ω!(LGR(), mesh.ngr-1)
            ξ2,ω2 = ξω2.ξ, ξω2.ω
            basis2 = build_Interpolation_basis!(ScaledLaguerreBasis(), ξ2, ξ2, TFloat)
            basis = (basis1, basis2)
            ω1 = ω
            ω = (ω1,ω2)
            metrics1 = build_metric_terms(SD, COVAR(), mesh, basis1, Nξ, Qξ, ξ, ω1, TFloat)
            metrics2 = build_metric_terms(SD, COVAR(), mesh, basis1, basis2, Nξ, Qξ, mesh.ngr, mesh.ngr, ξ, ω1, ω2, TFloat)
            metrics = (metrics1, metrics2)
    
            matrix = matrix_wrapper_laguerre(SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
        else
      
            basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
            ω1 = ω
            ω = ω1
            fx = init_filter(mesh.ngl-1,ξ,0.01)
            fy = init_filter(mesh.ngl-1,ξ,0.01) 
    #--------------------------------------------------------
    # Build metric terms
    #--------------------------------------------------------
           @info " metrics"
           @time metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, ω, TFloat)
    
           @info " periodicity_restructure!"
           @time periodicity_restructure!(mesh,inputs)
           #warp_mesh!(mesh,inputs)
            matrix = matrix_wrapper(SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
        end
    else 
        basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
        ω1 = ω
        ω = ω1
    #--------------------------------------------------------
    # Build metric terms
    #--------------------------------------------------------
        metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, TFloat)
    
        periodicity_restructure!(mesh,inputs)
     
        matrix = matrix_wrapper(SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
         
    end
    #--------------------------------------------------------
    # Build matrices
    #--------------------------------------------------------
    
    return (; QT, PT, mesh, metrics, basis, ω, matrix,fx,fy)
end
