include("../mesh/restructure_for_periodicity.jl")

function sem_setup(inputs::Dict)
    
    Nξ = inputs[:nop] 
    PT = inputs[:equations]
    lexact_integration = inputs[:lexact_integration]
    
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
    @info " --- BASES"
    @time basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
    
    #--------------------------------------------------------
    # Build metric terms
    #--------------------------------------------------------
    @info " --- METRICS"
    @time metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, TFloat)
    
    periodicity_restructure!(mesh,inputs)
    
    #--------------------------------------------------------
    # Build matrices
    #--------------------------------------------------------
    @info " --- MATRICES"
    @time matrix = matrix_wrapper(SD, QT, basis, ω, mesh, metrics, Nξ, Qξ, TFloat; ldss_laplace=inputs[:ldss_laplace], ldss_differentiation=inputs[:ldss_differentiation])
    
    return (; QT, PT, mesh, metrics, basis, ω, matrix)
end
