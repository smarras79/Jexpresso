const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../AbstractProblems.jl")
include("./rhs.jl")
include("./initialize.jl")
include("./timeLoop.jl")

include("../../io/mod_inputs.jl")
include("../../io/plotting/jeplots.jl")
include("../../io/print_matrix.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/basis/basis_structs.jl")
include("../../kernel/infrastructure/element_matrices.jl")
include("../../kernel/infrastructure/Kopriva_functions.jl")
include("../../kernel/infrastructure/2D_3D_structures.jl")
include("../../kernel/mesh/metric_terms.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/timeIntegration/TimeIntegrators.jl")  
include("../../kernel/boundaryconditions/BCs.jl")
#--------------------------------------------------------

function driver(DT::CG,       #Space discretization type
                PT::Adv2D,    #Equation subtype
                inputs::Dict, #input parameters from src/user_input.jl
                TFloat) 
    
    Nξ = inputs[:nop]
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
    if(mesh.nsd === 2)
        SD = NSD_2D()
    elseif(mesh.nsd === 3)
        SD = NSD_3D()
    else
        SD = NSD_1D()
    end
    
    #SM
    mesh.connijk[1,1,1] = 1
    mesh.connijk[2,1,1] = 2
    mesh.connijk[1,2,1] = 4
    mesh.connijk[2,2,1] = 5

    mesh.connijk[1,1,2] = 2
    mesh.connijk[2,1,2] = 3
    mesh.connijk[1,2,2] = 5
    mesh.connijk[2,2,2] = 6

    mesh.connijk[1,1,3] = 4
    mesh.connijk[2,1,3] = 5
    mesh.connijk[1,2,3] = 7
    mesh.connijk[2,2,3] = 8

    mesh.connijk[1,1,4] = 5
    mesh.connijk[2,1,4] = 6
    mesh.connijk[1,2,4] = 8
    mesh.connijk[2,2,4] = 9

    
    mesh.x[1], mesh.y[1] = -1.0, -1.0
    mesh.x[2], mesh.y[2] =  0.0, -1.0
    mesh.x[3], mesh.y[3] =  1.0, -1.0

    mesh.x[4], mesh.y[4] = -1.0, 0.0
    mesh.x[5], mesh.y[5] =  0.0, 0.0
    mesh.x[6], mesh.y[6] =  1.0, 0.0
    
    mesh.x[7], mesh.y[7] = -1.0, 1.0
    mesh.x[8], mesh.y[8] =  0.0, 1.0
    mesh.x[9], mesh.y[9] =  1.0, 1.0

    
    #END SM
    
    #--------------------------------------------------------
    println( " # Build INTERPOLATION points ξ and their weights ω ")
    ND = build_nodal_Storage([Nξ], LGL_1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
    ξ  = ND.ξ.ξ

    if lexact_integration
        #
        # Exact quadrature:
        # Quadrature order (Q = N+1) ≠ polynomial order (N)
        #
        QT  = Exact() #Quadrature Type
        QT_String = "Exact"
        Qξ  = Nξ + 1

        
        println( " # Build QUADRATURE points ξ and their weights ω ")
        NDQ = build_nodal_Storage([Qξ], LGL_1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
        ξq  = NDQ.ξ.ξ
        ω   = NDQ.ξ.ω
        
    else  
        #
        # Inexact quadrature:
        # Quadrature and interpolation orders coincide (Q = N)
        #
        QT  = Inexact() #Quadrature Type
        QT_String = "Inexact"
        Qξ  = Nξ
        NDQ = ND
        ξq  = ξ
        ω   = ND.ξ.ω
    end 
    SD = NSD_2D()
    
    #--------------------------------------------------------
    # Build Lagrange polynomials:
    #
    # Return:
    # ψ     = basis.ψ[N+1, Q+1]
    # dψ/dξ = basis.dψ[N+1, Q+1]
    #--------------------------------------------------------
    basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)

    #--------------------------------------------------------
    # Build metric terms
    #
    # Return:
    # dxdξ,dη[1:Q+1, 1:Q+1, 1:nelem]
    # dydξ,dη[1:Q+1, 1:Q+1, 1:nelem]
    # dzdξ,dη[1:Q+1, 1:Q+1, 1:nelem]
    # dξdx,dy[1:Q+1, 1:Q+1, 1:nelem]
    # dηdx,dy[1:Q+1, 1:Q+1, 1:nelem]
    #      Je[1:Q+1, 1:Q+1, 1:nelem]
    #--------------------------------------------------------
    metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, TFloat)
    
    #--------------------------------------------------------
    # Build element mass matrix
    #
    # Return:
    # M[1:N+1, 1:N+1, 1:N+1, 1:N+1, 1:nelem]
    #--------------------------------------------------------
    Me = build_mass_matrix(SD, QT, TensorProduct(), basis.ψ, ω, mesh, metrics, Nξ, Qξ, TFloat)
    #show(stdout, "text/plain", Me[:,:,1])
    
    M = DSSijk_mass(SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)
    #show(stdout, "text/plain", M)
    
    Le = build_laplace_matrix(SD, TensorProduct(), basis.ψ, basis.dψ, ω, mesh, metrics, Nξ, Qξ, TFloat)
    L = DSSijk_laplace(SD,  Le, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)
    #show(stdout, "text/plain", L)
    
    #Initialize q
    qp = initialize(Adv2D(), mesh, inputs, TFloat)
        
    Δt = inputs[:Δt]
    CFL = Δt/(abs(maximum(mesh.x) - minimum(mesh.x)/10/mesh.nop))
    println(" # CFL = ", CFL)    
    Nt = floor(Int64, (inputs[:tend] - inputs[:tinit])/Δt)
    
    # add a function to find the mesh mininum resolution
    TD = RK5()
    time_loop(TD, SD, QT, PT, mesh, metrics, basis, ω, qp, M, Le, Nt, Δt, inputs, TFloat)

    #Plot final solution
    jcontour(mesh.x, mesh.y, qp.qn[:,1], "Final solution at t=2π: tracer")
    
    return    
    
end
