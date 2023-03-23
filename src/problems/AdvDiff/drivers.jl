#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using PrettyTables
using Revise

#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../AbstractProblems.jl")

include("./initialize.jl")

include("../../kernel/operators/rhs.jl")


include("../../io/mod_inputs.jl")
include("../../io/plotting/jeplots.jl")
include("../../io/print_matrix.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/globalStructs.jl")
include("../../kernel/bases/basis_structs.jl")
include("../../kernel/infrastructure/element_matrices.jl")
include("../../kernel/infrastructure/Kopriva_functions.jl")
include("../../kernel/infrastructure/2D_3D_structures.jl")
include("../../kernel/mesh/metric_terms.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/solvers/TimeIntegrators.jl")  
include("../../kernel/mesh/restructure_for_periodicity.jl")
include("../../kernel/boundaryconditions/BCs.jl")
#--------------------------------------------------------
function driver(DT::ContGal,       #Space discretization type
                inputs::Dict,      #input parameters from src/user_input.jl
                OUTPUT_DIR::String,
                TFloat) 

    Nξ = inputs[:nop]
    lexact_integration = inputs[:lexact_integration]    
    PT    = inputs[:problem]
    neqns = inputs[:neqns]
    
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
    if (mesh.nsd == 1)
        SD = NSD_1D()
    elseif (mesh.nsd == 2)
        SD = NSD_2D()
    elseif (mesh.nsd == 3)
        SD = NSD_3D()
    else
        error(" Drivers.jl: Number of space dimnnsions unknow! CHECK Your grid!")
    end
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
    #--------------------------------------------------------
    metrics = build_metric_terms(SD, COVAR(), mesh, basis, Nξ, Qξ, ξ, TFloat)
        
    periodicity_restructure!(mesh,inputs)    
    #--------------------------------------------------------
    # Build element mass matrix
    #
    # Return:
    # M[1:N+1, 1:N+1, 1:N+1, 1:N+1, 1:nelem]
    #--------------------------------------------------------    
    Le =         build_laplace_matrix(SD,     basis.ψ, basis.dψ, ω, mesh, metrics, Nξ, Qξ, TFloat)
    De = build_differentiation_matrix(SD,     basis.ψ, basis.dψ, ω, mesh,          Nξ, Qξ, TFloat)
    Me =            build_mass_matrix(SD, QT, basis.ψ,           ω, mesh, metrics, Nξ, Qξ, TFloat)
    M  =                     DSS_mass(SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)
   
    #--------------------------------------------------------
    # Initialize q
    #--------------------------------------------------------
    qp = initialize(SD, PT, mesh, inputs, OUTPUT_DIR, TFloat)
    write_vtk(qp.qn[:,1], SD, mesh, OUTPUT_DIR,inputs)
    
    Δt = inputs[:Δt]
    CFL = Δt/(abs(maximum(mesh.x) - minimum(mesh.x)/10/mesh.nop))
    println(" # CFL = ", CFL)    
    Nt = floor(Int64, (inputs[:tend] - inputs[:tinit])/Δt)
    
    # NOTICE add a function to find the mesh mininum resolution
    solution = time_loop!(SD, QT, PT, mesh, metrics, basis, ω, qp, M, De, Le, Nt, Δt, neqns, inputs, AdvDiff_Circ(), OUTPUT_DIR, TFloat)

    #Out-to-file:
    write_output(solution, SD, mesh, OUTPUT_DIR, inputs, inputs[:outformat])
    
end
