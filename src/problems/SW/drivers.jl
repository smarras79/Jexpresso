#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using PrettyTables
using Revise
using WriteVTK

#Constants
const TInt   = Int64
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
include("../../kernel/bases/basis_structs.jl")
include("../../kernel/infrastructure/element_matrices.jl")
include("../../kernel/infrastructure/Kopriva_functions.jl")
include("../../kernel/infrastructure/2D_3D_structures.jl")
include("../../kernel/mesh/metric_terms.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/solver/mod_solution.jl")
include("../../kernel/timeIntegration/TimeIntegrators.jl")  
include("../../kernel/boundaryconditions/BCs.jl")
#--------------------------------------------------------

function driver(DT::CG,       #Space discretization type
                PT::SW,       #Equation subtype
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
    
    #--------------------------------------------------------
    ND = build_nodal_Storage([Nξ], LGL_1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
    ξ  = ND.ξ.ξ
    
    #
    # Inexact quadrature:
    # Quadrature and interpolation orders coincide (Q = N)
    #
    QT  = Inexact() #Quadrature Type
    Qξ  = Nξ
    NDQ = ND
    ξq  = ξ
    ω   = ND.ξ.ω
    SD  = NSD_2D()
    
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
    
    #--------------------------------------------------------
    # Build element mass matrix
    #--------------------------------------------------------    
    Me = build_mass_matrix!(SD, TensorProduct(), basis.ψ, ω, mesh, metrics, Nξ, Qξ, TFloat)
    M  = DSSijk_mass(SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)

    #--------------------------------------------------------
    #Initialize q
    #--------------------------------------------------------
    qp = initialize(PT, mesh, inputs, TFloat)
    
    Δt  = inputs[:Δt]
    CFL = Δt/(abs(maximum(mesh.x) - minimum(mesh.x)/10/mesh.nop))
    Nt  = floor(Int64, (inputs[:tend] - inputs[:tinit])/Δt)
    println(" # CFL = ", CFL)
        
    # add a function to find the mesh mininum resolution
    TD = RK5()
    time_loop!(TD, SD, QT, PT, mesh, metrics, basis, ω, qp, M, Nt, Δt, inputs, TFloat)

    #Plot final solution
    title = string("Final solution at t=inputs[:tend] for tracer")
    jcontour(mesh.x, mesh.y, qp.qn[:,1], title)
    
    return    
    
end
