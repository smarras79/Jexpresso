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
#include("../AbstractProblems.jl")

include("../../io/mod_inputs.jl")
include("../../io/write_output.jl")
include("../../io/print_matrix.jl")

include("../abstractTypes.jl")
include("../bases/basis_structs.jl")
include("../boundaryconditions/BCs.jl")
include("../globalStructs.jl")
include("../operators/rhs.jl")
include("../solvers/Axb.jl")
include("./element_matrices.jl")
include("./Kopriva_functions.jl")
include("./2D_3D_structures.jl")


#struct SEM{}

#--------------------------------------------------------
function SEMsetup(DT::ContGal,       #Space discretization type
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
    
    #Build L = DSS(∫∇ψᵢ∇ψⱼdΩₑ)
    Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh, metrics, Nξ, Qξ, TFloat)
    L  = DSS_laplace(SD, Le, mesh, TFloat)

    #Build De = ∫ψᵢ∇ψⱼdΩₑ
    De = build_differentiation_matrix(SD, basis.ψ, basis.dψ, ω, mesh, Nξ, Qξ, TFloat)
    
    #Build M = DSS(∫ψᵢψⱼdΩₑ)
    Me = build_mass_matrix(SD, QT, basis.ψ,   ω, mesh, metrics, Nξ, Qξ, TFloat)
    M  = DSS_mass(SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)

    #SEM = (;M, L, basis, metrics, mesh)
    return M, Me, De, Le, L, basis, ω, metrics, mesh, SD, QT
end
