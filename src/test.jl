#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using DifferentialEquations
using Revise

#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("./IO/mod_initialize.jl")
include("./IO/mod_inputs.jl")
include("./Mesh/mod_mesh.jl")
include("./solver/mod_solution.jl")
include("./basis/basis_structs.jl")
include("./Infrastructure/Kopriva_functions.jl")
include("./Infrastructure/2D_3D_structures.jl")
include("./element_matrices.jl")
include("./IO/plotting/jeplots.jl")
#--------------------------------------------------------

abstract type AbstractDiscretization end
struct CG <:  AbstractDiscretization end

abstract type AbstractProblem end
struct INTERPOLATION <: AbstractProblem end
struct INTEGRATION   <: AbstractProblem end


function test_driver(DT::CG,            #Discretization Type
                     SD::NSD_1D,        #Space Dimensions
                     PT::INTERPOLATION, #Problem Type
                     inputs::Dict,      #input parameters from src/user_input.jl
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
    ND = build_nodal_Storage([Nξ], LGL1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
    ξ  = ND.ξ.ξ
    qj = zeros(length(ξ))

    Nx = 101    
    x  = range(-1, 1, Nx)
    qe = zeros(Nx)
   
    
    #--------------------------------------------------------
    # Build Lagrange polynomials:
    #
    # Return:
    # ψ     = basis.ψ[N+1, Q+1]
    # dψ/dξ = basis.dψ[N+1, Q+1]
    #--------------------------------------------------------
    basis = build_Interpolation_basis!(LagrangeBasis(), SD, TFloat, ξ, x)

    #
    # Build exact function on a uniform grid of Nx points:
    #
    for ix = 1:length(x)
        x2 = x[ix]*x[ix]
        qe[ix] = 1.0/(1 + 50.0*x2)
        #qe[ix] = cos(π*x[ix]/2)
    end

    #
    # Calculate qⱼ at the LGL nodes:
    #
    for ilgl = 1:length(ξ)
        ξ2 = ξ[ilgl]*ξ[ilgl]
        qj[ilgl] = 1.0/(1 + 50.0*ξ2)
        #qj[ilgl] = cos(π*ξ[ilgl]/2)
    end

    #
    # Interpolate: q(x) = ∑{1,Nξ} ψⱼ(x)qⱼ
    #
    qh = zeros(Nx) #Interpolated data
    for ix = 1:length(qh)
        qh[ix] = 0
        for jlgl = 1:length(qj)
            qh[ix] = qh[ix] + basis.ψ[jlgl, ix]*qj[jlgl]
        end
    end

    #
    # Plot:
    #
    display(plot())
    plot!(x, qe, seriestype = :line,  title="Interpolation test", label="Exact")
    plot!(x, qh, seriestype = :line, label="qʰ(x)")
    display(plot!(ξ, qj, seriestype = :scatter, label="qⱼ(ξ)"))

    
end


function test_driver(DT::CG,            #Discretization Type
                     SD::NSD_1D,        #Space Dimensions
                     PT::INTEGRATION,   #Problem Type
                     inputs::Dict,      #input parameters from src/user_input.jl
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
    ND = build_nodal_Storage([Nξ], LGL1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
    ξ  = ND.ξ.ξ
    
    if lexact_integration
        #
        # Exact quadrature:
        # Quadrature order (Q = N+1) ≠ polynomial order (N)
        #
        QT  = Exact() #Quadrature Type
        Qξ  = Nξ + 1
        
        NDQ = build_nodal_Storage([Qξ], LGL1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
        ξq  = NDQ.ξ.ξ
        ω   = NDQ.ξ.ω
        
    else  
        #
        # Inexact quadrature:
        # Quadrature and interpolation orders coincide (Q = N)
        #
        QT  = Inexact() #Quadrature Type
        Qξ  = Nξ
        NDQ = ND
        ξq  = ξ
        ω   = ND.ξ.ω
    end
    
    #--------------------------------------------------------
    # Build Lagrange polynomials:
    #
    # Return:
    # ψ     = basis.ψ[N+1, Q+1]
    # dψ/dξ = basis.dψ[N+1, Q+1]
    #--------------------------------------------------------
    basis = build_Interpolation_basis!(LagrangeBasis(), SD, TFloat, ξ, ξq)

    
    
    display(plot())
    display(plot!(mesh.x, qp, seriestype = :scatter,  title="Final"))
    
end
