#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using DifferentialEquations
using Revise
using WriteVTK

#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../IO/mod_initialize.jl")
include("../IO/mod_inputs.jl")
include("../Mesh/mod_mesh.jl")
include("../solver/mod_solution.jl")
include("../basis/basis_structs.jl")
include("../Infrastructure/Kopriva_functions.jl")
include("../Infrastructure/2D_3D_structures.jl")
include("../element_matrices.jl")
include("../IO/plotting/jeplots.jl")
#--------------------------------------------------------


abstract type AbstractDiscretization end
struct CG <:  AbstractDiscretization end

abstract type AbstractProblem end
struct AD1D <: AbstractProblem end
struct NS1D <: AbstractProblem end
struct BURGERS1D <: AbstractProblem end

function driver(DT::CG,        #Space discretization type
                ET::AD1D,      #Equation subtype
                inputs::Dict,  #input parameters from src/user_input.jl
                TFloat) 
    
    N = inputs[:nop]
    lexact_integration = inputs[:lexact_integration]
    
    #--------------------------------------------------------
    # Create/read mesh
    # return mesh::St_mesh
    #--------------------------------------------------------
    mesh = mod_mesh_mesh_driver(inputs)
    
    #--------------------------------------------------------
    # Build interpolation nodes:
    #             the user decides among LGL, GL, etc. 
    # Return:
    # ξ = ND.ξ.ξ
    # ω = ND.ξ.ω
    #--------------------------------------------------------
    ND = build_nodal_Storage([N], LGL1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
    ξ  = ND.ξ.ξ
    
    if lexact_integration
        #
        # Exact quadrature:
        # Quadrature order (Q = N+1) ≠ polynomial order (N)
        #
        QT  = Exact() #Quadrature Type
        Q   = N + 1
        
        NDq = build_nodal_Storage([Q], LGL1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
        ξq  = NDq.ξ.ξ
        ω   = NDq.ξ.ω
        
    else  
        #
        # Inexact quadrature:
        # Quadrature and interpolation orders coincide (Q = N)
        #
        QT  = Inexact() #Quadrature Type
        Q   = N
        NDq = ND
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
    basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)

    
    #--------------------------------------------------------
    # Build element mass matrix
    #
    # Return:
    # el_mat.M[iel, i, j] <-- if exact (full)
    # el_mat.M[iel, i]    <-- if inexact (diagonal)
    # el_mat.D[iel, i, j] <-- either exact (full) OR inexact (sparse)
    #--------------------------------------------------------
    el_mat = build_element_matrices!(QT, basis.ψ, basis.dψ, ω, mesh, N, Q, TFloat)
    (M, Minv)= DSS(QT, el_mat.M, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)
    
    q = mod_initialize_initialize(mesh, inputs, TFloat)

    Δt = inputs[:Δt]
    Nt = floor((inputs[:tend] - inputs[:tinit])/Δt)
    #for it = 1:Nt

    rhs = build_rhs(AD1D(), mesh, el_mat, q.qn)
    RHS = DSSarray(rhs, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)
    
    #end
    
    
    
    
end

function build_rhs(PT::AD1D, mesh::St_mesh, el_mat, q)

    rhs = zeros(mesh.ngl^mesh.nsd, mesh.nelem)
    f   = zeros(mesh.ngl^mesh.nsd)
    u   = 2.0 #m/s

    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl
            
            ip = mesh.conn[i, iel]
            f[i] = u*q[ip]
            
        end
    end
    
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl

            #HERE IS WHERE the equation terms should come from a user defined tuple
            for j = 1:mesh.ngl
                rhs[i, iel] = -el_mat.D[i,j,iel]*f[j]
            end
        end
    end
    
    return rhs
    
end


