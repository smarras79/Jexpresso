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


abstract type AbstractBC end
struct PERIODIC1D_CG <: AbstractBC end

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
    C = 0.1
    u = 2.0
    Δt = 0.25 #C*u*minimum(mesh.Δx)
    Nt = floor((inputs[:tend] - inputs[:tinit])/Δt)

    plt = scatter() #Clear plot
    display(scatter(mesh.x, q.qn))
    rhs = zeros(mesh.ngl^mesh.nsd, mesh.nelem)

    #periodicity flag array
    periodicity = zeros(Int64, mesh.npoin)
    for ip=1:mesh.npoin
        periodicity[ip]=ip
    end
    periodicity[end]=1
    
    for it = 1:Nt
        @show it, Δt
        #rhs = drivers_build_rhs(AD1D(), mesh, el_mat, q.qn)
        #RHS = DSSarray(rhs, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)
        
        RHS = zeros(mesh.npoin)
        #S1
        for iel = 1:mesh.nelem
            for i = 1:mesh.ngl
                ip = mesh.conn[i, iel]
                for j = 1:mesh.ngl
                    #rhs[i, iel] = -el_mat.D[i,j,iel]*u*q.qn[ip]
                    RHS[ip] = RHS[ip] + el_mat.D[i,j,iel]*u*q.qn[ip]
                end
            end
        end
        q.qnp1 = q.qn + Δt*RHS
        
        display(scatter())
        display(scatter!(mesh.x, q.qn))
        
        q.qn = q.qnp1
        #RHS1 = DSSarray(rhs, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)
        #=qs1 = q.qn + Δt*RHS1
        

        #S2
        for iel = 1:mesh.nelem
            for i = 1:mesh.ngl
                ip = mesh.conn[i, iel]
                
                for j = 1:mesh.ngl
                    rhs[i, iel] = -el_mat.D[i,j,iel]*u*qs1[ip]
                end
            end
        end
        RHS2 = DSSarray(rhs, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)
        qs2 = (3/4)*q.qn + (1/4)*qs1 + (1/4)*Δt*RHS1

        #S3
        for iel = 1:mesh.nelem
            for i = 1:mesh.ngl
                ip = mesh.conn[i, iel]
                
                for j = 1:mesh.ngl
                    rhs[i, iel] = -el_mat.D[i,j,iel]*u*qs2[ip]
                end
            end
        end
        RHS2 = DSSarray(rhs, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)
        q.qnp1 = (1/4)*q.qn + (2/3)*qs2 + (2/3)*Δt*RHS2

        #Periodic b.c.
        q.qn[end] = q.qn[1]
        
        #=q1 = q.qn + Δt*RHS
        rhs1 =  drivers_build_rhs(AD1D(), mesh, el_mat, q1)
        RHS1 = DSSarray(rhs1, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)
        
        q2 = 3/4*q.qn + 1/4*q1 + 1/4*Δt*RHS1
        rhs2 = drivers_build_rhs(AD1D(), mesh, el_mat, q1)
        RHS2 = DSSarray(rhs1, mesh.conn, mesh.nelem, mesh.npoin, N, TFloat)

        qnp1 = 1/3*q.qn + 2/3*q2 + 2/3*Δt*RHS2=#
        =#
        #Periodic b.c.
        q.qn[end] = q.qn[1]
        
    end
    
    
end

function drivers_build_rhs(PT::AD1D, mesh::St_mesh, el_mat, q)

    rhs = zeros(mesh.ngl^mesh.nsd, mesh.nelem)
    RHS = zeros(mesh.npoin)
    f   = zeros(mesh.ngl^mesh.nsd)
    u   = 2.0 #m/s
    
    
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.conn[i, iel]
            f[i] = u*q[ip]
        end
        
        for i = 1:mesh.ngl
            ip = mesh.conn[i, iel]

            for j = 1:mesh.ngl
                #rhs[i, iel] = -el_mat.D[i,j,iel]*f[j]
                RHS[periodicity[ip]] = RHS[periodicity[ip]] + el_mat.D[i,j,iel]*f[j]
            end
        end
    end
    
    return rhs, RHS  
end

function drivers_apply_bc(PT::PERIODIC1D_CG, qn::Array)

    qn[end] = qn[1]
    
end
