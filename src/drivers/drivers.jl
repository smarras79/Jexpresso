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
include("../IO/print_matrix.jl")
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
    
    if lexact_integration
        #
        # Exact quadrature:
        # Quadrature order (Q = N+1) ≠ polynomial order (N)
        #
        QT  = Exact() #Quadrature Type
        Qξ  = Nξ + 1
        
        NDQ = build_nodal_Storage([Qξ], LGL_1D(), NodalGalerkin()) # --> ξ <- ND.ξ.ξ
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
    
    
    if (mesh.nsd == 1)
        SD = NSD_1D()
    elseif (mesh.nsd == 2)
        SD = NSD_2D()        
    elseif (mesh.nsd == 3)
        SD = NSD_3D()
    end
    
    
 
    #--------------------------------------------------------
    # Build Lagrange polynomials:
    #
    # Return:
    # ψ     = basis.ψ[N+1, Q+1]
    # dψ/dξ = basis.dψ[N+1, Q+1]
    #--------------------------------------------------------
    basis = build_Interpolation_basis!(LagrangeBasis(), SD, TFloat, ξ, ξq)
    
    #periodicity flag array
    periodicity = zeros(Int64, mesh.npoin)
    for iel = 1:mesh.nelem
        for i = 1:mesh.ngl
            ip = mesh.conn[i, iel]
            periodicity[ip]=ip
        end
    end
    periodicity[mesh.npoin_linear]=1
    
    #--------------------------------------------------------
    # Build element mass matrix
    #
    # Return:
    # el_mat.M[iel, i, j] <-- if exact (full)
    # el_mat.M[iel, i]    <-- if inexact (diagonal)
    # el_mat.D[iel, i, j] <-- either exact (full) OR inexact (sparse)
    #--------------------------------------------------------
    el_mat    = build_element_matrices!(QT, basis.ψ, basis.dψ, ω, mesh, Nξ, Qξ, TFloat)

    #show(stdout, "text/plain", mesh.conn)
    
    (M, Minv) = DSS(QT, el_mat.M, periodicity, mesh.conn, mesh.nelem, mesh.npoin, Nξ, TFloat)
    (D, Dinv) = DSS(Exact(), el_mat.D, periodicity, mesh.conn, mesh.nelem, mesh.npoin, Nξ, TFloat)

    #show(stdout, "text/plain", el_mat.D)
    #@show "-----"
    #show(stdout, "text/plain", D)
    
    #initial condition --> q.qn
    q         = mod_initialize_initialize(mesh, inputs, TFloat)
    
    display(plot())
    display(plot!(mesh.x, q.qn, seriestype = :scatter,  title="Initial"))
    
    Δt = inputs[:Δt]
    C = 0.1
    u = 2.0
    Δt = C*u*minimum(mesh.Δx)/mesh.nop
    Nt = floor((inputs[:tend] - inputs[:tinit])/Δt)
    
    plt = scatter() #Clear plot
    #display(scatter(mesh.x, q.qn))

    RKA = [(0), 
           (-567301805773) / (1357537059087), 
           (-2404267990393) / (2016746695238), 
           (-3550918686646) / (2091501179385), 
           (-1275806237668) / (842570457699 )];

    RKB = [(1432997174477) / (9575080441755 ),
           (5161836677717) / (13612068292357),
           (1720146321549) / (2090206949498 ),
           (3134564353537) / (4481467310338 ),
           (2277821191437) / (14882151754819)];

    RKC = [(0),
           (1432997174477) / (9575080441755),
           (2526269341429) / (6820363962896),
           (2006345519317) / (3224310063776),
           (2802321613138) / (2924317926251)];

    qel  = zeros(mesh.ngl, mesh.nelem);
    Rel  = zeros(mesh.ngl, mesh.nelem);
    R    = zeros(mesh.npoin);
    dq   = zeros(mesh.npoin);   
    qp   = copy(q.qn)
    
    for iel = 1:mesh.nelem
        for l=1:mesh.ngl
            ip         = mesh.conn[l, iel]                
            x          = mesh.x[ip]
            qel[l,iel] = exp(-64.0*x*x)
        end
    end
    
    t    = 0.0
    q0 = copy(q.qn)
    for it = 1:Nt
        #@show it, Δt
        t = t + Δt

        dq = zeros(mesh.npoin)
        
        #for s = 1:length(RKA)

        #Create RHS Matrix
        for iel=1:mesh.nelem
            Rel  = zeros(mesh.ngl, mesh.nelem)
            for i=1:mesh.ngl
                I = mesh.conn[i, iel]                
                for j=1:mesh.ngl
                    J = mesh.conn[j, iel]
                    #qel[j,iel] = q0[J]
                    qel[j,iel] = qp[J]
                    Rel[i,iel] = Rel[i,iel] - el_mat.D[j,i,iel]*(u*qel[j,iel])
                end
            end
        end
        
        #DSS R
        R = zeros(mesh.npoin);
        for iel=1:mesh.nelem
            for i=1:mesh.ngl
                I = mesh.conn[i,iel]
                R[I] = R[I] + Rel[i,iel]
            end
        end
        
        for iel=1:mesh.nelem
            for i=1:mesh.ngl
                I = mesh.conn[i, iel]
                #if (I != 1 && I != mesh.npoin_linear)
                #    R[I] = Minv[I]*R[I]
                #else
                #    R[I] = 0.0
                #end
            end
        end
        
        #=  for I=1:mesh.npoin
        if (I != 1 && I != mesh.npoin_linear)
        dq[I] = RKA[s]*dq[I] + Δt*R[I]
        qp[I] = qp[I] + RKB[s]*dq[I]
        end
        end=#
        #for iel=1:mesh.nelem
        #   for i=1:mesh.ngl
        for I=1:mesh.npoin
            #      I = mesh.conn[i, iel]
            qp[I] = qp[I] + Δt*R[I]
        end
        #end
        #qp[1] = 0
        #qp[mesh.npoin_linear] = 0
        #for I=1:mesh.npoin
        #    if (periodicity[mesh.npoin_linear] == periodicity[I])
        #        qp[I] = qp[1] #periodicity
        #    end
    #end
        
        q0 = copy(qp)
        
        display(plot())
        display(plot!(mesh.x, qp, seriestype = :scatter,  title="solution"))
        
#    end #s
    end
end

function drivers_build_rhs(PT::AD1D, mesh::St_mesh, el_mat, q, periodicity)
    
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
                RHS[ip] = RHS[ip] + el_mat.D[j,i,iel]*f[i]
            end
        end
    end

    #Zero-out the RHS row corresponding to the periodic node
    for ip=1:mesh.npoin
        #if (periodicity[ip] == 1)
        if( mesh.x[ip] > 0.999999)
            @show ip mesh.x[ip]
            RHS[ip] = 0
        end
    end
    
    return RHS  
end
