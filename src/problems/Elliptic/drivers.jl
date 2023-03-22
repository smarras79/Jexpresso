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

include("./rhs.jl")
include("./initialize.jl")

include("../../io/mod_inputs.jl")
include("../../io/write_output.jl")
include("../../io/print_matrix.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/globalStructs.jl")
include("../../kernel/bases/basis_structs.jl")
include("../../kernel/infrastructure/element_matrices.jl")
include("../../kernel/infrastructure/Kopriva_functions.jl")
include("../../kernel/infrastructure/2D_3D_structures.jl")
include("../../kernel/solvers/Axb.jl")
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
    
    #Build L = DSS(∫∇ψᵢ∇ψⱼdΩₑ)
    Le = build_laplace_matrix(SD, basis.ψ, basis.dψ, ω, mesh, metrics, Nξ, Qξ, TFloat)
    L  =          DSS_laplace(SD, Le, mesh, TFloat)

    #Build M = DSS(∫ψᵢψⱼdΩₑ)
    Me =    build_mass_matrix(SD, QT, basis.ψ,           ω, mesh, metrics, Nξ, Qξ, TFloat)
    M  =             DSS_mass(SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)
    
    #--------------------------------------------------------
    # Initialize q
    #--------------------------------------------------------
    qp = define_q(SD, mesh.nelem, mesh.npoin, mesh.ngl, neqns, TFloat)

    #Build ∫S(q)dΩ
    RHS = build_rhs_source(SD, QT, inputs[:problem], qp.qn, mesh, M, TFloat)

    # Dirichlet B.C.
    # NOTICE these will be replaced with tbe general way of building B.C.
    # Yassine is working on it.
    ϵ = eps(Float32)
    #=
    if (occursin(lowercase(inputs[:gmsh_filename]), "circle"))
        rmax = (maximum(mesh.x) - minimum(mesh.x))/2
        for ip=1:mesh.npoin
            x, y = mesh.x[ip], mesh.y[ip]
            r = sqrt(x*x + y*y)
            if( (r > rmax - ϵ) || (r < -rmax + ϵ))
                qp.qn[ip,1] = 0.0
                for jp=1:mesh.npoin
                    L[ip,jp] = 0.0
                end
                L[ip,ip] = 1.0
            end
        end
    else
        for ip=1:mesh.npoin
            x, y = mesh.x[ip], mesh.y[ip]
            if( (x > 1.0 - ϵ) || (x < -1.0 + ϵ))
                qp.qn[ip,1] = sinpi(2*y)
                for jp=1:mesh.npoin
                    L[ip,jp] = 0.0
                end
                L[ip,ip] = 1.0
            end
            if( (y > 1.0 - ϵ) || (y < -1.0 + ϵ))
                qp.qn[ip,1] = 0.0
                for jp=1:mesh.npoin
                    L[ip,jp] = 0.0
                end
                L[ip,ip] = 1.0
            end        
        end
    en=#
    #END Dirichlet B.C..
    BCT = AdvDiff_Circ()
    apply_boundary_conditions!(SD, zeros(mesh.ngl,mesh.ngl,mesh.nelem), qp.qn, mesh, inputs, QT, metrics, basis.ψ, basis.dψ, ω, 0.0, BCT, neqns; L=L)
    
    println(" # Solve Lq=RHS ................................")    
    solution = solveAx(L, RHS, inputs[:ode_solver])
    println(" # Solve Lq=RHS ................................ DONE")

    #Out-to-file:
    write_vtk(solution.u, SD, mesh, OUTPUT_DIR, inputs)
    write_output(solution, SD, mesh, OUTPUT_DIR, inputs, inputs[:outformat])
    
end
