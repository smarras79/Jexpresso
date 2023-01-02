#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
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

include("../../io/mod_inputs.jl")
include("../../io/plotting/jeplots.jl")
include("../../io/print_matrix.jl")

include("../../kernel/AbstractTypes.jl")
include("../../kernel/basis/basis_structs.jl")
include("../../kernel/infrastructure/element_matrices.jl")
include("../../kernel/infrastructure/Kopriva_functions.jl")
include("../../kernel/infrastructure/2D_3D_structures.jl")
include("../../kernel/mesh/metric_terms.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/solver/mod_solution.jl")
include("../../kernel/timeIntegration/TimeIntegrators.jl")  
#--------------------------------------------------------
function driver(DT::CG,        #Space discretization type
                ET::Wave1D,    #Equation subtype
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
        QT_String = "Exact"
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
        QT_String = "Inexact"
        Qξ  = Nξ
        NDQ = ND
        ξq  = ξ
        ω   = ND.ξ.ω
    end 
    
    SD = NSD_1D()
     
    #--------------------------------------------------------
    # Build Lagrange polynomials:
    #
    # Return:
    # ψ     = basis.ψ[N+1, Q+1]
    # dψ/dξ = basis.dψ[N+1, Q+1]
    #--------------------------------------------------------
    basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)

    MT = COVAR() #Metric type: COVAR or CNVAR
    mestrics = build_metric_terms(SD, MT, mesh, basis, Nξ, Qξ, ξ, TFloat)
   
    
    #--------------------------------------------------------
    # Build element mass matrix
    #
    # Return:
    # el_mat.M[iel, i, j] <-- if exact (full)
    # el_mat.M[iel, i]    <-- if inexact (diagonal)
    # el_mat.D[iel, i, j] <-- either exact (full) OR inexact (sparse)
    #--------------------------------------------------------
    el_mat    = build_element_matrices!(SD, QT, basis.ψ, basis.dψ, ω, mesh, Nξ, Qξ, TFloat)
    
    #show(stdout, "text/plain", mesh.conn)
    (M, Minv) = DSS(SD, QT,      el_mat.M, mesh.conn, mesh.nelem, mesh.npoin, Nξ, TFloat)
    (D, ~)    = DSS(SD, Exact(), el_mat.D, mesh.conn, mesh.nelem, mesh.npoin, Nξ, TFloat)
    
    #Initialize q
    q = initialize(Wave1D(), mesh, inputs, TFloat)

    dq   = zeros(mesh.npoin);   
    qp   = copy(q.qn)
    
    #Plot I.C.
    plt1 = plot(mesh.x, q.qn, seriestype = :scatter,  title="Initial", reuse = false)
    display(plt1)
    
    Δt = inputs[:Δt]
    C = 0.25
    u = 2.0
    Δt = C*u*minimum(mesh.Δx)/mesh.nop
    Nt = floor((inputs[:tend] - inputs[:tinit])/Δt)
    
    #
    # ALGO 5.6 FROM GIRALDO: GLOBAL VERSION WITH SOLID-WALL B.C. AS A FIRST TEST
    #
    plt2 = scatter() #Clear plot
    
    RK = RK_Integrator{TFloat}(zeros(TFloat,5),zeros(TFloat,5),zeros(TFloat,5))
    buildRK5Integrator!(RK)
    for it = 1:Nt
        
        dq = zeros(mesh.npoin);
        qe = zeros(mesh.ngl);
        for s = 1:length(RK.a)
            
            #
            # RHS
            #
            rhs = build_rhs(SD, QT, Wave1D(), mesh, metrics, M, el_mat, u*qp)

            for I=1:mesh.npoin
                dq[I] = RK.a[s]*dq[I] + Δt*rhs[I]
                qp[I] = qp[I] + RK.b[s]*dq[I]
            end

            #
            # B.C.: solid wall
            #
            qp[1] = 0.0
            qp[mesh.npoin_linear] = 0.0

        end #stages

        title = string("Solution for N=", Nξ, " & ", QT_String, " integration")
        plt2 = scatter(mesh.x, qp,  title=title)
        display(plt2)
    end

end

function driver(DT::CG,       #Space discretization type
                ET::Adv2D,    #Equation subtype
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
    
    if lexact_integration
        #
        # Exact quadrature:
        # Quadrature order (Q = N+1) ≠ polynomial order (N)
        #
        QT  = Exact() #Quadrature Type
        QT_String = "Exact"
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
    #show(stdout, "text/plain", mesh.conn)
    #@printf("\n")
    #show(stdout, "text/plain", mesh.connijk)
    #@printf("\n")
    
    Me = build_mass_matrix!(SD, QT, TensorProduct(), basis.ψ, ω, mesh, metrics, Nξ, Qξ, TFloat)
    #show(stdout, "text/plain", Me[:,:,1:mesh.nelem])
    
    M =              DSSijk_mass(SD, QT, Me, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)
    #show(stdout, "text/plain", M)
    
    
    #Initialize q
    qp = initialize(Adv2D(), mesh, inputs, TFloat)
    
    dq   = zeros(mesh.npoin);
    qpel = copy(qp.qnel)
    RHS  = zeros(mesh.npoin);
    
    Δt = inputs[:Δt]
    C = 0.25
    #@info Δt = C*u*minimum(mesh.Δx)/mesh.nop
    # add a function to find the mesh mininum resolution
    #
    Δt = 0.025
    Nt = floor(Int64, (inputs[:tend] - inputs[:tinit])/Δt)
    #Nt = 500
    
    #
    # ALGO 5.6 FROM GIRALDO: GLOBAL VERSION WITH SOLID-WALL B.C. AS A FIRST TEST
    #
    RK = RK_Integrator{TFloat}(zeros(TFloat,5),zeros(TFloat,5),zeros(TFloat,5))
    buildRK5Integrator!(RK)
    for it = 1:Nt
        
        dq = zeros(mesh.npoin);
        qe = zeros(mesh.ngl);
        for s = 1:length(RK.a)
            
            #
            # rhs[ngl,ngl,nelem]
            #
            rhs_el = build_rhs(SD, QT, Adv2D(), qp, basis.ψ, basis.dψ, ω, mesh, metrics)
            #
            # RHS[npoin] = DSS(rhs)
            #
            RHS = DSSijk_rhs(SD, QT, rhs_el, mesh.connijk, mesh.nelem, mesh.npoin, Nξ, TFloat)
            RHS .= RHS./M
            for I=1:mesh.npoin
                dq[I] = RK.a[s]*dq[I] + Δt*RHS[I]
                qp.qn[I,1] = qp.qn[I,1] + RK.b[s]*dq[I]
            end
            
            #
            # B.C.: solid wall
            #
            #qp[1] = 0.0
            #qp[mesh.npoin_linear] = 0.0

        end #stages
        #@info size(RHS) size(mesh.x) size(mesh.y)
                
        #clf()
        #title = string(" RHS for N=", Nξ, " & ", QT_String, " integration")        
        #frhs = PyPlot.tricontourf(mesh.x, mesh.y, qp.qn[:,1], levels=30)
        #PyPlot.colorbar(frhs)        
        #plt[:show]()

        title = string(" RHS for N=", Nξ, " & ", QT_String, " integration")        
        display(PyPlot.tricontourf(mesh.x, mesh.y, qp.qn[:,1], levels=30))
        #PyPlot.colorbar(frhs)        
        #plt[:show]()
    end
    
    error("QUI AdvDiff/drivers.jl")
    
    return    
    
end
