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
include("../basis/basis_structs.jl")
include("../IO/mod_initialize.jl")
include("../IO/mod_inputs.jl")
include("../IO/plotting/jeplots.jl")
include("../IO/print_matrix.jl")
include("../Infrastructure/element_matrices.jl")
include("../Infrastructure/Kopriva_functions.jl")
include("../Infrastructure/2D_3D_structures.jl")
include("../Mesh/metric_terms.jl")
include("../Mesh/mesh.jl")
include("../solver/mod_solution.jl")
include("../TimeIntegration/TimeIntegrators.jl")  
#--------------------------------------------------------


abstract type AbstractDiscretization end
struct CG <:  AbstractDiscretization end

abstract type AbstractProblem end
struct Wave1D <: AbstractProblem end
struct AD1D <: AbstractProblem end
struct NS1D <: AbstractProblem end
struct BURGERS1D <: AbstractProblem end
struct Adv2D <: AbstractProblem end


abstract type AbstractBC end
struct PERIODIC1D_CG <: AbstractBC end

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
    basis = build_Interpolation_basis!(LagrangeBasis(), ξ, ξq, TFloat)
    
    mestrics = build_metric_terms(SD, mesh, basis, Nξ, Qξ, ξ, TFloat)
    #if (mesh.nsd > 1)
    #    error("drivers.jl TEMPORARY STOP WHILE TESTING 2D/3D grids.")
    #end
    
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
    q = mod_initialize_initialize(mesh, inputs, TFloat)

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
            rhs = drivers_build_rhs(SD, QT, Wave1D(), mesh, M, el_mat, u*qp)

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

#=
#
# ALGO 5.5 FROM GIRALDO: GLOBAL VERSION WITH SOLID-WALL B.C. AS A FIRST TEST
#
Dstar = copy(D)
if QT == Exact()
Dstar = Minv*D
else
for I=1:mesh.npoin
Dstar[I,I] = Minv[I]*D[I,I]
end
end
    @info size(M) size(Minv) size(Dstar)
    for it = 1:Nt  
        dq = zeros(mesh.npoin)
        for s = 1:length(RK.a)
            rhs = zeros(mesh.npoin)
            for iel=1:mesh.nelem
                for i=1:mesh.ngl
                    I = mesh.conn[i,iel]
                    for j=1:mesh.ngl
                        J = mesh.conn[j,iel]
                        rhs[I] = rhs[I] - u*Dstar[I,J]*qp[J]
                    end
                end
            end
            rhs[1] = 0.0
            rhs[mesh.npoin_linear] = 0.0
            
            for I=1:mesh.npoin
                dq[I] = RK.a[s]*dq[I] + Δt*rhs[I]
                qp[I] = qp[I] + RK.b[s]*dq[I]
            end
 
            #
            # B.C. Solid wall
            #
            qp[1] = 0.0
            qp[mesh.npoin_linear] = 0.0

        end #s
        
        display(plot())
        display(plot!(mesh.x, qp, seriestype = :scatter,  title="solution"))
    end
=#

end

###SM
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
   error(" 2d qui") 
    #Initialize q
    q = mod_initialize_initialize(mesh, inputs, TFloat)
    return
    
    dq   = zeros(mesh.npoin);   
    qp   = copy(q.qn)
    
    
    #Plot I.C.
    plt1 = plot(mesh.x, q.qn, seriestype = :scatter,  title="Initial", reuse = false)
    display(plt1)

    return
    
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
            rhs = drivers_build_rhs(SD, QT, Wave1D(), mesh, M, el_mat, u*qp)

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

#=
#
# ALGO 5.5 FROM GIRALDO: GLOBAL VERSION WITH SOLID-WALL B.C. AS A FIRST TEST
#
Dstar = copy(D)
if QT == Exact()
Dstar = Minv*D
else
for I=1:mesh.npoin
Dstar[I,I] = Minv[I]*D[I,I]
end
end
    @info size(M) size(Minv) size(Dstar)
    for it = 1:Nt  
        dq = zeros(mesh.npoin)
        for s = 1:length(RK.a)
            rhs = zeros(mesh.npoin)
            for iel=1:mesh.nelem
                for i=1:mesh.ngl
                    I = mesh.conn[i,iel]
                    for j=1:mesh.ngl
                        J = mesh.conn[j,iel]
                        rhs[I] = rhs[I] - u*Dstar[I,J]*qp[J]
                    end
                end
            end
            rhs[1] = 0.0
            rhs[mesh.npoin_linear] = 0.0
            
            for I=1:mesh.npoin
                dq[I] = RK.a[s]*dq[I] + Δt*rhs[I]
                qp[I] = qp[I] + RK.b[s]*dq[I]
            end
 
            #
            # B.C. Solid wall
            #
            qp[1] = 0.0
            qp[mesh.npoin_linear] = 0.0

        end #s
        
        display(plot())
        display(plot!(mesh.x, qp, seriestype = :scatter,  title="solution"))
    end
=#

end

###END SM

function drivers_build_rhs(SD::NSD_1D, QT::Inexact, PT::Wave1D, mesh::St_mesh, M, el_mat, f)

    #
    # Linear RHS in flux form: f = u*q
    #
    
    rhs = zeros(mesh.npoin)
    fe  = zeros(mesh.ngl)
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            fe[i] = f[I]
        end
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            for j=1:mesh.ngl
                rhs[I] = rhs[I] - el_mat.D[i,j,iel]*fe[j]
            end
        end
    end

    # M⁻¹*rhs where M is diagonal
    rhs .= rhs./M
    
    return rhs
end

function drivers_build_rhs(SD::NSD_1D, QT::Exact, PT::Wave1D, mesh::St_mesh, M, el_mat, f)

    #
    # Linear RHS in flux form: f = u*q
    #
    
    rhs = zeros(mesh.npoin)
    fe  = zeros(mesh.ngl)
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            fe[i] = f[I]
        end
        for i=1:mesh.ngl
            I = mesh.conn[i,iel]
            for j=1:mesh.ngl
                rhs[I] = rhs[I] - el_mat.D[i,j,iel]*fe[j]
            end
        end
    end
    
    # M⁻¹*rhs where M is not full
    rhs = M\rhs
    
    return rhs
end
