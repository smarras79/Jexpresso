using Test

include("../AbstractProblems.jl")

include("../../kernel/abstractTypes.jl")
include("../../kernel/mesh/mesh.jl")
include("../../kernel/mesh/metric_terms.jl")
include("../../kernel/basis/basis_structs.jl")


function build_rhs(SD::NSD_2D, QT, AP::Adv2D, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics)

    qnel = zeros(mesh.ngl,mesh.ngl,mesh.nelem,3)
    
    rhs_el = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    #rhs_el = zeros(mesh.ngl*mesh.ngl,mesh.nelem)

     for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                m = mesh.connijk[i,j,iel]

                qnel[i,j,iel,1] = qp.qn[m,1]
                qnel[i,j,iel,2] = qp.qn[m,2]
                qnel[i,j,iel,3] = qp.qn[m,3]
                
            end
        end
     end
    
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                m = i + (j-1)*mesh.ngl
                                
                u  = qnel[i,j,iel,2]
                v  = qnel[i,j,iel,3]
                
                dqdξ = 0
                dqdη = 0
                for k = 1:mesh.ngl
                    dqdξ = dqdξ + dψ[k,i]*qnel[k,j,iel,1]
                    dqdη = dqdη + dψ[k,j]*qnel[i,k,iel,1]
                end
                dqdx = dqdξ*metrics.dξdx[i,j,iel] + dqdη*metrics.dηdx[i,j,iel]
                dqdy = dqdξ*metrics.dξdy[i,j,iel] + dqdη*metrics.dηdy[i,j,iel]

                #rhs_el[m, iel] = ω[i]*ω[j]*metrics.Je[i,j,iel]*(u*dqdx + v*dqdy)
                rhs_el[i,j,iel] = ω[i]*ω[j]*metrics.Je[i,j,iel]*(u*dqdx + v*dqdy)
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)

    return rhs_el
end



function build_rhs(SD::NSD_1D, QT::Inexact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, el_mat, f)

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

function build_rhs(SD::NSD_1D, QT::Exact, PT::Wave1D, mesh::St_mesh, metrics::St_metrics, M, el_mat, f)

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
