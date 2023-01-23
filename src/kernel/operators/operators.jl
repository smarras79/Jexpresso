function interpolate!(qh, qj, ψ)
    #
    # in/out: qh(x)  of size [1:Nx] where Nx is the number of interpolation points
    # in:     qj(ξ)  of size [1:Nξ] where Nξ = polynomial order + 1 (e.g. LGL points)
    # in:     ψ(ξ,x) Lagrange basis of size [1:Nξ, 1:Nx]
    #
    # Interpolation rule using Lagrange basis ψ with zeros at the ξ nodes:
    # q(x)     = ∑ⱼ{1,Nξ} ψⱼ(x)qⱼ
    # ∂q(x)/∂ξ = ∑ⱼ{1,Nξ} ∂ψⱼ(x)/∂x qⱼ
    #
    for ix = 1:length(qh)
        qh[ix] = 0
        for jlgl = 1:length(qj)
            qh[ix] = qh[ix] + ψ[jlgl, ix]*qj[jlgl]
        end
    end
    
    return qh
end

function build_gradient(SD::NSD_2D, QT::Inexact, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics)
    nvar = size(qp.qn,2)
    gradq = zeros(2,mesh.npoin,nvar)

    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                m = i + (j-1)*mesh.ngl
                for var = 1:nvar
                    dqdξ = 0
                    dqdη = 0
                    for k = 1:mesh.ngl
                        dqdξ = dqdξ + dψ[k,i]*qp.qn[ip,1]
                        dqdη = dqdη + dψ[k,j]*qp.qn[ip,1]
                    end
                    gradq[1,ip,var] = dqdξ*metrics.dξdx[i,j,iel] + dqdη*metrics.dηdx[i,j,iel]
                    gradq[2,ip,var] = dqdξ*metrics.dξdy[i,j,iel] + dqdη*metrics.dηdy[i,j,iel]
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)

    return gradq
end

function build_gradient(SD::NSD_3D, QT::Inexact, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics)
    nvar =size(qp.qn,2)
    gradq = zeros(3,mesh.npoin,nvar)
 
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                for k=1:mesh.ngl
                    ip = mesh.connijk[i,j,k,iel]
                    m = i + (j-1)*mesh.ngl
                    for var=1:nvar
                        dqdξ = 0.0
                        dqdη = 0.0
                        for l = 1:mesh.ngl
                            dqdξ = dqdξ + dψ[l,i]*qp.qn[ip,1]
                            dqdη = dqdη + dψ[l,j]*qp.qn[ip,1]
                            dqdζ = dqdζ + dψ[l,k]*qp.qn[ip,1]
                        end
                        gradq[1,ip,var] = dqdξ*metrics.dξdx[i,j,k,iel] + dqdη*metrics.dηdx[i,j,k,iel] + dqdζ*metrics.dζdx[i,j,k,iel]
                        gradq[2,ip,var] = dqdξ*metrics.dξdy[i,j,k,iel] + dqdη*metrics.dηdy[i,j,k,iel] + dqdζ*metrics.dζdy[i,j,k,iel]
                        gradq[3,ip,var] = dqdξ*metrics.dξdz[i,j,k,iel] + dqdη*metrics.dηdz[i,j,k,iel] + dqdζ*metrics.dζdz[i,j,k,iel]
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)

    return dqdx, dqdy, dqdz
end

function build_gradient(SD::NSD_1D, QT::Inexact, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics)
    nvar = size(qp.qn,2)
    gradq = zeros(1,mesh.npoin,nvar)

    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            ip = mesh.connijk[i,j,iel]
            m = i + (j-1)*mesh.ngl

            for var=1:nvar
                dqdξ = 0.0
                for k = 1:mesh.ngl
                    dqdξ = dqdξ + dψ[k,i]*qp.qn[ip,1]
                end
                gradq[ip,var] = dqdξ*metrics.dξdx[i,j,iel]
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)

    return dqdx
end
