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

function build_gradient(SD::NSD_2D, QT::Inexact, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics,gradq,nvars)
     
    for iel=1:mesh.nelem
        for i=1:mesh.ngl
            for j=1:mesh.ngl
                ip = mesh.connijk[i,j,iel]
                m = i + (j-1)*mesh.ngl
                for var = 1:nvars
                    dqdξ = 0
                    dqdη = 0
                    for k = 1:mesh.ngl
                        dqdξ = dqdξ + dψ[k,i]*qp[ip,1]
                        dqdη = dqdη + dψ[k,j]*qp[ip,1]
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

function build_laplacian(SD::NSD_2D, QT::Inexact, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics)
    nvar = size(qp.qn,2)
    L = zeros(mesh.ngl*mesh.ngl,mesh.ngl*mesh.ngl,mesh.nelem)
    Δq = zeros(mesh.ngl,mesh.ngl,mesh.nelem,nvar)    

    for iel=1:mesh.nelem
        for l=1:mesh.ngl
            for k=1:mesh.ngl #Integration points
                for j=1:mesh.ngl
                    for i=1:mesh.ngl
                        IJ = i + (j-1)*mesh.ngl
                        dψJKdx = dψ[i,k]*ψ[j,l]*metrics.dξdx[k,l,iel] + ψ[i,k]*dψ[j,l]*metrics.dηdx[k,l,iel]
                        dψJKdy = dψ[i,k]*ψ[j,l]*metrics.dξdy[k,l,iel] + ψ[i,k]*dψ[j,l]*metrics.dηdy[k,l,iel]
                        for n=1:mesh.ngl
                            for m=1:mesh.ngl
                                MN = m + (n-1)*mesh.ngl
                                dψIKdx = dψ[m,k]*ψ[n,l]*metrics.dξdx[k,l,iel] + ψ[m,k]*dψ[n,l]*metrics.dηdx[k,l,iel]
                                dψIKdy = dψ[m,k]*ψ[n,l]*metrics.dξdy[k,l,iel] + ψ[m,k]*dψ[n,l]*metrics.dηdy[k,l,iel]
                                L[MN,IJ,iel] += ω[k]*ω[l]*metrics.Je[k,l,iel]*(dψIKdx*dψJKdx + dψIKdy*dψJKdy)
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    for var=1:nvar
        for iel=1:mesh.nelem
            for j=1:mesh.ngl
                for i=1:mesh.ngl
                    IJ = i + (j-1)*mesh.ngl
                    ip = mesh.connijk[i,j,iel]
                    for n=1:mesh.ngl
                        for m=1:mesh.ngl
                            MN = m + (n-1)*mesh.ngl
                            Δq[i,j,iel,var] += qp.qn[ip,var]*L[MN,IJ,iel]
                        end
                    end
                end
            end
        end
    end
    
   return Δq
end

function build_laplacian(SD::NSD_3D, QT::Inexact, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics)
    nvar = size(qp.qn,2)
    L = zeros(mesh.ngl*mesh.ngl*mesh.ngl,mesh.ngl*mesh.ngl*mesh.ngl,mesh.nelem)
    Δq = zeros(mesh.ngl,mesh.ngl,mesh.ngl,mesh.nelem,nvar) 

    for iel=1:mesh.nelem
        for q=1:mesh.ngl
            for p=1:mesh.ngl #Integration points
                for o=1:mesh.ngl 
                    for n=1:mesh.ngl
                        for m=1:mesh.ngl
                            for l=1:mesh.ngl
                                IJ = n + mesh.ngl*(m-1 + (l-1)*mesh.ngl)
                                dψJKdx = dψ[n,q]*ψ[m,p]*ψ[l,o]*metrics.dξdx[q,p,o,iel] + ψ[n,q]*dψ[m,p]*ψ[l,o]metrics.dηdx[q,p,o,iel] + ψ[n,q]*ψ[m,p]*dψ[l,o]*metrics.dζdx[q,p,o,iel]
                                dψJKdy = dψ[n,q]*ψ[m,p]*ψ[l,o]*metrics.dξdy[q,p,o,iel] + ψ[n,q]*dψ[m,p]*ψ[l,o]metrics.dηdy[q,p,o,iel] + ψ[n,q]*ψ[m,p]*dψ[l,o]*metrics.dζdy[q,p,o,iel]
                                dψJKdz = dψ[n,q]*ψ[m,p]*ψ[l,o]*metrics.dξdz[q,p,o,iel] + ψ[n,q]*dψ[m,p]*ψ[l,o]metrics.dηdz[q,p,o,iel] + ψ[n,q]*ψ[m,p]*dψ[l,o]*metrics.dζdz[q,p,o,iel]
                                for i=1:mesh.ngl
                                    for j=1:mesh.ngl
                                        for k=1:mesh.ngl
                                            MN = i + mesh.ngl*(j-1 + (k-1)*mesh.ngl)
                                            dψIKdx = dψ[i,q]*ψ[j,p]*ψ[k,o]*metrics.dξdx[q,p,o,iel] + ψ[i,q]*dψ[j,p]*ψ[k,o]metrics.dηdx[q,p,o,iel] + ψ[i,q]*ψ[j,p]*dψ[k,o]*metrics.dζdx[q,p,o,iel]
                                            dψIKdy = dψ[i,q]*ψ[j,p]*ψ[k,o]*metrics.dξdy[q,p,o,iel] + ψ[i,q]*dψ[j,p]*ψ[k,o]metrics.dηdy[q,p,o,iel] + ψ[i,q]*ψ[j,p]*dψ[k,o]*metrics.dζdy[q,p,o,iel]
                                            dψIKdz = dψ[i,q]*ψ[j,p]*ψ[k,o]*metrics.dξdz[q,p,o,iel] + ψ[i,q]*dψ[j,p]*ψ[k,o]metrics.dηdz[q,p,o,iel] + ψ[i,q]*ψ[j,p]*dψ[k,o]*metrics.dζdz[q,p,o,iel]
                                            L[MN,IJ,iel] += ω[q]*ω[p]*ω[o]*metrics.Je[q,p,o,iel]*(dψIKdx*dψJKdx + dψIKdy*dψJKdy + dψIKdz*dψJKdz)
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    for var=1:nvar
        for iel=1:mesh.nelem
            for j=1:mesh.ngl
                for i=1:mesh.ngl
                    for k=1:mesh.ngl
                        IJ = i + mesh.ngl*(j-1+ (k-1)*mesh.ngl)
                        ip = mesh.connijk[i,j,k,iel]
                        for n=1:mesh.ngl
                            for m=1:mesh.ngl
                                for l=1:mesh.ngl
                                    MN = n + mesh.ngl*(m-1+(l-1)*mesh.ngl)
                                    Δq[i,j,k,iel,var] += qp.qn[ip,var]*L[MN,IJ,iel]
                                end
                            end
                        end
                    end
                end
            end
        end
    end

   return Δq
end

function build_laplacian(SD::NSD_1D, QT::Inexact, qp, ψ, dψ, ω, mesh::St_mesh, metrics::St_metrics)
    nvar = size(qp.qn,2)
    L = zeros(mesh.ngl,mesh.ngl,mesh.nelem)
    Δq = zeros(mesh.ngl,mesh.nelem,nvar)

    for iel=1:mesh.nelem
        for l=1:mesh.ngl
            for i=1:mesh.ngl
                dψJKdx = dψ[i,l]*metrics.dξdx[l,iel]
                for n=1:mesh.ngl
                    dψIKdx = dψ[n,l]*metrics.dξdx[k,l,iel]
                    L[n,i,iel] += ω[k]*metrics.Je[k,l,iel]*(dψIKdx*dψJKdx)
                end
            end
        end
    end
    #show(stdout, "text/plain", el_matrices.D)
    for var=1:nvar
        for iel=1:mesh.nelem
            for i=1:mesh.ngl
                ip = mesh.connijk[i,iel]
                for m=1:mesh.ngl
                    Δq[i,iel,var] += qp.qn[ip,var]*L[i,m,iel]
                end
            end
        end
    end

   return Δq
end
