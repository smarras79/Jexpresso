include("../abstractTypes.jl")
include("../mesh/mesh.jl")
include("../bases/basis_structs.jl")

Base.@kwdef mutable struct St_metrics{TFloat}
    dxdξ::Union{Array{TFloat}, Missing} = zeros(1)
    dxdη::Union{Array{TFloat}, Missing} = zeros(1)
    dxdζ::Union{Array{TFloat}, Missing} = zeros(1)
    
    dydξ::Union{Array{TFloat}, Missing} = zeros(1)
    dydη::Union{Array{TFloat}, Missing} = zeros(1)
    dydζ::Union{Array{TFloat}, Missing} = zeros(1)

    dzdξ::Union{Array{TFloat}, Missing} = zeros(1)
    dzdη::Union{Array{TFloat}, Missing} = zeros(1)
    dzdζ::Union{Array{TFloat}, Missing} = zeros(1)
    
    dξdx::Union{Array{TFloat}, Missing} = zeros(1)
    dξdy::Union{Array{TFloat}, Missing} = zeros(1)
    dξdz::Union{Array{TFloat}, Missing} = zeros(1)
    
    dηdx::Union{Array{TFloat}, Missing} = zeros(1)
    dηdy::Union{Array{TFloat}, Missing} = zeros(1)
    dηdz::Union{Array{TFloat}, Missing} = zeros(1)

    dζdx::Union{Array{TFloat}, Missing} = zeros(1)
    dζdy::Union{Array{TFloat}, Missing} = zeros(1)
    dζdz::Union{Array{TFloat}, Missing} = zeros(1)
    Jef ::Array{TFloat} = zeros(1)
    #
    # Contravariant arrays
    #
    vⁱ::Union{Array{TFloat}, Missing} = zeros(3) #contravariant unit vectors
        
    #
    # Element jacobian determinant
    #
    nx  ::Array{TFloat} = zeros(1)
    ny  ::Array{TFloat} = zeros(1)
    nz  ::Array{TFloat} = zeros(1)
    Je  ::Array{TFloat} = zeros(1)
end

function build_metric_terms(SD::NSD_1D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(1))
    
    return metrics
end


function build_metric_terms(SD::NSD_2D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(Q+1, Q+1, mesh.nelem), #∂x/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dxdη = zeros(Q+1, Q+1, mesh.nelem), #∂x/∂η[1:Nq, 1:Nq, 1:nelem]
                            dydξ = zeros(Q+1, Q+1, mesh.nelem), #∂y/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dydη = zeros(Q+1, Q+1, mesh.nelem), #∂y/∂η[1:Nq, 1:Nq, 1:nelem]

                            dξdx = zeros(Q+1, Q+1, mesh.nelem), #∂ξ/∂x[1:Nq, 1:Nq, 1:nelem]
                            dηdx = zeros(Q+1, Q+1, mesh.nelem), #∂η/∂x[1:Nq, 1:Nq, 1:nelem]
                            dξdy = zeros(Q+1, Q+1, mesh.nelem), #∂ξ/∂y[1:Nq, 1:Nq, 1:nelem]
                            dηdy = zeros(Q+1, Q+1, mesh.nelem), #∂η/∂y[1:Nq, 1:Nq, 1:nelem]
                            Jef  = zeros(Q+1, mesh.nedges_bdy),
                            nx   = zeros(Q+1, mesh.nedges_bdy),
                            ny   = zeros(Q+1, mesh.nedges_bdy),
                            Je   = zeros(Q+1, Q+1, mesh.nelem)) #   Je[1:Nq, 1:Nq, 1:nelem]

    ψ  = basis.ψ
    dψ = basis.dψ
    
    for iel = 1:mesh.nelem
        for l = 1:Q+1
            for k = 1:Q+1
                for j = 1:N+1
                    for i = 1:N+1

                        ip = mesh.connijk[i,j,iel]

                        xij = mesh.x[ip]
                        yij = mesh.y[ip]
                        
                        metrics.dxdξ[k, l, iel] += dψ[i,k]* ψ[j,l]*xij
                        metrics.dxdη[k, l, iel] +=  ψ[i,k]*dψ[j,l]*xij
                        
                        metrics.dydξ[k, l, iel] += dψ[i,k]* ψ[j,l]*yij
                        metrics.dydη[k, l, iel] +=  ψ[i,k]*dψ[j,l]*yij                        
                        #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[k, l, iel],  metrics.dxdη[k, l, iel], metrics.dydξ[k, l, iel],  metrics.dydη[k, l, iel] )
            end
        end
        
        for l = 1:Q+1
            for k = 1:Q+1
                metrics.Je[k, l, iel] = metrics.dxdξ[k, l, iel]*metrics.dydη[k, l, iel] - metrics.dydξ[k, l, iel]*metrics.dxdη[k, l, iel]
                
                metrics.dξdx[k, l, iel] =  metrics.dydη[k, l, iel]/metrics.Je[k, l, iel]
                metrics.dξdy[k, l, iel] = -metrics.dxdη[k, l, iel]/metrics.Je[k, l, iel]
                metrics.dηdx[k, l, iel] = -metrics.dydξ[k, l, iel]/metrics.Je[k, l, iel]
                metrics.dηdy[k, l, iel] =  metrics.dxdξ[k, l, iel]/metrics.Je[k, l, iel]
                
            end
        end
        #show(stdout, "text/plain", metrics.Je[:,:,iel])
    end
    #show(stdout, "text/plain", metrics.Je)
    for iedge =1:size(mesh.bdy_edge_comp,1)
        comp = mesh.bdy_edge_comp[iedge]
        for k=1:Q+1
            ip = mesh.poin_in_bdy_edge[iedge,k]
            if (k < Q+1)
              ip1 = mesh.poin_in_bdy_edge[iedge,k+1]
            else
              ip1 = mesh.poin_in_bdy_edge[iedge,k-1]
            end
            x1 = mesh.x[ip]
            x2 = mesh.x[ip1]
            y1 = mesh.y[ip]
            y2 = mesh.y[ip1]
            mag = sqrt((x1-x2)^2+(y1-y2)^2)
            metrics.Jef[k,iedge] = mag/2
            comp1 = (x1-x2)/mag
            comp2 = (y1-y2)/mag
            metrics.nx[k,iedge] = comp2
            metrics.ny[k,iedge] = -comp1
        end
    end

    return metrics
end

function build_metric_terms(SD::NSD_2D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Laguerre,N, Q, NGR, QGR, ξ, T;dir="x")
    
    metrics = St_metrics{T}(dxdξ = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂x/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dxdη = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂x/∂η[1:Nq, 1:Nq, 1:nelem]
                            dydξ = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂y/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dydη = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂y/∂η[1:Nq, 1:Nq, 1:nelem]

                            dξdx = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂ξ/∂x[1:Nq, 1:Nq, 1:nelem]
                            dηdx = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂η/∂x[1:Nq, 1:Nq, 1:nelem]
                            dξdy = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂ξ/∂y[1:Nq, 1:Nq, 1:nelem]
                            dηdy = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf), #∂η/∂y[1:Nq, 1:Nq, 1:nelem]
                            Jef  = zeros(mesh.ngl, size(mesh.bound_elem,1)+4),
                            Je   = zeros(mesh.ngl, mesh.ngr, mesh.nelem_semi_inf)) #   Je[1:Nq, 1:Nq, 1:nelem]
    
    ψ  = basis.ψ
    dψ = basis.dψ
    ψ1  = basisGR.ψ
    dψ1 = basisGR.dψ
    if ("Laguerre" in mesh.bdy_edge_type)
        for iel=1:mesh.nelem_semi_inf
            tan = mesh.bdy_tangents[iel,:]
            nor = mesh.bdy_normals[iel,:]
            for l=1:mesh.ngr
                for k =1:mesh.ngl
                    for j=1:mesh.ngr
                        for i=1:mesh.ngl
                            ip = mesh.connijk_lag[i,j,iel]
                            xij = mesh.x[ip]
                            yij = mesh.y[ip]
                            metrics.dxdξ[k, l, iel] += dψ[i,k]*ψ1[j,l]*xij 
                            metrics.dxdη[k, l, iel] +=  ψ[i,k]*dψ1[j,l]*xij
                            metrics.dydξ[k, l, iel] += dψ[i,k]* ψ1[j,l]*yij
                            metrics.dydη[k, l, iel] +=  ψ[i,k]*dψ1[j,l]*yij
                        end
                    end
                end
            end
            for l = 1:mesh.ngr
                for k = 1:mesh.ngl
                    metrics.Je[k, l, iel] = metrics.dxdξ[k, l, iel]*metrics.dydη[k, l, iel] - metrics.dydξ[k, l, iel]*metrics.dxdη[k, l, iel]

                    metrics.dξdx[k, l, iel] =  metrics.dydη[k, l, iel]/metrics.Je[k, l, iel]
                    metrics.dξdy[k, l, iel] = -metrics.dxdη[k, l, iel]/metrics.Je[k, l, iel]
                    metrics.dηdx[k, l, iel] = -metrics.dydξ[k, l, iel]/metrics.Je[k, l, iel]
                    metrics.dηdy[k, l, iel] =  metrics.dxdξ[k, l, iel]/metrics.Je[k, l, iel]

            end
        end
        #show(stdout, "text/plain", metrics.Je[:,:,iel])
    end

    return metrics
end

function build_metric_terms(SD::NSD_2D, MT::CNVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(Q+1, Q+1, mesh.nelem), #∂x/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dxdη = zeros(Q+1, Q+1, mesh.nelem), #∂x/∂η[1:Nq, 1:Nq, 1:nelem]
                            dydξ = zeros(Q+1, Q+1, mesh.nelem), #∂y/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dydη = zeros(Q+1, Q+1, mesh.nelem), #∂y/∂η[1:Nq, 1:Nq, 1:nelem]

                            dξdx = zeros(Q+1, Q+1, mesh.nelem), #∂ξ/∂x[1:Nq, 1:Nq, 1:nelem]
                            dηdx = zeros(Q+1, Q+1, mesh.nelem), #∂η/∂x[1:Nq, 1:Nq, 1:nelem]
                            dξdy = zeros(Q+1, Q+1, mesh.nelem), #∂ξ/∂y[1:Nq, 1:Nq, 1:nelem]
                            dηdy = zeros(Q+1, Q+1, mesh.nelem), #∂η/∂y[1:Nq, 1:Nq, 1:nelem]

                            vⁱ = zeros(3, Q+1, Q+1, mesh.nelem), #contravariant unit vectors
                            Jef  = zeros(Q+1, mesh.nedges_bdy),
                            nx   = zeros(Q+1, mesh.nedges_bdy),
                            ny   = zeros(Q+1, mesh.nedges_bdy),
                            Je   = zeros(Q+1, Q+1, mesh.nelem)) #   Je[1:Nq, 1:Nq, 1:nelem]

    ψ  = basis.ψ
    dψ = basis.dψ

    @info " WIP CONTRAVIARIANT metric terms"
    for iel = 1:mesh.nelem
        for l = 1:Q+1
            for k = 1:Q+1
                for j = 1:N+1
                    for i = 1:N+1

                        ip = mesh.connijk[i,j,iel]

                        xij = mesh.x[ip]
                        yij = mesh.y[ip]
                        
                        metrics.dxdξ[k, l, iel] = metrics.dxdξ[k, l, iel] + dψ[i,k]*ψ[j,l]*xij
                        metrics.dxdη[k, l, iel] = metrics.dxdη[k, l, iel] + ψ[i,k]*dψ[j,l]*xij
                        
                        metrics.dydξ[k, l, iel] = metrics.dydξ[k, l, iel] + dψ[i,k]*ψ[j,l]*yij
                        metrics.dydη[k, l, iel] = metrics.dydη[k, l, iel] + ψ[i,k]*dψ[j,l]*yij                        
                        #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[k, l, iel],  metrics.dxdη[k, l, iel], metrics.dydξ[k, l, iel],  metrics.dydη[k, l, iel] )
            end
        end
        
        for l = 1:Q+1
            for k = 1:Q+1
                metrics.Je[k, l, iel] = metrics.dxdξ[k, l, iel]*metrics.dydη[k, l, iel] - metrics.dydξ[k, l, iel]*metrics.dxdη[k, l, iel]
                
                metrics.dξdx[k, l, iel] =  metrics.dydη[k, l, iel]/metrics.Je[k, l, iel]
                metrics.dξdy[k, l, iel] = -metrics.dxdη[k, l, iel]/metrics.Je[k, l, iel]
                metrics.dηdx[k, l, iel] = -metrics.dydξ[k, l, iel]/metrics.Je[k, l, iel]
                metrics.dηdy[k, l, iel] =  metrics.dxdξ[k, l, iel]/metrics.Je[k, l, iel]

                vⁱ[1, k, l, iel] = metrics.dξdx[k, l, iel] + metrics.dξdy[k, l, iel]
                vⁱ[2, k, l, iel] = metrics.dηdx[k, l, iel] + metrics.dηdy[k, l, iel]
                vⁱ[3, k, l, iel] = 0
                
            end
        end
        #show(stdout, "text/plain", metrics.Je[:,:,iel])
    end
    #show(stdout, "text/plain", metrics.Je)
    #face jacobians for boundary faces
    for iedge =1:size(mesh.bdy_edge_comp,1)
        comp = mesh.bdy_edge_comp[iedge]
        for k=1:Q+1
            ip = mesh.poin_in_bdy_edge[iedge,k]
            if (k < Q+1)
              ip1 = mesh.poin_in_bdy_edge[iedge,k+1]
            else
              ip1 = mesh.poin_in_bdy_edge[iedge,k-1]
            end
            x1 = mesh.x[ip]
            x2 = mesh.x[ip1] 
            y1 = mesh.y[ip]
            y2 = mesh.y[ip2]
            mag = sqrt((x1-x2)^2+(y1-y2)^2)
            metrics.Jef[k,iedge] = mag/2
            comp1 = (x1-x2)/mag
            comp2 = (y1-y2)/mag
            metrics.nx[k,iedge] = comp2
            metrics.ny[k,iedge] = -comp1
        end
    end
    return metrics
end


function build_metric_terms(SD::NSD_3D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, T)

     metrics = St_metrics{T}(dxdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dydζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdz[2, 1:Nq, 1:Nq, 1:nelem]
                             Jef  = zeros(Q+1, Q+1, size(mesh.bound_elem,1)),
                             Je   = zeros(Q+1, Q+1, Q+1, mesh.nelem))  #   Je[1:Nq, 1:Nq, 1:nelem]
    ψ  = basis.ψ
    dψ = basis.dψ

    @info " metric terms WIP"
    for iel = 1:mesh.nelem
        for n = 1:Q+1
            for m = 1:Q+1
                for l = 1:Q+1
                    for k = 1:N+1
                        for j = 1:N+1
                            for i = 1:N+1

                                ip = mesh.connijk[i,j,k,iel]

                                xijk = mesh.x[ip]
                                yijk = mesh.y[ip]
                                zijk = mesh.z[ip]

                                metrics.dxdξ[l, m, n, iel] = metrics.dxdξ[l, m, n, iel] + dψ[i,l]*ψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdη[l, m, n, iel] = metrics.dxdη[l, m, n, iel] + ψ[i,l]*dψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdζ[l, m, n, iel] = metrics.dxdζ[l, m, n, iel] + ψ[i,l]*ψ[j,m]*dψ[k,n]*xijk

                                metrics.dydξ[l, m, n, iel] = metrics.dydξ[l, m, n, iel] + dψ[i,l]*ψ[j,m]*ψ[k,n]*yijk
                                metrics.dydη[l, m, n, iel] = metrics.dydη[l, m, n, iel] + ψ[i,l]*dψ[j,m]*ψ[k,n]*yijk
                                metrics.dydζ[l, m, n, iel] = metrics.dydζ[l, m, n, iel] + ψ[i,l]*ψ[j,m]*dψ[k,n]*yijk

                                metrics.dzdξ[l, m, n, iel] = metrics.dzdξ[l, m, n, iel] + dψ[i,l]*ψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdη[l, m, n, iel] = metrics.dzdη[l, m, n, iel] + ψ[i,l]*dψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdζ[l, m, n, iel] = metrics.dzdζ[l, m, n, iel] + ψ[i,l]*ψ[j,m]*dψ[k,n]*zijk
                                #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                            end
                        end
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[k, l, iel],  metrics.dxdη[k, l, iel], metrics.dydξ[k, l, iel],  metrics.dydη[k, l, iel] )
            end
        end

        for l = 1:Q+1
            for m = 1:Q+1
                for n =1:Q+1
                    metrics.Je[l, m, n, iel] = metrics.dxdξ[l, m, n, iel]*(metrics.dydη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])
                    metrics.Je[l, m, n, iel] += metrics.dydξ[l, m, n, iel]*(metrics.dxdζ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel])
                    metrics.Je[l, m, n, iel] += metrics.dzdξ[l, m, n, iel]*(metrics.dxdη[l, m, n, iel]*metrics.dydζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dydη[l, m, n, iel])
                    
                    metrics.dξdx[l, m, n, iel] =  (metrics.dydη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dξdy[l, m, n, iel] =  (metrics.dxdζ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dξdz[l, m, n, iel] =  (metrics.dxdη[l, m, n, iel]*metrics.dydζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dydη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdx[l, m, n, iel] =  (metrics.dydζ[l, m, n, iel]*metrics.dzdξ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdζ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdy[l, m, n, iel] =  (metrics.dxdξ[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dzdξ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdz[l, m, n, iel] =  (metrics.dxdζ[l, m, n, iel]*metrics.dydξ[l, m, n, iel] - metrics.dxdξ[l, m, n, iel]*metrics.dydζ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdx[l, m, n, iel] =  (metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dydη[l, m, n, iel]*metrics.dzdξ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdy[l, m, n, iel] =  (metrics.dxdη[l, m, n, iel]*metrics.dzdξ[l, m, n, iel] - metrics.dxdξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdz[l, m, n, iel] =  (metrics.dxdξ[l, m, n, iel]*metrics.dydη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dydξ[l, m, n, iel])/metrics.Je[l, m, n, iel] 
                
                 end
            end
        end
        show(stdout, "text/plain", metrics.Je[:,:,:,iel])
    end
    for iface = 1:size(mesh.xmin_faces,2) 
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.xmin_facetoelem[iface]
                metrics.Jef[l, m, iface] = metrics.dydη[1, l, m, iel]*metrics.dzdζ[1, l, m, iel] - metrics.dydζ[1, l, m, iel]*metrics.dzdη[1, l, m, iel]
            end
        end
    end
    disp = size(mesh.xmin_faces,2)
    for iface = 1:size(mesh.xmax_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.xmax_facetoelem[iface]
                metrics.Jef[l, m, iface+disp] = metrics.dydη[Q+1, l, m, iel]*metrics.dzdζ[Q+1, l, m, iel] - metrics.dydζ[Q+1, l, m, iel]*metrics.dzdη[Q+1, l, m, iel]
            end
        end
    end
    disp += size(mesh.xmax_faces,2)
    for iface = 1:size(mesh.ymin_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.ymin_facetoelem[iface]
                metrics.Jef[l, m, iface+disp] = metrics.dxdξ[l, 1, m, iel]*metrics.dzdζ[l, 1, m, iel] - metrics.dzdξ[l, 1, m, iel]*metrics.dxdζ[l, 1, m, iel]
            end
        end
    end
    disp += size(mesh.ymin_faces,2)
    for iface = 1:size(mesh.ymax_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.ymax_facetoelem[iface]
                metrics.Jef[l, m, iface+disp] = metrics.dxdξ[l, Q+1, m, iel]*metrics.dzdζ[l, Q+1, m, iel] - metrics.dzdξ[l, Q+1, m, iel]*metrics.dxdζ[l, Q+1, m, iel]
            end
        end
    end
    disp += size(mesh.ymax_faces,2)
    for iface = 1:size(mesh.zmin_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.zmin_facetoelem[iface]
                metrics.Jef[l, m, iface+disp] = metrics.dxdξ[l, m, 1, iel]*metrics.dydη[l, m, 1, iel] - metrics.dydξ[l, m, 1, iel]*metrics.dxdη[l, m, 1, iel]
            end
        end
    end
    disp += size(mesh.zmin_faces,2)
    for iface = 1:size(mesh.zmax_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.zmax_facetoelem[iface]
                metrics.Jef[l, m, iface+disp] = metrics.dxdξ[l, m, Q+1, iel]*metrics.dydη[l, m, Q+1, iel] - metrics.dydξ[l, m, Q+1, iel]*metrics.dxdη[l, m, Q+1, iel]
            end
        end
    end
     
    #show(stdout, "text/plain", metrics.Je)    
    return metrics
end

function build_metric_terms(SD::NSD_3D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Laguerre,N, Q, NGR, QGR, ξ, T;dir="x",side ="min")

     metrics = St_metrics{T}(dxdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dydζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdz[2, 1:Nq, 1:Nq, 1:nelem]
                             Jef  = zeros(Q+1, Q+1, size(mesh.bound_elem,1)),
                             Je   = zeros(Q+1, Q+1, Q+1, mesh.nelem))  #   Je[1:Nq, 1:Nq, 1:nelem]
    
    if (dir == "x")
      ψ1  = basisGR.ψ
      dψ1 = basisGR.dψ
      ψ  = basis.ψ
      dψ = basis.dψ
      ψ2  = basis.ψ
      dψ2 = basis.dψ
      N1 = NGR
      Q1 = QGR
      N2 = N
      Q2 = Q
      N3 = N
      Q3 = Q
      if (side == "min")
        nelem = size(mesh.xmin_faces,3)
        conn = mesh.conn_xminlag
        x = mesh.xmin_xtra
        y = mesh.ymin_xtra
        z = mesh.zmin_xtra
      else
        nelem = size(mesh.xmax_faces,3)
        conn = mesh.conn_xmaxlag
        x = mesh.xmax_xtra
        y = mesh.ymax_xtra
        z = mesh.zmax_xtra
      end
    elseif (dir == "y")
      ψ  = basis.ψ
      dψ = basis.dψ
      ψ  = basisGR.ψ
      dψ = basisGR.dψ
      ψ2  = basis.ψ
      dψ2 = basis.dψ
      N1 = N
      Q1 = Q
      N2 = NGR
      Q2 = QGR
      N3 = N
      Q3 = Q
      if (side == "min")
        nelem = size(mesh.ymin_faces,3)
        conn = mesh.conn_yminlag
        x = mesh.xmin_ytra
        y = mesh.ymin_ytra
        z = mesh.zmin_ytra
      else
        nelem = size(mesh.ymax_faces,3)
        conn = mesh.conn_ymaxlag
        x = mesh.xmax_ytra
        y = mesh.ymax_ytra
        z = mesh.zmax_ytra
      end
    else
      ψ  = basis.ψ
      dψ = basis.dψ
      ψ  = basis.ψ
      dψ = basis.dψ
      ψ2  = basisGR.ψ
      dψ2 = basisGR.dψ
      N1 = N
      Q1 = Q
      N2 = N
      Q2 = Q
      N3 = NGR
      Q3 = QGR
      if (side == "min")
        nelem = size(mesh.zmin_faces,3)
        conn = mesh.conn_zminlag
        x = mesh.xmin_ztra
        y = mesh.ymin_ztra
        z = mesh.zmin_ztra
      else
        nelem = size(mesh.zmax_faces,3)
        conn = mesh.conn_zmaxlag
        x = mesh.xmax_ztra
        y = mesh.ymax_ztra
        z = mesh.zmax_ztra
      end
    end

    @info " metric terms WIP"
    for iel = 1:nelem
        for n = 1:Q3+1
            for m = 1:Q2+1
                for l = 1:Q1+1
                    for k = 1:N1+1
                        for j = 1:N2+1
                            for i = 1:N3+1

                                ip = conn[i,j,k,iel]

                                xijk = x[ip]
                                yijk = y[ip]
                                zijk = z[ip]

                                metrics.dxdξ[l, m, n, iel] = metrics.dxdξ[l, m, n, iel] + dψ1[i,l]*ψ[j,m]*ψ2[k,n]*xijk
                                metrics.dxdη[l, m, n, iel] = metrics.dxdη[l, m, n, iel] + ψ1[i,l]*dψ[j,m]*ψ2[k,n]*xijk
                                metrics.dxdζ[l, m, n, iel] = metrics.dxdζ[l, m, n, iel] + ψ1[i,l]*ψ[j,m]*dψ2[k,n]*xijk

                                metrics.dydξ[l, m, n, iel] = metrics.dydξ[l, m, n, iel] + dψ1[i,l]*ψ[j,m]*ψ2[k,n]*yijk
                                metrics.dydη[l, m, n, iel] = metrics.dydη[l, m, n, iel] + ψ1[i,l]*dψ[j,m]*ψ2[k,n]*yijk
                                metrics.dydζ[l, m, n, iel] = metrics.dydζ[l, m, n, iel] + ψ1[i,l]*ψ[j,m]*dψ2[k,n]*yijk

                                metrics.dzdξ[l, m, n, iel] = metrics.dzdξ[l, m, n, iel] + dψ1[i,l]*ψ[j,m]*ψ2[k,n]*zijk
                                metrics.dzdη[l, m, n, iel] = metrics.dzdη[l, m, n, iel] + ψ1[i,l]*dψ[j,m]*ψ2[k,n]*zijk
                                metrics.dzdζ[l, m, n, iel] = metrics.dzdζ[l, m, n, iel] + ψ1[i,l]*ψ[j,m]*dψ2[k,n]*zijk
                                #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                            end
                        end
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[k, l, iel],  metrics.dxdη[k, l, iel], metrics.dydξ[k, l, iel],  metrics.dydη[k, l, iel] )
            end
        end
        for l = 1:Q3+1
            for m = 1:Q2+1
                for n =1:Q1+1
                    metrics.Je[l, m, n, iel] = metrics.dxdξ[l, m, n, iel]*(metrics.dydη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])
                    metrics.Je[l, m, n, iel] += metrics.dydξ[l, m, n, iel]*(metrics.dxdζ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel])
                    metrics.Je[l, m, n, iel] += metrics.dzdξ[l, m, n, iel]*(metrics.dxdη[l, m, n, iel]*metrics.dydζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dydη[l, m, n, iel])

                    metrics.dξdx[l, m, n, iel] =  (metrics.dydη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dξdy[l, m, n, iel] =  (metrics.dxdζ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dξdz[l, m, n, iel] =  (metrics.dxdη[l, m, n, iel]*metrics.dydζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dydη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdx[l, m, n, iel] =  (metrics.dydζ[l, m, n, iel]*metrics.dzdξ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdζ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdy[l, m, n, iel] =  (metrics.dxdξ[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dzdξ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdz[l, m, n, iel] =  (metrics.dxdζ[l, m, n, iel]*metrics.dydξ[l, m, n, iel] - metrics.dxdξ[l, m, n, iel]*metrics.dydζ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdx[l, m, n, iel] =  (metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dydη[l, m, n, iel]*metrics.dzdξ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdy[l, m, n, iel] =  (metrics.dxdη[l, m, n, iel]*metrics.dzdξ[l, m, n, iel] - metrics.dxdξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdz[l, m, n, iel] =  (metrics.dxdξ[l, m, n, iel]*metrics.dydη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dydξ[l, m, n, iel])/metrics.Je[l, m, n, iel]

                 end
            end
        end
        show(stdout, "text/plain", metrics.Je[:,:,:,iel])
    end
end

function build_metric_terms(SD::NSD_3D, MT::CNVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, T)

     metrics = St_metrics{T}(dxdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂x/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dydζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂y/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdξ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdη = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdζ = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂z/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdx = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdy = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ξdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ηdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdz = zeros(Q+1, Q+1, Q+1, mesh.nelem),  #∂ζdz[2, 1:Nq, 1:Nq, 1:nelem]

                             vⁱ = zeros(3, Q+1, Q+1, Q+1, mesh.nelem), #contravariant unit vectors

                             Je   = zeros(Q+1, Q+1, Q+1, mesh.nelem))  #   Je[1:Nq, 1:Nq, 1:nelem]
    ψ  = basis.ψ
    dψ = basis.dψ

    @info " metric terms WIP"
    for iel = 1:mesh.nelem
        for n = 1:Q+1
            for m = 1:Q+1
                for l = 1:Q+1
                    for k = 1:N+1
                        for j = 1:N+1
                            for i = 1:N+1

                                ip = mesh.connijk[i,j,k,iel]

                                xijk = mesh.x[ip]
                                yijk = mesh.y[ip]
                                zijk = mesh.z[ip]

                                metrics.dxdξ[l, m, n, iel] = metrics.dxdξ[l, m, n, iel] + dψ[i,l]*ψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdη[l, m, n, iel] = metrics.dxdη[l, m, n, iel] + ψ[i,l]*dψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdζ[l, m, n, iel] = metrics.dxdζ[l, m, n, iel] + ψ[i,l]*ψ[j,m]*dψ[k,n]*xijk

                                metrics.dydξ[l, m, n, iel] = metrics.dydξ[l, m, n, iel] + dψ[i,l]*ψ[j,m]*ψ[k,n]*yijk
                                metrics.dydη[l, m, n, iel] = metrics.dydη[l, m, n, iel] + ψ[i,l]*dψ[j,m]*ψ[k,n]*yijk
                                metrics.dydζ[l, m, n, iel] = metrics.dydζ[l, m, n, iel] + ψ[i,l]*ψ[j,m]*dψ[k,n]*yijk

                                metrics.dzdξ[l, m, n, iel] = metrics.dzdξ[l, m, n, iel] + dψ[i,l]*ψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdη[l, m, n, iel] = metrics.dzdη[l, m, n, iel] + ψ[i,l]*dψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdζ[l, m, n, iel] = metrics.dzdζ[l, m, n, iel] + ψ[i,l]*ψ[j,m]*dψ[k,n]*zijk
                        #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                            end
                        end
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[k, l, iel],  metrics.dxdη[k, l, iel], metrics.dydξ[k, l, iel],  metrics.dydη[k, l, iel] )
            end
        end

        for l = 1:Q+1
            for m = 1:Q+1
                for n =1:Q+1
                    metrics.Je[l, m, n, iel] = metrics.dxdξ[l, m, n, iel]*(metrics.dydη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])
                    metrics.Je[l, m, n, iel] += metrics.dydξ[l, m, n, iel]*(metrics.dxdζ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel])
                    metrics.Je[l, m, n, iel] += metrics.dzdξ[l, m, n, iel]*(metrics.dxdη[l, m, n, iel]*metrics.dydζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dydη[l, m, n, iel])
                    
                    metrics.dξdx[l, m, n, iel] =  (metrics.dydη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dξdy[l, m, n, iel] =  (metrics.dxdζ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dξdz[l, m, n, iel] =  (metrics.dxdη[l, m, n, iel]*metrics.dydζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dydη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdx[l, m, n, iel] =  (metrics.dydζ[l, m, n, iel]*metrics.dzdξ[l, m, n, iel] - metrics.dydξ[l, m, n, iel]*metrics.dzdζ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdy[l, m, n, iel] =  (metrics.dxdξ[l, m, n, iel]*metrics.dzdζ[l, m, n, iel] - metrics.dxdζ[l, m, n, iel]*metrics.dzdξ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dηdz[l, m, n, iel] =  (metrics.dxdζ[l, m, n, iel]*metrics.dydξ[l, m, n, iel] - metrics.dxdξ[l, m, n, iel]*metrics.dydζ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdx[l, m, n, iel] =  (metrics.dydξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dydη[l, m, n, iel]*metrics.dzdξ[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdy[l, m, n, iel] =  (metrics.dxdη[l, m, n, iel]*metrics.dzdξ[l, m, n, iel] - metrics.dxdξ[l, m, n, iel]*metrics.dzdη[l, m, n, iel])/metrics.Je[l, m, n, iel]
                    metrics.dζdz[l, m, n, iel] =  (metrics.dxdξ[l, m, n, iel]*metrics.dydη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dydξ[l, m, n, iel])/metrics.Je[l, m, n, iel] 

                    

                    vⁱ[1, l, m, n, iel] = metrics.dξdx[l, m, n, iel] + metrics.dξdy[l, m, n, iel]
                    vⁱ[2, l, m, n, iel] = metrics.dηdx[l, m, n, iel] + metrics.dηdy[l, m, n, iel]
                    vⁱ[3, l, m, n, iel] = metrics.dζdx[l, m, n, iel] + metrics.dζdy[l, m, n, iel]
                    
                 end
            end
        end
        show(stdout, "text/plain", metrics.Je[:,:,:,iel])
    end
    #show(stdout, "text/plain", metrics.Je)    
    return metrics
end
