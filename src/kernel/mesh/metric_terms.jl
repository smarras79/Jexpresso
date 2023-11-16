Base.@kwdef mutable struct St_metrics{TFloat}

    
    
    dxdξ::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dxdη::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dxdζ::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    
    dydξ::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dydη::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dydζ::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)

    dzdξ::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dzdη::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dzdζ::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    
    dξdx::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dξdy::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dξdz::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    
    dηdx::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dηdy::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dηdz::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)

    dζdx::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dζdy::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    dζdz::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    #
    # Contravariant arrays
    #
    vⁱ::Union{Array{TFloat}, Missing} = zeros(3) #contravariant unit vectors
        
    #
    # Element jacobian determinant
    #
    nx::Array{TFloat, 2} = zeros(TFloat, 0, 0)
    ny::Array{TFloat, 2} = zeros(TFloat, 0, 0)
    nz::Array{TFloat, 2} = zeros(TFloat, 0, 0)
    Je::Array{TFloat, 3} = zeros(TFloat, 0, 0, 0)
    Jef::Array{TFloat,2} = zeros(TFloat, 0, 0)
end

function build_metric_terms(SD::NSD_1D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(T, mesh.nelem, Q+1, 1), #∂x/∂ξ[1:Nq, 1:nelem]
                            dξdx = zeros(T, mesh.nelem, Q+1, 1), #∂ξ/∂x[1:Nq, 1:nelem]
                            Je   = zeros(T, mesh.nelem, Q+1, 1),
                            nx   = zeros(T, mesh.nedges_bdy, 1))

    for iel = 1:mesh.nelem
        for i = 1:N+1
            for k = 1:Q+1
                metrics.dxdξ[iel, k,1]  = mesh.Δx[iel]/2
                metrics.Je[iel, k, 1]   = metrics.dxdξ[iel, k, 1]
                metrics.dξdx[iel, k, 1] = 1.0/metrics.Je[iel, k, 1]
            end
        end        
    end     
    
    return metrics
end

function build_metric_terms_1D_Laguerre(SD::NSD_1D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, inputs,T)
   
    metrics = St_metrics{T}(dxdξ = zeros(T, mesh.nelem_semi_inf, Q+1, 1), #∂x/∂ξ[1:Nq, 1:nelem]
                            dξdx = zeros(T, mesh.nelem_semi_inf, Q+1, 1), #∂ξ/∂x[1:Nq, 1:nelem]
                            Je   = zeros(T, mesh.nelem_semi_inf, Q+1, 1),
                            nx   = zeros(T, mesh.nedges_bdy, 1))

    dψ = basis.dψ
    for iel = 1:mesh.nelem_semi_inf
        for i = 1:mesh.ngr
            ip = mesh.connijk_lag[iel,i,1]
            xij = mesh.x[ip]
            for k = 1:mesh.ngr
                metrics.dxdξ[iel, k,1]  += dψ[i,k] * xij
                metrics.Je[iel, k, 1]   = inputs[:yfac_laguerre]#abs(metrics.dxdξ[iel, k, 1])
                metrics.dξdx[iel, k, 1] = 1.0/metrics.Je[iel, k, 1]
            end
        end
    end    

    return metrics
end

function build_metric_terms(SD::NSD_2D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(T, mesh.nelem, Q+1, Q+1), #∂x/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dxdη = zeros(T, mesh.nelem, Q+1, Q+1), #∂x/∂η[1:Nq, 1:Nq, 1:nelem]
                            dydξ = zeros(T, mesh.nelem, Q+1, Q+1), #∂y/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dydη = zeros(T, mesh.nelem, Q+1, Q+1), #∂y/∂η[1:Nq, 1:Nq, 1:nelem]
                            dξdx = zeros(T, mesh.nelem, Q+1, Q+1), #∂ξ/∂x[1:Nq, 1:Nq, 1:nelem]
                            dηdx = zeros(T, mesh.nelem, Q+1, Q+1), #∂η/∂x[1:Nq, 1:Nq, 1:nelem]
                            dξdy = zeros(T, mesh.nelem, Q+1, Q+1), #∂ξ/∂y[1:Nq, 1:Nq, 1:nelem]
                            dηdy = zeros(T, mesh.nelem, Q+1, Q+1), #∂η/∂y[1:Nq, 1:Nq, 1:nelem]
                            Je   = zeros(T, mesh.nelem, Q+1, Q+1),
                            Jef  = zeros(T, mesh.nedges_bdy, Q+1),
                            nx   = zeros(T, mesh.nedges_bdy, Q+1),
                            ny   = zeros(T, mesh.nedges_bdy, Q+1)) #   Je[1:Nq, 1:Nq, 1:nelem]
    
    ψ  = @view(basis.ψ[:,:])
    dψ = @view(basis.dψ[:,:])
    
    for iel = 1:mesh.nelem
        for j = 1:N+1
            for i = 1:N+1
                ip = mesh.connijk[iel,i,j]
                
                xij = mesh.x[ip]; yij = mesh.y[ip]
                for l = 1:Q+1
                    for k = 1:Q+1
                        
                        metrics.dxdξ[iel, k, l] += dψ[i,k]*ψ[j,l] * xij
                        metrics.dxdη[iel, k, l] += ψ[i,k]*dψ[j,l] * xij

                        metrics.dydξ[iel, k, l] += dψ[i,k]*ψ[j,l] * yij
                        metrics.dydη[iel, k, l] += ψ[i,k]*dψ[j,l] * yij
                        
                        #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                    end
                end
            end
        end
             
        for l = 1:Q+1
            for k = 1:Q+1
              
                # Extract values from memory once per iteration
                dxdξ_val = metrics.dxdξ[iel, k, l]
                dydη_val = metrics.dydη[iel, k, l]
                dydξ_val = metrics.dydξ[iel, k, l]
                dxdη_val = metrics.dxdη[iel, k, l]
                # Compute Je once and reuse its value
                Je_val = dxdξ_val * dydη_val - dydξ_val * dxdη_val
                metrics.Je[iel, k, l] = Je_val

                # Use the precomputed Je value for the other calculations
                inv_Je = 1.0 / Je_val

                metrics.dξdx[iel, k, l] =  dydη_val * inv_Je
                metrics.dξdy[iel, k, l] = -dxdη_val * inv_Je
                metrics.dηdx[iel, k, l] = -dydξ_val * inv_Je
                metrics.dηdy[iel, k, l] =  dxdξ_val * inv_Je
                
            end
        end
        #show(stdout, "text/plain", metrics.Je[iel, :,:])
    end
    #show(stdout, "text/plain", metrics.Je)

    #=
    nbdy_edges = size(mesh.bdy_edge_comp,1)
    for iedge =1:nbdy_edges
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
            metrics.Jef[iedge, k] = mag/2
            comp1 = (x1-x2)/mag
            comp2 = (y1-y2)/mag
            metrics.nx[iedge, k] = comp2
            metrics.ny[iedge, k] = -comp1
        end
    end=#

    ##### CHECK BELOW WITH YASSINE WHETHER THIS
    ##### SHOULD BE done on the quadrature or grid points
    nbdy_edges = size(mesh.bdy_edge_comp,1)
    for iedge =1:nbdy_edges
        comp = mesh.bdy_edge_comp[iedge]
        for k=1:N+1
            ip = mesh.poin_in_bdy_edge[iedge,k]
            if (k < N+1)
              ip1 = mesh.poin_in_bdy_edge[iedge,k+1]
            else
              ip1 = mesh.poin_in_bdy_edge[iedge,k-1]
            end
            x1 = mesh.x[ip]
            x2 = mesh.x[ip1]
            y1 = mesh.y[ip]
            y2 = mesh.y[ip1]
            mag = sqrt((x1-x2)^2+(y1-y2)^2)
            metrics.Jef[iedge, k] = mag/2
            comp1 = (x1-x2)/mag
            comp2 = (y1-y2)/mag
            metrics.nx[iedge, k] = comp2
            metrics.ny[iedge, k] = -comp1
        end
    end

    return metrics
end

function build_metric_terms(SD::NSD_2D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange ,N, Q, NGR, QGR, ξ, ω1, ω2, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂x/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dxdη = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂x/∂η[1:Nq, 1:Nq, 1:nelem]
                            dydξ = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂y/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dydη = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂y/∂η[1:Nq, 1:Nq, 1:nelem]

                            dξdx = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂ξ/∂x[1:Nq, 1:Nq, 1:nelem]
                            dηdx = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂η/∂x[1:Nq, 1:Nq, 1:nelem]
                            dξdy = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂ξ/∂y[1:Nq, 1:Nq, 1:nelem]
                            dηdy = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr), #∂η/∂y[1:Nq, 1:Nq, 1:nelem]
                            Jef  = zeros(T, mesh.nedges_bdy, Q+1),
                            Je   = zeros(mesh.nelem_semi_inf, mesh.ngl, mesh.ngr))
    ψ  = basis.ψ
    dψ = basis.dψ
    ψ1  = basisGR.ψ
    dψ1 = basisGR.dψ
    if ("Laguerre" in mesh.bdy_edge_type)
        for iel=1:mesh.nelem_semi_inf
            for j=1:mesh.ngr
                for i =1:mesh.ngl
                    ip = mesh.connijk_lag[iel,i,j]
                    xij = mesh.x[ip]
                    yij = mesh.y[ip]
                    for l=1:mesh.ngr
                        for k=1:mesh.ngl
                            metrics.dxdξ[iel, k, l] += dψ[i,k]*ψ1[j,l]*xij 
                            metrics.dxdη[iel, k, l] +=  ψ[i,k]*inputs[:xfac_laguerre]#*dψ1[j,l]*xij
                            ###metrics.dxdη[iel, k, l] +=  ψ[i,k]*dψ1[j,l]*xij
                            metrics.dydξ[iel, k, l] += dψ[i,k]* ψ1[j,l]*yij
                            metrics.dydη[iel, k, l] +=  ψ[i,k]*inputs[:yfac_laguerre]/mesh.ngr#*dψ1[j,l]*yij
                            ###metrics.dydη[iel, k, l] +=  ψ[i,k]*dψ1[j,l]*yij
                        end
                    end
                    #@info metrics.dxdξ[iel, k, l],metrics.dxdη[iel, k, l],metrics.dydξ[iel, k, l],metrics.dydη[iel, k, l],xij,yij
                end
            end
            for l = 1:mesh.ngr
                for k = 1:mesh.ngl
                    #@info metrics.dxdξ[iel, k, l],metrics.dydη[iel, k, l], metrics.dydξ[iel, k, l],metrics.dxdη[iel, k, l]
                    metrics.Je[iel, k, l] = metrics.dxdξ[iel, k, l]*metrics.dydη[iel, k, l] - metrics.dydξ[iel, k, l]*metrics.dxdη[iel, k, l]
                    metrics.dξdx[iel, k, l] =  metrics.dydη[iel, k, l]/metrics.Je[iel, k, l]
                    metrics.dξdy[iel, k, l] = -metrics.dxdη[iel, k, l]/metrics.Je[iel, k, l]
                    metrics.dηdx[iel, k, l] = -metrics.dydξ[iel, k, l]/metrics.Je[iel, k, l]
                    metrics.dηdy[iel, k, l] =  metrics.dxdξ[iel, k, l]/metrics.Je[iel, k, l]
       
                end
            end
        #show(stdout, "text/plain", metrics.Je[:,:,iel])
        end
    end

    return metrics
end

function build_metric_terms(SD::NSD_2D, MT::CNVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T)
    
    metrics = St_metrics{T}(dxdξ = zeros(mesh.nelem, Q+1, Q+1), #∂x/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dxdη = zeros(mesh.nelem, Q+1, Q+1), #∂x/∂η[1:Nq, 1:Nq, 1:nelem]
                            dydξ = zeros(mesh.nelem, Q+1, Q+1), #∂y/∂ξ[1:Nq, 1:Nq, 1:nelem]
                            dydη = zeros(mesh.nelem, Q+1, Q+1), #∂y/∂η[1:Nq, 1:Nq, 1:nelem]

                            dξdx = zeros(mesh.nelem, Q+1, Q+1), #∂ξ/∂x[1:Nq, 1:Nq, 1:nelem]
                            dηdx = zeros(mesh.nelem, Q+1, Q+1), #∂η/∂x[1:Nq, 1:Nq, 1:nelem]
                            dξdy = zeros(mesh.nelem, Q+1, Q+1), #∂ξ/∂y[1:Nq, 1:Nq, 1:nelem]
                            dηdy = zeros(mesh.nelem, Q+1, Q+1), #∂η/∂y[1:Nq, 1:Nq, 1:nelem]

                            vⁱ = zeros(3, Q+1, Q+1, mesh.nelem), #contravariant unit vectors
                            Jef  = zeros(mesh.nedges_bdy, Q+1),
                            nx   = zeros(mesh.nedges_bdy, Q+1),
                            ny   = zeros(mesh.nedges_bdy, Q+1),
                            Je   = zeros(mesh.nelem, Q+1, Q+1)) #   Je[1:Nq, 1:Nq, 1:nelem]

    ψ  = basis.ψ
    dψ = basis.dψ

    @info " WIP CONTRAVIARIANT metric terms"
    for iel = 1:mesh.nelem
        for l = 1:Q+1
            for k = 1:Q+1
                for j = 1:N+1
                    for i = 1:N+1

                        ip = mesh.connijk[iel,i,j]

                        xij = mesh.x[ip]
                        yij = mesh.y[ip]
                        
                        metrics.dxdξ[iel, k, l] = metrics.dxdξ[iel, k, l] + dψ[i,k]*ψ[j,l]*xij
                        metrics.dxdη[iel, k, l] = metrics.dxdη[iel, k, l] + ψ[i,k]*dψ[j,l]*xij
                        
                        metrics.dydξ[iel, k, l] = metrics.dydξ[iel, k, l] + dψ[i,k]*ψ[j,l]*yij
                        metrics.dydη[iel, k, l] = metrics.dydη[iel, k, l] + ψ[i,k]*dψ[j,l]*yij                        
                        #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[iel, k, l],  metrics.dxdη[iel, k, l], metrics.dydξ[iel, k, l],  metrics.dydη[iel, k, l] )
            end
        end
        
        for l = 1:Q+1
            for k = 1:Q+1
                metrics.Je[iel, k, l] = metrics.dxdξ[iel, k, l]*metrics.dydη[iel, k, l] - metrics.dydξ[iel, k, l]*metrics.dxdη[iel, k, l]
                
                metrics.dξdx[iel, k, l] =  metrics.dydη[iel, k, l]/metrics.Je[iel, k, l]
                metrics.dξdy[iel, k, l] = -metrics.dxdη[iel, k, l]/metrics.Je[iel, k, l]
                metrics.dηdx[iel, k, l] = -metrics.dydξ[iel, k, l]/metrics.Je[iel, k, l]
                metrics.dηdy[iel, k, l] =  metrics.dxdξ[iel, k, l]/metrics.Je[iel, k, l]

                vⁱ[iel, 1, k, l] = metrics.dξdx[iel, k, l] + metrics.dξdy[iel, k, l]
                vⁱ[iel, 2, k, l] = metrics.dηdx[iel, k, l] + metrics.dηdy[iel, k, l]
                vⁱ[iel, 3, k, l] = 0
                
            end
        end
        #show(stdout, "text/plain", metrics.Je[iel, :,:])
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
            metrics.Jef[iedge,k] = mag/2
            comp1 = (x1-x2)/mag
            comp2 = (y1-y2)/mag
            metrics.nx[iedge, k] = comp2
            metrics.ny[iedge, k] = -comp1
        end
    end
    return metrics
end


function build_metric_terms(SD::NSD_3D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T)

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
                             Jef  = zeros(size(mesh.bound_elem,1), Q+1, Q+1),
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

                                ip = mesh.connijk[iel,i,j,k]

                                xijk = mesh.x[ip]
                                yijk = mesh.y[ip]
                                zijk = mesh.z[ip]

                                metrics.dxdξ[iel, l, m, n] = metrics.dxdξ[iel, l, m, n] + dψ[i,l]*ψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdη[iel, l, m, n] = metrics.dxdη[iel, l, m, n] + ψ[i,l]*dψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdζ[iel, l, m, n] = metrics.dxdζ[iel, l, m, n] + ψ[i,l]*ψ[j,m]*dψ[k,n]*xijk

                                metrics.dydξ[iel, l, m, n] = metrics.dydξ[iel, l, m, n] + dψ[i,l]*ψ[j,m]*ψ[k,n]*yijk
                                metrics.dydη[iel, l, m, n] = metrics.dydη[iel, l, m, n] + ψ[i,l]*dψ[j,m]*ψ[k,n]*yijk
                                metrics.dydζ[iel, l, m, n] = metrics.dydζ[iel, l, m, n] + ψ[i,l]*ψ[j,m]*dψ[k,n]*yijk

                                metrics.dzdξ[iel, l, m, n] = metrics.dzdξ[iel, l, m, n] + dψ[i,l]*ψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdη[iel, l, m, n] = metrics.dzdη[iel, l, m, n] + ψ[i,l]*dψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdζ[iel, l, m, n] = metrics.dzdζ[iel, l, m, n] + ψ[i,l]*ψ[j,m]*dψ[k,n]*zijk
                                #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                            end
                        end
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[iel, k, l],  metrics.dxdη[iel, k, l], metrics.dydξ[iel, k, l],  metrics.dydη[iel, k, l] )
            end
        end

        for l = 1:Q+1
            for m = 1:Q+1
                for n =1:Q+1
                    metrics.Je[iel, l, m, n] = metrics.dxdξ[iel, l, m, n]*(metrics.dydη[iel, l, m, n]*metrics.dzdζ[iel, l, m, n] - metrics.dydξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n])
                    metrics.Je[iel, l, m, n] += metrics.dydξ[iel, l, m, n]*(metrics.dxdζ[iel, l, m, n]*metrics.dzdη[iel, l, m, n] - metrics.dxdη[iel, l, m, n]*metrics.dzdζ[iel, l, m, n])
                    metrics.Je[iel, l, m, n] += metrics.dzdξ[iel, l, m, n]*(metrics.dxdη[iel, l, m, n]*metrics.dydζ[iel, l, m, n] - metrics.dxdζ[iel, l, m, n]*metrics.dydη[iel, l, m, n])
                    
                    metrics.dξdx[iel, l, m, n] =  (metrics.dydη[iel, l, m, n]*metrics.dzdζ[iel, l, m, n] - metrics.dydξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dξdy[iel, l, m, n] =  (metrics.dxdζ[iel, l, m, n]*metrics.dzdη[iel, l, m, n] - metrics.dxdη[iel, l, m, n]*metrics.dzdη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dξdz[iel, l, m, n] =  (metrics.dxdη[iel, l, m, n]*metrics.dydζ[iel, l, m, n] - metrics.dxdζ[iel, l, m, n]*metrics.dydη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dηdx[iel, l, m, n] =  (metrics.dydζ[iel, l, m, n]*metrics.dzdξ[iel, l, m, n] - metrics.dydξ[iel, l, m, n]*metrics.dzdζ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dηdy[iel, l, m, n] =  (metrics.dxdξ[iel, l, m, n]*metrics.dzdζ[iel, l, m, n] - metrics.dxdζ[iel, l, m, n]*metrics.dzdξ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dηdz[iel, l, m, n] =  (metrics.dxdζ[iel, l, m, n]*metrics.dydξ[iel, l, m, n] - metrics.dxdξ[iel, l, m, n]*metrics.dydζ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dζdx[iel, l, m, n] =  (metrics.dydξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n] - metrics.dydη[iel, l, m, n]*metrics.dzdξ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dζdy[iel, l, m, n] =  (metrics.dxdη[iel, l, m, n]*metrics.dzdξ[iel, l, m, n] - metrics.dxdξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dζdz[iel, l, m, n] =  (metrics.dxdξ[iel, l, m, n]*metrics.dydη[iel, l, m, n] - metrics.dxdη[iel, l, m, n]*metrics.dydξ[iel, l, m, n])/metrics.Je[iel, l, m, n] 
                
                 end
            end
        end
        show(stdout, "text/plain", metrics.Je[iel,:,:,:])
    end
    for iface = 1:size(mesh.xmin_faces,2) 
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.xmin_facetoelem[iface]
                metrics.Jef[iface, l, m] = metrics.dydη[iel, 1, l, m]*metrics.dzdζ[iel, 1, l, m] - metrics.dydζ[iel, 1, l, m]*metrics.dzdη[iel, 1, l, m]
            end
        end
    end
    disp = size(mesh.xmin_faces,2)
    for iface = 1:size(mesh.xmax_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.xmax_facetoelem[iface]
                metrics.Jef[iface+disp, l, m] = metrics.dydη[iel, Q+1, l, m]*metrics.dzdζ[iel, Q+1, l, m] - metrics.dydζ[iel, Q+1, l, m]*metrics.dzdη[iel, Q+1, l, m]
            end
        end
    end
    disp += size(mesh.xmax_faces,2)
    for iface = 1:size(mesh.ymin_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.ymin_facetoelem[iface]
                metrics.Jef[iface+disp, l, m] = metrics.dxdξ[iel, l, 1, m]*metrics.dzdζ[iel, l, 1, m] - metrics.dzdξ[iel, l, 1, m]*metrics.dxdζ[iel, l, 1, m]
            end
        end
    end
    disp += size(mesh.ymin_faces,2)
    for iface = 1:size(mesh.ymax_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.ymax_facetoelem[iface]
                metrics.Jef[iface+disp, l, m] = metrics.dxdξ[iel, l, Q+1, m]*metrics.dzdζ[iel, l, Q+1, m] - metrics.dzdξ[iel, l, Q+1, m]*metrics.dxdζ[iel, l, Q+1, m]
            end
        end
    end
    disp += size(mesh.ymax_faces,2)
    for iface = 1:size(mesh.zmin_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.zmin_facetoelem[iface]
                metrics.Jef[iface+disp, l, m] = metrics.dxdξ[iel, l, m, 1]*metrics.dydη[iel, l, m, 1] - metrics.dydξ[iel, l, m, 1]*metrics.dxdη[iel, l, m, 1]
            end
        end
    end
    disp += size(mesh.zmin_faces,2)
    for iface = 1:size(mesh.zmax_faces,2)
        for l = 1:Q+1
            for m = 1:Q+1
                iel = mesh.zmax_facetoelem[iface]
                metrics.Jef[iface+disp, l, m] = metrics.dxdξ[iel, l, m, Q+1]*metrics.dydη[iel, l, m, Q+1] - metrics.dydξ[iel, l, m, Q+1]*metrics.dxdη[iel, l, m, Q+1]
            end
        end
    end
     
    #show(stdout, "text/plain", metrics.Je)    
    return metrics
end

function build_metric_terms(SD::NSD_3D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange,N, Q, NGR, QGR, ξ, T;dir="x",side ="min")

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

function build_metric_terms(SD::NSD_3D, MT::CNVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T)
     metrics = St_metrics{T}(dxdξ = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂x/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdη = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂x/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dxdζ = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂x/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydξ = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂y/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dydη = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂y/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dydζ = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂y/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdξ = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂z/∂ξ[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdη = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂z/∂η[2, 1:Nq, 1:Nq, 1:nelem]
                             dzdζ = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂z/∂ζ[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdx = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ξdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdx = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ηdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdx = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ζdx[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdy = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ξdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdy = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ηdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdy = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ζdy[2, 1:Nq, 1:Nq, 1:nelem]
                             dξdz = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ξdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dηdz = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ηdz[2, 1:Nq, 1:Nq, 1:nelem]
                             dζdz = zeros(mesh.nelem, Q+1, Q+1, Q+1),  #∂ζdz[2, 1:Nq, 1:Nq, 1:nelem]

                             vⁱ = zeros(mesh.nelem, 3, Q+1, Q+1, Q+1), #contravariant unit vectors

                             Je   = zeros(mesh.nelem, Q+1, Q+1, Q+1))  #   Je[1:Nq, 1:Nq, 1:nelem]
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

                                ip = mesh.connijk[iel,i,j,k]

                                xijk = mesh.x[ip]
                                yijk = mesh.y[ip]
                                zijk = mesh.z[ip]

                                metrics.dxdξ[iel, l, m, n] = metrics.dxdξ[iel, l, m, n] + dψ[i,l]*ψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdη[iel, l, m, n] = metrics.dxdη[iel, l, m, n] + ψ[i,l]*dψ[j,m]*ψ[k,n]*xijk
                                metrics.dxdζ[iel, l, m, n] = metrics.dxdζ[iel, l, m, n] + ψ[i,l]*ψ[j,m]*dψ[k,n]*xijk

                                metrics.dydξ[iel, l, m, n] = metrics.dydξ[iel, l, m, n] + dψ[i,l]*ψ[j,m]*ψ[k,n]*yijk
                                metrics.dydη[iel, l, m, n] = metrics.dydη[iel, l, m, n] + ψ[i,l]*dψ[j,m]*ψ[k,n]*yijk
                                metrics.dydζ[iel, l, m, n] = metrics.dydζ[iel, l, m, n] + ψ[i,l]*ψ[j,m]*dψ[k,n]*yijk

                                metrics.dzdξ[iel, l, m, n] = metrics.dzdξ[iel, l, m, n] + dψ[i,l]*ψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdη[iel, l, m, n] = metrics.dzdη[iel, l, m, n] + ψ[i,l]*dψ[j,m]*ψ[k,n]*zijk
                                metrics.dzdζ[iel, l, m, n] = metrics.dzdζ[iel, l, m, n] + ψ[i,l]*ψ[j,m]*dψ[k,n]*zijk
                        #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                            end
                        end
                    end
                end
               # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[iel, k, l],  metrics.dxdη[iel, k, l], metrics.dydξ[iel, k, l],  metrics.dydη[iel, k, l]y )
            end
        end

        for l = 1:Q+1
            for m = 1:Q+1
                for n =1:Q+1
                    metrics.Je[iel, l, m, n] = metrics.dxdξ[iel, l, m, n]*(metrics.dydη[iel, l, m, n]*metrics.dzdζ[iel, l, m, n] - metrics.dydξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n])
                    metrics.Je[iel, l, m, n] += metrics.dydξ[iel, l, m, n]*(metrics.dxdζ[iel, l, m, n]*metrics.dzdη[iel, l, m, n] - metrics.dxdη[iel, l, m, n]*metrics.dzdζ[iel, l, m, n])
                    metrics.Je[iel, l, m, n] += metrics.dzdξ[iel, l, m, n]*(metrics.dxdη[iel, l, m, n]*metrics.dydζ[iel, l, m, n] - metrics.dxdζ[iel, l, m, n]*metrics.dydη[iel, l, m, n])
                    
                    metrics.dξdx[iel, l, m, n] =  (metrics.dydη[iel, l, m, n]*metrics.dzdζ[iel, l, m, n] - metrics.dydξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dξdy[iel, l, m, n] =  (metrics.dxdζ[iel, l, m, n]*metrics.dzdη[iel, l, m, n] - metrics.dxdη[iel, l, m, n]*metrics.dzdη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dξdz[iel, l, m, n] =  (metrics.dxdη[iel, l, m, n]*metrics.dydζ[iel, l, m, n] - metrics.dxdζ[iel, l, m, n]*metrics.dydη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dηdx[iel, l, m, n] =  (metrics.dydζ[iel, l, m, n]*metrics.dzdξ[iel, l, m, n] - metrics.dydξ[iel, l, m, n]*metrics.dzdζ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dηdy[iel, l, m, n] =  (metrics.dxdξ[iel, l, m, n]*metrics.dzdζ[iel, l, m, n] - metrics.dxdζ[iel, l, m, n]*metrics.dzdξ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dηdz[iel, l, m, n] =  (metrics.dxdζ[iel, l, m, n]*metrics.dydξ[iel, l, m, n] - metrics.dxdξ[iel, l, m, n]*metrics.dydζ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dζdx[iel, l, m, n] =  (metrics.dydξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n] - metrics.dydη[iel, l, m, n]*metrics.dzdξ[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dζdy[iel, l, m, n] =  (metrics.dxdη[iel, l, m, n]*metrics.dzdξ[iel, l, m, n] - metrics.dxdξ[iel, l, m, n]*metrics.dzdη[iel, l, m, n])/metrics.Je[iel, l, m, n]
                    metrics.dζdz[iel, l, m, n] =  (metrics.dxdξ[iel, l, m, n]*metrics.dydη[iel, l, m, n] - metrics.dxdη[iel, l, m, n]*metrics.dydξ[iel, l, m, n])/metrics.Je[iel, l, m, n] 

                    vⁱ[iel, 1, l, m, n] = metrics.dξdx[iel, l, m, n] + metrics.dξdy[iel, l, m, n]
                    vⁱ[iel, 2, l, m, n] = metrics.dηdx[iel, l, m, n] + metrics.dηdy[iel, l, m, n]
                    vⁱ[iel, 3, l, m, n] = metrics.dζdx[iel, l, m, n] + metrics.dζdy[iel, l, m, n]
                    
                 end
            end
        end
        show(stdout, "text/plain", metrics.Je[iel, :,:,:])
    end
    #show(stdout, "text/plain", metrics.Je)    
    return metrics
end
