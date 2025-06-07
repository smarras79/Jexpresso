Base.@kwdef mutable struct St_metrics{TFloat <: AbstractFloat, dims1, dims2, backend}

    dxdξ = KernelAbstractions.zeros(backend,TFloat, dims1)
    dxdη = KernelAbstractions.zeros(backend,TFloat, dims1)
    dxdζ = KernelAbstractions.zeros(backend,TFloat, dims1)
    
    dydξ = KernelAbstractions.zeros(backend,TFloat, dims1)
    dydη = KernelAbstractions.zeros(backend,TFloat, dims1)
    dydζ = KernelAbstractions.zeros(backend,TFloat, dims1)

    dzdξ = KernelAbstractions.zeros(backend,TFloat, dims1)
    dzdη = KernelAbstractions.zeros(backend,TFloat, dims1)
    dzdζ = KernelAbstractions.zeros(backend,TFloat, dims1)
    
    dξdx = KernelAbstractions.zeros(backend,TFloat, dims1)
    dξdy = KernelAbstractions.zeros(backend,TFloat, dims1)
    dξdz = KernelAbstractions.zeros(backend,TFloat, dims1)
    
    dηdx = KernelAbstractions.zeros(backend,TFloat, dims1)
    dηdy = KernelAbstractions.zeros(backend,TFloat, dims1)
    dηdz = KernelAbstractions.zeros(backend,TFloat, dims1)

    dζdx = KernelAbstractions.zeros(backend,TFloat, dims1)
    dζdy = KernelAbstractions.zeros(backend,TFloat, dims1)
    dζdz = KernelAbstractions.zeros(backend,TFloat, dims1)
    
    #
    # Element jacobian determinant
    #
    Je  = KernelAbstractions.zeros(backend,TFloat, dims1)
    Jef = KernelAbstractions.zeros(backend,TFloat, dims2)
    nx  = KernelAbstractions.zeros(backend,TFloat, dims2)
    ny  = KernelAbstractions.zeros(backend,TFloat, dims2)
    nz  = KernelAbstractions.zeros(backend,TFloat, dims2)
    
    dxdξ_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dxdη_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dxdζ_f = KernelAbstractions.zeros(backend,TFloat, dims2)

    dydξ_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dydη_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dydζ_f = KernelAbstractions.zeros(backend,TFloat, dims2)

    dzdξ_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dzdη_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dzdζ_f = KernelAbstractions.zeros(backend,TFloat, dims2)

    dξdx_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dξdy_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dξdz_f = KernelAbstractions.zeros(backend,TFloat, dims2)

    dηdx_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dηdy_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dηdz_f = KernelAbstractions.zeros(backend,TFloat, dims2)

    dζdx_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dζdy_f = KernelAbstractions.zeros(backend,TFloat, dims2)
    dζdz_f = KernelAbstractions.zeros(backend,TFloat, dims2)


    #
    # Contravariant arrays
    #
    vⁱ::Union{Array{TFloat}, Missing} = zeros(3) #contravariant unit vectors
    
end

function allocate_metrics(SD, nelem, nfaces_bdy, Q, T, backend)

    if SD == NSD_1D()
        dims1 = (nelem,      Q+1, 1, 1)
        dims2 = (nfaces_bdy, 1)
    elseif SD == NSD_2D()
        dims1 = (nelem,      Q+1, Q+1, 1)
        dims2 = (nfaces_bdy, Q+1, 1)
    elseif SD == NSD_3D()
        dims1 = (nelem,      Q+1, Q+1, Q+1)
        dims2 = (nfaces_bdy, Q+1, Q+1)
    end

    metrics = St_metrics{T, dims1, dims2, backend}()
    
    return metrics
end

function allocate_metrics_laguerre(SD, nelem, nfaces_bdy, Q, Qgr, T, backend)

    if SD == NSD_1D()
        dims1 = (nelem,      Q+1, 1, 1)
        dims2 = (nfaces_bdy, 1)
    elseif SD == NSD_2D()
        dims1 = (nelem,      Q+1, Qgr+1, 1)
        dims2 = (nfaces_bdy, Q+1, 1)
    elseif SD == NSD_3D()
        @mystop( " LAGUERRE ERROR: metric_terms.jl: Laguerre not implemented in 3D yet!")
        dims1 = (nelem,      Q+1, Q+1, Qgr+1)
        dims2 = (nfaces_bdy, Q+1, Q+1)
    end

    metrics = St_metrics{T, dims1, dims2, backend}()
    
    return metrics
end


function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T, MT::COVAR, SD::NSD_1D; backend = CPU())
#function build_metric_terms(SD::NSD_1D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T; backend = CPU())
    
    #metrics = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Q, T, backend)
    
    if (backend == CPU())
        for iel = 1:mesh.nelem
            for i = 1:N+1
                for k = 1:Q+1
                    metrics.dxdξ[iel, k, 1]  = mesh.Δx[iel]/2
                    metrics.Je[iel, k, 1]   = metrics.dxdξ[iel, k, 1]
                    metrics.dξdx[iel, k, 1] = 1.0/metrics.Je[iel, k, 1]
                end
            end        
        end
    else
        x = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
        connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem),N+1)
        Δx = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.nelem))
        KernelAbstractions.copyto!(backend, x, mesh.x)
        KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
        KernelAbstractions.copyto!(backend, Δx, mesh.Δx)
        k = build_1D_gpu_metrics!(backend,(N+1))
        k(metrics.dxdξ, metrics.Je, metrics.dξdx, basis.ψ, basis.dψ, x, connijk, Δx, Q; ndrange = (mesh.nelem*(N+1)), workgroupsize = (N+1))
    end
    
    return metrics
end

@kernel function build_1D_gpu_metrics!(dxdξ, Je, dξdx, ψ, dψ, x, connijk, Δx, Q)

    ie = @index(Group, Linear)
    il = @index(Local, Linear)
    T = eltype(x)
    for k=1:Q+1
        dxdξ[ie, k, 1] = Δx[ie]/2
        Je[ie, k, 1] = dxdξ[ie, k, 1]
        dξdx[ie, k, 1] = T(1.0)/Je[ie, k, 1]
    end

end

function build_metric_terms_1D_Laguerre!(metrics, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, inputs,T, MT::COVAR, SD::NSD_1D; backend = CPU())
#function build_metric_terms_1D_Laguerre(SD::NSD_1D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, inputs,T; backend = CPU())
    
    #metrics = allocate_metrics(SD, mesh.nelem_semi_inf, mesh.nedges_bdy, Q, T, backend)
    
    if (backend == CPU())
        dψ = basis.dψ
        for iel = 1:mesh.nelem_semi_inf
            for i = 1:mesh.ngr
                ip = mesh.connijk_lag[iel,i,1]
                xij = mesh.x[ip]
            
                for k = 1:mesh.ngr
                    metrics.dxdξ[iel, k,1]  += dψ[i,k] * (xij) * inputs[:yfac_laguerre]
                    metrics.Je[iel, k, 1]   = inputs[:yfac_laguerre]#abs(metrics.dxdξ[iel, k, 1])
                    if (xij > 0.1)
                        metrics.dξdx[iel, k, 1] = 1.0/metrics.Je[iel, k, 1]
                    else
                        metrics.dξdx[iel, k, 1] = -1.0/metrics.Je[iel, k, 1]
                    end
                end
            end
        end
    else
        x = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
        connijk_lag = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem_semi_inf),Int64(mesh.ngr))
        KernelAbstractions.copyto!(backend, x, mesh.x)
        KernelAbstractions.copyto!(backend, connijk_lag, mesh.connijk_lag)
        k = build_1D_gpu_metrics_laguerre!(backend,(Int64(mesh.ngr)))
        k(metrics.dxdξ, metrics.Je, metrics.dξdx, basis.ψ, basis.dψ, x, connijk_lag, TFloat(inputs[:yfac_laguerre]), Q; ndrange = (mesh.nelem_semi_inf*(mesh.ngr)), workgroupsize = (mesh.ngr))
    end

    return metrics
end

@kernel function build_1D_gpu_metrics_laguerre!(dxdξ, Je, dξdx, ψ, dψ, x, connijk_lag, yfac_laguerre, Q)

    ie = @index(Group, Linear)
    i = @index(Local, Linear)
    T = eltype(x)
    ip = connijk_lag[ie,i,1]
    xij = x[ip]
    for k=1:Q+1
        dxdξ[ie, k, 1] += dψ[i,k] * xij * yfac_laguerre
        Je[ie, k, 1] = yfac_laguerre
        if (xij > T(0.1))
            dξdx[ie, k, 1] = T(1.0)/Je[ie, k, 1]
        else
            dξdx[ie, k, 1] = -T(1.0)/Je[ie, k, 1]
        end
    end

end

#function build_metric_terms(SD::NSD_2D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T; backend = CPU())
function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T, MT::COVAR, SD::NSD_2D; backend = CPU())
    
    #metrics = allocate_metrics(SD, mesh.nelem, mesh.nedges_bdy, Q, T, backend)
    
    ψ  = @view(basis.ψ[:,:])
    dψ = @view(basis.dψ[:,:])
    if (backend == CPU())
        
        xij = 0.0
        yij = 0.0
        @inbounds for iel = 1:mesh.nelem
            for j = 1:N+1
                for i = 1:N+1

                    ip = mesh.connijk[iel, i, j]
                    xij = mesh.x[ip]
                    yij = mesh.y[ip]
                    
                    @turbo for l=1:Q+1
                        for k=1:Q+1

                            a = dψ[i,k]*ψ[j,l]
                            b = ψ[i,k]*dψ[j,l]
                            metrics.dxdξ[iel, k, l] += a * xij
                            metrics.dxdη[iel, k, l] += b * xij

                            metrics.dydξ[iel, k, l] += a * yij
                            metrics.dydη[iel, k, l] += b * yij
                            
                            #@printf(" i,j=%d, %d. x,y=%f,%f \n",i,j,xij, yij)
                        end
                    end
                end
            end
            
            @inbounds for l = 1:Q+1
                for k = 1:Q+1
                    
                    # Extract values from memory once per iteration
                    dxdξ_val = metrics.dxdξ[iel, k, l]
                    dydη_val = metrics.dydη[iel, k, l]
                    dydξ_val = metrics.dydξ[iel, k, l]
                    dxdη_val = metrics.dxdη[iel, k, l]
                    # Compute Je once and reuse its value
                    metrics.Je[iel, k, l] = dxdξ_val * dydη_val - dydξ_val * dxdη_val
                    
                    # Use the precomputed Je value for the other calculations
                    Jinv = 1.0/metrics.Je[iel, k, l]

                    metrics.dξdx[iel, k, l] =  dydη_val * Jinv
                    metrics.dξdy[iel, k, l] = -dxdη_val * Jinv
                    metrics.dηdx[iel, k, l] = -dydξ_val * Jinv
                    metrics.dηdy[iel, k, l] =  dxdξ_val * Jinv
                    
                end
            end
            #show(stdout, "text/plain", metrics.Je[iel, :,:])
        end
        
        nbdy_edges = size(mesh.poin_in_bdy_edge,1)
        @inbounds for iedge =1:nbdy_edges
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
                ip2 = mesh.poin_in_bdy_edge[iedge,1]
                ip3 = mesh.poin_in_bdy_edge[iedge,N+1]
                metrics.Jef[iedge, k] = sqrt((mesh.x[ip2]-mesh.x[ip3])^2+(mesh.y[ip2]-mesh.y[ip3])^2)/2
                comp1 = (x1-x2)/mag
                comp2 = (y1-y2)/mag
                metrics.nx[iedge, k] = comp2
                metrics.ny[iedge, k] = -comp1
            end
        end
    else
        x = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
        y = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
        connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem),N+1,N+1)
        KernelAbstractions.copyto!(backend, x, mesh.x)
        KernelAbstractions.copyto!(backend, y, mesh.y)
        KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
        k = build_2D_gpu_metrics!(backend,(N+1,N+1))
        k(metrics.dxdξ,metrics.dxdη,metrics.dydξ,metrics.dydη, ψ, dψ, x, y, connijk, Q; ndrange = (mesh.nelem*(N+1),mesh.ngl), workgroupsize = (N+1,N+1))
        metrics.Je   .= metrics.dxdξ.*metrics.dydη .- metrics.dydξ .* metrics.dxdη
        metrics.dξdx .= metrics.dydη ./ metrics.Je
        metrics.dξdy .= -metrics.dxdη ./ metrics.Je
        metrics.dηdx .= -metrics.dydξ ./ metrics.Je
        metrics.dηdy .= metrics.dxdξ ./ metrics.Je
        nbdy_edges    = size(mesh.poin_in_bdy_edge,1)
        poin_in_bdy_edge = KernelAbstractions.allocate(backend, TInt, Int64(nbdy_edges), N+1)
        KernelAbstractions.copyto!(backend, poin_in_bdy_edge,mesh.poin_in_bdy_edge)
        k = build_2D_gpu_bdy_metrics!(backend)
        k(metrics.Jef, metrics.nx, metrics.ny, x, y, poin_in_bdy_edge, N; ndrange = (nbdy_edges*(N+1)), workgroupsize = (N+1))
    end
    
    return metrics
    
end

@kernel function build_2D_gpu_metrics!(dxdξ, dxdη, dydξ, dydη, ψ, dψ, x, y, connijk, Q)
    s = Int32(@groupsize()[1])
    #n = div(@ndrange()[1],s)#div(length(A),s)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    ip = connijk[ie,i_x,i_y]
    xij = x[ip]
    yij = y[ip]
    for l=1:Q+1
        for k=1:Q+1
            KernelAbstractions.@atomic dxdξ[ie, k, l] += dψ[i_x,k]*ψ[i_y,l] * xij
            KernelAbstractions.@atomic dxdη[ie, k, l] += ψ[i_x,k]*dψ[i_y,l] * xij

            KernelAbstractions.@atomic dydξ[ie, k, l] += dψ[i_x,k]*ψ[i_y,l] * yij
            KernelAbstractions.@atomic dydη[ie, k, l] += ψ[i_x,k]*dψ[i_y,l] * yij
        end
    end
end

@kernel function build_3D_gpu_metrics!(dxdξ, dxdη, dxdζ, dydξ, dydη, dydζ, dzdξ, dzdη, dzdζ, ψ, dψ, x, y, z, connijk, Q)
    s = Int32(@groupsize()[1])
    #n = div(@ndrange()[1],s)#div(length(A),s)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    i_x = il[1]
    i_y = il[2]
    i_z = il[3]
    ip = connijk[ie,i_x,i_y,i_z]
    xijk = x[ip]
    yijk = y[ip]
    zijk = z[ip]
    for m=1:Q+1
        for l=1:Q+1
            for k=1:Q+1
                KernelAbstractions.@atomic dxdξ[ie, k, l, m] += dψ[i_x,k]*ψ[i_y,l] * ψ[i_z,m] * xijk
                KernelAbstractions.@atomic dxdη[ie, k, l, m] += ψ[i_x,k]*dψ[i_y,l] * ψ[i_z,m] * xijk
                KernelAbstractions.@atomic dxdζ[ie, k, l, m] += ψ[i_x,k]*ψ[i_y,l] * dψ[i_z,m] * xijk

                KernelAbstractions.@atomic dydξ[ie, k, l, m] += dψ[i_x,k]*ψ[i_y,l] * ψ[i_z,m] * yijk
                KernelAbstractions.@atomic dydη[ie, k, l, m] += ψ[i_x,k]*dψ[i_y,l] * ψ[i_z,m] * yijk
                KernelAbstractions.@atomic dydζ[ie, k, l, m] += ψ[i_x,k]*ψ[i_y,l] * dψ[i_z,m] * yijk

                KernelAbstractions.@atomic dzdξ[ie, k, l, m] += dψ[i_x,k]*ψ[i_y,l] * ψ[i_z,m] * zijk
                KernelAbstractions.@atomic dzdη[ie, k, l, m] += ψ[i_x,k]*dψ[i_y,l] * ψ[i_z,m] * zijk
                KernelAbstractions.@atomic dzdζ[ie, k, l, m] += ψ[i_x,k]*ψ[i_y,l] * dψ[i_z,m] * zijk
            end
        end
    end
end

@kernel function build_2D_gpu_bdy_metrics!(Jef, nx, ny, x, y, poin_in_bdy_edge,N)
    s = Int32(@groupsize()[1])
    #n = div(@ndrange()[1],s)#div(length(A),s)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    iedge = ie
    k = il[1]
    ip = poin_in_bdy_edge[iedge,k]
    if (k < N+1)
        ip1 = poin_in_bdy_edge[iedge,k+1]
    else
        ip1 = poin_in_bdy_edge[iedge,k-1]
    end
    x1 = x[ip]
    x2 = x[ip1]
    y1 = y[ip]
    y2 = y[ip1]
    mag = sqrt((x1-x2)^2+(y1-y2)^2)
    Jef[iedge, k] = mag/2
    comp1 = (x1-x2)/mag
    comp2 = (y1-y2)/mag
    nx[iedge, k] = comp2
    ny[iedge, k] = -comp1
end

@kernel function build_3D_gpu_bdy_metrics!(Jef, nx, ny, nz, x, y, z, poin_in_bdy_face,N)
    s = Int32(@groupsize()[1])
    #n = div(@ndrange()[1],s)#div(length(A),s)
    ie = @index(Group, Linear)
    il = @index(Local, NTuple)
    iface = ie
    i = il[1]
    j = il[2]
    ip = poin_in_bdy_face[iface,i,j]
    if (i < N+1)
        ip1 = poin_in_bdy_face[iface,i+1,j]
    else
        ip1 = poin_in_bdy_face[iface,i-1,j]
    end
    if (j < N+1)
        ip2 = poin_in_bdy_face[iface,i,j+1]
    else
        ip2 = poin_in_bdy_face[iface,i,j-1]
    end
    x1 = x[ip]
    x2 = x[ip1]
    x3 = x[ip2]
    y1 = y[ip]
    y2 = y[ip1]
    y3 = y[ip2]
    z1 = z[ip]
    z2 = z[ip1]
    z3 = z[ip2]
    a1 = x1 - x2
    a2 = y1 - y2
    a3 = z1 - z2
    b1 = x1 - x3
    b2 = y1 - y3
    b3 = z1 - z3
    comp1 = a2*b3 - a3*b2 
    comp2 = a3*b1 - a1*b3
    comp3 = a1*b2 - a2*b1
    mag = sqrt(comp1^2 + comp2^2 + comp3^2)
    if (mag < Float32(1e-6))
        mag = max(mag,abs(comp1),abs(comp2),abs(comp3),Float32(1e-7))
    end
    Jef[iface, i, j] = mag/2
    nx[iface, i, j] = comp1/mag
    ny[iface, i, j] = comp2/mag
    nz[iface, i, j] = comp3/mag
    if (abs(nx[iface,i,j]) < Float32(1e-2))
        nx[iface, i,j] = zero(Float32)
    end
    if (abs(ny[iface,i,j]) < Float32(1e-2))
        ny[iface, i,j] = zero(Float32)
    end
    if (abs(nz[iface,i,j]) < Float32(1e-2))
        nz[iface, i,j] = zero(Float32)
    end
end

function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange ,N, Q, NGR, QGR, ξ, ω1, ω2, T, MT::COVAR, SD::NSD_2D; backend = CPU())
#function build_metric_terms(SD::NSD_2D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange ,N, Q, NGR, QGR, ξ, ω1, ω2, T; backend = CPU())
    
    #metrics = allocate_metrics_laguerre(SD, mesh.nelem_semi_inf, mesh.nedges_bdy, Q, QGR, T, backend)
    
    ψ  = basis.ψ
    dψ = basis.dψ
    ψ1  = basisGR.ψ
    dψ1 = basisGR.dψ
    if ("Laguerre" in mesh.bdy_edge_type)
        if (backend == CPU())
            for iel=1:mesh.nelem_semi_inf
                for j=1:mesh.ngr
                    for i =1:mesh.ngl
                        ip = mesh.connijk_lag[iel,i,j]
                        xij = mesh.x[ip]
                        yij = mesh.y[ip]
                        for l=1:mesh.ngr
                            for k=1:mesh.ngl
                                if (inputs[:xfac_laguerre] == 0.0)
                                    metrics.dxdξ[iel, k, l] += dψ[i,k]*ψ1[j,l]*xij 
                                    metrics.dxdη[iel, k, l] +=  ψ[i,k]*inputs[:xfac_laguerre]/mesh.ngr#*dψ1[j,l]*xij
                                    metrics.dydξ[iel, k, l] += dψ[i,k]* ψ1[j,l]*yij
                                    metrics.dydη[iel, k, l] +=  ψ[i,k]*inputs[:yfac_laguerre]/mesh.ngr#*dψ1[j,l]*yij
                                else
                                    metrics.dxdξ[iel, k, l] += dψ[i,k]*ψ1[j,l]*xij
                                    metrics.dxdη[iel, k, l] +=  ψ[i,k]*inputs[:xfac_laguerre]/mesh.ngr#*dψ1[j,l]*xij
                                    metrics.dydξ[iel, k, l] += dψ[i,k]* ψ1[j,l]*yij
                                    metrics.dydη[iel, k, l] +=  ψ[i,k]*inputs[:yfac_laguerre]/mesh.ngr#*dψ1[j,l]*yij
                                end   
                            end
                        end
                        #@info metrics.dxdξ[iel, k, l],metrics.dxdη[iel, k, l],metrics.dydξ[iel, k, l],metrics.dydη[iel, k, l],xij,yij
                    end
                end
                for l = 1:mesh.ngr
                    for k = 1:mesh.ngl
                        ip = mesh.connijk_lag[iel,k,l]
                        #xij = mesh.x[ip]
                        #yij = mesh.y[ip]
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
        else
            x = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
            y = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
            connijk_lag = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem_semi_inf),mesh.ngl,mesh.ngr)

            KernelAbstractions.copyto!(backend, x, mesh.x)
            KernelAbstractions.copyto!(backend, y, mesh.y)
            KernelAbstractions.copyto!(backend, connijk_lag, mesh.connijk_lag)
            k = build_2D_gpu_metrics_lag!(backend)
            k(metrics.dxdξ, metrics.dxdη, metrics.dydξ, metrics.dydη, ψ, dψ, ψ1, dψ1, x, y, connijk_lag, mesh.ngl, mesh.ngr, TFloat(inputs[:xfac_laguerre]), TFloat(inputs[:yfac_laguerre]);
              ndrange = (mesh.nelem_semi_inf * mesh.ngl, mesh.ngr),
              workgroupsize = (mesh.ngl, mesh.ngr))
            metrics.Je .= metrics.dxdξ.*metrics.dydη .- metrics.dydξ .* metrics.dxdη
            metrics.dξdx .= metrics.dydη ./ metrics.Je
            metrics.dξdy .= -metrics.dxdη ./ metrics.Je
            metrics.dηdx .= -metrics.dydξ ./ metrics.Je
            metrics.dηdy .= metrics.dxdξ ./ metrics.Je
        end
    end

    return metrics
end

@kernel function build_2D_gpu_metrics_lag!(dxdξ, dxdη, dydξ, dydη, ψ, dψ, ψ1, dψ1, x, y, connijk_lag, ngl, ngr, xfac_lag, yfac_lag)
    s = Int32(@groupsize()[1])
    iel = @index(Group, Linear)
    il = @index(Local, NTuple)
    i = il[1]
    j = il[2]
    ip = connijk_lag[iel,i,j]
    xij = x[ip]
    yij = y[ip]
    for l=1:ngr
        for k=1:ngl
            KernelAbstractions.@atomic dxdξ[iel, k, l] += dψ[i,k]*ψ1[j,l]*xij 
            KernelAbstractions.@atomic dxdη[iel, k, l] +=  ψ[i,k]*xfac_lag/ngr#*dψ1[j,l]*xij
            KernelAbstractions.@atomic dydξ[iel, k, l] += dψ[i,k]* ψ1[j,l]*yij
            KernelAbstractions.@atomic dydη[iel, k, l] +=  ψ[i,k]*yfac_lag/ngr#*dψ1[j,l]*yij
        end   
    end
end

function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T, MT::COVAR, SD::NSD_3D; backend = CPU())
    
#function build_metric_terms(SD::NSD_3D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T; backend = CPU())

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)
    
    #metrics = allocate_metrics(SD, mesh.nelem, mesh.nfaces_bdy, Q, T, backend)
    
    ψ  = @view(basis.ψ[:,:])
    dψ = @view(basis.dψ[:,:])

    println_rank(" # 3D metric terms "; msg_rank = rank, suppress = mesh.msg_suppress)
    
    if (backend == CPU())
        
        xijk = 0.0
        yijk = 0.0
        zijk = 0.0
        @inbounds for iel = 1:mesh.nelem
            
            for k = 1:N+1
                for j = 1:N+1
                    for i = 1:N+1

                        ip = mesh.connijk[iel, i, j, k]
                        xijk = mesh.x[ip]
                        yijk = mesh.y[ip]
                        zijk = mesh.z[ip]
                        
                        @turbo for n = 1:Q+1
                            for m = 1:Q+1
                                for l = 1:Q+1
                                    
                                    a = dψ[i,l]* ψ[j,m]* ψ[k,n]
                                    b = ψ[i,l]*dψ[j,m]* ψ[k,n]
                                    c = ψ[i,l]* ψ[j,m]*dψ[k,n]
                                    
                                    metrics.dxdξ[iel, l, m, n] += a*xijk
                                    metrics.dxdη[iel, l, m, n] += b*xijk
                                    metrics.dxdζ[iel, l, m, n] += c*xijk

                                    metrics.dydξ[iel, l, m, n] += a*yijk
                                    metrics.dydη[iel, l, m, n] += b*yijk
                                    metrics.dydζ[iel, l, m, n] += c*yijk

                                    metrics.dzdξ[iel, l, m, n] += a*zijk
                                    metrics.dzdη[iel, l, m, n] += b*zijk
                                    metrics.dzdζ[iel, l, m, n] += c*zijk
                                end
                            end
                        end
                    end
                    # @printf(" dxdξ=%f, dxdη=%f, dydξ=%f dydη=%f \n",  metrics.dxdξ[iel, k, l],  metrics.dxdη[iel, k, l], metrics.dydξ[iel, k, l],  metrics.dydη[iel, k, l] )
                end
            end
            
            
            @inbounds for l = 1:Q+1
                for m = 1:Q+1
                    for n =1:Q+1

                        dxdξ = metrics.dxdξ[iel, l, m, n]
                        dydη = metrics.dydη[iel, l, m, n]
                        dzdζ = metrics.dzdζ[iel, l, m, n]
                        dydξ = metrics.dydξ[iel, l, m, n]
                        dzdη = metrics.dzdη[iel, l, m, n]
                        dxdζ = metrics.dxdζ[iel, l, m, n]
                        dxdη = metrics.dxdη[iel, l, m, n]
                        dydζ = metrics.dydζ[iel, l, m, n]
                        dzdξ = metrics.dzdξ[iel, l, m, n]

                        # Calculate metrics.Je[iel, l, m, n] with inlined expressions
                        metrics.Je[iel, l, m, n] = dxdξ * (dydη * dzdζ - dydξ * dzdη) +
                                                   dydξ * (dxdζ * dzdη - dxdη * dzdζ) +
                                                   dzdξ * (dxdη * dydζ - dxdζ * dydη)
                        
                        Jinv = 1.0/metrics.Je[iel, l, m, n]
                        
                        metrics.dξdx[iel, l, m, n] = (dydη*dzdζ- dydξ*dzdη)*Jinv
                        metrics.dξdy[iel, l, m, n] = (dxdζ*dzdη- dxdη*dzdη)*Jinv
                        metrics.dξdz[iel, l, m, n] = (dxdη*dydζ- dxdζ*dydη)*Jinv
                        metrics.dηdx[iel, l, m, n] = (dydζ*dzdξ- dydξ*dzdζ)*Jinv
                        metrics.dηdy[iel, l, m, n] = (dxdξ*dzdζ- dxdζ*dzdξ)*Jinv
                        metrics.dηdz[iel, l, m, n] = (dxdζ*dydξ- dxdξ*dydζ)*Jinv
                        metrics.dζdx[iel, l, m, n] = (dydξ*dzdη- dydη*dzdξ)*Jinv
                        metrics.dζdy[iel, l, m, n] = (dxdη*dzdξ- dxdξ*dzdη)*Jinv
                        metrics.dζdz[iel, l, m, n] = (dxdξ*dydη- dxdη*dydξ)*Jinv
                        
                    end
                end
            end

        end
        @inbounds for iface=1:mesh.nfaces_bdy
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    ip = mesh.poin_in_bdy_face[iface,i,j]
                    x1 = mesh.x[ip]
                    y1 = mesh.y[ip]
                    z1 = mesh.z[ip]
                    @turbo for k=1:mesh.ngl
                        for l=1:mesh.ngl


                            a = dψ[i,k]*ψ[j,l]
                            b = ψ[i,k]*dψ[j,l]

                            metrics.dxdξ_f[iface, k, l] += a * x1
                            metrics.dxdη_f[iface, k, l] += b * x1

                            metrics.dydξ_f[iface, k, l] += a * y1
                            metrics.dydη_f[iface, k, l] += b * y1

                            metrics.dzdξ_f[iface, k, l] += a * z1
                            metrics.dzdη_f[iface, k, l] += b * z1

                        end
                    end
                        
                end
            end
            for i=1:mesh.ngl
                for j=1:mesh.ngl
                    ip = mesh.poin_in_bdy_face[iface,i,j]
                    if (i < N+1)
                        ip1 = mesh.poin_in_bdy_face[iface,i+1,j]
                    else
                        ip1 = mesh.poin_in_bdy_face[iface,i-1,j]
                    end
                    if (j < N+1)
                        ip2 = mesh.poin_in_bdy_face[iface,i,j+1]
                    else
                        ip2 = mesh.poin_in_bdy_face[iface,i,j-1]
                    end
                    x1 = mesh.x[ip]
                    x2 = mesh.x[ip1]
                    x3 = mesh.x[ip2]
                    y1 = mesh.y[ip]
                    y2 = mesh.y[ip1]
                    y3 = mesh.y[ip2]
                    z1 = mesh.z[ip]
                    z2 = mesh.z[ip1]
                    z3 = mesh.z[ip2]
                    a1 = x1 - x2
                    a2 = y1 - y2
                    a3 = z1 - z2
                    b1 = x1 - x3
                    b2 = y1 - y3
                    b3 = z1 - z3
                    comp1 = a2*b3 - a3*b2
                    comp2 = a3*b1 - a1*b3
                    comp3 = a1*b2 - a2*b1
                    mag    = sqrt(comp1^2 + comp2^2 + comp3^2)
                    maginv = 1.0/mag
                    # Extract values from memory once per iteration
                    dxdξ = metrics.dxdξ_f[iface, i, j]
                    dydη = metrics.dydη_f[iface, i, j]
                    dydξ = metrics.dydξ_f[iface, i, j]
                    dxdη = metrics.dxdη_f[iface, i, j]
                    dzdξ = metrics.dzdξ_f[iface, i, j]
                    dzdη = metrics.dzdη_f[iface, i, j]
                    # Compute Je once and reuse its value
                    metrics.Jef[iface, i, j] = dxdξ * (dydη - dydξ * dzdη) +
                                             dydξ * (dzdη - dxdη) +
                                             dzdξ * (dxdη  - dydη) 

                    # Use the precomputed Je value for the other calculations
                    Jinv = 1.0/metrics.Jef[iface, i, j]

                    metrics.dξdx_f[iface, i, j] = (dydη - dydξ*dzdη)*Jinv
                    metrics.dξdy_f[iface, i, j] = (dzdη - dxdη*dzdη)*Jinv
                    metrics.dξdz_f[iface, i, j] = (dxdη - dydη)*Jinv
                    metrics.dηdx_f[iface, i, j] = (dzdξ - dydξ)*Jinv
                    metrics.dηdy_f[iface, i, j] = (dxdξ - dzdξ)*Jinv
                    metrics.dηdz_f[iface, i, j] = (dydξ - dxdξ)*Jinv

                    metrics.nx[iface, i, j] = comp1*maginv
                    metrics.ny[iface, i, j] = comp2*maginv
                    metrics.nz[iface, i, j] = comp3*maginv
                end
            end
        end

    else
        x = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
        y = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
        z = KernelAbstractions.allocate(backend, TFloat, Int64(mesh.npoin))
        connijk = KernelAbstractions.allocate(backend, TInt, Int64(mesh.nelem),N+1,N+1,N+1)
        KernelAbstractions.copyto!(backend, x, mesh.x)
        KernelAbstractions.copyto!(backend, y, mesh.y)
        KernelAbstractions.copyto!(backend, z, mesh.z)
        KernelAbstractions.copyto!(backend, connijk, mesh.connijk)
        k = build_3D_gpu_metrics!(backend,(N+1,N+1,N+1))
        k(metrics.dxdξ,metrics.dxdη,metrics.dxdζ,metrics.dydξ,metrics.dydη,metrics.dydζ,metrics.dzdξ,metrics.dzdη,metrics.dzdζ, ψ, dψ, x, y, z, connijk, Q;
          ndrange = (mesh.nelem*(N+1),mesh.ngl,mesh.ngl), workgroupsize = (N+1,N+1,N+1))
        metrics.Je .= metrics.dxdξ.*(metrics.dydη.*metrics.dzdζ .- metrics.dydξ.*metrics.dzdη)
        metrics.Je .+= metrics.dydξ.*(metrics.dxdζ.*metrics.dzdη .- metrics.dxdη.*metrics.dzdζ)
        metrics.Je .+= metrics.dzdξ.*(metrics.dxdη.*metrics.dydζ .- metrics.dxdζ.*metrics.dydη)
        metrics.dξdx .= (metrics.dydη.*metrics.dzdζ .- metrics.dydξ.*metrics.dzdη) ./ metrics.Je
        metrics.dξdy .= (metrics.dxdζ.*metrics.dzdη .- metrics.dxdη.*metrics.dzdη) ./ metrics.Je
        metrics.dξdz .= (metrics.dxdη.*metrics.dydζ .- metrics.dxdζ.*metrics.dydη) ./ metrics.Je
        metrics.dηdx .= (metrics.dydζ.*metrics.dzdξ .- metrics.dydξ.*metrics.dzdζ) ./ metrics.Je
        metrics.dηdy .= (metrics.dxdξ.*metrics.dzdζ .- metrics.dxdζ.*metrics.dzdξ) ./ metrics.Je
        metrics.dηdz .= (metrics.dxdζ.*metrics.dydξ .- metrics.dxdξ.*metrics.dydζ) ./ metrics.Je
        metrics.dζdx .= (metrics.dydξ.*metrics.dzdη .- metrics.dydη.*metrics.dzdξ) ./ metrics.Je
        metrics.dζdy .= (metrics.dxdη.*metrics.dzdξ .- metrics.dxdξ.*metrics.dzdη) ./ metrics.Je
        metrics.dζdz .= (metrics.dxdξ.*metrics.dydη .- metrics.dxdη.*metrics.dydξ) ./ metrics.Je
        nbdy_faces = size(mesh.poin_in_bdy_face,1)
        poin_in_bdy_face = KernelAbstractions.allocate(backend, TInt, Int64(nbdy_faces), N+1,N+1)
        KernelAbstractions.copyto!(backend, poin_in_bdy_face,mesh.poin_in_bdy_face)
        k = build_3D_gpu_bdy_metrics!(backend)
        k(metrics.Jef, metrics.nx, metrics.ny, metrics.nz, x, y, z, poin_in_bdy_face, N; ndrange = (nbdy_faces*(N+1),N+1), workgroupsize = (N+1,N+1))
        #show(stdout, "text/plain", metrics.Je[iel,:,:,:])
    end

    #show(stdout, "text/plain", metrics.Je)    
    return metrics
end



#function build_metric_terms(SD::NSD_3D, MT::COVAR, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange,N, Q, NGR, QGR, ξ, T;dir="x",side ="min")
   
function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange,N, Q, NGR, QGR, ξ, T, MT::COVAR, SD::NSD_3D; dir="x",side ="min")
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    
    #metrics = allocate_metrics(SD, mesh.nelem, mesh.nfaces_bdy, Q, T, backend)
    
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

    #@info " COVARIANT metric terms WIP"
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
        #show(stdout, "text/plain", metrics.Je[:,:,:,iel])
    end
end
