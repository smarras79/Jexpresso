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
    
    if (backend == CPU())
        dψ = basis.dψ
        @inbounds for iel = 1:mesh.nelem_semi_inf  # PERF: Added @inbounds
            for i = 1:mesh.ngr
                ip = mesh.connijk_lag[iel,i,1]
                xij = mesh.x[ip]

                for k = 1:mesh.ngr
                    metrics.dxdξ[iel, k,1]  += dψ[i,k] * (xij) * inputs[:yfac_laguerre]
                    metrics.Je[iel, k, 1]   = inputs[:yfac_laguerre]#abs(metrics.dxdξ[iel, k, 1])
                    if (xij > 0.1)
                        metrics.dξdx[iel, k, 1] = T(1.0)/metrics.Je[iel, k, 1]  # FIXED: use type parameter T
                    else
                        metrics.dξdx[iel, k, 1] = -T(1.0)/metrics.Je[iel, k, 1]  # FIXED: use type parameter T
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

function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T, MT::COVAR, SD::NSD_1D; backend = CPU())
    
    if (backend == CPU())
        @inbounds for iel = 1:mesh.nelem  # PERF: Added @inbounds
            for i = 1:N+1
                for k = 1:Q+1
                    metrics.dxdξ[iel, k, 1]  = mesh.Δx[iel]/2
                    metrics.Je[iel, k, 1]   = metrics.dxdξ[iel, k, 1]
                    metrics.dξdx[iel, k, 1] = T(1.0)/metrics.Je[iel, k, 1]  # FIXED: use type parameter T
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

function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T, MT::COVAR, SD::NSD_2D; backend = CPU())

    ψ  = @view(basis.ψ[:,:])
    dψ = @view(basis.dψ[:,:])

    if (backend == CPU())
        # Pre-allocate temporary variables outside loops to avoid repeated allocation
        xij = 0.0
        yij = 0.0
        
        @inbounds for iel = 1:mesh.nelem
            # Cache views to avoid repeated indexing overhead
            dxdξ_iel = @view metrics.dxdξ[iel, :, :]
            dxdη_iel = @view metrics.dxdη[iel, :, :]
            dydξ_iel = @view metrics.dydξ[iel, :, :]
            dydη_iel = @view metrics.dydη[iel, :, :]
            connijk_iel = @view mesh.connijk[iel, :, :]
            
            for j = 1:N+1
                for i = 1:N+1
                    ip = connijk_iel[i, j]
                    xij = mesh.x[ip]
                    yij = mesh.y[ip]
                    
                    # Unroll and optimize the inner loops
                    @turbo for l=1:Q+1
                        dψ_j_l = dψ[j,l]
                        ψ_j_l = ψ[j,l]
                        for k=1:Q+1
                            a = dψ[i,k] * ψ_j_l
                            b = ψ[i,k] * dψ_j_l
                            dxdξ_iel[k, l] += a * xij
                            dxdη_iel[k, l] += b * xij
                            dydξ_iel[k, l] += a * yij
                            dydη_iel[k, l] += b * yij
                        end
                    end
                end
            end
            
            # Second loop with cached views and optimized calculations
            Je_iel = @view metrics.Je[iel, :, :]
            dξdx_iel = @view metrics.dξdx[iel, :, :]
            dξdy_iel = @view metrics.dξdy[iel, :, :]
            dηdx_iel = @view metrics.dηdx[iel, :, :]
            dηdy_iel = @view metrics.dηdy[iel, :, :]
            
            @turbo for l = 1:Q+1
                for k = 1:Q+1
                    # Compute Jacobian determinant
                    Je_val = dxdξ_iel[k, l] * dydη_iel[k, l] - dydξ_iel[k, l] * dxdη_iel[k, l]
                    Je_iel[k, l] = Je_val

                    # Compute inverse Jacobian components using single division
                    Jinv = T(1.0) / Je_val  # FIXED: use type parameter T
                    dξdx_iel[k, l] =  dydη_iel[k, l] * Jinv
                    dξdy_iel[k, l] = -dxdη_iel[k, l] * Jinv
                    dηdx_iel[k, l] = -dydξ_iel[k, l] * Jinv
                    dηdy_iel[k, l] =  dxdξ_iel[k, l] * Jinv
                end
            end
        end
        
        # Optimize boundary edge calculations
        nbdy_edges = size(mesh.poin_in_bdy_edge, 1)
        @inbounds for iedge = 1:nbdy_edges
            poin_edge = @view mesh.poin_in_bdy_edge[iedge, :]
            
            # Pre-compute edge endpoints for Jef calculation
            ip_first = poin_edge[1]
            ip_last = poin_edge[N+1]
            edge_length = sqrt((mesh.x[ip_first] - mesh.x[ip_last])^2 + 
                (mesh.y[ip_first] - mesh.y[ip_last])^2)
            Jef_val = edge_length * 0.5  # Avoid division by 2
            
            for k = 1:N+1
                ip = poin_edge[k]
                
                # Determine next/previous point more efficiently
                ip1 = (k < N+1) ? poin_edge[k+1] : poin_edge[k-1]
                
                # Cache coordinates
                x1, y1 = mesh.x[ip], mesh.y[ip]
                x2, y2 = mesh.x[ip1], mesh.y[ip1]
                
                # Compute normal vector components
                dx, dy = x1 - x2, y1 - y2
                mag_inv = 1.0 / sqrt(dx*dx + dy*dy)  # Use single sqrt and invert
                
                # Store results
                metrics.Jef[iedge, k] = Jef_val
                metrics.nx[iedge, k] = dy * mag_inv
                metrics.ny[iedge, k] = -dx * mag_inv
                e = mesh.bdy_edge_in_elem[iedge]
                ip2 = mesh.connijk[e,2,2]
                idx1 = 0
                idx2 = 0
                #=for j=1:N+1
                    for i=1:N+1
                        if (mesh.connijk[e,i,j] == ip)
                            idx1 = i
                            idx2 = j
                        end
                    end
                end=#
                #if (idx1 + metrics.nx[iedge, k] < 1 || idx1 + metrics.nx[iedge, k] > N+1 || idx2 + metrics.ny[iedge, k] < 1 || idx2 + metrics.ny[iedge, k] > N+1)
                if (metrics.nx[iedge, k]*(mesh.x[ip2]-mesh.x[ip]) + metrics.ny[iedge, k]*(mesh.y[ip2] -mesh.y[ip]) > 0)
                    metrics.nx[iedge, k] = - metrics.nx[iedge, k]
                    metrics.ny[iedge, k] = - metrics.ny[iedge, k]
                end

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
end

function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, N, Q, ξ, ω, T, MT::COVAR, SD::NSD_3D; backend = CPU())
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)
    
    ψ  = @view(basis.ψ[:,:])
    dψ = @view(basis.dψ[:,:])

    println_rank(" # 3D metric terms "; msg_rank = rank, suppress = mesh.msg_suppress)
    
    if (backend == CPU())
        
        # Precompute frequently used values
        N1 = N + 1
        Q1 = Q + 1
        ngl = mesh.ngl
        
        # Pre-allocate temporary arrays for better memory access patterns
        temp_coords = Vector{NTuple{3,Float64}}(undef, N1*N1*N1)
        temp_basis = Matrix{Float64}(undef, Q1, 3)  # For storing ψ and dψ values
        
        @inbounds for iel = 1:mesh.nelem
            
            # Cache all coordinate data for current element upfront
            connijk_iel = @view mesh.connijk[iel, :, :, :]
            coord_idx = 1
            for k = 1:N1, j = 1:N1, i = 1:N1
                ip = connijk_iel[i, j, k]
                temp_coords[coord_idx] = (mesh.x[ip], mesh.y[ip], mesh.z[ip])
                coord_idx += 1
            end
            
            # Cache all metric views for current element
            dxdξ_iel = @view metrics.dxdξ[iel, :, :, :]
            dxdη_iel = @view metrics.dxdη[iel, :, :, :]
            dxdζ_iel = @view metrics.dxdζ[iel, :, :, :]
            dydξ_iel = @view metrics.dydξ[iel, :, :, :]
            dydη_iel = @view metrics.dydη[iel, :, :, :]
            dydζ_iel = @view metrics.dydζ[iel, :, :, :]
            dzdξ_iel = @view metrics.dzdξ[iel, :, :, :]
            dzdη_iel = @view metrics.dzdη[iel, :, :, :]
            dzdζ_iel = @view metrics.dzdζ[iel, :, :, :]
            
            # Optimized triple loop with better memory access
            coord_idx = 1
            for k = 1:N1, j = 1:N1, i = 1:N1
                xijk, yijk, zijk = temp_coords[coord_idx]
                coord_idx += 1
                
                # Precompute basis function values for current (i,j,k)
                @simd for idx = 1:Q1
                    temp_basis[idx, 1] = dψ[i, idx]  # dψ_i
                    temp_basis[idx, 2] =  ψ[j, idx]  # ψ_j  
                    temp_basis[idx, 3] = dψ[j, idx]  # dψ_j
                end
                
                # More cache-friendly nested loops
                @turbo for n = 1:Q1
                    ψ_k_n = ψ[k, n]
                    dψ_k_n = dψ[k, n]
                    for m = 1:Q1
                        ψ_j_m = temp_basis[m, 2]  # ψ[j, m]
                        dψ_j_m = temp_basis[m, 3]  # dψ[j, m]
                        for l = 1:Q1
                            dψ_i_l = temp_basis[l, 1]  # dψ[i, l]
                            ψ_i_l = ψ[i, l]
                            
                            # Compute coefficients once
                            a = dψ_i_l * ψ_j_m * ψ_k_n
                            b = ψ_i_l * dψ_j_m * ψ_k_n
                            c = ψ_i_l * ψ_j_m * dψ_k_n
                            
                            # Vectorized updates
                            dxdξ_iel[l, m, n] += a * xijk
                            dxdη_iel[l, m, n] += b * xijk
                            dxdζ_iel[l, m, n] += c * xijk

                            dydξ_iel[l, m, n] += a * yijk
                            dydη_iel[l, m, n] += b * yijk
                            dydζ_iel[l, m, n] += c * yijk

                            dzdξ_iel[l, m, n] += a * zijk
                            dzdη_iel[l, m, n] += b * zijk
                            dzdζ_iel[l, m, n] += c * zijk
                        end
                    end
                end
            end
            
            # Optimized Jacobian calculations with better memory access
            Je_iel   = @view metrics.Je[iel, :, :, :]
            dξdx_iel = @view metrics.dξdx[iel, :, :, :]
            dξdy_iel = @view metrics.dξdy[iel, :, :, :]
            dξdz_iel = @view metrics.dξdz[iel, :, :, :]
            dηdx_iel = @view metrics.dηdx[iel, :, :, :]
            dηdy_iel = @view metrics.dηdy[iel, :, :, :]
            dηdz_iel = @view metrics.dηdz[iel, :, :, :]
            dζdx_iel = @view metrics.dζdx[iel, :, :, :]
            dζdy_iel = @view metrics.dζdy[iel, :, :, :]
            dζdz_iel = @view metrics.dζdz[iel, :, :, :]
            
            @turbo for n = 1:Q1, m = 1:Q1, l = 1:Q1
                # Load derivatives once with better naming
                dxdξ = dxdξ_iel[l, m, n]
                dydη = dydη_iel[l, m, n]
                dzdζ = dzdζ_iel[l, m, n]
                dydξ = dydξ_iel[l, m, n]
                dzdη = dzdη_iel[l, m, n]
                dxdζ = dxdζ_iel[l, m, n]
                dxdη = dxdη_iel[l, m, n]
                dydζ = dydζ_iel[l, m, n]
                dzdξ = dzdξ_iel[l, m, n]

                # Compute cross products first
                cross1 = dydη * dzdζ - dydζ * dzdη
                cross2 = dxdζ * dzdη - dxdη * dzdζ
                cross3 = dxdη * dydζ - dxdζ * dydη
                cross4 = dydζ * dzdξ - dydξ * dzdζ
                cross5 = dxdξ * dzdζ - dxdζ * dzdξ
                cross6 = dxdζ * dydξ - dxdξ * dydζ
                cross7 = dydξ * dzdη - dydη * dzdξ
                cross8 = dxdη * dzdξ - dxdξ * dzdη
                cross9 = dxdξ * dydη - dxdη * dydξ

                # Calculate Jacobian determinant using precomputed cross products
                Je_val = dxdξ * cross1 + dydξ * cross2 + dzdξ * cross3
                
                Je_iel[l, m, n] = Je_val
                Jinv = 1.0 / Je_val
                
                # Calculate inverse Jacobian terms using precomputed values
                dξdx_iel[l, m, n] = cross1 * Jinv
                dξdy_iel[l, m, n] = cross2 * Jinv
                dξdz_iel[l, m, n] = cross3 * Jinv
                dηdx_iel[l, m, n] = cross4 * Jinv
                dηdy_iel[l, m, n] = cross5 * Jinv
                dηdz_iel[l, m, n] = cross6 * Jinv
                dζdx_iel[l, m, n] = cross7 * Jinv
                dζdy_iel[l, m, n] = cross8 * Jinv
                dζdz_iel[l, m, n] = cross9 * Jinv
            end
        end
        
        # Optimized boundary face calculations with better memory management
        temp_face_coords = Vector{NTuple{3,Float64}}(undef, ngl*ngl)
        
        @inbounds for iface = 1:mesh.nfaces_bdy
            # Cache all face coordinate data upfront
            poin_face = @view mesh.poin_in_bdy_face[iface, :, :]
            coord_idx = 1
            for j = 1:ngl, i = 1:ngl
                ip = poin_face[i, j]
                temp_face_coords[coord_idx] = (mesh.x[ip], mesh.y[ip], mesh.z[ip])
                coord_idx += 1
            end
            
            # Cache views for current face
            dxdξ_f_face = @view metrics.dxdξ_f[iface, :, :]
            dxdη_f_face = @view metrics.dxdη_f[iface, :, :]
            dydξ_f_face = @view metrics.dydξ_f[iface, :, :]
            dydη_f_face = @view metrics.dydη_f[iface, :, :]
            dzdξ_f_face = @view metrics.dzdξ_f[iface, :, :]
            dzdη_f_face = @view metrics.dzdη_f[iface, :, :]
            
            coord_idx = 1
            for j = 1:ngl, i = 1:ngl
                x1, y1, z1 = temp_face_coords[coord_idx]
                coord_idx += 1
                
                @turbo for l = 1:ngl, k = 1:ngl
                    dψ_i_k = dψ[i, k]
                    ψ_i_k  = ψ[i, k]
                    ψ_j_l  = ψ[j, l]
                    dψ_j_l = dψ[j, l]
                    
                    a = dψ_i_k * ψ_j_l
                    b = ψ_i_k * dψ_j_l

                    dxdξ_f_face[k, l] += a * x1
                    dxdη_f_face[k, l] += b * x1
                    dydξ_f_face[k, l] += a * y1
                    dydη_f_face[k, l] += b * y1
                    dzdξ_f_face[k, l] += a * z1
                    dzdη_f_face[k, l] += b * z1
                end
            end
            
            # Optimized second loop with precomputed values
            Jef_face = @view metrics.Jef[iface, :, :]
            dξdx_f_face = @view metrics.dξdx_f[iface, :, :]
            dξdy_f_face = @view metrics.dξdy_f[iface, :, :]
            dξdz_f_face = @view metrics.dξdz_f[iface, :, :]
            dηdx_f_face = @view metrics.dηdx_f[iface, :, :]
            dηdy_f_face = @view metrics.dηdy_f[iface, :, :]
            dηdz_f_face = @view metrics.dηdz_f[iface, :, :]
            nx_face = @view metrics.nx[iface, :, :]
            ny_face = @view metrics.ny[iface, :, :]
            nz_face = @view metrics.nz[iface, :, :]
            
            coord_idx = 1
            for j = 1:ngl, i = 1:ngl
                x1, y1, z1 = temp_face_coords[coord_idx]
                coord_idx += 1
                
                # More efficient neighbor point determination
                i_neighbor = (i < N1) ? i + 1 : i - 1
                j_neighbor = (j < N1) ? j + 1 : j - 1
                ip = poin_face[i,j]
                ip1 = poin_face[i_neighbor, j]
                ip2 = poin_face[i, j_neighbor]
                
                # Vectorized coordinate differences
                dx1, dy1, dz1 = x1 - mesh.x[ip1], y1 - mesh.y[ip1], z1 - mesh.z[ip1]
                dx2, dy2, dz2 = x1 - mesh.x[ip2], y1 - mesh.y[ip2], z1 - mesh.z[ip2]
                
                # Cross product for normal vector
                nx_comp = dy1 * dz2 - dz1 * dy2
                ny_comp = dz1 * dx2 - dx1 * dz2
                nz_comp = dx1 * dy2 - dy1 * dx2
                
                # Single reciprocal calculation
                norm_inv = 1.0 / sqrt(nx_comp*nx_comp + ny_comp*ny_comp + nz_comp*nz_comp)
                
                # Load derivative values
                dxdξ = dxdξ_f_face[i, j]
                dydη = dydη_f_face[i, j]
                dydξ = dydξ_f_face[i, j]
                dxdη = dxdη_f_face[i, j]
                dzdξ = dzdξ_f_face[i, j]
                dzdη = dzdη_f_face[i, j]
                
                # Corrected surface Jacobian calculation
                Je_val = sqrt((dydη * dzdξ - dydξ * dzdη)^2 + 
                             (dxdξ * dzdη - dxdη * dzdξ)^2 + 
                             (dxdη * dydξ - dxdξ * dydη)^2)
                
                Jef_face[i, j] = Je_val
                Jinv = 1.0 / Je_val
                
                # Calculate surface metric terms
                dξdx_f_face[i, j] = (dydη * dzdξ - dydξ * dzdη) * Jinv
                dξdy_f_face[i, j] = (dxdξ * dzdη - dxdη * dzdξ) * Jinv
                dξdz_f_face[i, j] = (dxdη * dydξ - dxdξ * dydη) * Jinv
                dηdx_f_face[i, j] = (dydξ * dzdη - dydη * dzdξ) * Jinv
                dηdy_f_face[i, j] = (dxdη * dzdξ - dxdξ * dzdη) * Jinv
                dηdz_f_face[i, j] = (dxdξ * dydη - dxdη * dydξ) * Jinv
                
                # Store normalized normal components
                nx_face[i, j] = nx_comp * norm_inv
                ny_face[i, j] = ny_comp * norm_inv
                nz_face[i, j] = nz_comp * norm_inv
                e = mesh.bdy_face_in_elem[iface]
                idx1 = 0
                idx2 = 0
                idx3 = 0
                #find arbitrary interior point
                ip3 = mesh.connijk[e,2,2,2]
                #=for o=1:mesh.ngl
                    for n=1:mesh.ngl
                        for m=1:mesh.ngl
                            if (mesh.connijk[e,m,n,o] == ip)
                                    idx1 = m
                                    idx2 = n
                                    idx3 = o
                            end
                        end
                    end
                end=#
                #if (idx1 + metrics.nx[iface, i, j] < 1 || idx1 + metrics.nx[iface, i, j] > N+1 || idx2 + metrics.ny[iface, i, j] < 1 || idx2 + metrics.ny[iface, i, j] > N+1 || idx3 + metrics.nz[iface, i, j] < 1 || idx3 + metrics.nz[iface, i, j] > N+1)
                if (metrics.nx[iface, i, j]*(mesh.x[ip3]-mesh.x[ip])+ metrics.ny[iface, i, j]*(mesh.y[ip3]-mesh.y[ip]) + metrics.nz[iface, i, j]*(mesh.z[ip3]-mesh.z[ip]) > 0)
                    metrics.nx[iface, i, j] = - metrics.nx[iface, i, j]
                    metrics.ny[iface, i, j] = - metrics.ny[iface, i, j]
                    metrics.nz[iface, i, j] = - metrics.nz[iface, i, j] 
                end
            end
        end
    else
        # GPU backend remains unchanged but could be optimized similarly
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
        metrics.dξdy .= (metrics.dxdζ.*metrics.dzdη .- metrics.dxdη.*metrics.dzdζ) ./ metrics.Je  # FIXED: was dzdη (typo)
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
    end
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
    # FIXED: Infer type from arrays instead of hardcoding Float32
    T = eltype(x)

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
    if (mag < T(1e-6))  # FIXED: use type parameter T
        mag = max(mag,abs(comp1),abs(comp2),abs(comp3),T(1e-7))  # FIXED: use type parameter T
    end
    Jef[iface, i, j] = mag/2
    nx[iface, i, j] = comp1/mag
    ny[iface, i, j] = comp2/mag
    nz[iface, i, j] = comp3/mag
    if (abs(nx[iface,i,j]) < T(1e-2))  # FIXED: use type parameter T
        nx[iface, i,j] = zero(T)  # FIXED: use type parameter T
    end
    if (abs(ny[iface,i,j]) < T(1e-2))  # FIXED: use type parameter T
        ny[iface, i,j] = zero(T)  # FIXED: use type parameter T
    end
    if (abs(nz[iface,i,j]) < T(1e-2))  # FIXED: use type parameter T
        nz[iface, i,j] = zero(T)  # FIXED: use type parameter T
    end
end

function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange ,N, Q, NGR, QGR, ξ, ω1, ω2, T, MT::COVAR, SD::NSD_2D; backend = CPU())
    
    ψ  = basis.ψ
    dψ = basis.dψ
    ψ1  = basisGR.ψ
    dψ1 = basisGR.dψ
    if ("Laguerre" in mesh.bdy_edge_type)
        if (backend == CPU())
            @inbounds for iel=1:mesh.nelem_semi_inf  # PERF: Added @inbounds for 4 nested loops
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

function build_metric_terms!(metrics, mesh::St_mesh, basis::St_Lagrange, basisGR::St_Lagrange,N, Q, NGR, QGR, ξ, T, MT::COVAR, SD::NSD_3D; dir="x",side ="min")
    
    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_sz = MPI.Comm_size(comm)
    
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
        # FIXED: Removed duplicate assignments (dead code)
        ψ1  = basis.ψ
        dψ1 = basis.dψ
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
        # FIXED: Removed duplicate assignments (dead code)
        ψ1  = basis.ψ
        dψ1 = basis.dψ
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
    @inbounds for iel = 1:nelem  # PERF: Added @inbounds for 6 nested loops
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
                    metrics.dξdy[l, m, n, iel] =  (metrics.dxdζ[l, m, n, iel]*metrics.dzdη[l, m, n, iel] - metrics.dxdη[l, m, n, iel]*metrics.dzdζ[l, m, n, iel])/metrics.Je[l, m, n, iel]  # FIXED: was dzdη (typo)
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


