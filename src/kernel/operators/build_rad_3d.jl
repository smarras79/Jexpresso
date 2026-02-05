using SparseArrays

function build_radiative_transfer_problem(mesh, inputs, neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz,
        nx, ny, nz, elem_to_face,
        extra_mesh, QT::Inexact, SD::NSD_3D, AD::ContGal)
    nc_mat = zeros(Float64,1)
    P = zeros(Float64,1)
    rest = zeros(Float64,1)
    nc_non_global_nodes = []
    n_non_global_nodes = 0
    n_spa = 0
    npoin = mesh.npoin
    nelem = mesh.nelem
    npoin_ang_total = 0

    if (inputs[:adaptive_extra_meshes])
        extra_meshes_coords = [Array{Float64}(undef, size(extra_mesh[e].extra_coords,1), size(extra_mesh[e].extra_coords,2)) for e in 1:nelem]
        extra_meshes_connijk = [Array{Int}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_Je = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dξdx = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dxdξ = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dξdy = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dydξ = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dηdx = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dxdη = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dηdy = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dydη = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_nops = [Array{Int}(undef, extra_mesh[e].extra_nelem) for e in 1:nelem]
        extra_meshes_extra_npoins = zeros(Int, nelem)
        extra_meshes_extra_nelems = zeros(Int, nelem)
        extra_meshes_ref_level = [Array{Int}(undef, extra_mesh[e].extra_nelem) for e in 1:nelem]
        npoin_ang_total = 0
        for e=1:nelem
            extra_meshes_coords[e] = extra_mesh[e].extra_coords[:,:]
            extra_meshes_connijk[e] = extra_mesh[e].extra_connijk
            extra_meshes_extra_Je[e] = extra_mesh[e].extra_metrics.Je[:,:,:]
            extra_meshes_extra_dξdx[e] = extra_mesh[e].extra_metrics.dξdx[:,:,:]
            extra_meshes_extra_dxdξ[e] = extra_mesh[e].extra_metrics.dxdξ[:,:,:]
            extra_meshes_extra_dξdy[e] = extra_mesh[e].extra_metrics.dξdy[:,:,:]
            extra_meshes_extra_dydξ[e] = extra_mesh[e].extra_metrics.dydξ[:,:,:]
            extra_meshes_extra_dηdx[e] = extra_mesh[e].extra_metrics.dηdx[:,:,:]
            extra_meshes_extra_dxdη[e] = extra_mesh[e].extra_metrics.dxdη[:,:,:]
            extra_meshes_extra_dηdy[e] = extra_mesh[e].extra_metrics.dηdy[:,:,:]
            extra_meshes_extra_dydη[e] = extra_mesh[e].extra_metrics.dydη[:,:,:] 
            extra_meshes_extra_npoins[e] = extra_mesh[e].extra_npoin
            extra_meshes_extra_nelems[e] = extra_mesh[e].extra_nelem
            extra_meshes_extra_nops[e] = extra_mesh[e].extra_nop
            extra_meshes_ref_level[e] = extra_mesh[e].ref_level
            npoin_ang_total += mesh.ngl*mesh.ngl*extra_mesh[e].extra_npoin
        end
        connijk_spa = [Array{Int}(undef, ngl, ngl, ngl, extra_meshes_extra_nelems[iel], extra_meshes_extra_nops[iel][1]+1, extra_meshes_extra_nops[iel][1]+1) for iel = 1:nelem]
        @info "building initial adaptive connectivity"
        neighbors = zeros(Int,nelem,26,2)
        adapted = false

        nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa  = adaptive_spatial_angular_numbering_3D_2D!(connijk_spa,nelem, ngl, mesh.connijk,
                                                                                extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                                                                                extra_meshes_coords, mesh.x, mesh.y, mesh.z, extra_meshes_ref_level, neighbors, adapted, extra_meshes_extra_Je)
        
        @info "built initial adaptive spatial angular connectivity"
        
        @info "construct global numbering pre adapting"
        ip2gip_spa, gip2ip, gip2owner_spa, gnpoin = setup_global_numbering_adaptive_angular_scalable(
            mesh.ip2gip, mesh.gip2owner, mesh, connijk_spa,
            extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops, extra_meshes_extra_nelems,
            n_spa, n_non_global_nodes, nc_non_global_nodes
            )
        @info maximum(ip2gip_spa), rank
        if rank == 0
            @info "Global spatial-angular numbering complete:"
            @info "  Total DOF: $gnpoin"
            @info "  Range: 1:$gnpoin (compact)"
            @info " next build ghost layer"
        end

        ghost_layer = build_nonconforming_ghost_layer_corrected(
            mesh, connijk_spa, mesh.ip2gip, ip2gip_spa, gip2owner_spa,
            extra_meshes_coords, extra_meshes_connijk,
            extra_meshes_extra_nops, extra_meshes_extra_nelems,
            extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dξdy,
            extra_meshes_extra_dηdx, extra_meshes_extra_dηdy,
            extra_meshes_ref_level,
            n_spa, neighbors
            )
        
        
        
        @time LHS = sparse_lhs_assembly_3Dby2D_adaptive(ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
                                                        mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                                        extra_meshes_extra_Je,
                                                        extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl, extra_meshes_extra_nelems,
                                                        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa)

        @time M = sparse_mass_assembly_3Dby2D_adaptive(ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl, extra_meshes_extra_nelems,
                                   extra_meshes_extra_npoins, connijk_spa)
        
        @info nnz(M), nnz(LHS)
        @info "built pre-adaptivity matrices"
        @info maximum(LHS), minimum(LHS)
        @info maximum(M), minimum(M)
        one_vec = Vector{Float64}(undef, size(LHS,1))
        fill!(one_vec,Float64(1))
        pointwise_interaction = abs.(LHS) * one_vec
        @info maximum(one_vec), minimum(one_vec), maximum(pointwise_interaction), minimum(pointwise_interaction)
        @time criterion = compute_adaptivity_criterion3D_2D(pointwise_interaction, nelem, ngl, mesh.connijk, 
                                                          extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems, extra_meshes_coords,
                                                          connijk_spa, extra_mesh[1].ψ, extra_mesh[1].dψ, extra_meshes_extra_dξdx, extra_meshes_extra_dηdx,
                                                          extra_meshes_extra_dξdy, extra_meshes_extra_dηdy)

        thresholds = [0.1]
        @info "criterion computed"
        @time adapt_angular_grid_3Dby2D!(criterion,thresholds, extra_meshes_ref_level,nelem,ngl,extra_meshes_extra_nelems, extra_meshes_extra_nops, neighbors, extra_meshes_extra_npoins,
                                            extra_meshes_connijk, extra_meshes_coords, extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dxdξ,
                                            extra_meshes_extra_dξdy, extra_meshes_extra_dydξ, extra_meshes_extra_dηdy, extra_meshes_extra_dydη, extra_meshes_extra_dηdx, extra_meshes_extra_dxdη,
                                            mesh.connijk,
                                            mesh.x, mesh.y, mesh.z, mesh.xmin, mesh.ymin, mesh.zmin, mesh.xmax, mesh.ymax, mesh.zmax, extra_mesh[1].ψ, extra_mesh[1].dψ)
        
        @info "angular mesh adapted"
        if !(maximum(extra_meshes_ref_level[:][:]) == 0)
            connijk_spa = [Array{Int}(undef, ngl, ngl, ngl, extra_meshes_extra_nelems[iel], extra_meshes_extra_nops[iel][1]+1, extra_meshes_extra_nops[iel][1]+1) for iel = 1:nelem]
            @time nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa  = adaptive_spatial_angular_numbering_3D_2D!(connijk_spa,nelem, ngl, mesh.connijk,
                                                                                        extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                                                                                        extra_meshes_coords, mesh.x, mesh.y, mesh.z, extra_meshes_ref_level, neighbors, adapted, extra_meshes_extra_Je)

            @info "adapted connectivity"
            adapted = true
            @info "number of hanging nodes", n_non_global_nodes
            
            @info "construct global numbering post adapted mesh"
            ip2gip_spa, gip2ip, gip2owner_spa, gnpoin = setup_global_numbering_adaptive_angular_scalable(
                mesh.ip2gip, mesh.gip2owner, mesh, connijk_spa,
                extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                n_spa, n_non_global_nodes, nc_non_global_nodes
                )
            @info maximum(ip2gip_spa), rank
            if rank == 0
                @info "Global spatial-angular numbering complete:"
                @info "  Total DOF: $gnpoin"
                @info "  Range: 1:$gnpoin (compact)"
                @info " next build ghost layer"
            end

            ghost_layer = build_nonconforming_ghost_layer_corrected(
                mesh, connijk_spa, mesh.ip2gip, ip2gip_spa, gip2owner_spa,
                extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dξdy,
                extra_meshes_extra_dηdx, extra_meshes_extra_dηdy,
                extra_meshes_ref_level,
                n_spa, neighbors
                )
            @info "building ghost layer index extension"
            gid_to_extended_local, extended_local_to_gid, n_total = build_extended_local_numbering(
                n_spa, ghost_layer,
                ip2gip_spa, rank
                )
            # build parallel nc_mat structures:
            nc_mat, P, ghost_constraint_data, all_hanging_nodes = build_restriction_matrices_local_and_ghost(
                connijk_spa, nc_non_global_nodes, n_spa, 
                ghost_layer,
                extra_meshes_coords, extra_meshes_connijk,
                extra_meshes_extra_nops, extra_meshes_extra_nelems,
                extra_meshes_extra_Je,
                mesh, ngl, nelem,
                neighbors,
                ip2gip_spa, gid_to_extended_local, extended_local_to_gid,
                rank
                )

            # Build reverse ghost constraint map for solution prolongation (reuse across timesteps)
            @info "[Rank $rank] Building reverse ghost constraint map for solution prolongation..."
            @time reverse_ghost_map = build_reverse_ghost_constraint_map(
                ghost_constraint_data, ip2gip_spa, gip2owner_spa, rank
            )

            @time LHS = sparse_lhs_assembly_3Dby2D_adaptive(ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ,
                                                        mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                        extra_meshes_extra_Je,
                                        extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl, extra_meshes_extra_nelems,
                                        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa)

            #Try this alternative assembly approach
            @info size(nc_mat), size(LHS), size(nc_mat')
            P = nc_mat'
            rest = nc_mat#spzeros(Float64,n_spa-n_non_global_nodes,n_spa)

            @time M = sparse_mass_assembly_3Dby2D_adaptive(ω, Je, mesh.connijk, extra_mesh[1].ωθ, extra_mesh[1].ωϕ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl, extra_meshes_extra_nelems,
                                   extra_meshes_extra_npoins, connijk_spa
                                   )

            Md = diag(M)
            pM = setup_assembler(SD, Md, ip2gip_spa, gip2owner_spa)
            if  pM != nothing 
                assemble_mpi!(Md,pM)
                M = Diagonal(Md)
                M = sparse(M)
        
                #assemble_mpi!(LHS,pM)
            end
        end
        


        
        M_inv = spdiagm(0 => 1 ./ diag(M))

        MLHS = sparse(M_inv * LHS)
        #do contruction for parallel nc_mat

        # Do parallel restriction Process
        
        # Compute row effects BEFORE restriction
    
        @info "[Rank $rank] Computing interface hanging row effects (before nc_mat)..."
        @time row_effects_to_send = compute_hanging_row_effects_before_restriction(
            ghost_constraint_data, MLHS, ip2gip_spa, gip2owner_spa, rank
        )
        #Apply local Restriction
        @info "[Rank $rank] Applying nc_mat from left..."
        @time A_left_restricted = nc_mat * MLHS
    
        @info "[Rank $rank] After nc_mat: $(size(A_left_restricted)), nnz=$(nnz(A_left_restricted))"
        #Apply restriction from processor interface
        @info "[Rank $rank] Exchanging row effects..."
        @time received_row_effects = exchange_hanging_effects(row_effects_to_send, rank, comm)
    
        @info "[Rank $rank] Adding received row effects..."
        @time A_with_row_effects = add_hanging_row_effects(
            A_left_restricted, received_row_effects, ip2gip_spa, n_spa, rank
            )
        @info "[Rank $rank] After row effects: $(size(A_with_row_effects)), nnz=$(nnz(A_with_row_effects))"
        # Compute prolongation to send to other processors
        @info "[Rank $rank] Computing interface hanging column effects (before P)..."
        @time col_effects_to_send = compute_hanging_col_effects_before_prolongation(
            ghost_constraint_data, A_with_row_effects, ip2gip_spa, gip2owner_spa, n_spa, rank
            )
        
        #Apply local prolongation
        @info "[Rank $rank] Applying P from right..."
        @time A_both_restricted = A_with_row_effects * P
    
        @info "[Rank $rank] After P: $(size(A_both_restricted)), nnz=$(nnz(A_both_restricted))"
        
        #Apply inter-processor prolongation
        @info "[Rank $rank] Exchanging column effects..."
        @time received_col_effects = exchange_hanging_effects(col_effects_to_send, rank, comm)
    
        @info "[Rank $rank] Adding received column effects..."
        @time A_with_col_effects = add_hanging_col_effects(
            A_both_restricted, received_col_effects, ip2gip_spa, n_spa, rank
        )
    
        @info "[Rank $rank] After column effects: $(size(A_with_col_effects)), nnz=$(nnz(A_with_col_effects))"
    
        n_free = n_spa - length(nc_non_global_nodes)
        @info "[Rank $rank] Extracting free node submatrix..."
        @time A_free = extract_free_submatrix_remove_all_hanging(
                A_with_col_effects, all_hanging_nodes, n_free, n_spa, rank)
    
        @info "[Rank $rank] Free submatrix: $(size(A_free)), nnz=$(nnz(A_free))"
        
        
        #A = sparse(rest * MLHS * P)#sparse(A_test)
        A = sparse(A_free)
        #@info "compare MLHS with A_free and A", rank, maximum(MLHS), minimum(MLHS), maximum(A_free), minimum(A_free), maximum(A), minimum(A)
        RHS = zeros(TFloat, n_spa)#npoin_ang_total)
        ref = zeros(TFloat, n_spa)
        BDY = zeros(TFloat, n_spa)
        @info size(RHS), size(A),n_spa-n_non_global_nodes
    else
        @info extra_mesh.extra_coords[1,:]
        @info extra_mesh.extra_coords[2,:]
        npoin_ang_total = npoin*extra_mesh.extra_npoin
        @info npoin_ang_total, extra_mesh.extra_npoin, npoin
        @time LHS = sparse_lhs_assembly_3Dby2D(ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ, 
                                           mesh.x, mesh.y, mesh.z, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk, 
                                        extra_mesh.extra_metrics.Je, 
                                        extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                        dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, extra_mesh.extra_npoin, inputs[:rad_HG_g])
        @info "assembled LHS"
        @time M = sparse_mass_assembly_3Dby2D(ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ, mesh.x, mesh.y, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk,
                                    extra_mesh.extra_metrics.Je,
                                    extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                   extra_mesh.extra_npoin)
        
        M2 = sparse_mass_assembly_3Dby2D_jacc(ω, Je, mesh.connijk, extra_mesh.ωθ, extra_mesh.ωϕ, ψ,
                                     extra_mesh.extra_connijk, extra_mesh.extra_metrics.Je, extra_mesh.extra_nop,
                                     npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem, extra_mesh.extra_npoin,
                                     backend = :CPU)
        @info "compare JACC vs standard" maximum(M-M2), minimum(M-M2)
        @info "assembled Mass matrix"
        @info nnz(M), nnz(LHS), npoin_ang_total^2, nnz(M)/npoin_ang_total^2, nnz(LHS)/npoin_ang_total^2
        @info maximum(LHS), minimum(LHS)
        @info maximum(M), minimum(M)
        # inexact integration makes M diagonal, build the sparse inverse to save space
        # inexact integration makes M diagonal, build the sparse inverse to save space
        ip2gip_spa, gip2owner_spa, gnpoin = setup_global_numbering_extra_dim(mesh.ip2gip, mesh.gip2owner, npoin, extra_mesh.extra_npoin, npoin_ang_total)
        Md = diag(M)
        pM = setup_assembler(SD, Md, ip2gip_spa, gip2owner_spa)
        if  pM != nothing 
            assemble_mpi!(Md,pM)
            M = Diagonal(Md)
            M = sparse(M)
        
            #assemble_mpi!(LHS,pM)
        end
        
        I_vec = Vector{Int}()
        J_vec = Vector{Int}()
        V_vec = Vector{Float64}()
        max_entries = npoin_ang_total^2
        sizehint!(I_vec, Int64(round(max_entries*0.0001)))
        sizehint!(J_vec, Int64(round(max_entries*0.0001)))
        sizehint!(V_vec, Int64(round(max_entries*0.0001)))
        for ip=1:npoin_ang_total
            val = 1/M[ip,ip]
            push!(I_vec, ip)
            push!(J_vec, ip)
            push!(V_vec, val)
        end
        M_inv = sparse(I_vec, J_vec, V_vec)
        #@time M_inv = M \ Matrix(I, size(M)) #M\Diagonal(ones(npoin_ang_total))
        #M_inv = sparse(M_inv)
        M = nothing
        LHS_sp = nothing
        LHS_el_ang = nothing
        M_rad_el = nothing
        LHS_ang_spat = nothing
        Mass_ang_spat = nothing
        E_rad_el = nothing
        P_rad_el = nothing
        S_rad_el = nothing
        M_sp = nothing
        GC.gc()
    
        A = M_inv * LHS
        @info maximum(A), minimum(A)
        M_inv = nothing
        LHS = nothing
        GC.gc()
        BDY = zeros(TFloat, npoin_ang_total)
        RHS = zeros(TFloat, npoin_ang_total)
        ref = zeros(TFloat, npoin_ang_total)
    end

    nc_rows = zeros(Int,1,1)
    rows = rowvals(A)
    vals = nonzeros(A)
    if (size(nc_mat,1) > 2)
        nc_rows = rowvals(nc_mat)
    end
    n_free = 0
    if (inputs[:adaptive_extra_meshes])
        n_free = n_spa-n_non_global_nodes
    else
        n_free = npoin_ang_total
    end
    bdy_nodes = []
    bdy_values = []
    @time for iel=1:nelem
        for i = 1:ngl
            for j = 1:ngl
                for k = 1:ngl
                    ip = mesh.connijk[iel,i,j,k]
                    x = mesh.x[ip]
                    y = mesh.y[ip]
                    z = mesh.z[ip]
                    gip = exp(-((1. / 3) * (x - (3 / 3.)))^2)#exp(-((x-3/3)/3)^2)
                    hip = exp(-4. * (2 - y) / 2)#exp(-4*(2-y)/2)
                    fip = 1.0
                    is_boundary = ip in mesh.poin_in_bdy_face
                    if (is_boundary)
                        iface = elem_to_face[iel,i,j,k,1]
                        face_i = elem_to_face[iel,i,j,k,2]
                        face_j = elem_to_face[iel,i,j,k,3]
                        nx_new = nx[iface,face_i,face_j]
                        ny_new = ny[iface,face_i,face_j] 
                        nz_new = nz[iface,face_i,face_j]
                        matchx = (x == mesh.xmax || x == mesh.xmin)
                        matchy = (y == mesh.ymax || y == mesh.ymin)
                        matchz = (z == mesh.zmax || z == mesh.zmin)
                        nmatches = 0
                        if (matchx) nmatches +=1 end
                        if (matchy) nmatches +=1 end
                        if (matchz) nmatches +=1 end
                        if (nmatches == 2) ##Do 2 boundary corner case
                            iface1 = 1
                            found = false
                            iface_found = 0
                            face_found_j = 0
                            face_found_i = 0
                            while (iface1 <= mesh.nfaces_bdy && found == false)
                                if (mesh.bdy_face_in_elem[iface1] == iel)
                                    for iter_i = 1:ngl
                                        for iter_j = 1:ngl
                                            ip1 = mesh.poin_in_bdy_face[iface1,iter_i, iter_j]
                                            if (ip1 == ip && iface != iface1)
                                                if (abs(nx[iface,face_i,face_j]-nx[iface1,iter_i,iter_j]) > 1e-12)|| (abs(ny[iface,face_i,face_j]-ny[iface1,iter_i,iter_j]) > 1e-12) || (abs(nz[iface,face_i,face_j]-nz[iface1,iter_i,iter_j]) > 1e-12)
                                                    found = true
                                                    iface_found = iface1
                                                    face_found_i = iter_i
                                                    face_found_j = iter_j
                                                end
                                            end
                                        end
                                    end
                                end
                                iface1 +=1
                            end
                            nx_new = nx[iface,face_i,face_j] + nx[iface_found,face_found_i,face_found_j]
                            ny_new = ny[iface,face_i,face_j] + ny[iface_found,face_found_i,face_found_j]
                            nz_new = nz[iface,face_i,face_j] + nz[iface_found,face_found_i,face_found_j]
                                            
                        elseif (nmatches == 3)
                            iface1 = 1
                            found = false
                            iface_found = 0
                            face_found_j = 0
                            face_found_i = 0
                            while (iface1 <= mesh.nfaces_bdy && found == false)
                                if (mesh.bdy_face_in_elem[iface1] == iel)
                                    for iter_i = 1:ngl
                                        for iter_j = 1:ngl
                                            ip1 = mesh.poin_in_bdy_face[iface1,iter_i, iter_j]
                                            if (ip1 == ip && iface != iface1)
                                                found = true
                                                iface_found = iface1
                                                face_found_i = iter_i
                                                face_found_j = iter_j
                                            end
                                        end
                                    end
                                end
                                iface1 +=1
                            end
                            iface1 = 1
                            found = false
                            iface_found_2 = 0
                            face_found_j_2 = 0
                            face_found_i_2 = 0
                            while (iface1 <= mesh.nfaces_bdy && found == false)
                                if (mesh.bdy_face_in_elem[iface1] == iel)
                                    for iter_i = 1:ngl
                                        for iter_j = 1:ngl
                                            ip1 = mesh.poin_in_bdy_face[iface1,iter_i, iter_j]
                                            if (ip1 == ip && iface != iface1 && iface_found != iface1)
                                                found = true
                                                iface_found_2 = iface1
                                                face_found_i_2 = iter_i
                                                face_found_j_2 = iter_j
                                            end
                                        end
                                    end
                                end
                                iface1 +=1
                            end
                            nx_new = nx[iface,face_i,face_j] + nx[iface_found,face_found_i,face_found_j] + nx[iface_found_2,face_found_i_2,face_found_j_2]
                            ny_new = ny[iface,face_i,face_j] + ny[iface_found,face_found_i,face_found_j] + ny[iface_found_2,face_found_i_2,face_found_j_2]
                            nz_new = nz[iface,face_i,face_j] + nz[iface_found,face_found_i,face_found_j] + nz[iface_found_2,face_found_i_2,face_found_j_2]
                        end
                    end
                    if (inputs[:adaptive_extra_meshes])
                        for e_ext = 1:extra_meshes_extra_nelems[iel]
                            for jθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                                for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                                    ip_ext = extra_meshes_connijk[iel][e_ext,iθ,jθ]
                                    θ = extra_meshes_coords[iel][1,ip_ext]
                                    ϕ = extra_meshes_coords[iel][2,ip_ext]
                                    ip_g = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]

                                    sip = exp(-((6/ (2. * π)) * (θ - (3. * π / 5.)))^2)#exp(-((96/(2*π))*(θ-7*π/5))^2)
                                    bip = exp(-((6/ (2. * π)) * (ϕ - (2. * π / 3.)))^2)
                                    uip = gip*hip*fip*sip*bip
                                    ref[ip_g] = uip
                
                                    if (is_boundary)
                                        applied = false
                                        prodx = nx_new*sin(θ)*cos(ϕ)
                                        prody = ny_new*sin(θ)*sin(ϕ)
                                        prodz = nz_new*cos(θ)
                                       
                                        if (prodx + prody + prodz < 0)
                                            
                                            if (ip_g <= n_free)
                                                if (gip2owner_spa[ip_g] == rank)
                                                    BDY[ip_g] = user_rad_bc(x,y,z,θ,ϕ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                                    push!(bdy_nodes, ip_g)
                                                    push!(bdy_values, BDY[ip_g])
                                                else
                                                    BDY[ip_g] = user_rad_bc(x,y,z,θ,ϕ)
                                                    push!(bdy_nodes, ip_g)
                                                    push!(bdy_values, 0.0)
                                                end
                                                applied = true
                                            
                                            end
                                           

                                        end
                                        if (applied == false)
                                            if (gip2owner_spa[ip_g] == rank)
                                                RHS[ip_g] = user_rhs(x,y,z,θ,ϕ)
                                            end
                                        end
                                    else
                                        if (gip2owner_spa[ip_g] == rank)
                                            RHS[ip_g] = user_rhs(x,y,z,θ,ϕ)
                                        end
                                    end
                                end
                            end
                        end
                    else

                        for e_ext = 1:extra_mesh.extra_nelem
                            for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                                for iϕ = 1:extra_mesh.extra_nop[e_ext]+1
                                    ip_ext = extra_mesh.extra_connijk[e_ext,iϕ,iθ]
                                    θ = extra_mesh.extra_coords[1,ip_ext]
                                    ϕ = extra_mesh.extra_coords[2,ip_ext]
                                    ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                                    sip = exp(-((6/ (2. * π)) * (θ - (3. * π / 5.)))^2)#exp(-((96/(2*π))*(θ-7*π/5))^2)
                                    bip = exp(-((6/ (2. * π)) * (ϕ - (2. * π / 3.)))^2)
                                    uip = gip*hip*fip*sip*bip
                                    ref[ip_g] = uip

                                    if (is_boundary)
                                        
                                        applied = false
                                        prodx = nx_new*sin(θ)*cos(ϕ)
                                        prody = ny_new*sin(θ)*sin(ϕ)
                                        prodz = nz_new*cos(θ)
                                        #@info nx[iface,face_i,face_j], ny[iface,face_i,face_j], nz[iface,face_i,face_j], x, y, z
                                        if (prodx + prody + prodz < 0)
                                            if (gip2owner_spa[ip_g] == rank)
                                                BDY[ip_g] = user_rad_bc(x,y,z,θ,ϕ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                                push!(bdy_nodes, ip_g)
                                                push!(bdy_values, BDY[ip_g])
                                            else
                                                BDY[ip_g] = user_rad_bc(x,y,z,θ,ϕ)
                                                push!(bdy_nodes, ip_g)
                                                push!(bdy_values, 0.0)
                                            end
                                            applied = true
                                            
                                        end
                                        if (applied == false)
                                            if (gip2owner_spa[ip_g] == rank)
                                                RHS[ip_g] = user_rhs(x,y,z,θ,ϕ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip)
                                            end
                                        end
                                    else
                                        if (gip2owner_spa[ip_g] == rank)
                                            RHS[ip_g] = user_rhs(x,y,z,θ,ϕ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip) 
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
    B = RHS

    I_orig, J_orig, V_orig = findnz(A)
    
    # Build modifications in COO format
    I_mod = Int[]
    J_mod = Int[]
    V_mod = Float64[]
    boundary_set = Set(bdy_nodes)
    boundary_dict = Dict(zip(bdy_nodes, bdy_values))
    for (i, j, v) in zip(I_orig, J_orig, V_orig)
        if i in boundary_set
            # This row should become identity
            if i == j
                # Keep diagonal, but set to 1

                if (gip2owner_spa[i] == rank)
                    push!(I_mod, i)
                    push!(J_mod, j)
                    push!(V_mod, 1.0 - v)# Correction to make it 1
                else
                    push!(I_mod, i)
                    push!(J_mod, j)
                    push!(V_mod, -v) 
                    
                end
            else
                # Zero out off-diagonal
                
                push!(I_mod, i)
                push!(J_mod, j)
                push!(V_mod, -v)  # Correction to make it 0
            end
        elseif j in boundary_set
            # For symmetric elimination: zero out column and transfer to RHS
            # This is the key addition for symmetry
            if i != j  # Don't touch diagonal elements of boundary nodes
               # Modify RHS to account for known boundary value
                
                RHS[i] -= v * BDY[j]
                
                # Zero out this column entry
                push!(I_mod, i)
                push!(J_mod, j)
                push!(V_mod, -v)
            end
        end
    end

    
    @info size(I_mod), size(J_mod), size(V_mod)
    @info maximum(J_mod), minimum(J_mod)
    
    #=existing_diag = Set{Int}()
    for (i, j) in zip(I_orig, J_orig)
        if i == j && i in boundary_set
            push!(existing_diag, i)
        end
    end
    
    # Add missing diagonals
    for boundary_node in bdy_nodes
        if !(boundary_node in existing_diag)
            push!(I_mod, boundary_node)
            push!(J_mod, boundary_node)
            push!(V_mod, 1.0)
        end
    end=#

    @info size(I_mod), size(J_mod), size(V_mod)
    @info maximum(J_mod), minimum(J_mod)
    
    # Build correction matrix and add
    C = sparse(I_mod, J_mod, V_mod, n_free, n_free)
    A = sparse(A + C)
    
    if (inputs[:adaptive_extra_meshes])

        #RHS_red = nc_mat * RHS
        # Parallel RHS restriction
        @info "[Rank $rank] Computing interface hanging RHS effects..."
        @time rhs_effects_to_send = compute_hanging_rhs_effects_before_restriction(
            ghost_constraint_data, RHS, ip2gip_spa, gip2owner_spa, rank
        )

        @info "[Rank $rank] Applying nc_mat to RHS..."
        @time RHS_restricted = nc_mat * RHS

        @info "[Rank $rank] Exchanging RHS effects..."
        @time received_rhs_effects = exchange_hanging_effects_vector(rhs_effects_to_send, rank, comm)

        @info "[Rank $rank] Adding received RHS effects..."
        @time RHS_with_effects = add_hanging_rhs_effects(
            RHS_restricted, received_rhs_effects, ip2gip_spa, n_spa, rank
        )

        @info "[Rank $rank] Extracting free RHS subvector..."
        @time RHS_red = extract_free_rhs_subvector(
            RHS_with_effects, all_hanging_nodes, n_free, n_spa, rank
        )

        @time for boundary_node in bdy_nodes
            RHS_red[boundary_node] = boundary_dict[boundary_node]
        end
        #@info maximum(abs.(RHS_red-U_red_proj))
        B = RHS_red
    else
        @time for boundary_node in bdy_nodes
            RHS[boundary_node] = boundary_dict[boundary_node]
        end
        B = RHS
    end
    

    #A_inv = inv(A)
    @info "built RHS"
    #@info RHS
    @info "solving system"
    As = sparse(A)
    #A = nothing
    GC.gc()
    #@time solution = As \ B
    @info "sizes of A and B", size(As), size(B)
    npoin_ang_total = size(B,1)
    @time solution = solve_parallel_lsqr(ip2gip_spa, gip2owner_spa, As, B, gnpoin, npoin_ang_total, pM)
    #=@time solution, stats = Krylov.cgs(As, B;
                   atol = 1e-7,
                   rtol = 1e-7,
                   #btol = 1e-13,
                   #etol = 1e-13,
                   #axtol = 1e-13,
                   itmax = n_spa,
                   verbose = 1)=#
    @info "done radiation solved"
    @info maximum(solution), minimum(solution)
    @info "dof", npoin_ang_total
    A = nothing
    RHS = nothing
    GC.gc()
    @info "integrating solution and reference in angle"
    int_sol = zeros(TFloat, npoin,1)
    int_ref = zeros(TFloat, npoin,1)
    L2_err = 0.0
    L2_ref = 0.0
    solution_new = zeros(Float64,n_spa)
    if inputs[:adaptive_extra_meshes]
        # solution_new = nc_mat' * solution
        # Parallel solution prolongation

        # Step 1: Compute contributions from free nodes that are ghost parents
        @info "[Rank $rank] Computing solution prolongation contributions..."
        @time solution_contributions_to_send = compute_solution_prolongation_contributions(
            reverse_ghost_map, solution, ip2gip_spa, n_free, rank
        )

        # Step 2: Apply local prolongation
        @info "[Rank $rank] Applying local prolongation P * solution..."
        # Extend solution to n_spa with zeros for hanging nodes
        solution_extended = zeros(eltype(solution), n_spa)
        solution_extended[1:n_free] .= solution

        @time solution_local = nc_mat' * solution_extended

        # Step 3: Exchange and add contributions from ghost parents
        @info "[Rank $rank] Exchanging solution prolongation contributions..."
        @time received_solution_contributions = exchange_hanging_effects_vector(
            solution_contributions_to_send, rank, comm
        )

        @info "[Rank $rank] Adding received solution contributions..."
        @time solution_new = add_solution_prolongation_contributions(
            solution_local, received_solution_contributions, ip2gip_spa, n_spa, rank
        )

        @info "[Rank $rank] Solution prolongation complete"
    end
    if inputs[:adaptive_extra_meshes]
        for iel=1:nelem
            for i=1:ngl
                for j =1:ngl
                    for k=1:ngl
                        ip = mesh.connijk[iel,i,j,k]
                        x = mesh.x[ip]
                        y = mesh.y[ip]
                        z = mesh.z[ip]
                        div = 1
                        corner_match = 0
                        matchx = (abs(x-mesh.xmin)<1e-10 || abs(x- mesh.xmax)< 1e-10)
                        matchy = (abs(y-mesh.ymin)<1e-10 || abs(y-mesh.ymax) < 1e-10)
                        matchz = (abs(z-mesh.zmin)<1e-10 || abs(z-mesh.zmax) < 1e-10)
                        nmatches = 0
                        if (matchx) nmatches +=1 end
                        if (matchy) nmatches +=1 end
                        if (matchz) nmatches +=1 end
                        if (i==1 && j == 1 && k == 1) || (i==1 && j == 1 && k == ngl) || (i==1 && j == ngl && k == 1) || (i==ngl && j == 1 && k == 1) || (i==1 && j == ngl && k == ngl)
                            corner_match = 1
                        end
                        if  (i==ngl && j == 1 && k == ngl) || (i==ngl && j == ngl && k == 1) || (i==ngl && j == ngl && k == ngl)
                            corner_match = 1
                        end
                        if (corner_match == 1)
                            if !(ip in  mesh.poin_in_bdy_face)
                                div = 8
                            elseif (nmatches == 1)
                                div = 4
                            elseif (nmatches == 2)
                                div = 2
                            end
                        elseif (i==1 && j == 1) || (i==ngl && j==ngl) || (i==1 && j == ngl) || (i==ngl && j == 1)||(i==1 && k == 1) || (i == 1 && k== ngl) || (i==ngl && k==1) || (i==ngl && k== ngl) || (j==1 && k==1) || (j==1 && k==ngl) || (j==ngl && k==1) || (j==ngl && k==ngl)


                            if !(ip in  mesh.poin_in_bdy_face)
                                div = 4
                            else
                                if ( nmatches == 1)
                                    div = 2
                                end
                            end

                        elseif (i==1 || j==1 || i==ngl || j==ngl) || (k==1 || k==ngl)
                            if !(ip in  mesh.poin_in_bdy_face)
                                div = 2
                            end
                        end
                        for e_ext = 1:extra_meshes_extra_nelems[iel]
                            for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                                for iϕ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                                    if (iθ == 1 || iθ == extra_meshes_extra_nops[iel][e_ext]+1)
                                        div1 = div *2
                                    else
                                        div1 = div
                                    end
                                    ip_ext = extra_meshes_connijk[iel][e_ext,iθ,iϕ]
                                    θ = extra_meshes_coords[iel][1,ip_ext]
                                    ϕ = extra_meshes_coords[iel][2,ip_ext]
                                    ip_g = connijk_spa[iel][i,j,k,e_ext,iθ,iϕ]#(ip-1) * extra_mesh.extra_npoin + ip_ext
                                    int_sol[ip] += solution_new[ip_g]*extra_meshes_extra_Je[iel][e_ext,iθ,iϕ]*extra_mesh[iel].ωθ[iθ]*extra_mesh[iel].ωθ[iϕ]/div
                                    int_ref[ip] += (ref[ip_g])*extra_meshes_extra_Je[iel][e_ext,iθ,iϕ]*extra_mesh[iel].ωθ[iθ]*extra_mesh[iel].ωθ[iϕ]/div
                                    #if (abs(y-1.0)<1e-5 && abs(x-1.0) < 1e-5) @info int_sol[ip], int_ref[ip], x, y, z, θ, ϕ, solution_new[ip_g], ref[ip_g], e_ext, ip_g, ip, ip_ext end
                                    L2_ref += (ref[ip_g])^2*extra_meshes_extra_Je[iel][e_ext,iθ,iϕ]*extra_mesh[iel].ωθ[iθ]*extra_mesh[iel].ωθ[iϕ]*ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]#/div1
                                    L2_err += (ref[ip_g]-solution_new[ip_g])^2*extra_meshes_extra_Je[iel][e_ext,iθ,iϕ]*extra_mesh[iel].ωθ[iθ]*extra_mesh[iel].ωθ[iϕ]*ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]#/div1
                                end
                            end
                        end
                    end
                end
            end
        end
  
    else
        
        for iel=1:nelem
            for i=1:ngl
                for j =1:ngl
                    for k=1:ngl
                        ip = mesh.connijk[iel,i,j,k]
                        x = mesh.x[ip]
                        y = mesh.y[ip]
                        z = mesh.z[ip]
                        div = 1
                        corner_match = 0
                        matchx = (abs(x-mesh.xmin)<1e-10 || abs(x- mesh.xmax)< 1e-10)
                        matchy = (abs(y-mesh.ymin)<1e-10 || abs(y-mesh.ymax) < 1e-10)
                        matchz = (abs(z-mesh.zmin)<1e-10 || abs(z-mesh.zmax) < 1e-10)
                        nmatches = 0
                        if (matchx) nmatches +=1 end
                        if (matchy) nmatches +=1 end
                        if (matchz) nmatches +=1 end
                        if (i==1 && j == 1 && k == 1) || (i==1 && j == 1 && k == ngl) || (i==1 && j == ngl && k == 1) || (i==ngl && j == 1 && k == 1) || (i==1 && j == ngl && k == ngl)
                            corner_match = 1
                        end
                        if  (i==ngl && j == 1 && k == ngl) || (i==ngl && j == ngl && k == 1) || (i==ngl && j == ngl && k == ngl)
                            corner_match = 1
                        end
                        if (corner_match == 1)
                            if !(ip in  mesh.poin_in_bdy_face)
                                div = 8
                            elseif (nmatches == 1)
                                div = 4
                            elseif (nmatches == 2)
                                div = 2
                            end
                        elseif (i==1 && j == 1) || (i==ngl && j==ngl) || (i==1 && j == ngl) || (i==ngl && j == 1)||(i==1 && k == 1) || (i == 1 && k== ngl) || (i==ngl && k==1) || (i==ngl && k== ngl) || (j==1 && k==1) || (j==1 && k==ngl) || (j==ngl && k==1) || (j==ngl && k==ngl)

                    
                            if !(ip in  mesh.poin_in_bdy_face)
                                div = 4
                            else
                                if ( nmatches == 1)
                                    div = 2
                                end
                            end
                
                        elseif (i==1 || j==1 || i==ngl || j==ngl) || (k==1 || k==ngl)
                            if !(ip in  mesh.poin_in_bdy_face)
                                div = 2
                            end 
                        end
                        for e_ext = 1:extra_mesh.extra_nelem
                            for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                                for iϕ = 1:extra_mesh.extra_nop[e_ext]+1
                                    if (iθ == 1 || iθ == extra_mesh.extra_nop[e_ext]+1)
                                        div1 = div *2
                                    else
                                        div1 = div
                                    end
                                    ip_ext = extra_mesh.extra_connijk[e_ext,iθ,iϕ]
                                    θ = extra_mesh.extra_coords[1,ip_ext]
                                    ϕ = extra_mesh.extra_coords[2,ip_ext]
                                    ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                                    int_sol[ip] += solution[ip_g]*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]*extra_mesh.ωθ[iθ]*extra_mesh.ωθ[iϕ]/div
                                    int_ref[ip] += (ref[ip_g]-solution[ip_g])*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]*extra_mesh.ωθ[iθ]*extra_mesh.ωθ[iϕ]/div
                                    L2_ref += (ref[ip_g])^2*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]*extra_mesh.ωθ[iθ]*extra_mesh.ωθ[iϕ]*ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]#/div1
                                    L2_err += (ref[ip_g]-solution[ip_g])^2*extra_mesh.extra_metrics.Je[e_ext,iθ,iϕ]*extra_mesh.ωθ[iθ]*extra_mesh.ωθ[iϕ]*ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]#/div1
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    @info "new L2 norms", sqrt(L2_ref), sqrt(L2_err), sqrt(L2_err/L2_ref)
    title = @sprintf "Solution-Radiation"
    @info maximum(int_ref), minimum(int_ref)
    write_vtk(SD, mesh, int_sol, int_sol, nothing,
                  nothing, nothing,
                  0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
                  ["Ang_int"], ["Ang_int"];
                  iout=1,nvar = 1)
end

function sparse_lhs_assembly_3Dby2D(ω, Je, connijk, ωθ, ωϕ, x, y, z,
                                    ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
                                   dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, npoin_ang, rad_HG_g)

        
    entries_per_row = nelem_ang * ngl/2  # Spatial stencil × angular nodes
    max_entries = npoin_ang_total * entries_per_row
    @info max_entries
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries)))
    sizehint!(J_vec, Int64(round(max_entries)))
    sizehint!(V_vec, Int64(round(max_entries)))
    HG, error = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)), 0, 2*π, rtol=1e-13, atol = 1e-13)    
    for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ip = connijk[iel,i,j,k]
                    ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                    dξdx_ij = dξdx[iel,i,j,k]
                    dξdy_ij = dξdy[iel,i,j,k]
                    dξdz_ij = dξdz[iel,i,j,k]
                    dηdx_ij = dηdx[iel,i,j,k]
                    dηdy_ij = dηdy[iel,i,j,k]
                    dηdz_ij = dηdz[iel,i,j,k]
                    dζdx_ij = dζdx[iel,i,j,k]
                    dζdy_ij = dζdy[iel,i,j,k]
                    dζdz_ij = dζdz[iel,i,j,k]
                    κ = user_extinction(x[ip],y[ip],z[ip])
                    σ = user_scattering_coef(x[ip],y[ip],z[ip])
                    for e_ext = 1:nelem_ang
                        for jθ = 1:nop_ang[e_ext]+1
                            for iθ = 1:nop_ang[e_ext]+1
                                ip_ext = connijk_ang[e_ext,iθ,jθ]
                                
                                ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext,iθ,jθ]
                                ωJac_full = ωJac * ωJac_rad
                                #@info coords_ang[1,ip_ext], coords_ang[2,ip_ext], e_ext, iθ, jθ
                                
                                θ = coords_ang[1,ip_ext]
                                ϕ = coords_ang[2,ip_ext]
                                Ωx = sin(θ) * cos(ϕ)
                                Ωy = sin(θ) * sin(ϕ)
                                Ωz = cos(θ)
                                intΦ = 0.0
                                for e_ext_scatter = 1:nelem_ang
                                    for nθ = 1:nop_ang[e_ext]+1
                                        for mθ = 1:nop_ang[e_ext]+1
                                            div = 1
                                            ipθ = connijk_ang[e_ext_scatter,mθ,nθ]
                                            θ1 = coords_ang[1,ipθ]
                                            ϕ1 = coords_ang[2,ipθ]

                                            Φ = user_scattering_functions(θ,θ1,ϕ,ϕ1,HG)
                                            ωJac_rad_scatter = ωθ[mθ]*ωϕ[nθ]*Je_ang[e_ext_scatter,mθ,nθ]
                                            intΦ +=   ωJac_rad_scatter*Φ/div
                                        end
                                    end
                                end
                                
                                idx_ip = (ip-1)*(npoin_ang) + ip_ext
                                val_diagonal = κ * ωJac_full - σ * intΦ * ωJac_full
                                push!(I_vec, idx_ip)
                                push!(J_vec, idx_ip)
                                push!(V_vec, val_diagonal)
                                prop_coeff_i = (dξdx_ij * Ωx + dξdy_ij * Ωy + dξdz_ij * Ωz) * ωJac_full
                                for m = 1:ngl
                                    # Skip diagonal (already added)
                        
                                    dψ_mi = dψ[m, i]
                                    if abs(dψ_mi) < eps(Float64) continue end
                        
                                    jp = connijk[iel, m, j, k]
                                    idx_m = (jp-1)*(npoin_ang) + ip_ext
                        
                                    val = dψ_mi * prop_coeff_i
                        
                                    push!(I_vec, idx_ip)
                                    push!(J_vec, idx_m)
                                    push!(V_vec, val)
                                end 

                                prop_coeff_j = (dηdx_ij * Ωx + dηdy_ij * Ωy + dηdz_ij * Ωz) * ωJac_full

                                for n = 1:ngl
                                    # Skip diagonal
                        
                                    dψ_nj = dψ[n, j]
                                    if abs(dψ_nj) < eps(Float64) continue end
                        
                                    jp = connijk[iel, i, n, k]
                                    idx_n = (jp-1)*(npoin_ang) + ip_ext
                        
                                    val = dψ_nj * prop_coeff_j
                        
                                    push!(I_vec, idx_ip)
                                    push!(J_vec, idx_n)
                                    push!(V_vec, val)
                                end
                                prop_coeff_k = (dζdx_ij * Ωx + dζdy_ij * Ωy + dζdz_ij * Ωz) * ωJac_full
                                for o = 1:ngl
                                    # Skip diagonal
                        
                                    dψ_ok = dψ[o, k]
                                    if abs(dψ_ok) < eps(Float64) continue end
                        
                                    jp = connijk[iel, i, j, o]
                                    idx_o = (jp-1)*(npoin_ang) + ip_ext
                        
                                    val = dψ_ok * prop_coeff_k
                        
                                    push!(I_vec, idx_ip)
                                    push!(J_vec, idx_o)
                                    push!(V_vec, val)
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_lhs_assembly_3Dby2D_adaptive(ω, Je, connijk, ωθ, ωϕ, x, y, z,
                                    ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
                                   dξdx, dξdy, dξdz, dηdx, dηdy, dηdz, dζdx, dζdy, dζdz, npoin_ang, rad_HG_g, connijk_spa)
    nelem_ang_avg = 0
    for iel=1:nelem
        nelem_ang_avg += nelem_ang[iel]/nelem
    end
    entries_per_row = nelem_ang_avg * ngl/2  # Spatial stencil × angular nodes
    max_entries = npoin_ang_total * entries_per_row
    @info max_entries
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    sizehint!(I_vec, Int64(round(max_entries)))
    sizehint!(J_vec, Int64(round(max_entries)))
    sizehint!(V_vec, Int64(round(max_entries)))
    HG, error = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)), 0, 2*π, rtol=1e-13, atol = 1e-13)
    for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    eq1 = i == j
                    eq2 = i == k
                    eq3 = j == k
                    ip = connijk[iel,i,j,k]
                    ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                    dξdx_ij = dξdx[iel,i,j,k]
                    dξdy_ij = dξdy[iel,i,j,k]
                    dξdz_ij = dξdz[iel,i,j,k]
                    dηdx_ij = dηdx[iel,i,j,k]
                    dηdy_ij = dηdy[iel,i,j,k]
                    dηdz_ij = dηdz[iel,i,j,k]
                    dζdx_ij = dζdx[iel,i,j,k]
                    dζdy_ij = dζdy[iel,i,j,k]
                    dζdz_ij = dζdz[iel,i,j,k]
                    κ = user_extinction(x[ip],y[ip],z[ip])
                    σ = user_scattering_coef(x[ip],y[ip],z[ip])
                    for e_ext = 1:nelem_ang[iel]
                        for jθ = 1:nop_ang[iel][e_ext]+1
                            for iθ = 1:nop_ang[iel][e_ext]+1
                                ip_ext = connijk_ang[iel][e_ext,iθ,jθ]
                                sum = 0.0
                                ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[iel][e_ext,iθ,jθ]
                                ωJac_full = ωJac * ωJac_rad
                                #@info coords_ang[1,ip_ext], coords_ang[2,ip_ext], e_ext, iθ, jθ
                                extinction = κ*ωJac*ωJac_rad
                                θ = coords_ang[iel][1,ip_ext]
                                ϕ = coords_ang[iel][2,ip_ext]
                                Ωx = sin(θ) * cos(ϕ)
                                Ωy = sin(θ) * sin(ϕ)
                                Ωz = cos(θ)
                                intΦ = 0.0
                                for e_ext_scatter = 1:nelem_ang[iel]
                                    for nθ = 1:nop_ang[iel][e_ext]+1
                                        for mθ = 1:nop_ang[iel][e_ext]+1
                                            div = 1
                                            ipθ = connijk_ang[iel][e_ext_scatter,mθ,nθ]
                                            θ1 = coords_ang[iel][1,ipθ]
                                            ϕ1 = coords_ang[iel][2,ipθ]

                                            Φ = user_scattering_functions(θ,θ1,ϕ,ϕ1,HG)
                                            ωJac_rad_scatter = ωθ[mθ]*ωϕ[nθ]*Je_ang[iel][e_ext_scatter,mθ,nθ]
                                            intΦ +=   ωJac_rad_scatter*Φ/div
                                        end
                                    end
                                end
                                scattering = intΦ *ωJac*ωJac_rad*σ
                                idx_ip = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                                val_diagonal = κ * ωJac_full - σ * intΦ * ωJac_full
                                push!(I_vec, idx_ip)
                                push!(J_vec, idx_ip)
                                push!(V_vec, val_diagonal)
                                prop_coeff_i = (dξdx_ij * Ωx + dξdy_ij * Ωy + dξdz_ij * Ωz) * ωJac_full
                                for m = 1:ngl
                                    # Skip diagonal (already added)
                        
                                    dψ_mi = dψ[m, i]
                                    if abs(dψ_mi) < eps(Float64) continue end
                        
                                    jp = connijk[iel, m, j, k]
                                    idx_m = connijk_spa[iel][m, j, k, e_ext, iθ, jθ]
                        
                                    val = dψ_mi * prop_coeff_i
                        
                                    push!(I_vec, idx_ip)
                                    push!(J_vec, idx_m)
                                    push!(V_vec, val)
                                end 

                                prop_coeff_j = (dηdx_ij * Ωx + dηdy_ij * Ωy + dηdz_ij * Ωz) * ωJac_full

                                for n = 1:ngl
                                    # Skip diagonal
                        
                                    dψ_nj = dψ[n, j]
                                    if abs(dψ_nj) < eps(Float64) continue end
                        
                                    jp = connijk[iel, i, n, k]
                                    idx_n = connijk_spa[iel][i, n, k, e_ext, iθ, jθ]
                        
                                    val = dψ_nj * prop_coeff_j
                        
                                    push!(I_vec, idx_ip)
                                    push!(J_vec, idx_n)
                                    push!(V_vec, val)
                                end
                                prop_coeff_k = (dζdx_ij * Ωx + dζdy_ij * Ωy + dζdz_ij * Ωz) * ωJac_full
                                for o = 1:ngl
                                    # Skip diagonal
                        
                                    dψ_ok = dψ[o, k]
                                    if abs(dψ_ok) < eps(Float64) continue end
                        
                                    jp = connijk[iel, i, j, o]
                                    idx_o = connijk_spa[iel][i, j, o, e_ext, iθ, jθ]
                        
                                    val = dψ_ok * prop_coeff_k
                        
                                    push!(I_vec, idx_ip)
                                    push!(J_vec, idx_o)
                                    push!(V_vec, val)
                                end
                                
                            end
                        end
                    end
                end
            end
        end
    end
    return sparse(I_vec, J_vec, V_vec)
end


function sparse_mass_assembly_3Dby2D(ω, Je, connijk, ωθ, ωϕ, x, y, ψ, dψ, ψ_ang, 
        connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang)
  
    max_entries = npoin_ang_total
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    
    sizehint!(I_vec, Int64(round(max_entries)))
    sizehint!(J_vec, Int64(round(max_entries)))
    sizehint!(V_vec, Int64(round(max_entries)))

    for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ip = connijk[iel,i,j,k]
                    ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                    for e_ext = 1:nelem_ang
                        for jθ = 1:nop_ang[e_ext]+1
                            for iθ = 1:nop_ang[e_ext]+1
                                ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[e_ext,iθ,jθ]
                                ip_ext = connijk_ang[e_ext,iθ,jθ]
                                for o=1:ngl
                                    for n=1:ngl
                                        for m=1:ngl
                                            jp = connijk[iel,m,n,o]
                                            val = ωJac*ωJac_rad*ψ[k,o]*ψ[j,n]*ψ[i,m]#*ψ_ang[iθ,kθ]*ψ_ang[jθ,lθ]
                                            idx_ip = (ip-1)*(npoin_ang) + ip_ext
                                            idx_jp = (jp-1)*(npoin_ang) + ip_ext
                                            if abs(val) > eps(Float64)  # Skip near-zero entries
                                                push!(I_vec, idx_ip)
                                                push!(J_vec, idx_jp)
                                                push!(V_vec, val)
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
    end
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_mass_assembly_3Dby2D_adaptive(ω, Je, connijk, ωθ, ωϕ, x, y, ψ, dψ, ψ_ang,
        connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang, connijk_spa)

    max_entries = npoin_ang_total
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    
    sizehint!(I_vec, Int64(round(max_entries)))
    sizehint!(J_vec, Int64(round(max_entries)))
    sizehint!(V_vec, Int64(round(max_entries)))

    @inbounds for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ip = connijk[iel,i,j,k]
                    ωJac = ω[i]*ω[j]*ω[k]*Je[iel,i,j,k]
                    for e_ext = 1:nelem_ang[iel]
                        for jθ = 1:nop_ang[iel][e_ext]+1
                            for iθ = 1:nop_ang[iel][e_ext]+1
                                ωJac_rad = ωθ[iθ]*ωϕ[jθ]*Je_ang[iel][e_ext,iθ,jθ]
                                ip_ext = connijk_ang[iel][e_ext,iθ,jθ]
                                for o=1:ngl
                                    for n=1:ngl
                                        for m=1:ngl
                                            jp = connijk[iel,m,n,o]
                                            val = ωJac*ωJac_rad*ψ[k,o]*ψ[j,n]*ψ[i,m]#*ψ_ang[iθ,kθ]*ψ_ang[jθ,lθ]
                                            idx_ip = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]#(ip-1)*(npoin_ang) + ip_ext
                                            idx_jp = connijk_spa[iel][m,n,o,e_ext,iθ,jθ]#(jp-1)*(npoin_ang) + ip_ext
                                            if (idx_ip == 14122 && idx_jp == 14122) @info iel, i,j,k,e_ext, iθ, jθ, m,n,o, val, ωJac, ωJac_rad, ψ[k,o], ψ[j,n], ψ[i,m] end
                                            if abs(val) > eps(Float64)  # Skip near-zero entries
                                                push!(I_vec, idx_ip)
                                                push!(J_vec, idx_jp)
                                                push!(V_vec, val)
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
    end
    return sparse(I_vec, J_vec, V_vec)
end

function compute_adaptivity_criterion3D_2D(pointwise_interaction, nelem, ngl, connijk, connijk_ang, nop_ang, nelem_ang,
        coords_ang, connijk_spa, ψ_ang, dψ_ang,dξdθ, dηdθ, dξdϕ, dηdϕ)
    criterion = [Vector{Float64}(undef, nelem_ang[e]) for e=1:nelem]
    for iel=1:nelem
        for k=1:ngl
            for j=1:ngl
                for i=1:ngl
                    ip = connijk[iel,i,j,k]
                    for e_ext = 1:nelem_ang[iel]
                        criterion[iel][e_ext] = 0.0
                        for jθ = 1:nop_ang[iel][e_ext]+1
                            for iθ = 1:nop_ang[iel][e_ext]+1
                    
                                ip_ext = connijk_ang[iel][e_ext,iθ,jθ]
                                idx_ip = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                                gradξ = 0.0
                                gradη = 0.0
                                for kθ = 1:nop_ang[iel][e_ext]+1
                                    gradξ += pointwise_interaction[idx_ip]*dψ_ang[iθ,kθ]*ψ_ang[jθ,kθ]
                                    gradη += pointwise_interaction[idx_ip]*ψ_ang[iθ,kθ]*dψ_ang[jθ,kθ]
                                end
                                gradθ = gradξ * dξdθ[iel][e_ext,iθ,jθ] + gradη * dηdθ[iel][e_ext,iθ,jθ]
                                gradϕ = gradξ * dξdϕ[iel][e_ext,iθ,jθ] + gradη * dηdϕ[iel][e_ext,iθ,jθ]
                                
                                criterion[iel][e_ext] = max(criterion[iel][e_ext], abs(gradθ + gradϕ))
                            end        
                        end
                    end
                end
            end
        end
    end
    return criterion
end


function adapt_angular_grid_3Dby2D!(criterion,thresholds,ref_level,nelem,ngl,nelem_ang, nop_ang, neighbors, npoin_ang,
        connijk_ang, coords_ang, Je_ang, dξdx_ang, dxdξ_ang, dξdy_ang, dydξ_ang, dηdy_ang, dydη_ang, dηdx_ang, dxdη_ang, connijk, 
        x, y, z, xmin_grid, ymin_grid, zmin_grid, xmax_grid, ymax_grid, zmax_grid, ψ, dψ)

    lgl = basis_structs_ξ_ω!(LGL(), ngl-1, CPU())
    adapted_ang = zeros(Int, nelem)
    ang_adapted = zeros(Int, nelem)
    #loop through all spatial elements
    for iel = 1:nelem
        #loop through angular elements
        original_e_ext = nelem_ang[iel] #save original number of angular elements
        ang_adapt = false
        e_ext = 1
        while (e_ext <= nelem_ang[iel])
            #determine if angular element is to be adapted
            level = min(ref_level[iel][e_ext]+1, size(thresholds,1))
            #if (abs(criterion[iel][e_ext]) > 0.0001)
            #if criterion[iel][e_ext] > thresholds[level] && level < size(thresholds,1)
            if (iel == 1 && (e_ext == 5))#criterion[iel][e_ext] > 0.02)
                @info iel, e_ext, "to be adapted"
                adapted_ang[iel] = 1
                ref_level[iel][e_ext] += 1
                ang_adapted[iel] = 1
                # Make new angular elements and reconstruct extra_mesh arrays
                ang_adapt = true
                nelem_ang[iel] += 3
                npoin_ang[iel] += (4*(nop_ang[iel][e_ext]+1)^2-4*(nop_ang[iel][e_ext]+1))
                θmax = coords_ang[iel][1,connijk_ang[iel][e_ext,nop_ang[iel][e_ext]+1,nop_ang[iel][e_ext]+1]]
                θmin = coords_ang[iel][1,connijk_ang[iel][e_ext,1,1]]
                ϕmax = coords_ang[iel][2,connijk_ang[iel][e_ext,nop_ang[iel][e_ext]+1,nop_ang[iel][e_ext]+1]]
                ϕmin = coords_ang[iel][2,connijk_ang[iel][e_ext,1,1]]
                if (ϕmax == 0)
                    ϕmax = 2*π
                end
                
                θ12 = (θmax + θmin)/2
                ϕ12 = (ϕmax + ϕmin)/2
                connijk_ang_new = zeros(Int, nelem_ang[iel], nop_ang[iel][1]+1, nop_ang[iel][1]+1)
                coords_new = zeros(Float64, 2, npoin_ang[iel])
                metrics = allocate_metrics(NSD_2D(), nelem_ang[iel], 0, nop_ang[iel][e_ext], TFloat, CPU())
                nop_ang_new = zeros(Int,nelem_ang[iel])
                nop_ang_new[1:nelem_ang[iel]-1] .= nop_ang[iel][e_ext]
                nop_ang_new[nelem_ang[iel]] = nop_ang[iel][1]
                iter = 1
                criterion_new = zeros(Float64,nelem_ang[iel])
                ref_level_new = zeros(Int,nelem_ang[iel])
                points = []
                points_c = [5.5678199*10^23, 5.5678199*10^45] .+ Array{Float64}(undef,2)
                #populate the elements coming before
                if (e_ext > 1)
                    for e_ext1=1:e_ext-1
                        for i=1:nop_ang[iel][e_ext1]+1
                            for j=1:nop_ang[iel][e_ext1]+1
                                #connijk_ang_new[e_ext1,i,j] = iter
                                θ = coords_ang[iel][1,connijk_ang[iel][e_ext1,i,j]]
                                ϕ = coords_ang[iel][2,connijk_ang[iel][e_ext1,i,j]]
                                metrics.dxdξ[e_ext1, i, j]  = dxdξ_ang[iel][e_ext1, i, j]
                                metrics.Je[e_ext1, i, j]  = Je_ang[iel][e_ext1, i, j]
                                metrics.dξdx[e_ext1, i, j]  = dξdx_ang[iel][e_ext1, i, j]
                                metrics.dxdη[e_ext1, i, j]  = dxdη_ang[iel][e_ext1, i, j]
                                metrics.dηdx[e_ext1, i, j]  = dηdx_ang[iel][e_ext1, i, j]
                                metrics.dηdy[e_ext1, i, j]  = dηdy_ang[iel][e_ext1, i, j]
                                metrics.dξdy[e_ext1, i, j]  = dξdy_ang[iel][e_ext1, i, j]
                                metrics.dydη[e_ext1, i, j]  = dydη_ang[iel][e_ext1, i, j]
                                metrics.dydξ[e_ext1, i, j]  = dydξ_ang[iel][e_ext1, i, j]
                                criterion_new[e_ext1] = criterion[iel][e_ext1]
                                ref_level_new[e_ext1] = ref_level[iel][e_ext1]
                                if (vec([θ,ϕ]) in eachrow(points_c'))
                                    found = false
                                    iter1 = 1
                                    while (found == false && iter1 < iter)
                                        if (points_c[1,iter1] == θ) && (points_c[2,iter1] == ϕ)
                                            connijk_ang_new[e_ext1,i,j] = iter-iter1
                                            found = true
                                        else
                                            iter1 +=1
                                        end
                                    end
                                else
                                    connijk_ang_new[e_ext1,i,j] = iter
                                    coords_new[1,iter] = θ
                                    coords_new[2,iter] = ϕ
                                    points_c = hcat([θ,ϕ],points_c)
                                    push!(points,iter)
                                    iter += 1
                                end
                                if (connijk_ang_new[e_ext1,i,j] == 0) @info e_ext1,i,j, "1",vec([θ,ϕ]), points_c' end
                                #if (i != nop_ang[iel][e_ext1]+1) && (j != nop_ang[iel][e_ext1]+1) iter +=1 end
                            end
                        end
                    end
                end
                #populate for the new elements
                lgl = basis_structs_ξ_ω!(LGL(), nop_ang[iel][1], CPU())
                for i=1:nop_ang_new[e_ext]+1
                    ξi = lgl.ξ[i]
                    for j=1:nop_ang_new[e_ext]+1
                        ξj = lgl.ξ[j]
                        θ = θmin*(1.0-ξi)*0.5+θ12*(1.0 + ξi)*0.5
                        ϕ = ϕmin*(1.0-ξj)*0.5+ϕ12*(1.0 + ξj)*0.5
                        #connijk_ang_new[e_ext,i,j]    = iter
                        criterion_new[e_ext] = 0.0#criterion[iel][e_ext]
                        ref_level_new[e_ext] = ref_level[iel][e_ext]
                        if (vec([θ,ϕ]) in eachrow(points_c'))
                            found = false
                            iter1 = 1
                            while (found == false && iter1 < iter)
                                if (points_c[1,iter1] == θ) && (points_c[2,iter1] == ϕ)
                                    connijk_ang_new[e_ext,i,j] = iter-iter1
                                    found = true
                                else
                                    iter1 +=1
                                end
                            end
                        else
                            connijk_ang_new[e_ext,i,j] = iter
                            coords_new[1,iter] = θ
                            coords_new[2,iter] = ϕ
                            points_c = hcat([θ,ϕ],points_c)
                            push!(points,iter)
                            iter += 1
                        end
                        if (connijk_ang_new[e_ext,i,j] == 0) @info e_ext,i,j, "2" end
                        #if (i != nop_ang[iel][e_ext]+1) && (j != nop_ang[iel][e_ext]+1) iter +=1 end
                    end
                end
                for i=1:nop_ang_new[e_ext+1]+1
                    ξi = lgl.ξ[i]
                    for j=1:nop_ang_new[e_ext+1]+1
                        ξj = lgl.ξ[j]
                        θ = θ12*(1.0-ξi)*0.5+θmax*(1.0 + ξi)*0.5
                        ϕ = ϕmin*(1.0-ξj)*0.5+ϕ12*(1.0 + ξj)*0.5
                        #connijk_ang_new[e_ext+1,i,j] = iter
                        criterion_new[e_ext+1] = 0.0#criterion[iel][e_ext]
                        ref_level_new[e_ext+1] = ref_level[iel][e_ext]
                        if (vec([θ,ϕ]) in eachrow(points_c'))
                            found = false
                            iter1 = 1 
                            while (found == false && iter1 < iter)
                                if (points_c[1,iter1] == θ) && (points_c[2,iter1] == ϕ)
                                    connijk_ang_new[e_ext+1,i,j] = iter-iter1
                                    found = true
                                else
                                    iter1 +=1
                                end 
                            end 
                        else
                            connijk_ang_new[e_ext+1,i,j] = iter
                            coords_new[1,iter] = θ 
                            coords_new[2,iter] = ϕ 
                            points_c = hcat([θ,ϕ],points_c)
                            push!(points,iter)
                            iter += 1
                        end
                        if (connijk_ang_new[e_ext+1,i,j] == 0) @info e_ext+1,i,j, "3" end
                        #if (i != nop_ang_new[e_ext+1]+1) && (j != nop_ang[iel][e_ext]+1) iter +=1 end
                    end
                end

                for i=1:nop_ang_new[e_ext+2]+1
                    ξi = lgl.ξ[i]
                    for j=1:nop_ang_new[e_ext+2]+1
                        ξj = lgl.ξ[j]
                        θ = θmin*(1.0-ξi)*0.5+θ12*(1.0 + ξi)*0.5
                        ϕ = ϕ12*(1.0-ξj)*0.5+ϕmax*(1.0 + ξj)*0.5
                        #connijk_ang_new[e_ext+2,i,j]    = iter
                        criterion_new[e_ext+2] = 0.0#criterion[iel][e_ext]
                        ref_level_new[e_ext+2] = ref_level[iel][e_ext]
                        if (vec([θ,ϕ]) in eachrow(points_c'))
                            found = false
                            iter1 = 1 
                            while (found == false && iter1 < iter)
                                if (points_c[1,iter1] == θ) && (points_c[2,iter1] == ϕ)
                                    connijk_ang_new[e_ext+2,i,j] = iter-iter1
                                    found = true
                                else
                                    iter1 +=1
                                end 
                            end 
                        else
                            connijk_ang_new[e_ext+2,i,j] = iter
                            coords_new[1,iter] = θ 
                            coords_new[2,iter] = ϕ 
                            points_c = hcat([θ,ϕ],points_c)
                            push!(points,iter)
                            iter += 1
                        end
                        if (connijk_ang_new[e_ext+2,i,j] == 0) @info e_ext+2,i,j, "3" end
                        #if (i != nop_ang[iel][e_ext+2]+1) && (j != nop_ang[iel][e_ext]+1) iter +=1 end
                    end
                end
                for i=1:nop_ang_new[e_ext+3]+1
                    ξi = lgl.ξ[i]
                    for j=1:nop_ang_new[e_ext+3]+1
                        ξj = lgl.ξ[j]
                        θ = θ12*(1.0-ξi)*0.5+θmax*(1.0 + ξi)*0.5
                        ϕ = ϕ12*(1.0-ξj)*0.5+ϕmax*(1.0 + ξj)*0.5
                        #connijk_ang_new[e_ext+3,i,j] = iter
                        criterion_new[e_ext+3] = 0.0#criterion[iel][e_ext]
                        ref_level_new[e_ext+3] = ref_level[iel][e_ext]
                        if (vec([θ,ϕ]) in eachrow(points_c'))
                            found = false
                            iter1 = 1 
                            while (found == false && iter1 < iter)
                                if (points_c[1,iter1] == θ) && (points_c[2,iter1] == ϕ)
                                    connijk_ang_new[e_ext+3,i,j] = iter-iter1
                                    found = true
                                else
                                    iter1 +=1
                                end 
                            end 
                        else
                            connijk_ang_new[e_ext+3,i,j] = iter
                            coords_new[1,iter] = θ 
                            coords_new[2,iter] = ϕ 
                            points_c = hcat([θ,ϕ],points_c)
                            push!(points,iter)
                            iter += 1
                        end
                        if (connijk_ang_new[e_ext+3,i,j] == 0) @info e_ext+3,i,j, "4" end
                        #if (i != nop_ang_new[e_ext+3]+1) && (j != nop_ang[iel][e_ext]+1) iter +=1 end
                    end
                end
                if (e_ext < nelem_ang[iel]-1)
                    for e_ext1=e_ext+4:nelem_ang[iel]
                        for i=1:nop_ang_new[e_ext1]+1
                            for j=1:nop_ang_new[e_ext1]+1
                                #connijk_ang_new[e_ext1,i,j] = iter
                                θ = coords_ang[iel][1,connijk_ang[iel][e_ext1-3,i,j]]
                                ϕ = coords_ang[iel][2,connijk_ang[iel][e_ext1-3,i,j]]
                                metrics.dxdξ[e_ext1, i, j]  = dxdξ_ang[iel][e_ext1-3, i, j]
                                metrics.Je[e_ext1, i, j]  = Je_ang[iel][e_ext1-3, i, j]
                                metrics.dξdx[e_ext1, i, j]  = dξdx_ang[iel][e_ext1-3, i, j]
                                metrics.dxdη[e_ext1, i, j]  = dxdη_ang[iel][e_ext1-3, i, j] 
                                metrics.dηdx[e_ext1, i, j]  = dηdx_ang[iel][e_ext1-3, i, j]
                                metrics.dηdy[e_ext1, i, j]  = dηdy_ang[iel][e_ext1-3, i, j]
                                metrics.dξdy[e_ext1, i, j]  = dξdy_ang[iel][e_ext1-3, i, j]
                                metrics.dydη[e_ext1, i, j]  = dydη_ang[iel][e_ext1-3, i, j] 
                                metrics.dydξ[e_ext1, i, j]  = dydξ_ang[iel][e_ext1-3, i, j]
                                criterion_new[e_ext1] = criterion[iel][e_ext1-3]
                                ref_level_new[e_ext1] = ref_level[iel][e_ext1-3]
                                if (vec([θ,ϕ]) in eachrow(points_c'))
                                    found = false
                                    iter1 = 1
                                    while (found == false && iter1 < iter)
                                        if (points_c[1,iter1] == θ) && (points_c[2,iter1] == ϕ)
                                            connijk_ang_new[e_ext1,i,j] = iter-iter1
                                            found = true
                                        else
                                            iter1 +=1
                                        end
                                    end
                                else
                                    connijk_ang_new[e_ext1,i,j] = iter
                                    coords_new[1,iter] = θ
                                    coords_new[2,iter] = ϕ
                                    points_c = hcat([θ,ϕ],points_c)
                                    push!(points,iter)
                                    iter += 1
                                end
                                if (connijk_ang_new[e_ext1,i,j] == 0) @info e_ext1,i,j, "5" end
                                #if (i != nop_ang_new[e_ext1]+1) && (j != nop_ang[iel][e_ext]+1) iter +=1 end
                            end
                        end
                    end
                end
                npoin_ang[iel] = iter-1
                aux_coords = zeros(Float64, 2, npoin_ang[iel])
                aux_coords[:,:] .= coords_new[:,1:npoin_ang[iel]]
                coords_new = aux_coords
                #build metrics for new elements
                for e_ext1 = e_ext:e_ext+3
                    for i = 1:nop_ang_new[e_ext1]+1
                        for j=1:nop_ang_new[e_ext1] + 1
                            ip = connijk_ang_new[e_ext1, i, j]
                            θij = coords_new[1,ip]
                            ϕij = coords_new[2,ip]
                
                            @turbo for l=1:nop_ang_new[e_ext1]+1
                                for k=1:nop_ang_new[e_ext1]+1
        
                                    a = dψ[i,k]*ψ[j,l]
                                    b = ψ[i,k]*dψ[j,l]
                                    metrics.dxdξ[e_ext1, k, l] += a * θij
                                    metrics.dxdη[e_ext1, k, l] += b * θij
        
                                    metrics.dydξ[e_ext1, k, l] += a * ϕij
                                    metrics.dydη[e_ext1, k, l] += b * ϕij
                                end
                            end
                        end
                    end
                    @inbounds for l = 1:nop_ang_new[e_ext1]+1
                        for k = 1:nop_ang_new[e_ext1]+1

                            # Extract values from memory once per iteration
                            dxdξ_val = metrics.dxdξ[e_ext1, k, l]
                            dydη_val = metrics.dydη[e_ext1, k, l]
                            dydξ_val = metrics.dydξ[e_ext1, k, l]
                            dxdη_val = metrics.dxdη[e_ext1, k, l]
                            # Compute Je once and reuse its value
                            metrics.Je[e_ext1, k, l] = dxdξ_val * dydη_val - dydξ_val * dxdη_val
                            # Use the precomputed Je value for the other calculations
                            Jinv = 1.0/metrics.Je[e_ext1, k, l]

                            metrics.dξdx[e_ext1, k, l] =  dydη_val * Jinv
                            metrics.dξdy[e_ext1, k, l] = -dxdη_val * Jinv
                            metrics.dηdx[e_ext1, k, l] = -dydξ_val * Jinv
                            metrics.dηdy[e_ext1, k, l] =  dxdξ_val * Jinv

                        end
                    end

                end
                # enforce periodicity for new elements if applicable
                    for rep = 1:2
                        for iper=1:npoin_ang[iel]
                            θ = coords_new[1,iper]
                            ϕ = coords_new[2,iper]
                            if (abs(ϕ/π - 2.0) <= eps(Float64))
                                #found a periodic point
                                iper1 = 1
                                found = false
                                while (iper1 <= npoin_ang[iel] && found == false)
                                    θ1 = coords_new[1,iper1]
                                    ϕ1 = coords_new[2,iper1]
                                    if (ϕ1 <= eps(Float64) && abs(θ-θ1) <= eps(Float64))
                                        found = true
                                    end
                                    iper1 += 1
                                end
                                if (found)
                                    ip_old = iper
                                    ip_new = iper1-1
                                    for e=1:nelem_ang[iel]
                                        for i=1:nop_ang_new[e]+1
                                            for j=1:nop_ang_new[e]+1
                                                ip = connijk_ang_new[e,i,j]
                                                if (ip == ip_old)
                                                    connijk_ang_new[e,i,j] = ip_new
                                                    coords_new[1,ip] = coords_new[1,ip_new]
                                                    coords_new[2,ip] = coords_new[2,ip_new]
                                                end
                                            end
                                        end
                                    end
                                    for e=1:nelem_ang[iel]
                                        for i=1:nop_ang_new[e]+1
                                            for j=1:nop_ang_new[e]+1
                                                ip = connijk_ang_new[e,i,j]
                                                if (ip >= ip_old)
                                                    connijk_ang_new[e,i,j] -= 1
                                                end
                                            end
                                        end
                                    end
                                    for i = ip_old+1: npoin_ang[iel]
                                        coords_new[1,i-1] = coords_new[1,i]
                                        coords_new[2,i-1] = coords_new[2,i]
                                    end
                                    npoin_ang[iel] -= 1
                                end
                            end
                        end
                    end
 
                connijk_ang[iel] = connijk_ang_new
                coords_ang[iel] = coords_new
                dxdξ_ang[iel] = metrics.dxdξ[:,:,:]
                dxdη_ang[iel] = metrics.dxdη[:,:,:]
                dydξ_ang[iel] = metrics.dydξ[:,:,:]
                dydη_ang[iel] = metrics.dydη[:,:,:]
                Je_ang[iel] = metrics.Je[:,:,:]
                dξdx_ang[iel] = metrics.dξdx[:,:,:]
                dξdy_ang[iel] = metrics.dξdy[:,:,:]
                dηdx_ang[iel] = metrics.dηdx[:,:,:]
                dηdy_ang[iel] = metrics.dηdy[:,:,:]
                nop_ang[iel] = nop_ang_new
                criterion[iel] = criterion_new
                ref_level[iel] = ref_level_new
            end
            e_ext += 1
        end
    end
    for iel = 1:nelem
        # Find and save spatial neighbors for non-conforming assembly
            #First find element corners
            xmin = 10^10
            xmax = -10^10
            ymin = 10^10
            ymax = -10^10
            zmin = 10^10
            zmax = -10^10

            for i=1:ngl
                for j=1:ngl
                    for k=1:ngl
                        ip = connijk[iel,i,j,k]
                        if (x[ip] < xmin) xmin = x[ip] end
                        if (y[ip] < ymin) ymin = y[ip] end
                        if (z[ip] < zmin) zmin = z[ip] end
                        if (x[ip] > xmax) xmax = x[ip] end
                        if (y[ip] > ymax) ymax = y[ip] end
                        if (z[ip] > zmax) zmax = z[ip] end
                    end
                end
            end
            match_bdy = 0
            if (xmin == xmin_grid) match_bdy += 1 end
            if (xmax == xmax_grid) match_bdy += 1 end
            if (ymin == ymin_grid) match_bdy += 1 end
            if (ymax == ymax_grid) match_bdy += 1 end
            if (zmin == zmin_grid) match_bdy += 1 end
            if (zmax == zmax_grid) match_bdy += 1 end
            iter = 1
            found_neighbors = 0
            while (iter <= nelem && found_neighbors <26)
                #find corners for comparison
                xmin_i = 10^10
                xmax_i = -10^10
                ymin_i = 10^10
                ymax_i = -10^10
                zmin_i = 10^10
                zmax_i = -10^10
                for i=1:ngl
                    for j=1:ngl
                        for k=1:ngl
                            ip = connijk[iter,i,j,k]
                            if (x[ip] < xmin_i) xmin_i = x[ip] end
                            if (y[ip] < ymin_i) ymin_i = y[ip] end
                            if (x[ip] > xmax_i) xmax_i = x[ip] end
                            if (y[ip] > ymax_i) ymax_i = y[ip] end
                            if (z[ip] < zmin_i) zmin_i = z[ip] end
                            if (z[ip] > zmax_i) zmax_i = z[ip] end
                        end
                    end
                end
                match_neighbor1 = 0
                match_neighbor2 = 0
                if (abs(xmin - xmin_i) < 1e-5) match_neighbor1 += 1 end
                if (abs(xmax - xmax_i) < 1e-5) match_neighbor1 += 1 end
                if (abs(ymin - ymin_i) < 1e-5) match_neighbor1 += 1 end
                if (abs(ymax - ymax_i) < 1e-5) match_neighbor1 += 1 end
                if (abs(xmin - xmax_i) < 1e-5) match_neighbor2 += 1 end
                if (abs(xmax - xmin_i) < 1e-5) match_neighbor2 += 1 end
                if (abs(ymin - ymax_i) < 1e-5) match_neighbor2 += 1 end
                if (abs(ymax - ymin_i) < 1e-5) match_neighbor2 += 1 end
                if (abs(zmin - zmin_i) < 1e-5) match_neighbor1 += 1 end
                if (abs(zmax - zmax_i) < 1e-5) match_neighbor1 += 1 end
                if (abs(zmin - zmax_i) < 1e-5) match_neighbor2 += 1 end
                if (abs(zmax - zmin_i) < 1e-5) match_neighbor2 += 1 end

                if (match_neighbor1 > 2 && match_neighbor2 > 0 ) || (match_neighbor1 > 0 && match_neighbor2 > 1) || match_neighbor2 > 2
                    found_neighbors += 1
                    neighbors[iel,found_neighbors,1] = iter
                    #check for conformity here
                    
                    if !(adapted_ang[iel] == 0 && adapted_ang[iter] == 0) #if no angular refinement has taken place no need to check conformity
                        if (nelem_ang[iel] != nelem_ang[iter] || coords_ang[iel] != coords_ang[iter]) #non conforming
                            neighbors[iel,found_neighbors,2] = 1 # Save information that these neighbors are non conforming
                            
                        end
                    end
                end
                if (match_bdy == 1 && found_neighbors == 17)
                    found_neighbors = 26
                elseif (match_bdy == 2 && found_neighbors == 11)
                        found_neighbors = 26
                elseif (match_bdy == 3 && found_neighbors == 7)
                    found_neighbors = 26
                end
                iter += 1
            end
    end
end

function adaptive_spatial_angular_numbering_3D_2D!(connijk_spa,nelem, ngl, connijk, connijk_ang, nop_ang, nelem_ang, coords_ang, x, y, z, ref_level, neighbors, adapted, extra_meshes_extra_Je)
    points = []
    points_c = [5.5678199*10^23, 5.5678199*10^45, 5.5678199*10^23, 5.5678199*10^45, 5.5678199*10^23] .+ Array{Float64}(undef,5)
    iter = 1
    interp_sourcesθ = zeros(Float64,nop_ang[1][1]+1)
    interp_targetsθ = zeros(Float64,nop_ang[1][1]+1)
    interp_sourcesϕ = zeros(Float64,nop_ang[1][1]+1)
    interp_targetsϕ = zeros(Float64,nop_ang[1][1]+1)
    ωθ = zeros(Float64,nop_ang[1][1]+1)
    ωϕ = zeros(Float64,nop_ang[1][1]+1)
    Lθ = zeros(Float64,nop_ang[1][1]+1,nop_ang[1][1]+1)
    Lϕ = zeros(Float64,nop_ang[1][1]+1,nop_ang[1][1]+1)

    nc_non_global_nodes = []
    n_non_global_nodes = 0

    point_dict = Dict{NTuple{5,Float64}, Int}()
    @time @inbounds for iel = 1:nelem
        for k = 1:ngl
            for j = 1:ngl
                for i = 1:ngl
                    ip = connijk[iel, i, j, k]
                    x_p = x[ip]
                    y_p = y[ip]
                    z_p = z[ip]
                    
                    for e_ext = 1:nelem_ang[iel]
                        for jθ = 1:nop_ang[iel][e_ext]+1
                            for iθ = 1:nop_ang[iel][e_ext]+1
                                ip_ang = connijk_ang[iel][e_ext, iθ, jθ]
                                θ_p = coords_ang[iel][1, ip_ang]
                                ϕ_p = coords_ang[iel][2, ip_ang]
                                
                                # Create unique key - round to avoid floating point issues
                                key = (round(x_p, digits=12), round(y_p, digits=12), 
                                       round(z_p, digits=12), round(θ_p, digits=12), 
                                       round(ϕ_p, digits=12))
                                
                                # O(1) lookup instead of O(n) search
                                idx = get(point_dict, key, 0)
                                if idx == 0
                                    point_dict[key] = iter
                                    connijk_spa[iel][i, j, k, e_ext, iθ, jθ] = iter
                                    iter += 1
                                else
                                    connijk_spa[iel][i, j, k, e_ext, iθ, jθ] = idx
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    @info "finished non-adaptive connectivity"
    @info "total number of independent points", iter-1
    n_spa = iter-1
    rep_ip = zeros(Int,ngl*ngl)
    interpolation_cache = Dict{NTuple{4,Int}, Tuple{Matrix{Float64}, Matrix{Float64}, Int, Int}}()
    nc_non_global_set = Set{Int}()
    iter_nc = 1
    @info "Identifying non-conforming nodes..."
    for iel = 1:nelem
        if !(1 in neighbors[iel, :, 2])
            continue
        end
        
        for ineighbor = 1:26
            if neighbors[iel, ineighbor, 2] != 1
                continue
            end
            
            adapted = true
            iel1 = neighbors[iel, ineighbor, 1]
            
            # Find matching face nodes efficiently
            matching_nodes = find_matching_face_nodes_optimized(ngl, iel, iel1, connijk)
            
            for (i, j, k, i1, j1, k1) in matching_nodes
                for e_ext = 1:nelem_ang[iel]
                    # Check for exact angular match
                    if has_exact_angular_match(iel, iel1, e_ext, nelem_ang, coords_ang, connijk_ang, nop_ang)
                        continue
                    end
                    
                    # Get bounds
                    θmin, θmax, ϕmin, ϕmax = get_element_bounds_fast(
                        iel, e_ext, coords_ang, connijk_ang, nop_ang
                    )
                    
                    # Find child elements
                    for e_ext1 = 1:nelem_ang[iel1]
                        θmin1, θmax1, ϕmin1, ϕmax1 = get_element_bounds_fast(
                            iel1, e_ext1, coords_ang, connijk_ang, nop_ang
                        )
                        
                        if is_child_element(θmin, θmax, ϕmin, ϕmax, θmin1, θmax1, ϕmin1, ϕmax1)
                            # Build or retrieve interpolation matrices
                            cache_key = (iel, e_ext, iel1, e_ext1)
                            
                            if !haskey(interpolation_cache, cache_key)
                                Lθ_cached, Lϕ_cached = build_interpolation_matrices!(
                                    iel, e_ext, iel1, e_ext1, coords_ang, connijk_ang, nop_ang,
                                    θmin, θmax, ϕmin, ϕmax, ϕmin1, ϕmax1,
                                    interp_sourcesθ, interp_targetsθ, interp_sourcesϕ, interp_targetsϕ,
                                    ωθ, ωϕ, Lθ, Lϕ
                                )
                                nop_parent = nop_ang[iel][e_ext]
                                nop_child = nop_ang[iel1][e_ext1]
                                interpolation_cache[cache_key] = (
                                    copy(Lθ_cached), copy(Lϕ_cached), nop_parent, nop_child
                                )
                            end
                            
                            Lθ_use, Lϕ_use, nop_parent, nop_child = interpolation_cache[cache_key]
                            
                            # Identify hanging nodes
                            for jθ = 1:(nop_child+1), iθ = 1:(nop_child+1)
                                is_vertex_θ = (1 in Lθ_use[iθ,:])
                                is_vertex_ϕ = (1 in Lϕ_use[jθ,:])
                                
                                if !is_vertex_θ || !is_vertex_ϕ
                                    jp_spa = connijk_spa[iel1][i1, j1, k1, e_ext1, iθ, jθ]
                                    if !(jp_spa in nc_non_global_nodes)
                                        push!(nc_non_global_set, jp_spa)
                                        n_non_global_nodes += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    nc_non_global_nodes = sort(collect(nc_non_global_set))
    n_non_global_nodes = length(nc_non_global_nodes)
    @info "Number of NC nodes: $n_non_global_nodes"
    if n_non_global_nodes > 0
        @info "Renumbering nodes..."
        node_map = Dict{Int, Int}()
        
        free_counter = 1
        for old_idx = 1:n_spa
            if !(old_idx in nc_non_global_set)
                node_map[old_idx] = free_counter
                free_counter += 1
            end
        end
        
        hanging_counter = n_spa - n_non_global_nodes + 1
        for old_idx in nc_non_global_nodes
            node_map[old_idx] = hanging_counter
            hanging_counter += 1
        end
        
        # Apply renumbering
        @inbounds for iel = 1:nelem
            for k = 1:ngl, j = 1:ngl, i = 1:ngl
                for e_ext = 1:nelem_ang[iel]
                    for jθ = 1:nop_ang[iel][e_ext]+1, iθ = 1:nop_ang[iel][e_ext]+1
                        old_idx = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                        connijk_spa[iel][i, j, k, e_ext, iθ, jθ] = node_map[old_idx]
                    end
                end
            end
        end
        
        nc_non_global_nodes = [node_map[old_idx] for old_idx in nc_non_global_nodes]
        sort!(nc_non_global_nodes)
    end

    @info "Building constraint matrix..."
    
    n_free = n_spa - n_non_global_nodes
    
    # Estimate size for triplet format
    estimated_entries = n_free + n_non_global_nodes * 20  # Rough estimate
    I_vec = sizehint!(Int[], estimated_entries)
    J_vec = sizehint!(Int[], estimated_entries)
    V_vec = sizehint!(Float64[], estimated_entries)
    
    # Add identity for free nodes
    for ip_g = 1:n_free
        push!(I_vec, ip_g)
        push!(J_vec, ip_g)
        push!(V_vec, 1.0)
    end
    
    # Build constraint weights for hanging nodes
    @info "Processing hanging node constraints..."
    @inbounds for iel = 1:nelem
        if !(1 in neighbors[iel, :, 2])
            continue
        end
        
        for ineighbor = 1:26
            if neighbors[iel, ineighbor, 2] != 1
                continue
            end
            
            iel1 = neighbors[iel, ineighbor, 1]
            matching_nodes = find_matching_face_nodes_optimized(ngl, iel, iel1, connijk)
            
            for (i, j, k, i1, j1, k1) in matching_nodes
                for e_ext = 1:nelem_ang[iel]
                    if has_exact_angular_match(iel, iel1, e_ext, nelem_ang, coords_ang, connijk_ang, nop_ang)
                        continue
                    end
                    
                    θmin, θmax, ϕmin, ϕmax = get_element_bounds_fast(
                        iel, e_ext, coords_ang, connijk_ang, nop_ang
                    )
                    
                    for e_ext1 = 1:nelem_ang[iel1]
                        θmin1, θmax1, ϕmin1, ϕmax1 = get_element_bounds_fast(
                            iel1, e_ext1, coords_ang, connijk_ang, nop_ang
                        )
                        
                        if !is_child_element(θmin, θmax, ϕmin, ϕmax, θmin1, θmax1, ϕmin1, ϕmax1)
                            continue
                        end
                        
                        # Retrieve cached interpolation matrices
                        cache_key = (iel, e_ext, iel1, e_ext1)
                        Lθ_use, Lϕ_use, nop_parent, nop_child = interpolation_cache[cache_key]
                        
                        # Fill constraint matrix entries
                        for jθ = 1:(nop_parent+1), iθ = 1:(nop_parent+1)
                            ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                            
                            if ip_spa > n_free
                                continue  # Skip hanging nodes as rows
                            end
                            
                            Je_parent = extra_meshes_extra_Je[iel][e_ext, iθ, jθ]
                            
                            for lθ = 1:(nop_child+1), kθ = 1:(nop_child+1)
                                jp_spa = connijk_spa[iel1][i1, j1, k1, e_ext1, kθ, lθ]
                                
                                if jp_spa <= n_free
                                    continue  # Only process hanging nodes as columns
                                end
                                
                                Lθ_val = Lθ_use[kθ, iθ]
                                Lϕ_val = Lϕ_use[lθ, jθ]
                                
                                if abs(Lθ_val) < 1e-14 || abs(Lϕ_val) < 1e-14
                                    continue
                                end
                                
                                # Skip if this is identity (vertex node)
                                if abs(Lθ_val - 1.0) < 1e-14 && abs(Lϕ_val - 1.0) < 1e-14
                                    continue
                                end
                                
                                Je_child = extra_meshes_extra_Je[iel1][e_ext1, kθ, lθ]
                                weight = Lθ_val * Lϕ_val * (Je_parent / Je_child)
                                
                                push!(I_vec, ip_spa)
                                push!(J_vec, jp_spa)
                                push!(V_vec, weight)
                            end
                        end
                    end
                end
            end
        end
    end
    
    # Build sparse matrix
    nc_mat = sparse(I_vec, J_vec, V_vec, n_free, n_spa)
    
    # Normalize columns
    @info "Normalizing constraint matrix..."
    @inbounds for j = 1:n_spa
        col_sum = 0.0
        for idx in nzrange(nc_mat, j)
            col_sum += nonzeros(nc_mat)[idx]
        end
        
        if abs(col_sum - 1.0) > 1e-13 && abs(col_sum) > 1e-14
            for idx in nzrange(nc_mat, j)
                nonzeros(nc_mat)[idx] /= col_sum
            end
        end
    end
    
    @info nnz(nc_mat),estimated_entries, nnz(nc_mat)/estimated_entries
    return nc_mat, nc_non_global_nodes, n_non_global_nodes, n_spa
    
end

function find_face_node_match(ngl,iel,ip,connijk)
    found = false           
    ip1 = 0                 
    iter = 0                   
    i = 0
    j = 0                       
    k = 0                           
    while (found == false && iter <= ngl^3-(ngl-2)^3) 
        if (iter < ngl^2)            
            i = 1
            j = (iter ÷ ngl)+1
            k = iter % ngl+1
        elseif (iter >= ngl^2 && iter < (2*ngl^2)-ngl) 
            i = (iter ÷ ngl)%(ngl-1)+1
            if (i==1) i = ngl end
            j = ngl
            k = iter % ngl+1
        elseif (iter >= (2*ngl^2)-ngl && iter < (3*ngl^2)-(2*ngl))
            i = ngl
            j = (iter ÷ ngl)%(ngl-1)+1
            k = iter % ngl+1
        elseif (iter >= (3*ngl^2)-(2*ngl) && iter < (4*ngl^2)-(4*ngl))                             
            i = (iter ÷ ngl)%(ngl-1)+1
            j = 1
            k = iter % ngl+1
        elseif (iter >= (4*ngl^2)-(4*ngl) && iter < (4*ngl^2)-(4*ngl)+(ngl-2)^2)
            d = ((4*ngl^2)-(4*ngl)) ÷ (ngl-2)
            r = ((4*ngl^2)-(4*ngl)) - d*(ngl-2)
            i = ((iter-r) ÷ (ngl-2))-d + 2  
            j = iter % (ngl-2)+2
            k = 1
        else
            d = ((4*ngl^2)-(4*ngl)+(ngl-2)^2) ÷ (ngl-2)
            r = ((4*ngl^2)-(4*ngl)+(ngl-2)^2) - d*(ngl-2)
            i = ((iter-r) ÷ (ngl-2))-d + 2 
            j = iter % (ngl-2)+2
            k = ngl
        end                         
        if (connijk[iel,i,j,k] == ip)
            ip1 = connijk[iel,i,j,k]
            found = true
        end         
        iter += 1
            
    end 
    if (found == false)
        @info "failed to find", iel, iter, i, j, k 
    end
    return i, j, k, ip1
end     

function find_matching_face_nodes_optimized(ngl, iel, iel1, connijk)
    """Find all matching face nodes between two elements"""
    matching = Tuple{Int,Int,Int,Int,Int,Int}[]
    
    # Create set for O(1) lookup
    iel1_nodes = Set(connijk[iel1, :, :, :])
    
    # Check all face nodes
    @inbounds for k = 1:ngl, j = 1:ngl, i = 1:ngl
        is_face = (i == 1 || i == ngl || j == 1 || j == ngl || k == 1 || k == ngl)
        if !is_face
            continue
        end
        
        ip = connijk[iel, i, j, k]
        if ip in iel1_nodes
            i1, j1, k1, ip1 = find_face_node_match(ngl, iel1, ip, connijk)
            if ip1 != 0
                push!(matching, (i, j, k, i1, j1, k1))
            end
        end
    end
    
    return matching
end

function has_exact_angular_match(iel, iel1, e_ext, nelem_ang, coords_ang, connijk_ang, nop_ang)
    """Check if two angular elements are identical"""
    for e_check = 1:nelem_ang[iel1]
        if (coords_ang[iel][1, connijk_ang[iel][e_ext, :, :]] == 
            coords_ang[iel1][1, connijk_ang[iel1][e_check, :, :]] &&
            coords_ang[iel][2, connijk_ang[iel][e_ext, :, :]] == 
            coords_ang[iel1][2, connijk_ang[iel1][e_check, :, :]])
            return true
        end
    end
    return false
end

function get_element_bounds_fast(iel, e_ext, coords_ang, connijk_ang, nop_ang)
    """Get angular element bounds"""
    nop = nop_ang[iel][e_ext]
    
    θmin = coords_ang[iel][1, connijk_ang[iel][e_ext, 1, 1]]
    θmax = coords_ang[iel][1, connijk_ang[iel][e_ext, nop+1, nop+1]]
    ϕmin = coords_ang[iel][2, connijk_ang[iel][e_ext, 1, 1]]
    ϕmax = coords_ang[iel][2, connijk_ang[iel][e_ext, nop+1, nop+1]]
    
    ϕmax = ϕmax == 0.0 ? 2π : ϕmax
    
    return minmax(θmin, θmax)..., minmax(ϕmin, ϕmax)...
end

function is_child_element(θmin, θmax, ϕmin, ϕmax, θmin1, θmax1, ϕmin1, ϕmax1; tol=1e-14)
    """Check if element 1 is a child of element 0"""
    return (θmin1 >= θmin - tol && θmax1 <= θmax + tol &&
            ϕmin1 >= ϕmin - tol && ϕmax1 <= ϕmax + tol)
end

function build_interpolation_matrices!(iel, e_ext, iel1, e_ext1, coords_ang, connijk_ang, nop_ang,
                                      θmin, θmax, ϕmin, ϕmax, ϕmin1, ϕmax1,
                                      interp_sourcesθ, interp_targetsθ, interp_sourcesϕ, interp_targetsϕ,
                                      ωθ, ωϕ, Lθ, Lϕ)
    """Build and normalize interpolation matrices"""
    nop_parent = nop_ang[iel][e_ext]
    nop_child = nop_ang[iel1][e_ext1]
    
    # Setup source points
    for iθ = 1:(nop_parent+1)
        interp_sourcesθ[iθ] = coords_ang[iel][1, connijk_ang[iel][e_ext, iθ, iθ]]
        interp_sourcesϕ[iθ] = coords_ang[iel][2, connijk_ang[iel][e_ext, iθ, iθ]]
        if ϕmax == 2π && iθ == nop_parent + 1
            interp_sourcesϕ[iθ] = ϕmax
        end
    end
    
    # Setup target points
    for iθ = 1:(nop_child+1)
        interp_targetsθ[iθ] = coords_ang[iel1][1, connijk_ang[iel1][e_ext1, iθ, iθ]]
        interp_targetsϕ[iθ] = coords_ang[iel1][2, connijk_ang[iel1][e_ext1, iθ, iθ]]
        if ϕmax1 == 2π && iθ == nop_child + 1
            interp_targetsϕ[iθ] = ϕmax1
        end
    end
    
    # Build interpolation matrices
    fill!(Lθ, 0.0)
    fill!(Lϕ, 0.0)
    fill!(ωθ, 0.0)
    fill!(ωϕ, 0.0)
    
    BarycentricWeights!(view(interp_sourcesθ, 1:nop_parent+1), view(ωθ, 1:nop_parent+1))
    BarycentricWeights!(view(interp_sourcesϕ, 1:nop_parent+1), view(ωϕ, 1:nop_parent+1))
    
    PolynomialInterpolationMatrix!(
        view(interp_sourcesθ, 1:nop_parent+1), view(ωθ, 1:nop_parent+1),
        view(interp_targetsθ, 1:nop_child+1), view(Lθ, 1:nop_child+1, 1:nop_parent+1)
    )
    PolynomialInterpolationMatrix!(
        view(interp_sourcesϕ, 1:nop_parent+1), view(ωϕ, 1:nop_parent+1),
        view(interp_targetsϕ, 1:nop_child+1), view(Lϕ, 1:nop_child+1, 1:nop_parent+1)
    )
    
    # Normalize rows
    for i = 1:(nop_child+1)
        row_sum_θ = sum(Lθ[i, 1:nop_parent+1])
        if abs(row_sum_θ) > 1e-14
            Lθ[i, 1:nop_parent+1] ./= row_sum_θ
        end
        
        row_sum_ϕ = sum(Lϕ[i, 1:nop_parent+1])
        if abs(row_sum_ϕ) > 1e-14
            Lϕ[i, 1:nop_parent+1] ./= row_sum_ϕ
        end
    end
    
    return view(Lθ, 1:nop_child+1, 1:nop_parent+1), view(Lϕ, 1:nop_child+1, 1:nop_parent+1)
end