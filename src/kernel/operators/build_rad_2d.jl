using SparseArrays
using KrylovPreconditioners
function build_radiative_transfer_problem(mesh, inputs, neqs, ngl, dψ, ψ, ω, Je, dξdx, dξdy, dηdx, dηdy, nx, ny, elem_to_edge, 
        extra_mesh, QT::Inexact, SD::NSD_2D, AD::ContGal)
    comm = MPI.COMM_WORLD
    npoin = mesh.npoin
    nelem = mesh.nelem
    nc_mat = zeros(Float64,1)
    nc_non_global_nodes = []
    n_non_global_nodes = 0 
    n_spa = 0
    begin_time = time()
    mesh.xmax = MPI.Allreduce(maximum(mesh.x), MPI.MAX, comm)
    mesh.xmin = MPI.Allreduce(minimum(mesh.x), MPI.MIN, comm)
    mesh.ymax = MPI.Allreduce(maximum(mesh.y), MPI.MAX, comm)
    mesh.ymin = MPI.Allreduce(minimum(mesh.y), MPI.MIN, comm)
    if (inputs[:adaptive_extra_meshes])
        extra_meshes_coords = [Vector{Float64}(undef, size(extra_mesh[e].extra_coords,1)) for e in 1:nelem]
        extra_meshes_connijk = [Array{Int}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_Je = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dξdx = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_dxdξ = [Array{Float64}(undef, extra_mesh[e].extra_nelem, extra_mesh[e].extra_nop[1]+1) for e in 1:nelem]
        extra_meshes_extra_nops = [Array{Int}(undef, extra_mesh[e].extra_nelem) for e in 1:nelem]
        extra_meshes_extra_npoins = zeros(Int, nelem)
        extra_meshes_extra_nelems = zeros(Int, nelem)
        extra_meshes_ref_level = [Array{Int}(undef, extra_mesh[e].extra_nelem) for e in 1:nelem]
        npoin_ang_total = 0
        for e=1:nelem
            extra_meshes_coords[e] = extra_mesh[e].extra_coords[:]
            extra_meshes_connijk[e] = extra_mesh[e].extra_connijk
            extra_meshes_extra_Je[e] = extra_mesh[e].extra_metrics.Je[:,:]
            extra_meshes_extra_dξdx[e] = extra_mesh[e].extra_metrics.dξdx[:,:]
            extra_meshes_extra_dxdξ[e] = extra_mesh[e].extra_metrics.dxdξ[:,:]
            extra_meshes_extra_npoins[e] = extra_mesh[e].extra_npoin
            extra_meshes_extra_nelems[e] = extra_mesh[e].extra_nelem
            extra_meshes_extra_nops[e] = extra_mesh[e].extra_nop
            extra_meshes_ref_level[e] = extra_mesh[e].ref_level
            npoin_ang_total += mesh.ngl*mesh.ngl*extra_mesh[e].extra_npoin
        end
        connijk_spa = [Array{Int}(undef, ngl, ngl, extra_meshes_extra_nelems[iel], extra_meshes_extra_nops[iel][1]+1) for iel = 1:nelem]
        @info "building initial adaptive connectivity"
        neighbors = zeros(Int,nelem,8,2)
        adapted = false
        
        nc_mat, nc_mat_div, nc_non_global_nodes, n_non_global_nodes, n_spa  = adaptive_spatial_angular_numbering_2D_1D!(connijk_spa,nelem, ngl, mesh.connijk, 
                                                                                                     extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                                                  extra_meshes_coords, mesh.x, mesh.y,extra_meshes_ref_level, neighbors, adapted)
        @info "built initial adaptive spatial angular connectivity"
        @time LHS = sparse_lhs_assembly_2Dby1D_adaptive(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl, extra_meshes_extra_nelems,
                                   dξdx, dξdy, dηdx, dηdy, extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa, nc_mat, nc_mat_div, adapted, nc_non_global_nodes, n_non_global_nodes,
                                  n_spa)
        @time M = sparse_mass_assembly_2Dby1D_adaptive(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelems,
                                   extra_meshes_extra_npoins, connijk_spa, nc_mat, nc_mat_div, adapted, nc_non_global_nodes, n_non_global_nodes, n_spa)
        total_ip = size(LHS,1)
        @info "built pre-adaptivity matrices"
        @info maximum(LHS), minimum(LHS)
        @info maximum(M), minimum(M)
        one_vec = Vector{Float64}(undef, size(LHS,1))
        fill!(one_vec,Float64(1))
        pointwise_interaction = abs.(LHS) * one_vec
        @info maximum(one_vec), minimum(one_vec), maximum(pointwise_interaction), minimum(pointwise_interaction)
        @time criterion = compute_adaptivity_criterion(pointwise_interaction, nelem, ngl, mesh.connijk, extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems, extra_meshes_coords,
                                                connijk_spa)
        
        @info "criterion computed"
        @time adapt_angular_grid_2Dby1D!(criterion,inputs[:RT_amr_threshold], extra_meshes_ref_level,nelem,ngl,extra_meshes_extra_nelems, extra_meshes_extra_nops, neighbors, extra_meshes_extra_npoins,
                                  extra_meshes_connijk, extra_meshes_coords, extra_meshes_extra_Je, extra_meshes_extra_dξdx, extra_meshes_extra_dxdξ, mesh.connijk,
                                  mesh.x, mesh.y, mesh.xmin, mesh.ymin, mesh.xmax, mesh.ymax) 
        @info "angular mesh adapted"
        if !(maximum(extra_meshes_ref_level[:][:]) == 0)
            connijk_spa = [Array{Int}(undef, ngl, ngl, extra_meshes_extra_nelems[iel], extra_meshes_extra_nops[iel][1]+1) for iel = 1:nelem]
            @time nc_mat, nc_mat_div, nc_non_global_nodes, n_non_global_nodes, n_spa  = adaptive_spatial_angular_numbering_2D_1D!(connijk_spa,nelem, ngl, mesh.connijk, 
                                                               extra_meshes_connijk, extra_meshes_extra_nops, extra_meshes_extra_nelems,
                                                                extra_meshes_coords, mesh.x, mesh.y, extra_meshes_ref_level, neighbors, adapted)
            
            @info "adapted connectivity"
            adapted = true
            @info "number of hanging nodes", n_non_global_nodes
            @time LHS = sparse_lhs_assembly_2Dby1D_adaptive(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, n_spa, nelem, ngl, extra_meshes_extra_nelems,
                                   dξdx, dξdy, dηdx, dηdy, extra_meshes_extra_npoins, inputs[:rad_HG_g], connijk_spa, nc_mat, nc_mat_div, adapted, nc_non_global_nodes, n_non_global_nodes, n_spa)

            #Try this alternative assembly approach
            @info size(nc_mat), size(LHS), size(nc_mat')
            A_test = nc_mat*LHS*nc_mat'

            @time M = sparse_mass_assembly_2Dby1D_adaptive(extra_meshes_ref_level, ω, Je, mesh.connijk, extra_mesh[1].ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh[1].ψ, extra_meshes_connijk,
                                    extra_meshes_extra_Je,
                                    extra_meshes_coords, extra_meshes_extra_nops, npoin_ang_total, nelem, ngl, extra_meshes_extra_nelems,
                                   extra_meshes_extra_npoins, connijk_spa, nc_mat, nc_mat_div, adapted, nc_non_global_nodes, n_non_global_nodes, n_spa)
            @info "built adapted matrices"
            M_test = nc_mat*M*nc_mat'

            @info maximum(A_test), minimum(A_test)
            @info maximum(M_test), minimum(M_test)
        end
        npoin_ang_total = n_spa#maximum(connijk_spa[1])
        counter = 0
        for count=1:size(nc_non_global_nodes,1)
            if (nc_non_global_nodes[count] <= n_spa)
                counter += 1
            end
        end
        @info npoin_ang_total, counter
        npoin_ang_total -= counter
        #invert mass matrix
        I_vec = Vector{Int}()
        J_vec = Vector{Int}()
        V_vec = Vector{Float64}()
        max_entries = npoin_ang_total^2
        sizehint!(I_vec, Int64(round(max_entries*0.0001)))
        sizehint!(J_vec, Int64(round(max_entries*0.0001)))
        sizehint!(V_vec, Int64(round(max_entries*0.0001)))

        for ip=1:npoin_ang_total
            val = 1/M_test[ip,ip]
            push!(I_vec, ip)
            push!(J_vec, ip)
            push!(V_vec, val)
        end
        M_inv = sparse(I_vec, J_vec, V_vec)
        
        #@time M_inv = M \ Matrix(I, size(M)) #M\Diagonal(ones(npoin_ang_total))
        M_inv = sparse(M_inv)
        #M = nothing
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
        A = sparse(M_inv*A_test)
        #M_inv = nothing
        #LHS = nothing
        GC.gc()
        RHS = zeros(TFloat, n_spa)#npoin_ang_total)
        ref = zeros(TFloat, n_spa)
        BDY = zeros(TFloat, n_spa)
        @info size(RHS), size(A),n_spa-n_non_global_nodes
    else
        npoin_ang_total = npoin*extra_mesh.extra_npoin
        @time LHS = sparse_lhs_assembly_2Dby1D(ω, Je, mesh.connijk, extra_mesh.ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk, 
                                    extra_mesh.extra_metrics.Je, 
                                    extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                   dξdx, dξdy, dηdx, dηdy, extra_mesh.extra_npoin, inputs[:rad_HG_g])
        @info "assembled LHS"
        @time M = sparse_mass_assembly_2Dby1D(ω, Je, mesh.connijk, extra_mesh.ωθ, mesh.x, mesh.y, ψ, dψ, extra_mesh.ψ, extra_mesh.extra_connijk,
                                    extra_mesh.extra_metrics.Je,
                                    extra_mesh.extra_coords, extra_mesh.extra_nop, npoin_ang_total, nelem, ngl, extra_mesh.extra_nelem,
                                   extra_mesh.extra_npoin)
        @info "assembled Mass matrix"
        @info nnz(M), nnz(LHS), npoin_ang_total^2, nnz(M)/npoin_ang_total^2, nnz(LHS)/npoin_ang_total^2
        
        # inexact integration makes M diagonal, build the sparse inverse to save space
        ip2gip_extra, gip2owner_extra, gnpoin = setup_global_numbering_extra_dim(mesh.ip2gip, mesh.gip2owner, npoin, extra_mesh.extra_npoin, npoin_ang_total)
        Md = diag(M)
        pM = setup_assembler(SD, Md, ip2gip_extra, gip2owner_extra)
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
        @info size(M_inv), size(LHS)    
        #@time M_inv = M \ Matrix(I, size(M)) #M\Diagonal(ones(npoin_ang_total))
        #M_inv = sparse(M_inv)
        #M = nothing
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
    
        @info "LHS max/min"
        @info maximum(LHS), minimum(LHS)
        @info maximum(M), minimum(M)
        A = sparse(M_inv*LHS)
        #=x = real.(eigvals(Array(LHS)))
        y = imag.(eigvals(Array(LHS)))
        display(Makie.scatter(x, y, label="e-values"))=#
        
        #M_inv = nothing
        #LHS = nothing
        GC.gc()
        ref = zeros(TFloat, npoin_ang_total)
        RHS = zeros(TFloat, npoin_ang_total)
        BDY = zeros(TFloat, npoin_ang_total)
    end
    nc_rows = zeros(Int,1,1)
    A_rows = rowvals(A)
    if (size(nc_mat,1) > 2)
        nc_rows = rowvals(nc_mat)
    end

    @time for iel=1:nelem
        for i=1:ngl
            for j =1:ngl
                ip = mesh.connijk[iel,i,j]
                x = mesh.x[ip]
                y = mesh.y[ip]
                if (inputs[:adaptive_extra_meshes])
                    for e_ext = 1:extra_meshes_extra_nelems[iel]
                        for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                            ip_ext = extra_meshes_connijk[iel][e_ext,iθ]

                            θ = extra_meshes_coords[iel][ip_ext]
                            ip_g = connijk_spa[iel][i,j,e_ext,iθ]
                            κip = 10*exp(-((x-3/2)/3)^2)*exp(-y/2)
                            σip = 0.1*κip
                            gip = exp(-((1. / 3) * (x - (3 / 3.)))^2)#exp(-((x-3/3)/3)^2)
                            dgip = -(2. / 3^2) * (x - (3 / 3.)) * gip
                            hip = exp(-4. * (2 - y) / 2)#exp(-4*(2-y)/2)
                            dhip = (4. / 2) * hip
                            sip = exp(-((96 / (2. * π)) * (θ - (7. * π / 5.)))^2)#exp(-((96/(2*π))*(θ-7*π/5))^2)
                            uip = gip*hip*sip
                            ref[ip_g] = uip
                            propip = (cos(θ)*dgip*hip+sin(θ)*gip*dhip)*sip#cos(θ)*hip*sip*gip*((-2*x)/9 + 2/9) + sin(θ)*gip*sip*2*hip
                            if (ip in mesh.poin_in_bdy_edge)
                                applied = false
                                iedge = elem_to_edge[iel,i,j,1]
                                edge_i = elem_to_edge[iel,i,j,2]
                                matchx = (x == mesh.xmax || x == mesh.xmin)
                                matchy = (y == mesh.ymax || y == mesh.ymin)
                                if (matchx && matchy) ##Do corner boundary case
                                    iedge1 = 1
                                    found = false
                                    iedge_found = 0
                                    edge_found_i = 0
                                    while (iedge1 <= mesh.nedges_bdy && found == false)
                                        for iter = 1:ngl
                                            ip1 = mesh.poin_in_bdy_edge[iedge1,iter]
                                            if (ip1 == ip && iedge != iedge1)
                                                found = true
                                                iedge_found = iedge1
                                                edge_found_i = iter
                                            end
                                        end
                                        iedge1 +=1
                                    end
                                    if (cos(θ) * (nx[iedge,edge_i]+nx[iedge_found,edge_found_i]) + sin(θ)*(ny[iedge,edge_i]+ny[iedge_found,edge_found_i]) < -1e-13) #&& (cos(θ) * nx[iedge_found,edge_found_i] + sin(θ)*ny[iedge_found,edge_found_i] < 0.0)
                                        BDY[ip_g] = user_rad_bc(x,y,θ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                        if (ip_g <= n_spa - n_non_global_nodes)
                                            A[ip_g,:] .= 0.0
                                            A[ip_g,ip_g] = 1.0
                                        end
                                        applied = true
                                    end
                                elseif (cos(θ) * nx[iedge,edge_i] + sin(θ)*ny[iedge,edge_i] < -1e-13)
                                    #BDY[ip_g] = user_rad_bc(x,y,θ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                    if (ip_g <= n_spa - n_non_global_nodes)
                                        BDY[ip_g] = user_rad_bc(x,y,θ)
                                        A[ip_g,:] .= 0.0
                                        A[ip_g,ip_g] = 1.0
                                        applied = true
                                    end
                                    #applied = true
                                end
                                if (applied == false)
                                    RHS[ip_g] = user_rhs(x,y,θ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip)
                                end
                            else
                                RHS[ip_g] = user_rhs(x,y,θ)#(-gip*hip*(user_f!(x,y,θ))*σip + κip*uip +  propip)
                            end
                        end
                    end
                else

                    for e_ext = 1:extra_mesh.extra_nelem
                        for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                            ip_ext = extra_mesh.extra_connijk[e_ext,iθ]
                            θ = extra_mesh.extra_coords[ip_ext]
                            ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                            gip_g = ip2gip_extra[ip_g]
                            κip = 10*exp(-((x-3/2)/3)^2)*exp(-y/2)
                            σip = 0.1*κip
                            gip = exp(-((1. / 3) * (x - (3 / 3.)))^2)#exp(-((x-3/3)/3)^2)
                            dgip = -(2. / 3^2) * (x - (3 / 3.)) * gip
                            hip = exp(-4. * (2 - y) / 2)#exp(-4*(2-y)/2)
                            dhip = (4. / 2) * hip
                            sip = exp(-((96 / (2. * π)) * (θ - (7. * π / 5.)))^2)#exp(-((96/(2*π))*(θ-7*π/5))^2)
                            uip = gip*hip*sip
                            ref[ip_g] = uip
                            propip = (cos(θ)*dgip*hip+sin(θ)*gip*dhip)*sip#cos(θ)*hip*sip*gip*((-2*x)/9 + 2/9) + sin(θ)*gip*sip*2*hip
                            if (ip in mesh.poin_in_bdy_edge)
                                applied = false
                                iedge = elem_to_edge[iel,i,j,1]
                                edge_i = elem_to_edge[iel,i,j,2]
                                matchx = (x == mesh.xmax || x == mesh.xmin)
                                matchy = (y == mesh.ymax || y == mesh.ymin)
                                if (matchx && matchy) ##Do corner boundary case
                                    iedge1 = 1
                                    found = false
                                    iedge_found = 0
                                    edge_found_i = 0
                                    while (iedge1 <= mesh.nedges_bdy && found == false) 
                                        for iter = 1:ngl
                                            ip1 = mesh.poin_in_bdy_edge[iedge1,iter]
                                            if (ip1 == ip && iedge != iedge1)
                                                found = true
                                                iedge_found = iedge1
                                                edge_found_i = iter
                                            end
                                        end
                                        iedge1 +=1
                                    end
                                    if (cos(θ) * (nx[iedge,edge_i]+nx[iedge_found,edge_found_i]) + sin(θ)*(ny[iedge,edge_i]+ny[iedge_found,edge_found_i]) < -1e-13) 
                                        if (gip2owner_extra[ip_g] == rank)
                                            bc_value = user_rad_bc(x,y,θ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                            BDY .-= A[:, ip_g] .* bc_value
                                            A[:, ip_g] .= 0.0
                                            BDY[ip_g] = bc_value#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                            A[ip_g,:] .= 0.0
                                            A[ip_g,ip_g] = 1.0
                                            dropzeros!(A)
                                            applied = true
                                        else
                                            bc_value = user_rad_bc(x,y,θ)
                                            BDY .-= A[:, ip_g] .* bc_value
                                            A[:, ip_g] .= 0.0
                                            BDY[ip_g] = 0.0
                                            A[ip_g,:] .= 0.0
                                            dropzeros!(A)
                                            applied = true
                                        end
                                        
                                    end
                                elseif (cos(θ) * nx[iedge,edge_i] + sin(θ)*ny[iedge,edge_i] < -1e-13)
                                    if (gip2owner_extra[ip_g] == rank)
                                        bc_value = user_rad_bc(x,y,θ)#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                        BDY .-= A[:, ip_g] .* bc_value
                                        A[:, ip_g] .= 0.0
                                        BDY[ip_g] = bc_value
                                        A[ip_g,:] .= 0.0
                                        A[ip_g,ip_g] = 1.0
                                        dropzeros!(A)
                                        applied = true
                                    else
                                        bc_value = user_rad_bc(x,y,θ)
                                        BDY .-= A[:, ip_g] .* bc_value
                                        A[:, ip_g] .= 0.0
                                        BDY[ip_g] = 0.0
                                        A[ip_g,:] .= 0.0
                                        dropzeros!(A)
                                        applied = true
                                    end
                                end
                                if (applied == false)
                                    if (gip2owner_extra[ip_g] == rank)
                                        RHS[ip_g] = user_rhs(x,y,θ)
                                    end
                                end
                            else
                                if (gip2owner_extra[ip_g] == rank)
                                    RHS[ip_g] = user_rhs(x,y,θ)
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    B = RHS + BDY
    
    if (inputs[:adaptive_extra_meshes]) 
        BDY_free = BDY[1:(n_spa - n_non_global_nodes)]
        
        U_red_proj = BDY_free#M_test \ BDY_proj
        
        RHS_red = rest * RHS
        
        for iel=1:nelem
            for i=1:ngl
                for j =1:ngl
                    ip = mesh.connijk[iel,i,j]
                    x = mesh.x[ip]
                    y = mesh.y[ip]
                    for e_ext = 1:extra_meshes_extra_nelems[iel]
                        for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                            ip_ext = extra_meshes_connijk[iel][e_ext,iθ]
                            θ = extra_meshes_coords[iel][ip_ext]
                            ip_g = connijk_spa[iel][i,j,e_ext,iθ]
                            if (ip_g <= n_spa - n_non_global_nodes)
                                if (ip in mesh.poin_in_bdy_edge)
                                    applied = false
                                    iedge = elem_to_edge[iel,i,j,1]
                                    edge_i = elem_to_edge[iel,i,j,2]
                                    matchx = (x == mesh.xmax || x == mesh.xmin)
                                    matchy = (y == mesh.ymax || y == mesh.ymin)
                                    if (matchx && matchy) ##Do corner boundary case
                                        iedge1 = 1
                                        found = false
                                        iedge_found = 0
                                        edge_found_i = 0
                                        while (iedge1 <= mesh.nedges_bdy && found == false)
                                            for iter = 1:ngl
                                                ip1 = mesh.poin_in_bdy_edge[iedge1,iter]
                                                if (ip1 == ip && iedge != iedge1)
                                                    found = true
                                                    iedge_found = iedge1
                                                    edge_found_i = iter
                                                end
                                            end
                                            iedge1 +=1
                                        end
                                        if !(cos(θ) * (nx[iedge,edge_i]+nx[iedge_found,edge_found_i]) + sin(θ)*(ny[iedge,edge_i]+ny[iedge_found,edge_found_i]) <= 0.0)
                                            U_red_proj[ip_g] = RHS_red[ip_g]#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                        end
                                    elseif !(cos(θ) * nx[iedge,edge_i] + sin(θ)*ny[iedge,edge_i] <= 0.0)
                                        U_red_proj[ip_g] = RHS_red[ip_g]#exp(-((48/(2*π))*(θ-7*π/4))^2)#uip
                                    end
                                else
                                    U_red_proj[ip_g] = RHS_red[ip_g]
                                end
                            end
                        end
                    end
                end
            end
        end
        B = U_red_proj
    end
    @info "RHS max/min"
    @info maximum(B), minimum(B)
    #@info maximum(U_red_proj), minimum(U_red_proj)

    #x = real.(eigvals(Array(A)))
    #y = imag.(eigvals(Array(A)))
    #display(Makie.scatter(x, y, label="e-values"))

    #A_inv = inv(A)
    @info "built RHS"
    #@info RHS
    @info "solving system"
    As = sparse(A)
    A = nothing
    GC.gc()
    @info typeof(As)
    @info size(As), size(B)

    #@time solution = As \ B#RHS
    #=@time solution, stats = Krylov.fgmres(As, B;
                   atol = 1e-13,
                   rtol = 1e-13,
                   #btol = 1e-13,
                   #etol = 1e-13,
                   #axtol = 1e-13,
                   itmax = n_spa,
                   verbose = 1)=#
    @time solution = solve_parallel_lsqr(ip2gip_extra, gip2owner_extra, As, B, gnpoin, npoin_ang_total, pM; 
    npoin_g = npoin_ang_total)
   
    @info maximum(solution), minimum(solution)
    @info "done radiation solved"
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
        ### Do scatter operation before outputting solution
        for iel=1:nelem
            for i=1:ngl
                for j =1:ngl
                    ip = mesh.connijk[iel,i,j]
                    for e_ext = 1:extra_meshes_extra_nelems[iel]
                        for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                            ip_ext = extra_meshes_connijk[iel][e_ext,iθ]
                            ip_g = connijk_spa[iel][i,j,e_ext,iθ]
                            if (ip_g > n_spa-n_non_global_nodes)
                                #scatter here
                                for nc_i in nzrange(nc_mat,ip_g)
                                    count = 0
                                    for counter = 1:size(nc_non_global_nodes,1)
                                        if (nc_rows[nc_i] > nc_non_global_nodes[counter])
                                            count +=1
                                        end
                                    end
                                    #connijk_spa[iel][i,j,e_ext,iθ] = ip_g-n_spa
                                    solution_new[ip_g] += (nc_mat[nc_rows[nc_i],ip_g]/nc_mat_div[nc_rows[nc_i], ip_g])*solution[nc_rows[nc_i]]
                                end
                            else
                                #For global nodes adjust numbering
                                count = 0
                                for counter = 1:size(nc_non_global_nodes,1)
                                    if (ip_g > nc_non_global_nodes[counter])
                                        count +=1
                                    end
                                end
                                solution_new[ip_g] = solution[ip_g]
                            end
                        end
                    end
                end
            end
        end
    end
    
    int_sol      = zeros(TFloat, npoin, 1)
    int_ref      = zeros(TFloat, npoin, 1)
    int_sol_accum = zeros(npoin)
    int_ref_accum = zeros(npoin)
    L2_err = 0.0
    L2_ref = 0.0

    node_div = ones(Int, npoin)
    
        for iel = 1:nelem
            for j = 1:ngl, i = 1:ngl
                ip    = mesh.connijk[iel, i, j]
                x_ip  = mesh.x[ip]; y_ip = mesh.y[ip]
                matchx  = abs(x_ip - mesh.xmin) < 1e-10 || abs(x_ip - mesh.xmax) < 1e-10
                matchy  = abs(y_ip - mesh.ymin) < 1e-10 || abs(y_ip - mesh.ymax) < 1e-10
                nmatches = matchx + matchy
                on_bdy   = ip in mesh.poin_in_bdy_edge

                is_corner = (i ∈ (1,ngl)) && (j ∈ (1,ngl))
                is_edge  = !is_corner &&
                            (i ∈ (1,ngl) || j ∈ (1,ngl))

                div = 1
                if is_corner
                    if !on_bdy; div = 4
                    elseif nmatches == 1; div = 2
                    elseif nmatches == 2; div = 1
                    end
                elseif is_edge
                    if !on_bdy; div = 2
                    elseif nmatches == 1; div = 1
                    end
                end
                node_div[ip] = max(node_div[ip], div)
            end
        end


    for iel=1:nelem
        for i=1:ngl
            for j =1:ngl
                if (inputs[:adaptive_extra_meshes])
                    ip = mesh.connijk[iel,i,j]
                    
                    for e_ext = 1:extra_meshes_extra_nelems[iel]
                        for iθ = 1:extra_meshes_extra_nops[iel][e_ext]+1
                            
                            ip_ext = extra_meshes_connijk[iel][e_ext,iθ]
                            θ = extra_meshes_coords[iel][ip_ext]
                            ip_g = connijk_spa[iel][i,j,e_ext,iθ]
                            int_sol_accum[ip] += solution_new[ip_g]*extra_meshes_extra_Je[iel][e_ext,iθ]*extra_mesh[1].ωθ[iθ] / node_div[ip]
                            #if (iel == 11) @info int_sol[ip], solution_new[ip_g], solution_new[ip_g]*extra_meshes_extra_Je[iel][e_ext,iθ]*extra_mesh[1].ωθ[iθ]/div, e_ext, iθ, ip, x, y end
                            if (inputs[:lmanufactured_solution])
                                int_ref_accum[ip] += (ref[ip_g])*extra_meshes_extra_Je[iel][e_ext,iθ]*extra_mesh[1].ωθ[iθ] / node_div[ip]
                                L2_ref += (ref[ip_g])^2*extra_meshes_extra_Je[iel][e_ext,iθ]*extra_mesh[1].ωθ[iθ]*ω[i]*ω[j]*Je[iel,i,j]
                                L2_err += (ref[ip_g]-solution_new[ip_g])^2*extra_meshes_extra_Je[iel][e_ext,iθ]*extra_mesh[1].ωθ[iθ]*ω[i]*ω[j]*Je[iel,i,j]
                            end
                        end
                    end
                else
                    ip = mesh.connijk[iel,i,j]
                    
                    for e_ext = 1:extra_mesh.extra_nelem
                        for iθ = 1:extra_mesh.extra_nop[e_ext]+1
                            
                            ip_ext = extra_mesh.extra_connijk[e_ext,iθ]
                            θ = extra_mesh.extra_coords[ip_ext]
                            ip_g = (ip-1) * extra_mesh.extra_npoin + ip_ext
                            
                            int_sol_accum[ip] += solution[ip_g]*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ] / node_div[ip]
                            if (inputs[:lmanufactured_solution])
                                int_ref_accum[ip] += (ref[ip_g])*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ] / node_div[ip]
                                L2_ref += (ref[ip_g])^2*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ]*ω[i]*ω[j]*Je[iel,i,j]#/div1
                                L2_err += (ref[ip_g]-solution[ip_g])^2*extra_mesh.extra_metrics.Je[e_ext,iθ]*extra_mesh.ωθ[iθ]*ω[i]*ω[j]*Je[iel,i,j]#/div1
                            end
                        end
                    end
                end
            end
        end
    end


    # ── Global reduction of integrals and norms ───────────────────────────────
    
        gnpoin_spa  = mesh.gnpoin
        g_int_sol   = zeros(Float64, gnpoin_spa)
        g_int_ref   = zeros(Float64, gnpoin_spa)

        for ip = 1:npoin
            if int_sol_accum[ip] != 0.0 || int_ref_accum[ip] != 0.0
                gip = mesh.ip2gip[ip]
                g_int_sol[gip] += int_sol_accum[ip]
                g_int_ref[gip] += int_ref_accum[ip]
            end
        end
        MPI.Allreduce!(g_int_sol, MPI.SUM, comm)
        MPI.Allreduce!(g_int_ref, MPI.SUM, comm)
        for ip = 1:npoin
            gip = mesh.ip2gip[ip]
            int_sol[ip] = g_int_sol[gip]
            int_ref[ip] = g_int_ref[gip]
        end
        if (inputs[:lmanufactured_solution])
            L2_ref_g = MPI.Allreduce(L2_ref, MPI.SUM, comm)
            L2_err_g = MPI.Allreduce(L2_err, MPI.SUM, comm)
        end


    if (inputs[:lmanufactured_solution])
        @info "new L2 norms", sqrt(L2_ref), sqrt(L2_err), sqrt(L2_err/L2_ref)
        @info "infinity norms", maximum(abs.(solution-ref)), maximum(abs.(solution-ref))/maximum(ref)
        
    end
    title = @sprintf "Solution-Radiation"
    write_vtk(SD, mesh, int_sol, int_sol, nothing, nothing, nothing,
              0.0, 0.0, 0.0, 0.0, title, inputs[:output_dir], inputs,
              ["Ang_int"], ["Ang_int"]; iout=1, nvar=1)
end

function sparse_lhs_assembly_2Dby1D(ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
                                   dξdx, dξdy, dηdx, dηdy, npoin_ang, rad_HG_g)

    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries*0.0009)))
    sizehint!(J_vec, Int64(round(max_entries*0.0009)))
    sizehint!(V_vec, Int64(round(max_entries*0.0009)))
    HG, error = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)), 0, 2*π, rtol=1e-13, atol = 1e-13)    
    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                dξdx_ij = dξdx[iel,i,j]
                dξdy_ij = dξdy[iel,i,j]
                dηdx_ij = dηdx[iel,i,j]
                dηdy_ij = dηdy[iel,i,j]
                κ = user_extinction(x[ip],y[ip])
                σ = user_scattering_coef(x[ip],y[ip])
                for e_ext = 1:nelem_ang
                    for iθ = 1:nop_ang[e_ext]+1
                        ip_ext = connijk_ang[e_ext,iθ]
                        ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                        
                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]
                                extinction = κ*ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]#*ψ_ang[iθ,jθ]

                                θ = coords_ang[ip_ext]
                                propagation = (ωJac*ωJac_rad)*(ψ[n,j]*dψ[m,i]*dξdx_ij*cos(θ) + ψ[m,i]*dψ[n,j]*dηdy_ij*sin(θ))
                                intϕ = 0.0 
                                for e_ext_scatter = 1:nelem_ang
                                    for kθ = 1:nop_ang[e_ext]+1
                                        div = 1
                                        #if (kθ == nop_ang[e_ext]+1 || kθ == 1)
                                        #    div =2
                                        #end
                                        ipθ = connijk_ang[e_ext_scatter,kθ]
                                        θ1 = coords_ang[ipθ]
                                        Φ = user_scattering_functions(θ,θ1,HG)
                                        ωJac_rad_scatter = ωθ[kθ]*Je_ang[e_ext_scatter,kθ]
                                        intϕ +=   ωJac_rad_scatter*Φ/div
                                    end
                                end
                                scattering = intϕ * ψ[i,m] * ψ[j,n] * ωJac*ωJac_rad*σ               
                                val = extinction + propagation - scattering
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
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_lhs_assembly_2Dby1D_adaptive(ref_level, ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang, connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang,
                                   dξdx, dξdy, dηdx, dηdy, npoin_ang, rad_HG_g, connijk_spa, nc_mat, nc_mat_div, adapted, nc_non_global_nodes, n_non_global_nodes, n_spa)
    @info adapted
    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    max_col = size(nc_mat,2)
    nc_rows = rowvals(nc_mat)
    sizehint!(I_vec, Int64(round(max_entries*0.0009)))
    sizehint!(J_vec, Int64(round(max_entries*0.0009)))
    sizehint!(V_vec, Int64(round(max_entries*0.0009)))
    HG, error = quadgk(v -> (1-rad_HG_g^2)/((1+rad_HG_g^2-2*rad_HG_g*cos(v))^(3/2)), 0, 2*π, rtol=1e-13, atol = 1e-13)
    temp_sum = 0
    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                dξdx_ij = dξdx[iel,i,j]
                dξdy_ij = dξdy[iel,i,j]
                dηdx_ij = dηdx[iel,i,j]
                dηdy_ij = dηdy[iel,i,j]
                κ = user_extinction(x[ip],y[ip])
                σ = user_scattering_coef(x[ip],y[ip])
                for e_ext = 1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel][e_ext]+1
                        ip_ext = connijk_ang[iel][e_ext,iθ]
                        ωJac_rad = ωθ[iθ]*Je_ang[iel][e_ext,iθ]

                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]
                                extinction = κ*ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]#*ψ_ang[iθ,jθ]

                                θ = coords_ang[iel][ip_ext]
                                propagation = (ωJac*ωJac_rad)*(ψ[n,j]*dψ[m,i]*dξdx_ij*cos(θ) + ψ[m,i]*dψ[n,j]*dηdy_ij*sin(θ))
                                intϕ = 0.0
                                for e_ext_scatter = 1:nelem_ang[iel]
                                    for kθ = 1:nop_ang[iel][e_ext]+1
                                        div = 1
                                        if (kθ == nop_ang[iel][e_ext]+1 || kθ == 1)
                                            div =2
                                        end
                                        ipθ = connijk_ang[iel][e_ext_scatter,kθ]
                                        θ1 = coords_ang[iel][ipθ]
                                        Φ = user_scattering_functions(θ,θ1,HG)
                                        ωJac_rad_scatter = ωθ[kθ]*Je_ang[iel][e_ext_scatter,kθ]
                                        intϕ +=   ωJac_rad_scatter*Φ/div
                                    end
                                end
                                scattering = intϕ * ψ[i,m] * ψ[j,n] * ωJac*ωJac_rad*σ
                                val = extinction + propagation - scattering
                                idx_ip = connijk_spa[iel][i,j,e_ext,iθ]
                                idx_jp = connijk_spa[iel][m,n,e_ext,iθ]
                                if (abs(val) > eps(Float64)) #&& (idx_ip <= n_spa) && (idx_jp <= n_spa)   # Skip near-zero entries, make sure node is global
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
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_mass_assembly_2Dby1D(ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang, 
        connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang)
  
    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()

    sizehint!(I_vec, Int64(round(max_entries*0.0001)))
    sizehint!(J_vec, Int64(round(max_entries*0.0001)))
    sizehint!(V_vec, Int64(round(max_entries*0.0001)))

    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                for e_ext = 1:nelem_ang
                    for iθ = 1:nop_ang[e_ext]+1
                        ωJac_rad = ωθ[iθ]*Je_ang[e_ext,iθ]
                        ip_ext = connijk_ang[e_ext,iθ]
                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]

                                val = ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]#*ψ_ang[iθ,jθ]
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
    return sparse(I_vec, J_vec, V_vec)
end

function sparse_mass_assembly_2Dby1D_adaptive(ref_level, ω, Je, connijk, ωθ, x, y, ψ, dψ, ψ_ang,
        connijk_ang, Je_ang, coords_ang, nop_ang, npoin_ang_total, nelem, ngl, nelem_ang, npoin_ang, connijk_spa, nc_mat, nc_mat_div, adapted, nc_non_global_nodes, n_non_global_nodes, n_spa)

    max_entries = npoin_ang_total^2
    I_vec = Vector{Int}()
    J_vec = Vector{Int}()
    V_vec = Vector{Float64}()
    max_col = size(nc_mat,2)
    nc_rows = rowvals(nc_mat)
    sizehint!(I_vec, Int64(round(max_entries*0.0001)))
    sizehint!(J_vec, Int64(round(max_entries*0.0001)))
    sizehint!(V_vec, Int64(round(max_entries*0.0001)))

    for iel=1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                ωJac = ω[i]*ω[j]*Je[iel,i,j]
                for e_ext = 1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel][e_ext]+1
                        ωJac_rad = ωθ[iθ]*Je_ang[iel][e_ext,iθ]
                        ip_ext = connijk_ang[iel][e_ext,iθ]
                        for n=1:ngl
                            for m=1:ngl
                                jp = connijk[iel,m,n]

                                val = ωJac*ωJac_rad*ψ[j,n]*ψ[i,m]#*ψ_ang[iθ,jθ]
                                idx_ip = connijk_spa[iel][i,j,e_ext,iθ]
                                idx_jp = connijk_spa[iel][m,n,e_ext,iθ]
                                if (abs(val) > eps(Float64)) #&& (idx_ip <= n_spa) && (idx_jp <= n_spa)   # Skip near-zero entries, make sure node is global
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
    return sparse(I_vec, J_vec, V_vec)
end

function compute_adaptivity_criterion(pointwise_interaction, nelem, ngl, connijk, connijk_ang, nop_ang, nelem_ang, coords_ang, connijk_spa)
    criterion = [Vector{Float64}(undef, nelem_ang[e]) for e=1:nelem]
    for iel=1:nelem 
        for j=1:ngl
            for i=1:ngl 
                ip = connijk[iel,i,j]
                for e_ext = 1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel][e_ext]+1

                        ip_ext = connijk_ang[iel][e_ext,iθ]
                        idx_ip = connijk_spa[iel][i,j,e_ext,iθ]
                        θ =  coords_ang[iel][ip_ext]
                        e_ext1 = 0
                        e_ext2 = 0
                        iθ1 = 0
                        iθ2 = 0
                        if(iθ > 1 && iθ < nop_ang[iel][e_ext]+1)
                            e_ext1 = e_ext
                            e_ext2 = e_ext
                            iθ1 = iθ-1
                            iθ2 = iθ+1
                        elseif (iθ == 1)
                            if(e_ext == 1)
                                e_ext1 = nelem_ang[iel]
                                e_ext2 = e_ext
                                iθ1 = nop_ang[iel][e_ext1]+1
                                iθ2 = iθ+1
                            else
                                e_ext1 = e_ext-1
                                e_ext2 = e_ext
                                iθ1 = nop_ang[iel][e_ext1]+1
                                iθ2 = iθ+1
                            end
                        elseif (iθ == nop_ang[iel][e_ext]+1)
                            if(e_ext == nelem_ang[iel])
                                e_ext1 = e_ext
                                e_ext2 = 1
                                iθ1 = iθ-1
                                iθ2 = 1
                            else
                                e_ext1 = e_ext
                                e_ext2 = e_ext+1
                                iθ1 = iθ-1
                                iθ2 = 1
                            end
                        end
                        ip_ext1 = connijk_ang[iel][e_ext1,iθ1]
                        ip_ext2 = connijk_ang[iel][e_ext2,iθ2]
                        idx_ip1 = connijk_spa[iel][i,j,e_ext1,iθ1]
                        idx_ip2 = connijk_spa[iel][i,j,e_ext2,iθ2]
                        θ1 = coords_ang[iel][ip_ext1]
                        θ2 = coords_ang[iel][ip_ext2]
                        Δθ = θ2-θ1
                        if (Δθ < 0)
                            Δθ += 2*π
                        end
                        criterion[iel][e_ext] += (pointwise_interaction[idx_ip2]-pointwise_interaction[idx_ip1])/Δθ
                    end
                    criterion[iel][e_ext] = criterion[iel][e_ext]/(nop_ang[iel][e_ext]+1)
                end

            end
        end
    end
    return criterion
end

function adapt_angular_grid_2Dby1D!(criterion,thresholds,ref_level,nelem,ngl,nelem_ang, nop_ang, neighbors, npoin_ang,
        connijk_ang, coords_ang, Je_ang, dξdx_ang, dxdξ_ang, connijk, x, y, xmin_grid, ymin_grid, xmax_grid, ymax_grid)
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
            if (abs(criterion[iel][e_ext]) > thresholds[1])
                #@info e_ext, "to be adapted"
                adapted_ang[iel] = 1
                ref_level[iel][e_ext] += 1
                ang_adapted[iel] = 1
                # Make new angular elements and reconstruct extra_mesh arrays
                ang_adapt = true
                nelem_ang[iel] += 1
                npoin_ang[iel] += nop_ang[iel][e_ext]
                θmax = coords_ang[iel][connijk_ang[iel][e_ext,nop_ang[iel][e_ext]+1]]
                θmin = coords_ang[iel][connijk_ang[iel][e_ext,1]]
                if (θmax == 0)
                    θmax = 2*π
                end
                θ12 = (θmax + θmin)/2
                connijk_ang_new = zeros(Int, nelem_ang[iel], nop_ang[iel][1]+1)
                coords_new = zeros(Float64, npoin_ang[iel]) 
                metrics = allocate_metrics(NSD_1D(), nelem_ang[iel], 0, nop_ang[iel][e_ext], TFloat, CPU())
                nop_ang_new = zeros(Int,nelem_ang[iel])
                nop_ang_new[1:nelem_ang[iel]-1] .= nop_ang[iel][e_ext]
                nop_ang_new[nelem_ang[iel]] = nop_ang[iel][1]
                iter = 1
                criterion_new = zeros(Float64,nelem_ang[iel])
                ref_level_new = zeros(Int,nelem_ang[iel])
                #populate the elements coming before
                if (e_ext > 1)
                    for e_ext1=1:e_ext-1
                        for i=1:nop_ang[iel][e_ext1]+1
                            connijk_ang_new[e_ext1,i] = iter
                            coords_new[iter] = coords_ang[iel][connijk_ang[iel][e_ext1,i]]
                            metrics.dxdξ[e_ext1, i, 1]  = dxdξ_ang[iel][e_ext1, i, 1]
                            metrics.Je[e_ext1, i, 1]  = Je_ang[iel][e_ext1, i, 1]
                            metrics.dξdx[e_ext1, i, 1]  = dξdx_ang[iel][e_ext1, i, 1]
                            criterion_new[e_ext] = criterion[iel][e_ext1]
                            ref_level_new[e_ext] = ref_level[iel][e_ext1]
                            if (i != nop_ang[iel][e_ext1]+1) iter +=1 end
                        end
                    end
                end
                #populate for the new elements
                lgl = basis_structs_ξ_ω!(LGL(), nop_ang[iel][1], CPU())
                for i=1:nop_ang[iel][e_ext]+1
                    ξ = lgl.ξ[i]
                    coords_new[iter] = θmin*(1.0-ξ)*0.5+θ12*(1.0 + ξ)*0.5
                    connijk_ang_new[e_ext,i]    = iter
                    metrics.dxdξ[e_ext, i, 1]   = (θ12-θmin)/2
                    metrics.Je[e_ext, i, 1]     = metrics.dxdξ[e_ext, i, 1]
                    metrics.dξdx[e_ext, i, 1]  = 1.0/metrics.Je[e_ext, i, 1]
                    criterion_new[e_ext] = 0.0#criterion[iel][e_ext]
                    ref_level_new[e_ext] = ref_level[iel][e_ext]
                    if (i != nop_ang[iel][e_ext]+1) iter +=1 end
                end
                for i=1:nop_ang_new[e_ext+1]+1
                    ξ = lgl.ξ[i]
                    if (θmax == 2*π) && (i == nop_ang_new[e_ext+1]+1)
                        iter = 1
                        coords_new[iter] = 0.0
                        connijk_ang_new[e_ext+1,i] = iter
                    else
                        coords_new[iter] = θ12*(1.0-ξ)*0.5+θmax*(1.0 + ξ)*0.5
                        connijk_ang_new[e_ext+1,i] = iter
                    end
                    metrics.dxdξ[e_ext+1, i, 1]   = (θmax-θ12)/2
                    metrics.Je[e_ext+1, i, 1]     = metrics.dxdξ[e_ext+1, i, 1]
                    metrics.dξdx[e_ext+1, i, 1]  = 1.0/metrics.Je[e_ext+1, i, 1]
                    criterion_new[e_ext+1] = 0.0#criterion[iel][e_ext]
                    ref_level_new[e_ext+1] = ref_level[iel][e_ext]
                    if (i != nop_ang_new[e_ext+1]+1) iter +=1 end
                end
                if (e_ext < nelem_ang[iel]-1)
                    for e_ext1=e_ext+2:nelem_ang[iel]
                        for i=1:nop_ang_new[e_ext1]+1
                            if (i == nop_ang_new[e_ext1]+1) && (e_ext1 == nelem_ang[iel])
                                iter = 1
                            end
                            connijk_ang_new[e_ext1,i] = iter
                            coords_new[iter] = coords_ang[iel][connijk_ang[iel][e_ext1-1,i]]
                            metrics.dxdξ[e_ext1, i, 1]  = dxdξ_ang[iel][e_ext1-1, i, 1]
                            metrics.Je[e_ext1, i, 1]  = Je_ang[iel][e_ext1-1, i, 1]
                            metrics.dξdx[e_ext1, i, 1]  = dξdx_ang[iel][e_ext1-1, i, 1]
                            criterion_new[e_ext1] = criterion[iel][e_ext1-1]
                            ref_level_new[e_ext1] = ref_level[iel][e_ext1-1]
                            if (i != nop_ang_new[e_ext1]+1) iter +=1 end
                        end 
                    end
                end

                #=for e_ext1 = 1:nelem_ang[iel]
                    for i = 1:nop_ang[iel][1]+1
                        ip_ext1 = connijk_ang_new[e_ext1,1]
                        ip_ext2 = connijk_ang_new[e_ext1,nop_ang[iel][1]+1]
                        Δθe = abs(coords_new[ip_ext2]-coords_new[ip_ext1])
                        for k = 1:nop_ang[iel][1]+1
                            metrics.dxdξ[e_ext1, k, 1]  = Δθe/2
                            metrics.Je[e_ext1, k, 1]   = metrics.dxdξ[e_ext1, k, 1]
                            metrics.dξdx[e_ext1, k, 1] = 1.0/metrics.Je[e_ext1, k, 1]
                        end
                    end
                end=#
                connijk_ang[iel] = connijk_ang_new
                coords_ang[iel] = coords_new
                dxdξ_ang[iel] = metrics.dxdξ[:,:,1]
                Je_ang[iel] = metrics.Je[:,:,1]
                dξdx_ang[iel] = metrics.dξdx[:,:,1]
                nop_ang[iel] = nop_ang_new
                criterion[iel] = criterion_new
                ref_level[iel] = ref_level_new
            end
            e_ext += 1
        end
    end 
    for iel = 1:nelem
        # Find and save spatial neighbords for non-conforming assembly
            #First find element corners
            xmin = 10^10
            xmax = -10^10
            ymin = 10^10
            ymax = -10^10
        
            for i=1:ngl
                for j=1:ngl
                    ip = connijk[iel,i,j]
                    if (x[ip] < xmin) xmin = x[ip] end
                    if (y[ip] < ymin) ymin = y[ip] end
                    if (x[ip] > xmax) xmax = x[ip] end
                    if (y[ip] > ymax) ymax = y[ip] end
                end
            end
            match_bdy = 0
            if (xmin == xmin_grid) match_bdy += 1 end
            if (xmax == xmax_grid) match_bdy += 1 end
            if (ymin == ymin_grid) match_bdy += 1 end
            if (ymax == ymax_grid) match_bdy += 1 end

            iter = 1
            found_neighbors = 0
            while (iter <= nelem && found_neighbors <8)
                #find corners for comparison
                xmin_i = 10^10
                xmax_i = -10^10
                ymin_i = 10^10
                ymax_i = -10^10
                for i=1:ngl
                    for j=1:ngl
                        ip = connijk[iter,i,j]
                        if (x[ip] < xmin_i) xmin_i = x[ip] end
                        if (y[ip] < ymin_i) ymin_i = y[ip] end
                        if (x[ip] > xmax_i) xmax_i = x[ip] end
                        if (y[ip] > ymax_i) ymax_i = y[ip] end
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
                if (match_neighbor1 > 0 && match_neighbor2 > 0 ) || match_neighbor2 > 1
                    found_neighbors += 1
                    neighbors[iel,found_neighbors,1] = iter
                    #check for conformity here
                    if !(adapted_ang[iel] == 0 && adapted_ang[iter] == 0) #if no angular refinement has taken place no need to check conformity
                        if (nelem_ang[iel] != nelem_ang[iter] || coords_ang[iel] != coords_ang[iter]) #non conforming
                            neighbors[iel,found_neighbors,2] = 1 # Save information that these neighbors are non conforming
                            #@info neighbors[iel,found_neighbors,2], iel, iter 
                        end
                    end
                end
                if (match_bdy == 1 && found_neighbors == 5)
                    found_neighbors = 8
                elseif (match_bdy == 2 && found_neighbors == 3)
                        found_neighbors = 8
                end
                iter += 1
            end
    end
end

function adaptive_spatial_angular_numbering_2D_1D!(connijk_spa,nelem, ngl, connijk, connijk_ang, nop_ang, nelem_ang, coords_ang, x, y,ref_level, neighbors, adapted)
    points = []
    points_c = Array{Float64}(undef,3)
    iter = 1
    interp_sources = zeros(Float64,nop_ang[1][1]+1)
    interp_targets = zeros(Float64,nop_ang[1][1]+1)
    ω = zeros(Float64,nop_ang[1][1]+1)
    L = zeros(Float64,nop_ang[1][1]+1,nop_ang[1][1]+1)
     
    nc_non_global_nodes = []
    n_non_global_nodes = 0
    for iel = 1:nelem
        for j=1:ngl
            for i=1:ngl
                ip = connijk[iel,i,j]
                x_p = x[ip]
                y_p = y[ip]
                for e_ext=1:nelem_ang[iel]
                    for iθ = 1:nop_ang[iel][e_ext]+1
                        ip_ang = connijk_ang[iel][e_ext,iθ]
                        θ_p = coords_ang[iel][ip_ang]
                        if (vec([x_p, y_p, θ_p]) in eachrow(points_c'))
                            found = false
                            iter1 = 1
                            while (found == false && iter1 < iter)
                                if (points_c[1,iter1] == x_p) && (points_c[2,iter1] == y_p) && (points_c[3,iter1] == θ_p)
                                   connijk_spa[iel][i,j,e_ext,iθ] = iter-iter1  
                                   #@info connijk_spa[iel][i,j,ip_ang], x_p, y_p, θ_p 
                                   found = true
                               else
                                    iter1 +=1
                                end
                            end
                        else
                            connijk_spa[iel][i,j,e_ext,iθ] = iter
                            #@info connijk_spa[iel][i,j,ip_ang], x_p, y_p, θ_p
                            points_c = hcat([x_p, y_p, θ_p],points_c)
                            push!(points,iter)
                            iter += 1
                        end

                    end
                end
            end
        end
    end
    @info "finished non-adaptive connectivity"
    @info "total number of independent points", iter-1
    n_spa = iter-1
    max_entries = iter^2
    In_vec = Vector{Int}()
    Jn_vec = Vector{Int}()
    Vn_vec = Vector{Float64}()
    sizehint!(In_vec, Int64(round(max_entries*0.0001)))
    sizehint!(Jn_vec, Int64(round(max_entries*0.0001)))
    sizehint!(Vn_vec, Int64(round(max_entries*0.0001)))
    rep_ip = zeros(Int,ngl)
    iter_nc = 1
    ### Done with conforming connectivity
    #handle non-conformity
    for iel = 1:nelem
            if (1 in neighbors[iel,:,2]) #element has non-conforming neighbors
                for ineighbor = 1:8
                    if (neighbors[iel,ineighbor,2] == 1)
                        adapted = true
                        iel1 = neighbors[iel,ineighbor,1] # identify neighbor spatial element number
                        i = 0
                        j = 0
                        i1 = 0
                        j1 = 0
                        ip1 = 0
                        ip = 0
                        #loop through element edges to find matching spatial nodes
                        fill!(rep_ip,zero(Int))
                        for igl=1:ngl
                            match = false
                            if (connijk[iel,1,igl] in connijk[iel1,:,:]) && !(connijk[iel,1,igl] in rep_ip) # a matching spatial node exists
                                ip = connijk[iel,1,igl]
                                i = 1
                                j = igl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            elseif (connijk[iel,igl,1] in connijk[iel1,:,:]) && !(connijk[iel,igl,1] in rep_ip) # a matching spatial node exists
                                ip = connijk[iel,igl,1]
                                i = igl
                                j = 1
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            elseif (connijk[iel,ngl,igl] in connijk[iel1,:,:]) && !(connijk[iel,ngl,igl] in rep_ip) # a matching spatial node exists
                                ip = connijk[iel,ngl,igl]
                                i = ngl
                                j = igl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            elseif (connijk[iel,igl,ngl] in connijk[iel1,:,:]) && !(connijk[iel,igl,ngl] in rep_ip)
                                ip = connijk[iel,igl,ngl]
                                i = igl
                                j = ngl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            end
                            #@info "matching neighbor nodes", iel, iel1, match, ip, ip1, rep_ip 
                            if (match == true)
                                rep_ip[igl] = ip
                                #@info "matching nodes", ip, ip1, i, j, i1, j1, iel, iel1 

                                # matched spatial node found
                                # matching spatial nodes must communicate non-conforming angular nodes
                                # loop through angular elements
                            
                                for e_ext=1:nelem_ang[iel]
                                    #angular elements are ordered
                                    #check for an exact angular element match on neighboring element
                                    e_check = 1
                                    found = false
                                    while (found == false && e_check <= nelem_ang[iel1])
                                        matched = true
                                        iter = 1
                                        while (matched && iter <= nop_ang[iel][e_ext]+1)
                                            ip_ang = connijk_ang[iel][e_ext,iter]
                                            ip_ang1 = connijk_ang[iel1][e_check,iter]
                                            if !(coords_ang[iel][ip_ang] == coords_ang[iel1][ip_ang1])
                                                matched = false
                                            end
                                            iter +=1
                                        end
                                        if (matched == true)
                                            found = true
                                        end
                                        e_check += 1
                                    end
                                    # There is no need for nc-treatment for angular elements that are exact matches
                                    if (found == false) # no exact angular element matches were found for this element, it's necessary to find NC-DSS target elements
                                        #find this element's end nodes
                                        θmin = coords_ang[iel][connijk_ang[iel][e_ext,1]]
                                        θmax = coords_ang[iel][connijk_ang[iel][e_ext,nop_ang[iel][e_ext]+1]]
                                        if (θmax == 0.0) θmax = 2*π end
                                        #look for target elements
                                        for e_ext1 = 1:nelem_ang[iel1]
                                            θmin1 = coords_ang[iel1][connijk_ang[iel1][e_ext1,1]]
                                            θmax1 = coords_ang[iel1][connijk_ang[iel1][e_ext1,nop_ang[iel1][e_ext1]+1]]
                                            if (θmax1 == 0.0) θmax1 = 2*π end
                                            if (θmin1 >= θmin && θmax1 <= θmax)
                                                #First situation e_ext the "parent" element and e_ext1 is the "child" element
                                                #Find interpolating points
                                                #All parent nodes send information to all child nodes
                                                for iθ=1:nop_ang[iel][e_ext]+1
                                                    interp_sources[iθ] = coords_ang[iel][connijk_ang[iel][e_ext,iθ]]
                                                    if (θmax == 2*π && iθ == nop_ang[iel][e_ext]+1) interp_sources[iθ] = θmax end
                                                end
                                                #Find target points
                                                resize!(interp_targets,nop_ang[iel1][e_ext1]+1)
                                                L = zeros(Float64, nop_ang[iel][e_ext]+1, nop_ang[iel1][e_ext1]+1)
                                                for iθ=1:nop_ang[iel1][e_ext1]+1
                                                    interp_targets[iθ] = coords_ang[iel1][connijk_ang[iel1][e_ext1,iθ]]
                                                    if (θmax1 == 2*π && iθ == nop_ang[iel][e_ext]+1) interp_targets[iθ] = θmax1 end
                                                end
                                                #contruct lagrange interpolator
                                                #find barycentric weights
                                                BarycentricWeights!(interp_sources,ω)
                                                #build interpolation matrix
                                                PolynomialInterpolationMatrix!(interp_sources,ω,interp_targets,L)
                                                #@info L, interp_sources, interp_targets, "parent to child"
                                                #Store data for assembly
                                                for jθ=1:nop_ang[iel1][e_ext1]+1
                                                    if !(1 in L[jθ,:])
                                                        jp_spa = connijk_spa[iel1][i1,j1,e_ext1,jθ]
                                                        if !(jp_spa in nc_non_global_nodes)
                                                            push!(nc_non_global_nodes,jp_spa)
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
                    end
                end
            end
    end
    #Redo numbering to put N_C nodes at end order-wise.
    if (n_non_global_nodes > 1)
        for iel = 1:nelem
            for j=1:ngl
                for i=1:ngl
                    for e_ext=1:nelem_ang[iel]
                        for iθ = 1:nop_ang[iel][e_ext]+1 
                            ip_g = connijk_spa[iel][i,j,e_ext,iθ]
                            if (ip_g <= n_spa) && !(ip_g in nc_non_global_nodes)
                                count = 0
                                for counter =1:size(nc_non_global_nodes,1)
                                    if (ip_g > nc_non_global_nodes[counter])
                                        count += 1
                                    end
                                end
                                connijk_spa[iel][i,j,e_ext,iθ] = ip_g - count
                            elseif (ip_g in nc_non_global_nodes)
                                count = 1
                                found = false
                                while (count <= size(nc_non_global_nodes,1) && found == false)
                                    if (ip_g == nc_non_global_nodes[count])
                                        connijk_spa[iel][i,j,e_ext,iθ] = n_spa + count
                                        found = true
                                    else
                                        count += 1
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    # adjust numbering for nc nodes at the end
    for iel = 1:nelem
        if (1 in neighbors[iel,:,2])
            for j=1:ngl
                for i=1:ngl
                    for e_ext=1:nelem_ang[iel]
                        for iθ = 1:nop_ang[iel][e_ext]+1
                            ip_g = connijk_spa[iel][i,j,e_ext,iθ]
                            if (ip_g > n_spa)
                                dist = ip_g - n_spa
                                connijk_spa[iel][i,j,e_ext,iθ] = n_spa - (dist-1)
                            end
                        end
                    end
                end
            end
        end
    end
    nc_mat = spzeros(Float64, n_spa-n_non_global_nodes, n_spa)
    for ip_g = 1:n_spa-n_non_global_nodes
        nc_mat[ip_g,ip_g] = 1.0
        push!(In_vec, ip_g)
        push!(Jn_vec, ip_g)
        push!(Vn_vec, 1)
    end
    for iel = 1:nelem
            if (1 in neighbors[iel,:,2]) #element has non-conforming neighbors
                for ineighbor = 1:8
                    if (neighbors[iel,ineighbor,2] == 1)
                        adapted = true
                        iel1 = neighbors[iel,ineighbor,1] # identify neighbor spatial element number
                        #@info "non-conforming neighbors", iel, iel1
                        i = 0
                        j = 0
                        i1 = 0
                        j1 = 0
                        ip1 = 0
                        ip = 0
                        #loop through element edges to find matching spatial nodes
                        fill!(rep_ip,zero(Int))
                        for igl=1:ngl
                            match = false
                            if (connijk[iel,1,igl] in connijk[iel1,:,:]) && !(connijk[iel,1,igl] in rep_ip) # a matching spatial node exists
                                ip = connijk[iel,1,igl]
                                i = 1
                                j = igl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            elseif (connijk[iel,igl,1] in connijk[iel1,:,:]) && !(connijk[iel,igl,1] in rep_ip) # a matching spatial node exists
                                ip = connijk[iel,igl,1]
                                i = igl
                                j = 1
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            elseif (connijk[iel,ngl,igl] in connijk[iel1,:,:]) && !(connijk[iel,ngl,igl] in rep_ip) # a matching spatial node exists
                                ip = connijk[iel,ngl,igl]
                                i = ngl
                                j = igl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            elseif (connijk[iel,igl,ngl] in connijk[iel1,:,:]) && !(connijk[iel,igl,ngl] in rep_ip)
                                ip = connijk[iel,igl,ngl]
                                i = igl
                                j = ngl
                                i1, j1, ip1 = find_edge_node_match(ngl,iel1,ip,connijk)
                                match = true
                            end
                            #@info "matching neighbor nodes", iel, iel1, match, ip, ip1, rep_ip 
                            if (match == true)
                                rep_ip[igl] = ip
                                #@info "matching nodes", ip, ip1, i, j, i1, j1, iel, iel1 

                                # matched spatial node found
                                # matching spatial nodes must communicate non-conforming angular nodes
                                # loop through angular elements

                                for e_ext=1:nelem_ang[iel]
                                    #angular elements are ordered
                                    #check for an exact angular element match on neighboring element
                                    e_check = 1
                                    found = false
                                    while (found == false && e_check <= nelem_ang[iel1])
                                        matched = true
                                        iter = 1
                                        while (matched && iter <= nop_ang[iel][e_ext]+1)
                                            ip_ang = connijk_ang[iel][e_ext,iter]
                                            ip_ang1 = connijk_ang[iel1][e_check,iter]
                                            if !(coords_ang[iel][ip_ang] == coords_ang[iel1][ip_ang1])
                                                matched = false
                                            end
                                            iter +=1
                                        end
                                        if (matched == true)
                                            found = true
                                        end
                                        e_check += 1
                                    end
                                    # There is no need for nc-treatment for angular elements that are exact matches
                                    if (found == false)
                                    #find this element's end nodes
                                        θmin = coords_ang[iel][connijk_ang[iel][e_ext,1]]
                                        θmax = coords_ang[iel][connijk_ang[iel][e_ext,nop_ang[iel][e_ext]+1]]
                                        if (θmax == 0.0) θmax = 2*π end
                                        #look for target elements
                                        for e_ext1 = 1:nelem_ang[iel1]
                                            θmin1 = coords_ang[iel1][connijk_ang[iel1][e_ext1,1]]
                                            θmax1 = coords_ang[iel1][connijk_ang[iel1][e_ext1,nop_ang[iel1][e_ext1]+1]]
                                            if (θmax1 == 0.0) θmax1 = 2*π end
                                            if (θmin1 >= θmin && θmax1 <= θmax)
                                                #First situation e_ext the "parent" element and e_ext1 is the "child" element
                                                #Find interpolating points
                                                #All parent nodes send information to all child nodes
                                                for iθ=1:nop_ang[iel][e_ext]+1
                                                    interp_sources[iθ] = coords_ang[iel][connijk_ang[iel][e_ext,iθ]]
                                                    if (θmax == 2*π && iθ == nop_ang[iel][e_ext]+1) interp_sources[iθ] = θmax end
                                                end
                                                #Find target points
                                                resize!(interp_targets,nop_ang[iel1][e_ext1]+1)
                                                L = zeros(Float64, nop_ang[iel][e_ext]+1, nop_ang[iel1][e_ext1]+1)
                                                for iθ=1:nop_ang[iel1][e_ext1]+1
                                                    interp_targets[iθ] = coords_ang[iel1][connijk_ang[iel1][e_ext1,iθ]]
                                                    if (θmax1 == 2*π && iθ == nop_ang[iel][e_ext]+1) interp_targets[iθ] = θmax1 end
                                                end
                                                #contruct lagrange interpolator
                                                #find barycentric weights
                                                BarycentricWeights!(interp_sources,ω)
                                                #build interpolation matrix
                                                PolynomialInterpolationMatrix!(interp_sources,ω,interp_targets,L)
                                                #@info L, interp_sources, interp_targets, "parent to child"
                                                #Store data for assembly
                                                    #=if (connijk_spa[iel1][i1,j1,e_ext1,jθ] <= n_spa)
                                                        connijk_spa[iel1][i1,j1,e_ext1,jθ] += n_spa 
                                                        iter_nc += 1
                                                    end=#
                                                for iθ=1:nop_ang[iel][e_ext]+1
                                                    ip_spa = connijk_spa[iel][i,j,e_ext,iθ]
                                                    for jθ=1:nop_ang[iel1][e_ext1]+1
                                                        jp_spa = connijk_spa[iel1][i1,j1,e_ext1,jθ]
                                                        if (L[jθ,iθ]!=0)
                                                            nc_mat[ip_spa,jp_spa] = L[jθ,iθ]
                                                            push!(In_vec, ip_spa)
                                                            push!(Jn_vec, jp_spa)
                                                            push!(Vn_vec, 1)
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
            end
    end
    return nc_mat, sparse(In_vec, Jn_vec, Vn_vec), nc_non_global_nodes, n_non_global_nodes, n_spa
end

function find_edge_node_match(ngl,iel,ip,connijk)
    found = false
    ip1 = 0
    iter = 1
    j = 0
    k = 0
    while (found == false && iter <= 4*ngl-3)
        if (iter <= ngl)
            k=1
            j=iter
        elseif (iter > ngl && iter < 2*ngl)
            j = ngl
            k = iter % (ngl)+1
        elseif (iter >= 2*ngl && iter < (3*ngl-1))
            k= ngl
            j = (ngl-1) - (iter % (ngl))
        else
            j = 1
            k = (ngl-1) - ((iter + 1) % ngl)
        end
        if (connijk[iel,k,j] == ip)
            ip1 = connijk[iel,k,j]
            found = true
        end
        iter += 1

    end
    if (found == false) 
        @info "failed to find" 
    end
    return k, j, ip1
end
