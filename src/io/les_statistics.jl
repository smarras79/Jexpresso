struct LESStatCache
    z_levels::Vector{Float64}        # sorted unique z values (global)
    z_groups::Vector{Vector{Int64}}  # local indices of points at each z-level
    npts_per_z::Vector{Int64}        # global count of points at each z-level
    local_sum_u::Vector{Float64}     # work buffer for local sums
    global_sum_u::Vector{Float64}    # work buffer for global sums
end

function build_les_stat_cache(mesh)
    comm = MPI.COMM_WORLD

    z = Array(mesh.z)  # ensure CPU array

    # Gather all unique z values across all ranks
    local_z_unique = sort(unique(round.(z; digits=8)))
    rank = MPI.Comm_rank(comm)

    # Gather variable-length arrays to rank 0, merge, then broadcast back
    all_z_arrays = MPI.gather(local_z_unique, comm)
    if rank == 0
        z_levels = sort(unique(round.(vcat(all_z_arrays...); digits=8)))
    else
        z_levels = Float64[]
    end
    nz_global = MPI.bcast(length(z_levels), 0, comm)
    if rank != 0
        resize!(z_levels, nz_global)
    end
    MPI.Bcast!(z_levels, 0, comm)

    nz = length(z_levels)

    # Build local index groups for each global z-level
    z_groups = Vector{Vector{Int64}}(undef, nz)
    local_counts = zeros(Int64, nz)
    for iz in 1:nz
        ids = findall(x -> abs(x - z_levels[iz]) < 1e-6, z)
        z_groups[iz] = ids
        local_counts[iz] = length(ids)
    end

    # Get global point counts per z-level
    npts_per_z = zeros(Int64, nz)
    MPI.Allreduce!(local_counts, npts_per_z, MPI.SUM, comm)

    return LESStatCache(z_levels, z_groups, npts_per_z, zeros(nz), zeros(nz))
end

function les_statistics(u, params, time)
    mesh   = params.mesh
    npoin  = mesh.npoin
    neqs   = params.neqs
    cache  = params.les_stat_cache
    nz     = length(cache.z_levels)
    ngl    = params.mesh.ngl

    wθ     = params.WM.wθ
    wqv    = params.WM.wqv
    
    nfaces_bdy       = params.mesh.nfaces_bdy
    bdy_face_type    = params.mesh.bdy_face_type
    poin_in_bdy_face = params.mesh.poin_in_bdy_face

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)

    # Convert u (flat vector) to uaux (npoin × neqs)
    uaux = params.uaux
    qn   = params.mp.qn
    u2uaux!(@view(uaux[:,:]), u, neqs, npoin)

    # Compute local sum of u-velocity at each z-level
    local_sum_u  = cache.local_sum_u
    global_sum_u = cache.global_sum_u
    fill!(local_sum_u, 0.0)
    for iz in 1:nz
        s = 0.0
        for ip in cache.z_groups[iz]
            # ρ      = uaux[ip, 1]
            # ρu_val = uaux[ip, 2]
            # s += ρu_val / ρ
            qn_val   = qn[ip]
            s += qn_val
        end
        local_sum_u[iz] = s
    end

    PhysConst = PhysicalConst{Float64}() 
    wθ_sum  = 0.0
    wqv_sum = 0.0
    cnt     = 0.0
    for iface = 1:nfaces_bdy
        if (bdy_face_type[iface] == "MOST")
            for i = 1:ngl
                for j = 1:ngl
                    ip       = poin_in_bdy_face[iface,i,j]
                    ρ        = uaux[ip, 1]
                    wθ_sum  += ρ*PhysConst.cp*wθ[iface,i,j,1]
                    wqv_sum += ρ*PhysConst.Lc*wqv[iface,i,j,1]
                    cnt     += 1.0
                end
            end
        end
    end
    most_sum_local  = [wθ_sum, wqv_sum, cnt]
    most_sum_global = zeros(3)

    # Global reduction
    MPI.Allreduce!(local_sum_u, global_sum_u, MPI.SUM, comm)
    MPI.Allreduce!(most_sum_local, most_sum_global, MPI.SUM, comm)

    # Only rank 0 writes the file
    if rank == 0
        outfile_les = joinpath(params.inputs[:output_dir], "les_statistics.dat")
        if !isfile(outfile_les)
            open(outfile_les, "w") do io
                println(io, "# time  z  qn")
            end
        end
        open(outfile_les, "a") do io
            for iz in 1:nz
                @printf(io, "%.6e  %.6e  %.6e\n", time, cache.z_levels[iz], global_sum_u[iz] / cache.npts_per_z[iz])
            end
            println(io)
        end


        outfile_most = joinpath(params.inputs[:output_dir], "MOST_statistics.dat")
        if !isfile(outfile_most)
            open(outfile_most, "w") do io
                println(io, "# time  LHF  SHF")
            end
        end
        open(outfile_most, "a") do io
            @printf(io, "%.6e  %.6e  %.6e\n", time, most_sum_global[2]/most_sum_global[3], most_sum_global[1]/most_sum_global[3])
        end
    end
end
