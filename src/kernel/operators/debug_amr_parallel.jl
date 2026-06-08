# =========================================================================
# Debug utilities for comparing serial vs. parallel AMR results
# =========================================================================

using JLD2, Serialization, MPI

"""
    save_serial_rhs(RHS, connijk_spa, coords_spa, coords_ang, filename)

Save the serial RHS vector along with spatial-angular coordinates for later comparison.
Call this in serial mode (or with 1 processor) after computing the final RHS.

# Arguments:
- `RHS`: The final restricted RHS vector (n_free)
- `connijk_spa`: Spatial-angular connectivity (n_spa × 3 or similar)
- `coords_spa`: Spatial coordinates (3 × n_spatial or n_spatial × 3)
- `coords_ang`: Angular coordinates (2 × n_angular or n_angular × 2)
- `filename`: Where to save (e.g., "serial_rhs_reference.jld2")
"""
function save_serial_rhs(RHS, connijk_spa, connijk, x, y, z, n_spa, coords_ang, nelem, ngl, nelem_ang, nop_ang, connijk_ang, filename)

    # Convert to dense if sparse
    rhs_vec = vec(collect(RHS))
    n_free = length(rhs_vec)
    coords_spa = zeros(n_free,5)
    #collect data into n_free size vectors
    for iel =1:nelem
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
                                ip_spa = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                                if (ip_spa <= n_free)
                                    coords_spa[ip_spa,1] = x_p
                                    coords_spa[ip_spa,2] = y_p
                                    coords_spa[ip_spa,3] = z_p
                                    coords_spa[ip_spa,4] = θ_p
                                    coords_spa[ip_spa,5] = ϕ_p
                                end

                            end
                        end
                    end
                end
            end
        end
    end
    # Save data with coordinates
    JLD2.jldopen(filename, "w") do file
        file["RHS"] = rhs_vec
        file["coords_spa"] = coords_spa# Free node connectivity
        file["n_free"] = n_free
    end

    @info "Saved serial RHS to $filename: n_free=$n_free, nnz=$(count(!iszero, rhs_vec))"
end

"""
    get_spatial_angular_coords(local_idx, connijk_spa, coords_spa, coords_ang)

Get the (x, y, z, θ, ϕ) coordinates for a spatial-angular node.

# Arguments:
- `local_idx`: Local node index
- `connijk_spa`: Spatial-angular connectivity
- `coords_spa`: Spatial coordinates
- `coords_ang`: Angular coordinates
"""
function get_spatial_angular_coords(local_idx, connijk_spa, coords_spa, coords_ang)
    # Get spatial and angular indices from connectivity
    if size(connijk_spa, 2) >= 2
        ip_spatial = connijk_spa[local_idx, 1]
        ip_angular = connijk_spa[local_idx, 2]
    else
        error("connijk_spa should have at least 2 columns")
    end

    # Get spatial coordinates
    if size(coords_spa, 1) == 3
        x = coords_spa[1, ip_spatial]
        y = coords_spa[2, ip_spatial]
        z = coords_spa[3, ip_spatial]
    else
        x = coords_spa[ip_spatial, 1]
        y = coords_spa[ip_spatial, 2]
        z = coords_spa[ip_spatial, 3]
    end

    # Get angular coordinates
    if size(coords_ang, 1) == 2
        theta = coords_ang[1, ip_angular]
        phi = coords_ang[2, ip_angular]
    else
        theta = coords_ang[ip_angular, 1]
        phi = coords_ang[ip_angular, 2]
    end

    return (x, y, z, theta, phi)
end

"""
    reduce_and_compare_parallel_rhs(RHS_local, connijk_spa, coords_spa, coords_ang,
                                     n_free, n_spa, comm, rank, reference_filename)

Gather the parallel RHS from all processors using spatial-angular coordinates as matching key,
and compare with the serial reference.

# Arguments:
- `RHS_local`: Local RHS vector after restriction (n_free)
- `connijk_spa`: Local spatial-angular connectivity
- `coords_spa`: Spatial coordinates
- `coords_ang`: Angular coordinates
- `n_free`: Number of free nodes locally
- `n_spa`: Total local spatial nodes
- `comm`: MPI communicator
- `rank`: Rank of this processor
- `reference_filename`: Path to serial reference file
"""
function reduce_and_compare_parallel_rhs(RHS_local, connijk_spa, connijk, x, y, z,
                                        n_free, n_spa,
                                        all_hanging_nodes,
                                        coords_ang, nelem, ngl, nelem_ang, nop_ang, connijk_ang, 
                                        comm, rank, reference_filename)

    nproc = MPI.Comm_size(comm)

    # =====================================================================
    # Step 1: Extract local free node RHS entries with coordinates
    # =====================================================================

    local_rhs_vec = vec(collect(RHS_local))

    # Create mapping of coordinates to RHS values
    local_entries = Tuple{NTuple{5, Float64}, Float64}[]  # ((x,y,z,θ,ϕ), value)

    #Get coordinate data
    coords_spa = zeros(n_free,5)
    for iel =1:nelem
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
                                ip_spa = connijk_spa[iel][i,j,k,e_ext,iθ,jθ]
                                if (ip_spa <= n_free)
                                    coords_spa[ip_spa,1] = x_p
                                    coords_spa[ip_spa,2] = y_p
                                    coords_spa[ip_spa,3] = z_p
                                    coords_spa[ip_spa,4] = θ_p
                                    coords_spa[ip_spa,5] = ϕ_p
                                end

                            end
                        end
                    end
                end
            end
        end
    end

    for local_idx = 1:n_free
        if !(local_idx in all_hanging_nodes)
            val = local_rhs_vec[local_idx]
            #coord = get_spatial_angular_coords(local_idx, connijk_spa, coords_spa, coords_ang)
            coord = (coords_spa[local_idx,1], coords_spa[local_idx,2], coords_spa[local_idx,3], coords_spa[local_idx,4], coords_spa[local_idx,5])
            push!(local_entries, (coord, val))
        end
    end

    @info "[Rank $rank] Local free nodes: $(length(local_entries)) entries"

    # =====================================================================
    # Step 2: Gather all entries from all processors
    # =====================================================================

    # Use binary serialization to handle variable-length data
    local_buffer = IOBuffer()
    serialize(local_buffer, local_entries)
    local_data = take!(local_buffer)

    # Gather sizes
    local_size = Int32(length(local_data))
    all_sizes = MPI.Allgather([local_size], comm)

    # Gather data
    if sum(all_sizes) > 0
        all_data_flat = MPI.Allgatherv(local_data, Int32.(all_sizes), comm)

        # Deserialize all data
        global_entries = Dict{NTuple{5, Float64}, Float64}()

        offset = 1
        for rank_idx = 0:(nproc-1)
            size = all_sizes[rank_idx+1]
            if size > 0
                chunk = all_data_flat[offset:offset+size-1]
                rank_entries = deserialize(IOBuffer(chunk))
                for (coord, val) in rank_entries
                    # Sum contributions from all processors
                    # (Ideally each coordinate should appear only once)
                    if haskey(global_entries, coord)
                        global_entries[coord] += val
                    else
                        global_entries[coord] = val
                    end
                end
            end
            offset += size
        end
    else
        global_entries = Dict{NTuple{5, Float64}, Float64}()
    end

    # =====================================================================
    # Step 3: Load reference and compare (on rank 0 only)
    # =====================================================================

    if rank == 0
            ref_data = JLD2.jldopen(reference_filename, "r") do file
                (
                    RHS = file["RHS"],
                    coords_spa = file["coords_spa"],
                    n_free = file["n_free"]
                )
            end

            ref_rhs = ref_data.RHS
            ref_coords_spa = ref_data.coords_spa
            ref_n_free = ref_data.n_free

            @info "Loaded reference: n_free=$ref_n_free, nnz=$(count(!iszero, ref_rhs))"
            @info "checking coords size", size(ref_coords_spa)
            # Build reference dictionary for comparison
            ref_entries = Dict{NTuple{5, Float64}, Float64}()

            for idx = 1:ref_n_free
                #coord = get_spatial_angular_coords(idx, ref_connijk, ref_coords_spa, ref_coords_ang)
                coord = (ref_coords_spa[idx,1], ref_coords_spa[idx,2], ref_coords_spa[idx,3], ref_coords_spa[idx,4], ref_coords_spa[idx,5])
                ref_entries[coord] = ref_rhs[idx]
            end

            # ================================================================
            # Compare
            # ================================================================

            @info "=== RHS COMPARISON ==="

            # Check all coordinates that appear in either
            all_coords = union(Set(keys(global_entries)), Set(keys(ref_entries)))

            max_rel_error = 0.0
            max_abs_error = 0.0
            n_different = 0
            n_missing_parallel = 0
            n_missing_serial = 0

            for coord in all_coords
                parallel_val = get(global_entries, coord, 0.0)
                serial_val = get(ref_entries, coord, 0.0)

                abs_err = abs(parallel_val - serial_val)

                if abs_err > 1e-7
                    n_different += 1
                    max_abs_error = max(max_abs_error, abs_err)

                    if abs(serial_val) > 1e-15
                        rel_err = abs_err / abs(serial_val)
                        max_rel_error = max(max_rel_error, rel_err)
                    end

                    if n_different <= 20  # Print first 20 differences
                        @info "  Difference at coord=$coord: serial=$serial_val, parallel=$parallel_val, abs_err=$abs_err"
                    end
                end

                if !haskey(global_entries, coord) && haskey(ref_entries, coord)
                    n_missing_parallel += 1
                end
                if !haskey(ref_entries, coord) && haskey(global_entries, coord)
                    n_missing_serial += 1
                end
            end

            @info "Total entries: serial=$(length(ref_entries)), parallel=$(length(global_entries))"
            @info "Missing from parallel: $n_missing_parallel"
            @info "Extra in parallel: $n_missing_serial"
            @info "Value mismatches: $n_different entries with |error| > 1e-14"
            @info "Max absolute error: $max_abs_error"
            @info "Max relative error: $max_rel_error"

            if n_different == 0 && n_missing_parallel == 0 && n_missing_serial == 0
                @info "✓ RHS VECTORS MATCH PERFECTLY!"
                return true
            else
                @warn "✗ RHS VECTORS DO NOT MATCH"
                return false
            end

        
    end

    return nothing
end

# =========================================================================
# Matrix Save and Compare Functions
# =========================================================================

"""
    save_serial_matrix(A, connijk_spa, connijk, x, y, z, n_spa, coords_ang,
                       nelem, ngl, nelem_ang, nop_ang, connijk_ang, filename)

Save the serial A matrix along with spatial-angular coordinates for later comparison.
Call this in serial mode (or with 1 processor) after constructing the final A matrix.

# Arguments:
- `A`: The final restricted and prolonged A matrix (n_free × n_free)
- `connijk_spa`: Spatial-angular connectivity
- `connijk`: Spatial connectivity
- `x, y, z`: Spatial coordinates
- `n_spa`: Number of spatial nodes
- `coords_ang`: Angular coordinates
- `nelem`: Number of spatial elements
- `ngl`: Spectral order + 1 in spatial domain
- `nelem_ang`: Angular elements per spatial element
- `nop_ang`: Angular order per element
- `connijk_ang`: Angular connectivity
- `filename`: Where to save (e.g., "serial_matrix_reference.jld2")
"""
function save_serial_matrix(A, connijk_spa, connijk, x, y, z, n_spa, coords_ang,
                            nelem, ngl, nelem_ang, nop_ang, connijk_ang, filename)

    # Convert to dense if sparse for storage
    A_dense = sparse(A)
    n_free = size(A, 1)
    coords_spa = zeros(n_free, 5)

    # Collect coordinate data
    for iel = 1:nelem
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
                                ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                                if (ip_spa <= n_free)
                                    coords_spa[ip_spa, 1] = x_p
                                    coords_spa[ip_spa, 2] = y_p
                                    coords_spa[ip_spa, 3] = z_p
                                    coords_spa[ip_spa, 4] = θ_p
                                    coords_spa[ip_spa, 5] = ϕ_p
                                end
                            end
                        end
                    end
                end
            end
        end
    end

    # Save data with coordinates
    JLD2.jldopen(filename, "w") do file
        file["A"] = A_dense
        file["coords_spa"] = coords_spa
        file["n_free"] = n_free
    end

    @info "Saved serial matrix to $filename: n_free=$n_free, nnz=$(nnz(A_dense))"
end

"""
    reduce_and_compare_parallel_matrix(A_local, connijk_spa, connijk, x, y, z,
                                       n_free, n_spa, all_hanging_nodes, coords_ang,
                                       nelem, ngl, nelem_ang, nop_ang, connijk_ang,
                                       comm, rank, reference_filename)

Gather the parallel A matrix from all processors using spatial-angular coordinates as matching key,
and compare with the serial reference.

# Arguments:
- `A_local`: Local A matrix after restriction and prolongation (n_free × n_free)
- `connijk_spa`: Local spatial-angular connectivity
- `connijk`: Spatial connectivity
- `x, y, z`: Spatial coordinates
- `n_free`: Number of free nodes locally
- `n_spa`: Total local spatial nodes
- `all_hanging_nodes`: Indices of hanging nodes
- `coords_ang`: Angular coordinates
- `nelem`: Number of spatial elements
- `ngl`: Spectral order + 1
- `nelem_ang`: Angular elements per spatial element
- `nop_ang`: Angular order per element
- `connijk_ang`: Angular connectivity
- `comm`: MPI communicator
- `rank`: Rank of this processor
- `reference_filename`: Path to serial reference file
"""
function reduce_and_compare_parallel_matrix(A_local, connijk_spa, connijk, x, y, z,
                                            n_free, n_spa, all_hanging_nodes, coords_ang,
                                            nelem, ngl, nelem_ang, nop_ang, connijk_ang,
                                            extended_parents_to_gid, n_spa_g,
                                            extended_parents_x, extended_parents_y, extended_parents_z, extended_parents_θ, extended_parents_ϕ,
                                            comm, rank, reference_filename)

    nproc = MPI.Comm_size(comm)

    # =====================================================================
    # Step 1: Extract local free node matrix entries with coordinates
    # =====================================================================

    # Create mapping of (row_coord, col_coord) to matrix value
    local_matrix_entries = Tuple{NTuple{5, Float64}, NTuple{5, Float64}, Float64}[]

    # Get coordinate data
    n_extended = n_spa_g - n_spa
    coords_spa = zeros(n_free+n_extended, 5)
    #get regular node data
    for iel = 1:nelem
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
                                ip_spa = connijk_spa[iel][i, j, k, e_ext, iθ, jθ]
                                if (ip_spa <= n_free)
                                    coords_spa[ip_spa, 1] = x_p
                                    coords_spa[ip_spa, 2] = y_p
                                    coords_spa[ip_spa, 3] = z_p
                                    coords_spa[ip_spa, 4] = θ_p
                                    coords_spa[ip_spa, 5] = ϕ_p
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    #get extended node coordinates
    for i=1:n_extended
        ip_spa = n_free + i
        coords_spa[ip_spa, 1] = extended_parents_x[i]
        coords_spa[ip_spa, 2] = extended_parents_y[i]
        coords_spa[ip_spa, 3] = extended_parents_z[i]
        coords_spa[ip_spa, 4] = extended_parents_θ[i]
        coords_spa[ip_spa, 5] = extended_parents_ϕ[i]
    end



    # Extract non-zero entries from local sparse matrix (efficient)
    hanging_nodes_set = Set(all_hanging_nodes)

    # Iterate only over stored non-zero entries
    if issparse(A_local)
        rows = rowvals(A_local)
        vals = nonzeros(A_local)

        for col_idx = 1:n_spa_g
            if !(col_idx in hanging_nodes_set)
                c_idx = 0
                if (col_idx > n_spa)
                    c_idx = col_idx-n_spa+n_free
                else
                    c_idx = col_idx
                end
                col_coord = (coords_spa[c_idx, 1], coords_spa[c_idx, 2], coords_spa[c_idx, 3],
                            coords_spa[c_idx, 4], coords_spa[c_idx, 5])

                for idx in nzrange(A_local, col_idx)
                    row_idx = rows[idx]
                    val = vals[idx]

                    if !(row_idx in hanging_nodes_set)
                        r_idx = 0
                        if (row_idx > n_spa)
                            r_idx = row_idx-n_spa+n_free
                        else
                            r_idx = row_idx
                        end
                        row_coord = (coords_spa[r_idx, 1], coords_spa[r_idx, 2], coords_spa[r_idx, 3],
                                    coords_spa[r_idx, 4], coords_spa[r_idx, 5])
                        push!(local_matrix_entries, (row_coord, col_coord, val))
                    end
                end
            end
        end
        # do extended parent nodes
    else
        # Fallback for dense matrices
        for row_idx = 1:n_free
            if !(row_idx in hanging_nodes_set)
                row_coord = (coords_spa[row_idx, 1], coords_spa[row_idx, 2], coords_spa[row_idx, 3],
                            coords_spa[row_idx, 4], coords_spa[row_idx, 5])

                for col_idx = 1:n_free
                    val = A_local[row_idx, col_idx]
                    if abs(val) > 1e-15
                        if !(col_idx in hanging_nodes_set)
                            col_coord = (coords_spa[col_idx, 1], coords_spa[col_idx, 2], coords_spa[col_idx, 3],
                                        coords_spa[col_idx, 4], coords_spa[col_idx, 5])
                            push!(local_matrix_entries, (row_coord, col_coord, val))
                        end
                    end
                end
            end
        end
    end

    @info "[Rank $rank] Local matrix entries: $(length(local_matrix_entries)) non-zero entries"

    # =====================================================================
    # Step 2: Gather all entries from all processors
    # =====================================================================

    local_buffer = IOBuffer()
    serialize(local_buffer, local_matrix_entries)
    local_data = take!(local_buffer)

    # Gather sizes
    local_size = Int32(length(local_data))
    all_sizes = MPI.Allgather([local_size], comm)

    # Gather data
    if sum(all_sizes) > 0
        all_data_flat = MPI.Allgatherv(local_data, Int32.(all_sizes), comm)

        # Deserialize all data
        global_matrix_entries = Dict{Tuple{NTuple{5, Float64}, NTuple{5, Float64}}, Float64}()

        offset = 1
        for rank_idx = 0:(nproc-1)
            size = all_sizes[rank_idx+1]
            if size > 0
                chunk = all_data_flat[offset:offset+size-1]
                rank_entries = deserialize(IOBuffer(chunk))
                for (row_coord, col_coord, val) in rank_entries
                    key = (row_coord, col_coord)
                    # Sum contributions from all processors
                    if haskey(global_matrix_entries, key)
                        global_matrix_entries[key] += val
                    else
                        global_matrix_entries[key] = val
                    end
                end
            end
            offset += size
        end
    else
        global_matrix_entries = Dict{Tuple{NTuple{5, Float64}, NTuple{5, Float64}}, Float64}()
    end

    # =====================================================================
    # Step 3: Load reference and compare (on rank 0 only)
    # =====================================================================

    if rank == 0
        ref_data = JLD2.jldopen(reference_filename, "r") do file
            (
                A = file["A"],
                coords_spa = file["coords_spa"],
                n_free = file["n_free"]
            )
        end

        ref_A = ref_data.A
        ref_coords_spa = ref_data.coords_spa
        ref_n_free = ref_data.n_free

        @info "Loaded reference matrix: n_free=$ref_n_free, nnz=$(nnz(ref_A))"

        # Build reference dictionary for comparison (efficient sparse iteration)
        ref_matrix_entries = Dict{Tuple{NTuple{5, Float64}, NTuple{5, Float64}}, Float64}()

        if issparse(ref_A)
            rows = rowvals(ref_A)
            vals = nonzeros(ref_A)

            for col_idx = 1:ref_n_free
                col_coord = (ref_coords_spa[col_idx, 1], ref_coords_spa[col_idx, 2], ref_coords_spa[col_idx, 3],
                            ref_coords_spa[col_idx, 4], ref_coords_spa[col_idx, 5])

                for idx in nzrange(ref_A, col_idx)
                    row_idx = rows[idx]
                    val = vals[idx]

                    row_coord = (ref_coords_spa[row_idx, 1], ref_coords_spa[row_idx, 2], ref_coords_spa[row_idx, 3],
                                ref_coords_spa[row_idx, 4], ref_coords_spa[row_idx, 5])
                    key = (row_coord, col_coord)
                    ref_matrix_entries[key] = val
                end
            end
        else
            # Fallback for dense matrices
            for row_idx = 1:ref_n_free
                row_coord = (ref_coords_spa[row_idx, 1], ref_coords_spa[row_idx, 2], ref_coords_spa[row_idx, 3],
                            ref_coords_spa[row_idx, 4], ref_coords_spa[row_idx, 5])

                for col_idx = 1:ref_n_free
                    val = ref_A[row_idx, col_idx]
                    if abs(val) > 1e-15
                        col_coord = (ref_coords_spa[col_idx, 1], ref_coords_spa[col_idx, 2], ref_coords_spa[col_idx, 3],
                                    ref_coords_spa[col_idx, 4], ref_coords_spa[col_idx, 5])
                        key = (row_coord, col_coord)
                        ref_matrix_entries[key] = val
                    end
                end
            end
        end

        # ================================================================
        # Compare
        # ================================================================

        @info "=== MATRIX COMPARISON ==="

        all_keys = union(Set(keys(global_matrix_entries)), Set(keys(ref_matrix_entries)))

        max_rel_error = 0.0
        max_abs_error = 0.0
        n_different = 0
        n_missing_parallel = 0
        n_missing_serial = 0

        for key in all_keys
            parallel_val = get(global_matrix_entries, key, 0.0)
            serial_val = get(ref_matrix_entries, key, 0.0)

            abs_err = abs(parallel_val - serial_val)

            if abs_err > 1e-7
                n_different += 1
                max_abs_error = max(max_abs_error, abs_err)

                if abs(serial_val) > 1e-15
                    rel_err = abs_err / abs(serial_val)
                    max_rel_error = max(max_rel_error, rel_err)
                end

                if n_different <= 20  # Print first 20 differences
                    row_coord, col_coord = key
                    @info "  Difference at entry (row=$row_coord, col=$col_coord): serial=$serial_val, parallel=$parallel_val, abs_err=$abs_err"
                end
            end

            if !haskey(global_matrix_entries, key) && haskey(ref_matrix_entries, key)
                n_missing_parallel += 1
                #@info "missing parallel", key, ref_matrix_entries[key]
            end
            if !haskey(ref_matrix_entries, key) && haskey(global_matrix_entries, key)
                n_missing_serial += 1
                #@info "missing serial", key, global_matrix_entries[key]
            end
        end

        @info "Total entries: serial=$(length(ref_matrix_entries)), parallel=$(length(global_matrix_entries))"
        @info "Missing from parallel: $n_missing_parallel"
        @info "Extra in parallel: $n_missing_serial"
        @info "Value mismatches: $n_different entries with |error| > 1e-7"
        @info "Max absolute error: $max_abs_error"
        @info "Max relative error: $max_rel_error"

        if n_different == 0 && n_missing_parallel == 0 && n_missing_serial == 0
            @info "✓ MATRIX ENTRIES MATCH PERFECTLY!"
            return true
        else
            @warn "✗ MATRIX ENTRIES DO NOT MATCH"
            return false
        end
    end

    return nothing
end

# =========================================================================
# Spatial AMR + Uniform Angular Mesh Debug Utilities
# =========================================================================
# These mirror the adaptive-angular functions above but use the much simpler
# DOF index formula for the uniform angular case:
#   DOF i  →  ip_spa = div(i-1, n_ang) + 1,  ip_ang = mod(i-1, n_ang) + 1
# Coordinates come straight from mesh.x/y/z and extra_mesh.extra_coords.
# =========================================================================

"""
    dof_coords_spatial_amr(i, mesh, extra_mesh, n_ang) -> NTuple{5,Float64}

Return (x, y, z, θ, ϕ) coordinates for the i-th DOF in the uniform-angular
spatial-AMR system.
"""
function dof_coords_spatial_amr(i::Int, mesh, extra_mesh, n_ang::Int)
    ip_spa = div(i - 1, n_ang) + 1
    ip_ang = mod(i - 1, n_ang) + 1
    return (Float64(mesh.x[ip_spa]),
            Float64(mesh.y[ip_spa]),
            Float64(mesh.z[ip_spa]),
            Float64(extra_mesh.extra_coords[1, ip_ang]),
            Float64(extra_mesh.extra_coords[2, ip_ang]))
end

"""
    save_serial_spatial_amr(RHS, As, A_with_rows, mesh, extra_mesh, n_ang,
                            npoin_ang_total, hanging_nodes, filename)

Save the serial RHS vector and row-nnz counts for A_with_rows (post-restriction,
pre-prolongation) as the reference for later parallel comparison.  Call with 1 rank.
"""
function save_serial_spatial_amr(RHS, As, A_with_rows, mesh, extra_mesh, n_ang::Int,
                                  npoin_ang_total::Int, hanging_nodes,
                                  filename::String;
                                  dedup_coords::Union{Matrix{Float64},Nothing}=nothing)
    n = npoin_ang_total

    # Build coordinate array (one row per DOF).
    # With dedup numbering, the old (ip-1)*n_ang+ip_a formula is wrong — use the
    # pre-built dedup_coords matrix (n × 5) mapping ip_spa → (x,y,z,θ,ϕ) instead.
    coords = Matrix{Float64}(undef, n, 5)
    if dedup_coords !== nothing
        coords .= dedup_coords[1:n, :]
    else
        for i = 1:n
            c = dof_coords_spatial_amr(i, mesh, extra_mesh, n_ang)
            coords[i, 1] = c[1]; coords[i, 2] = c[2]; coords[i, 3] = c[3]
            coords[i, 4] = c[4]; coords[i, 5] = c[5]
        end
    end

    rhs_vec = vec(collect(Float64, RHS))

    # A_with_rows: store every non-zero as (row_ip, col_ip, val) so the parallel
    # run can do a value-based, coordinate-keyed comparison.
    # In serial all ip_spa indices are in [1,n] with valid dedup_coords entries.
    ar_row_ip = Int32[]
    ar_col_ip = Int32[]
    ar_val    = Float64[]
    if A_with_rows !== nothing
        Ar   = sparse(A_with_rows)
        Arvl = nonzeros(Ar)
        Arrv = rowvals(Ar)
        for col_ip = 1:min(size(Ar, 2), n)
            for ptr in nzrange(Ar, col_ip)
                row_ip = Arrv[ptr]
                row_ip > n && continue
                v = Arvl[ptr]
                abs(v) < 1e-20 && continue
                push!(ar_row_ip, Int32(row_ip))
                push!(ar_col_ip, Int32(col_ip))
                push!(ar_val,    v)
            end
        end
        @info "  A_with_rows: $(length(ar_val)) entries stored"
    end

    JLD2.jldopen(filename, "w") do f
        f["RHS"]           = rhs_vec
        f["coords"]        = coords
        f["n"]             = n
        f["hanging_nodes"] = collect(Set(hanging_nodes))
        f["ar_row_ip"]     = ar_row_ip
        f["ar_col_ip"]     = ar_col_ip
        f["ar_val"]        = ar_val
    end
    @info "Saved serial spatial-AMR reference to $filename: n=$n, A_with_rows entries=$(length(ar_val))"
end

"""
    reduce_and_compare_parallel_spatial_amr(RHS, As, A_with_rows, mesh, extra_mesh,
                                            n_ang, npoin_ang_total, hanging_nodes,
                                            ip2gip_spa, gnpoin,
                                            extended_parents_to_gid_spa, comm, rank, filename)

Compare the parallel assembly against the serial reference:
- RHS: coordinate-matched value comparison
- A_with_rows (post-restriction, pre-prolongation): coordinate-matched row-nnz comparison.
  Reports how many rows have fewer/more nonzero entries in parallel vs serial and a
  histogram of the differences, making it easy to identify whether the parallel matrix
  is structurally missing or gaining entries after restriction.
"""
function reduce_and_compare_parallel_spatial_amr(
    RHS, As, A_with_rows, mesh, extra_mesh, n_ang::Int, npoin_ang_total::Int,
    hanging_nodes, ip2gip_spa, gnpoin::Int,
    extended_parents_to_gid_spa::Vector{Int},
    comm, rank::Int, filename::String;
    dedup_coords::Union{Matrix{Float64},Nothing}=nothing,
    gip2owner_extra::Union{Vector{Int},Nothing}=nothing,
    ext_parent_coords::Union{Vector{NTuple{5,Float64}},Nothing}=nothing
)
    nproc       = MPI.Comm_size(comm)
    hanging_set = Set(hanging_nodes)
    n           = npoin_ang_total

    # With dedup numbering use the pre-built coordinate map; fall back to old formula.
    get_coord = if dedup_coords !== nothing
        i -> (dedup_coords[i,1], dedup_coords[i,2], dedup_coords[i,3],
              dedup_coords[i,4], dedup_coords[i,5])
    else
        i -> dof_coords_spatial_amr(i, mesh, extra_mesh, n_ang)
    end

    # Ownership predicate: use gip2owner_extra when available (dedup path),
    # otherwise fall back to mesh.gip2owner on spatial node index.
    owns_dof = if gip2owner_extra !== nothing
        i -> (i <= length(gip2owner_extra) && gip2owner_extra[i] == rank)
    else
        i -> begin
            npoin_spa = n ÷ n_ang
            spa_ip = div(i - 1, n_ang) + 1
            spa_ip <= npoin_spa && mesh.gip2owner[spa_ip] == rank
        end
    end

    # ── Step 1: RHS entries (coord → value) ─────────────────────────────────
    rhs_vec   = vec(collect(Float64, RHS))
    local_rhs = Tuple{NTuple{5,Float64}, Float64}[]
    for i = 1:n
        i in hanging_set && continue
        push!(local_rhs, (get_coord(i), rhs_vec[i]))
    end

    # ── Step 1b: Row-nnz entries for A_with_rows (post-restriction, pre-prolongation)
    # Collect from ALL local DOFs (not just owned). The gather sums contributions from
    # all ranks — equivalent to AllReducing the distributed matrix and comparing the
    # global nnz structure with serial. Shared DOFs (in both ranks' local meshes via
    # coordinate dedup) are counted on every rank that has them; the sum may therefore
    # overcount shared column entries, but if par_sum < serial any rank is missing entries.
    # ── Step 1b: A_with_rows entries (value-based, coordinate-keyed) ────────────
    # Pack each non-zero as 11 Float64: [row_x,y,z,θ,ϕ, col_x,y,z,θ,ϕ, val].
    # Regular indices (≤ n) use dedup_coords; extended parent indices (> n) use
    # ext_parent_coords[i-n].  Entries with NaN coords (unknown GID) are counted.
    local_ar_flat = Float64[]
    local_ar_nan  = 0
    if A_with_rows !== nothing && dedup_coords !== nothing
        function _coord_ar(ip)
            if ip <= n
                return (dedup_coords[ip,1], dedup_coords[ip,2],
                        dedup_coords[ip,3], dedup_coords[ip,4], dedup_coords[ip,5])
            elseif ext_parent_coords !== nothing && (ip - n) <= length(ext_parent_coords)
                return ext_parent_coords[ip - n]
            else
                return (NaN, NaN, NaN, NaN, NaN)
            end
        end
        Ar   = sparse(A_with_rows)
        Arvl = nonzeros(Ar)
        Arrv = rowvals(Ar)
        for col_ip = 1:size(Ar, 2)   # include extended cols (> n)
            col_coord = _coord_ar(col_ip)
            col_ok = !any(isnan, col_coord)
            for ptr in nzrange(Ar, col_ip)
                row_ip = Arrv[ptr]
                v = Arvl[ptr]
                abs(v) < 1e-20 && continue
                row_coord = _coord_ar(row_ip)
                row_ok = !any(isnan, row_coord)
                if !row_ok || !col_ok
                    local_ar_nan += 1
                    continue
                end
                push!(local_ar_flat,
                    row_coord[1], row_coord[2], row_coord[3], row_coord[4], row_coord[5],
                    col_coord[1], col_coord[2], col_coord[3], col_coord[4], col_coord[5], v)
            end
        end
    end
    n_ar_nan_g = MPI.Allreduce(local_ar_nan, +, comm)
    @info "[$rank] A_with_rows: $(length(local_ar_flat) ÷ 11) coord-matched entries (incl. extended), $local_ar_nan NaN-coord entries skipped (global=$n_ar_nan_g)"

    # ── Step 2: Gather to rank 0 ─────────────────────────────────────────────
    function gather_to_rank0(local_entries, tag::Int)
        buf  = IOBuffer(); serialize(buf, local_entries); data = take!(buf)
        sz   = Int64(length(data))
        all_sz = MPI.Allgather([sz], comm)
        if rank == 0
            result = collect(deserialize(IOBuffer(data)))
            for r = 1:nproc-1
                s = all_sz[r+1]
                if s > 0
                    rbuf = Vector{UInt8}(undef, s)
                    MPI.Recv!(rbuf, r, tag, comm)
                    append!(result, deserialize(IOBuffer(rbuf)))
                end
            end
            return result
        else
            MPI.Send(data, 0, tag, comm)
            return []
        end
    end

    all_rhs_entries = gather_to_rank0(local_rhs, 42)

    # Gather A_with_rows entry counts (tiny collective) then use point-to-point
    # so rank 0 processes one rank's data at a time and frees it immediately,
    # avoiding the peak memory of holding all ranks' data simultaneously.
    n_ar_counts = MPI.Allgather([Int32(length(local_ar_flat))], comm)
    if rank != 0
        MPI.Send(local_ar_flat, 0, 55, comm)
        local_ar_flat = nothing; GC.gc()
    end

    # ── Step 3: Compare on rank 0 ────────────────────────────────────────────
    if rank == 0
        ref = JLD2.jldopen(filename, "r") do f
            (coords    = f["coords"],
             n         = f["n"],
             RHS       = f["RHS"],
             ar_row_ip = haskey(f, "ar_row_ip") ? f["ar_row_ip"] : nothing,
             ar_col_ip = haskey(f, "ar_col_ip") ? f["ar_col_ip"] : nothing,
             ar_val    = haskey(f, "ar_val")    ? f["ar_val"]    : nothing)
        end
        ref_n = ref.n
        @info "Loaded serial reference: n=$ref_n"

        # ── Build serial RHS dict ─────────────────────────────────────────────
        ref_rhs_dict  = Dict{NTuple{5,Float64}, Float64}()
        n_serial_coll = 0
        for i = 1:ref_n
            c = (ref.coords[i,1], ref.coords[i,2], ref.coords[i,3],
                 ref.coords[i,4], ref.coords[i,5])
            if haskey(ref_rhs_dict, c)
                n_serial_coll += 1
            else
                ref_rhs_dict[c] = ref.RHS[i]
            end
        end
        n_serial_coll > 0 && @warn "serial reference: $n_serial_coll coord collisions (skipped)"
        # ref kept alive until after A_with_rows comparison (needs ref.ar_row_ip etc.)

        # ── RHS comparison ────────────────────────────────────────────────────
        par_rhs_dict  = Dict{NTuple{5,Float64}, Float64}()
        par_rhs_count = Dict{NTuple{5,Float64}, Int}()
        for (coord, val) in all_rhs_entries
            par_rhs_dict[coord]  = get(par_rhs_dict, coord, 0.0) + val
            par_rhs_count[coord] = get(par_rhs_count, coord, 0) + 1
        end
        n_dup = count(v -> v > 1, values(par_rhs_count))
        n_dup > 0 && @warn "parallel RHS dict: $n_dup coords with >1 contribution"

        ref_rhs_free = Dict(k => v for (k,v) in ref_rhs_dict if abs(v) > 1e-20)
        @info "=== SPATIAL-AMR RHS COMPARISON ==="
        all_rhs_coords = union(keys(ref_rhs_free), keys(par_rhs_dict))
        rhs_ndiff = 0; rhs_max_abs = 0.0; rhs_nmiss = 0; rhs_nextra = 0
        for coord in all_rhs_coords
            sv = get(ref_rhs_free,  coord, 0.0)
            pv = get(par_rhs_dict, coord, 0.0)
            ae = abs(pv - sv)
            ae > 1e-9 && (rhs_ndiff += 1; rhs_max_abs = max(rhs_max_abs, ae);
                          rhs_ndiff <= 20 && @info "  RHS diff $coord: ser=$sv par=$pv")
            !haskey(par_rhs_dict,  coord) && haskey(ref_rhs_free, coord) && (rhs_nmiss  += 1)
            !haskey(ref_rhs_free, coord) && haskey(par_rhs_dict,  coord) && (rhs_nextra += 1)
        end
        @info "RHS: mismatches=$rhs_ndiff, missing_par=$rhs_nmiss, extra_par=$rhs_nextra, max_err=$rhs_max_abs"
        rhs_ok = (rhs_ndiff == 0 && rhs_nmiss == 0 && rhs_nextra == 0)
        rhs_ok ? @info("✓ RHS MATCHES SERIAL") : @warn("✗ RHS DOES NOT MATCH SERIAL")

        ref_rhs_dict = nothing; ref_rhs_free = nothing
        par_rhs_dict = nothing; par_rhs_count = nothing
        all_rhs_entries = nothing; GC.gc()

        # ── A_with_rows value comparison (coordinate-keyed) ──────────────────
        mat_ok = true
        if ref.ar_row_ip === nothing
            @warn "=== A_WITH_ROWS COMPARISON SKIPPED: no ar_row_ip in reference — delete .jld2 and re-run with 1 rank ==="
        elseif sum(n_ar_counts) == 0
            @warn "=== A_WITH_ROWS COMPARISON SKIPPED: parallel sent no entries ==="
        else
            @info "=== A_WITH_ROWS VALUE COMPARISON (all-reduced post-restriction, coord-keyed) ==="
            if n_ar_nan_g > 0
                @warn "  $n_ar_nan_g parallel entries skipped due to NaN dedup_coords — these represent a parallel assembly error"
            end

            # Build serial dict: key = NTuple{10,Float16} (row+col coords, 20 bytes),
            # value = Float32 (4 bytes). Total ~24 bytes/entry vs 88 bytes for Float64.
            ser_coords = ref.coords
            ser_rip = ref.ar_row_ip; ser_cip = ref.ar_col_ip; ser_av = ref.ar_val
            ser_ar = Dict{NTuple{10,Float16}, Float32}()
            sizehint!(ser_ar, length(ser_av))
            for k in eachindex(ser_av)
                r = Int(ser_rip[k]); c = Int(ser_cip[k])
                key = (Float16(ser_coords[r,1]), Float16(ser_coords[r,2]),
                       Float16(ser_coords[r,3]), Float16(ser_coords[r,4]),
                       Float16(ser_coords[r,5]),
                       Float16(ser_coords[c,1]), Float16(ser_coords[c,2]),
                       Float16(ser_coords[c,3]), Float16(ser_coords[c,4]),
                       Float16(ser_coords[c,5]))
                ser_ar[key] = Float32(ser_av[k])
            end
            ser_rip = nothing; ser_cip = nothing; ser_av = nothing
            ser_coords = nothing; ref = nothing; GC.gc()

            # Build parallel dict by receiving one rank's data at a time.
            function _add_flat_to_dict!(data, d)
                n = length(data) ÷ 11
                sizehint!(d, length(d) + n)
                for k = 1:n
                    s = 11*(k-1)
                    key = ntuple(j -> Float16(data[s+j]), Val(10))
                    d[key] = get(d, key, Float32(0)) + Float32(data[s+11])
                end
            end
            par_ar = Dict{NTuple{10,Float16}, Float32}()
            _add_flat_to_dict!(local_ar_flat, par_ar)
            local_ar_flat = nothing; GC.gc()
            for r = 1:nproc-1
                n_recv = Int(n_ar_counts[r+1])
                if n_recv > 0
                    recv_buf = Vector{Float64}(undef, n_recv)
                    MPI.Recv!(recv_buf, r, 55, comm)
                    _add_flat_to_dict!(recv_buf, par_ar)
                    recv_buf = nothing; GC.gc()
                end
            end

            # Compare — key is NTuple{10,Float16} so coords are printed directly.
            tol = 1e-3   # relaxed to account for Float16 coord precision
            all_ar_keys = union(keys(ser_ar), keys(par_ar))
            ar_ndiff = 0; ar_nmiss = 0; ar_nextra = 0; ar_max_err = 0.0f0
            for key in all_ar_keys
                in_ser = haskey(ser_ar, key)
                in_par = haskey(par_ar, key)
                if !in_par
                    ar_nmiss += 1
                    if ar_nmiss <= 20
                        rk = key[1:5]; ck = key[6:10]
                        sv = ser_ar[key]
                        @info "  MISSING_PAR row=$rk col=$ck ser=$sv"
                    end
                    continue
                end
                if !in_ser; ar_nextra += 1; continue; end
                sv = ser_ar[key]; pv = par_ar[key]
                err = abs(pv - sv) / max(abs(sv), Float32(1e-10))
                if err > tol
                    ar_ndiff += 1
                    ar_max_err = max(ar_max_err, err)
                    if ar_ndiff <= 20
                        rk = key[1:5]; ck = key[6:10]
                        @info "  DIFF row=$rk col=$ck ser=$sv par=$pv rel_err=$(round(err,sigdigits=3))"
                    end
                end
            end
            @info "A_with_rows: entries_ser=$(length(ser_ar)) entries_par=$(length(par_ar))"
            @info "A_with_rows: value_diffs=$ar_ndiff missing_par=$ar_nmiss extra_par=$ar_nextra max_rel_err=$ar_max_err"
            mat_ok = (ar_ndiff == 0 && ar_nmiss == 0 && ar_nextra == 0)
            mat_ok ? @info("✓ A_WITH_ROWS MATCHES SERIAL") : @warn("✗ A_WITH_ROWS DOES NOT MATCH SERIAL")
            ser_ar = nothing; par_ar = nothing; GC.gc()
        end
        ref = nothing; GC.gc()

        return rhs_ok && mat_ok
    end

    return nothing
end

"""
    diagnose_row(As, target_coord, dedup_coords, n, ext_parent_coords,
                 ip2gip_spa, rank, comm, label; tol)

Gather and print all non-zero column entries of `As` for the row whose
(x,y,z,θ,ϕ) coordinates match `target_coord`.  Works in both serial and
parallel: each rank searches its local DOFs (and extended parent rows), collects
entries, sends to rank 0 which prints the merged list.

Usage: call once in serial and once in parallel with the same target_coord
(taken from the nnz-comparison output), then compare the printed column lists
to find exactly which relationships are missing.
"""
function diagnose_row(
    As, target_coord::NTuple{5,Float64},
    dedup_coords::Matrix{Float64}, n::Int,
    ext_parent_coords::Union{Vector{NTuple{5,Float64}},Nothing},
    ip2gip_spa::Vector{Int},
    rank::Int, comm, label::String;
    tol::Float64 = 1e-2
)
    nproc = MPI.Comm_size(comm)
    n_ext = ext_parent_coords !== nothing ? length(ext_parent_coords) : 0

    # ── Locate target row ────────────────────────────────────────────────────
    target_row = 0
    for i = 1:n
        if abs(dedup_coords[i,1]-target_coord[1]) < tol &&
           abs(dedup_coords[i,2]-target_coord[2]) < tol &&
           abs(dedup_coords[i,3]-target_coord[3]) < tol &&
           abs(dedup_coords[i,4]-target_coord[4]) < tol &&
           abs(dedup_coords[i,5]-target_coord[5]) < tol
            target_row = i; break
        end
    end
    if target_row == 0 && ext_parent_coords !== nothing
        for k = 1:n_ext
            c = ext_parent_coords[k]
            if abs(c[1]-target_coord[1]) < tol && abs(c[2]-target_coord[2]) < tol &&
               abs(c[3]-target_coord[3]) < tol && abs(c[4]-target_coord[4]) < tol &&
               abs(c[5]-target_coord[5]) < tol
                target_row = n + k; break
            end
        end
    end

    # ── Collect local non-zeros for target row ───────────────────────────────
    local_entries = Tuple{NTuple{5,Float64}, Int, Float64}[]  # (col_coord, col_gid, value)
    if target_row > 0
        As_rvals  = rowvals(As)
        As_nzvals = nonzeros(As)
        for col = 1:size(As, 2)
            for ptr in nzrange(As, col)
                As_rvals[ptr] == target_row || continue
                val = As_nzvals[ptr]
                abs(val) < 1e-20 && continue
                col_coord = if col <= n
                    (dedup_coords[col,1], dedup_coords[col,2], dedup_coords[col,3],
                     dedup_coords[col,4], dedup_coords[col,5])
                elseif ext_parent_coords !== nothing && col-n <= n_ext
                    ext_parent_coords[col-n]
                else
                    continue
                end
                if any(isnan, col_coord)
                    col_gid_nan = col <= length(ip2gip_spa) ? ip2gip_spa[col] : -1
                    @info "[$rank] $label: SKIPPED col=$col gid=$col_gid_nan val=$val (NaN in dedup_coords — ip_spa exists in As but has no entry in _dedup_coords)"
                    continue
                end
                col_gid = col <= n ? ip2gip_spa[col] : -1
                push!(local_entries, (col_coord, col_gid, val))
            end
        end
        row_gid = target_row <= n ? ip2gip_spa[target_row] : -1
        @info "[$rank] $label: target found local_idx=$target_row gid=$row_gid $(length(local_entries)) local non-zeros"
    else
        @info "[$rank] $label: target coord not in local dedup_coords"
    end

    # ── Gather all entries to rank 0 ─────────────────────────────────────────
    buf = IOBuffer(); serialize(buf, local_entries); data = take!(buf)
    sz  = Int32(length(data))
    all_sz = MPI.Allgather([sz], comm)

    if rank == 0
        all_entries = collect(local_entries)
        for r = 1:nproc-1
            s = all_sz[r+1]
            if s > 0
                rbuf = Vector{UInt8}(undef, s)
                MPI.Recv!(rbuf, r, 99, comm)
                append!(all_entries, deserialize(IOBuffer(rbuf)))
            end
        end
        # Sum values for same col_coord (handles extended+regular rows for same parent)
        combined = Dict{NTuple{5,Float64}, Tuple{Int,Float64}}()
        for (cc, cg, v) in all_entries
            combined[cc] = haskey(combined, cc) ? (cg, combined[cc][2]+v) : (cg, v)
        end
        @info "[0] $label: $(length(combined)) unique column entries for target=$target_coord"
        for (cc, (cg, v)) in sort(collect(combined), by=x->x[1])
            @info "[0]   col=$cc  gid=$cg  val=$v"
        end
    else
        MPI.Send(data, 0, 99, comm)
    end
end
