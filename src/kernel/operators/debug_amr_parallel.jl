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
