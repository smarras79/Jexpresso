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
    save_serial_spatial_amr(RHS, As, mesh, extra_mesh, n_ang, npoin_ang_total,
                            hanging_nodes, filename)

Save the serial RHS vector and A matrix with 5-tuple coordinates as the
reference for later parallel comparison.  Call with 1 MPI rank.
"""
function save_serial_spatial_amr(RHS, As, mesh, extra_mesh, n_ang::Int,
                                  npoin_ang_total::Int, hanging_nodes,
                                  filename::String)
    n = npoin_ang_total

    # Build coordinate array (one row per DOF)
    coords = Matrix{Float64}(undef, n, 5)
    for i = 1:n
        c = dof_coords_spatial_amr(i, mesh, extra_mesh, n_ang)
        coords[i, 1] = c[1]; coords[i, 2] = c[2]; coords[i, 3] = c[3]
        coords[i, 4] = c[4]; coords[i, 5] = c[5]
    end

    rhs_vec = vec(collect(Float64, RHS))
    A_spa   = sparse(As)

    # Matvec reference: x_test[ip] = coord_func(x,y,z,θ,ϕ)
    # Using a linear combination of physical coordinates so that the test vector
    # value at each DOF is the SAME whether running with 1 or N MPI ranks
    # (i.e. independent of local/global index assignment).
    coord_func_val(c1,c2,c3,c4,c5) = c1 + 1.7*c2 + 2.3*c3 + 3.1*c4 + 4.7*c5
    x_test = Vector{Float64}(undef, n)
    for i = 1:n
        c = dof_coords_spatial_amr(i, mesh, extra_mesh, n_ang)
        x_test[i] = coord_func_val(c[1], c[2], c[3], c[4], c[5])
    end
    y_matvec = A_spa * x_test

    JLD2.jldopen(filename, "w") do f
        f["RHS"]      = rhs_vec
        f["A"]        = A_spa
        f["y_matvec"] = y_matvec
        f["coords"]   = coords
        f["n"]        = n
        f["hanging_nodes"] = collect(Set(hanging_nodes))
    end
    @info "Saved serial spatial-AMR reference to $filename: n=$n, nnz_A=$(nnz(A_spa))"
end

"""
    reduce_and_compare_parallel_spatial_amr(RHS, As, mesh, extra_mesh, n_ang,
                                            npoin_ang_total, hanging_nodes,
                                            comm, rank, filename)

Gather the parallel RHS and A matrix from all ranks using (x,y,z,θ,ϕ)
as the matching key, then compare with the serial reference saved by
`save_serial_spatial_amr`.  Reports mismatches, missing entries, and
maximum errors.
"""
function reduce_and_compare_parallel_spatial_amr(
    RHS, As, mesh, extra_mesh, n_ang::Int, npoin_ang_total::Int,
    hanging_nodes, ip2gip_spa, gnpoin::Int,
    extended_parents_to_gid_spa::Vector{Int},
    comm, rank::Int, filename::String
)
    nproc       = MPI.Comm_size(comm)
    hanging_set = Set(hanging_nodes)
    n           = npoin_ang_total
    n_ext       = length(extended_parents_to_gid_spa)
    n_ext_spa   = n + n_ext

    # ── Step 1: build local (coord, value) RHS entries ───────────────────────
    get_coord = i -> dof_coords_spatial_amr(i, mesh, extra_mesh, n_ang)

    rhs_vec = vec(collect(Float64, RHS))
    local_rhs = Tuple{NTuple{5,Float64}, Float64}[]
    for i = 1:n
        i in hanging_set && continue
        push!(local_rhs, (get_coord(i), rhs_vec[i]))
    end

    # ── Step 1b: coordinate-based matvec comparison ──────────────────────────
    # Test vector: x[dof] = coord_func(x,y,z,θ,ϕ) — value is determined by the
    # physical location, not by the local/global DOF index, so serial and parallel
    # assign the same x value to the same physical point regardless of numbering.
    #
    # Extended parents (not in local ghost layer): their spatial coordinates are
    # fetched via an AllGather of owned spatial node coordinates so every rank
    # can resolve any global spatial IP.
    coord_func_val(c1,c2,c3,c4,c5) = c1 + 1.7*c2 + 2.3*c3 + 3.1*c4 + 4.7*c5

    # AllGather owned spatial node (spatial_GIP, x, y, z) pairs
    npoin_spa = n ÷ n_ang
    local_spa_xyz = [(Int(mesh.ip2gip[i]), Float64(mesh.x[i]),
                      Float64(mesh.y[i]),   Float64(mesh.z[i]))
                     for i = 1:npoin_spa if mesh.gip2owner[i] == rank]
    sp_buf  = IOBuffer(); serialize(sp_buf, local_spa_xyz); sp_data = take!(sp_buf)
    sp_sz   = Int64(length(sp_data))
    all_sp_sz   = MPI.Allgather([sp_sz], comm)
    all_sp_data = MPI.Allgatherv(sp_data, Int32.(all_sp_sz), comm)
    spa_gip_to_xyz = Dict{Int, NTuple{3,Float64}}()
    off_sp = 1
    for r = 0:nproc-1
        s = all_sp_sz[r+1]
        if s > 0
            for (gip, xv, yv, zv) in deserialize(IOBuffer(all_sp_data[off_sp:off_sp+s-1]))
                spa_gip_to_xyz[gip] = (xv, yv, zv)
            end
        end
        off_sp += s
    end

    @info "[$rank] Computing coord-based matvec: size(As)=$(size(As)), n=$n, n_ext=$n_ext"
    x_local_coord = zeros(Float64, n_ext_spa)
    for ip = 1:n
        c = dof_coords_spatial_amr(ip, mesh, extra_mesh, n_ang)
        x_local_coord[ip] = coord_func_val(c[1], c[2], c[3], c[4], c[5])
    end
    for k = 1:n_ext
        gid_ext  = extended_parents_to_gid_spa[k]
        spa_gip  = (gid_ext - 1) ÷ n_ang + 1
        ang_idx  = (gid_ext - 1) % n_ang + 1
        xyz      = get(spa_gip_to_xyz, spa_gip, (0.0, 0.0, 0.0))
        θ_ext    = Float64(extra_mesh.extra_coords[1, ang_idx])
        ϕ_ext    = Float64(extra_mesh.extra_coords[2, ang_idx])
        x_local_coord[n + k] = coord_func_val(xyz[1], xyz[2], xyz[3], θ_ext, ϕ_ext)
    end

    y_temp_coord = As * x_local_coord
    y_global_coord = zeros(Float64, gnpoin)
    for ip = 1:n
        y_global_coord[ip2gip_spa[ip]] += y_temp_coord[ip]
    end
    for k = 1:n_ext
        y_global_coord[extended_parents_to_gid_spa[k]] += y_temp_coord[n + k]
    end
    MPI.Allreduce!(y_global_coord, +, comm)
    @info "[$rank] Coord-based matvec AllReduce done"

    # Emit (coord, y_value) for OWNED DOFs only (avoids ghost duplicates)
    local_mat_y = Tuple{NTuple{5,Float64}, Float64}[]
    for ip = 1:n
        spa_ip = div(ip - 1, n_ang) + 1
        mesh.gip2owner[spa_ip] == rank || continue
        c = dof_coords_spatial_amr(ip, mesh, extra_mesh, n_ang)
        push!(local_mat_y, (c, y_global_coord[ip2gip_spa[ip]]))
    end



    # ── Step 2: Gather to rank 0 only (avoids Allgatherv memory explosion) ──────
    # Each non-0 rank sends its serialized entries to rank 0.
    # Rank 0 receives sequentially — no large per-rank receive buffers elsewhere.
    # MPI.Allgather([sz], comm) is called unconditionally so ALL ranks stay in sync.
    # Uses Int64 for sizes to avoid Int32 overflow on large datasets.
    function gather_to_rank0(local_entries, tag::Int; skip::Bool=false)
        buf  = IOBuffer(); serialize(buf, local_entries); data = take!(buf)
        sz   = Int64(length(data))
        all_sz = MPI.Allgather([sz], comm)   # Vector{Int64}
        

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
    all_mat_y       = gather_to_rank0(local_mat_y, 43)

    # ── Step 3: compare on rank 0 ────────────────────────────────────────────
    if rank == 0
        ref = JLD2.jldopen(filename, "r") do f
            (RHS      = f["RHS"],
             y_matvec = haskey(f, "y_matvec") ? f["y_matvec"] : nothing,
             coords   = f["coords"],
             n        = f["n"])
        end
        ref_n = ref.n
        @info "Loaded serial reference: n=$ref_n"

        # Build serial RHS dict
        ref_rhs_dict = Dict{NTuple{5,Float64}, Float64}()
        n_serial_overwrites = 0
        for i = 1:ref_n
            c = (ref.coords[i,1], ref.coords[i,2], ref.coords[i,3],
                 ref.coords[i,4], ref.coords[i,5])
            if haskey(ref_rhs_dict, c)
                n_serial_overwrites += 1
                n_serial_overwrites <= 5 &&
                    @warn "serial RHS dict overwrite at coord $c: old=$(ref_rhs_dict[c]) new=$(ref.RHS[i])"
            end
            ref_rhs_dict[c] = ref.RHS[i]
        end
        n_serial_overwrites > 0 &&
            @warn "serial RHS dict: $n_serial_overwrites coord collisions"
        ref_y_matvec   = ref.y_matvec
        ref_mat_coords = ref_y_matvec !== nothing ? ref.coords : nothing
        ref = nothing; GC.gc()

        # Aggregate parallel RHS
        par_rhs_dict  = Dict{NTuple{5,Float64}, Float64}()
        par_rhs_count = Dict{NTuple{5,Float64}, Int}()
        if all_rhs_entries !== nothing
            for (coord, val) in all_rhs_entries
                par_rhs_dict[coord]  = get(par_rhs_dict, coord, 0.0) + val
                par_rhs_count[coord] = get(par_rhs_count, coord, 0) + 1
            end
        end
        n_dup = count(v -> v > 1, values(par_rhs_count))
        if n_dup > 0
            @warn "parallel RHS dict: $n_dup coords with >1 contribution"
            shown = 0
            for (coord, cnt) in par_rhs_count
                cnt > 1 || continue; shown += 1; shown > 10 && break
                @info "  dup coord $coord: count=$cnt val=$(par_rhs_dict[coord])"
            end
        else
            @info "parallel RHS dict: no duplicate coord keys (good)"
        end

        # ── Compare RHS ──────────────────────────────────────────────────────
        ref_rhs_dict_free = Dict{NTuple{5,Float64}, Float64}()
        for (coord, val) in ref_rhs_dict
            abs(val) > 1e-20 && (ref_rhs_dict_free[coord] = val)
        end

        @info "=== SPATIAL-AMR RHS COMPARISON ==="
        all_rhs_coords = union(keys(ref_rhs_dict_free), keys(par_rhs_dict))
        rhs_ndiff = 0; rhs_max_abs = 0.0; rhs_max_rel = 0.0
        rhs_nmiss_par = 0; rhs_nextra_par = 0
        for coord in all_rhs_coords
            sv = get(ref_rhs_dict_free, coord, 0.0)
            pv = get(par_rhs_dict, coord, 0.0)
            ae = abs(pv - sv)
            if ae > 1e-9
                rhs_ndiff += 1
                rhs_max_abs = max(rhs_max_abs, ae)
                abs(sv) > 1e-15 && (rhs_max_rel = max(rhs_max_rel, ae / abs(sv)))
                rhs_ndiff <= 20 && @info "  RHS diff at $coord: serial=$sv parallel=$pv err=$ae"
            end
            !haskey(par_rhs_dict, coord) && haskey(ref_rhs_dict_free, coord) && (rhs_nmiss_par += 1)
            !haskey(ref_rhs_dict_free, coord) && haskey(par_rhs_dict, coord) && (rhs_nextra_par += 1)
        end
        @info "RHS totals: serial_nonzero=$(length(ref_rhs_dict_free)), parallel=$(length(par_rhs_dict))"
        @info "RHS: missing_from_parallel=$rhs_nmiss_par, extra_in_parallel=$rhs_nextra_par, mismatches=$rhs_ndiff"
        @info "RHS: max_abs_err=$rhs_max_abs, max_rel_err=$rhs_max_rel"
        rhs_ok = (rhs_ndiff == 0 && rhs_nmiss_par == 0 && rhs_nextra_par == 0)
        rhs_ok ? @info("✓ RHS MATCHES SERIAL") : @warn("✗ RHS DOES NOT MATCH SERIAL")

        ref_rhs_dict = nothing; par_rhs_dict = nothing; par_rhs_count = nothing
        ref_rhs_dict_free = nothing; all_rhs_coords = nothing
        all_rhs_entries = nothing; GC.gc()

        # ── Matvec-based matrix comparison (coordinate-based) ────────────────
        mat_ok = true
        if ref_y_matvec === nothing || ref_mat_coords === nothing
            @warn "=== MATRIX COMPARISON SKIPPED: no y_matvec/coords in reference file (re-run serial to regenerate) ==="
        else
            @info "=== SPATIAL-AMR MATRIX COMPARISON (coord-based matvec) ==="

            # Build serial coord → y_matvec dict
            ser_mat_dict = Dict{NTuple{5,Float64}, Float64}()
            n_ser_mat = length(ref_y_matvec)
            n_ser_coll = 0
            for i = 1:n_ser_mat
                c = (ref_mat_coords[i,1], ref_mat_coords[i,2], ref_mat_coords[i,3],
                     ref_mat_coords[i,4], ref_mat_coords[i,5])
                if haskey(ser_mat_dict, c)
                    n_ser_coll += 1
                else
                    ser_mat_dict[c] = ref_y_matvec[i]
                end
            end
            n_ser_coll > 0 && @warn "serial mat dict: $n_ser_coll coord collisions (skipped)"

            # Build parallel coord → y_matvec dict (sum contributions from all owned DOFs)
            par_mat_dict  = Dict{NTuple{5,Float64}, Float64}()
            par_mat_count = Dict{NTuple{5,Float64}, Int}()
            for (coord, val) in all_mat_y
                par_mat_dict[coord]  = get(par_mat_dict, coord, 0.0) + val
                par_mat_count[coord] = get(par_mat_count, coord, 0) + 1
            end
            n_dup_mat = count(v -> v > 1, values(par_mat_count))
            n_dup_mat > 0 && @warn "parallel mat dict: $n_dup_mat coords with >1 contribution (summed)"

            # Compare by coordinate key
            mat_ndiff = 0; mat_max_abs = 0.0; mat_max_rel = 0.0
            mat_nmiss_par = 0; mat_nextra_par = 0
            mat_ndiff_hang = 0; mat_ndiff_par_lt = 0; mat_ndiff_par_gt = 0
            n_nonzero_ser = 0

            # Build hanging node coord set for classification
            hang_coords = Set{NTuple{5,Float64}}()
            for hn_ip in hanging_nodes
                c = dof_coords_spatial_amr(hn_ip, mesh, extra_mesh, n_ang)
                push!(hang_coords, c)
            end

            all_mat_coords = union(keys(ser_mat_dict), keys(par_mat_dict))
            for coord in all_mat_coords
                sv = get(ser_mat_dict, coord, 0.0)
                pv = get(par_mat_dict, coord, 0.0)
                abs(sv) > 1e-20 && (n_nonzero_ser += 1)
                ae = abs(pv - sv)
                if ae > 1e-6
                    mat_ndiff += 1
                    mat_max_abs = max(mat_max_abs, ae)
                    abs(sv) > 1e-15 && (mat_max_rel = max(mat_max_rel, ae / abs(sv)))
                    is_hang = coord in hang_coords
                    if pv < sv; mat_ndiff_par_lt += 1
                    else;        mat_ndiff_par_gt += 1; end
                    is_hang && (mat_ndiff_hang += 1)
                    if mat_ndiff <= 20
                        @info "  A*x diff at coord=$coord: serial=$sv parallel=$pv err=$ae hang=$is_hang"
                    end
                end
                !haskey(par_mat_dict, coord) && haskey(ser_mat_dict, coord) && (mat_nmiss_par += 1)
                !haskey(ser_mat_dict, coord) && haskey(par_mat_dict, coord) && (mat_nextra_par += 1)
            end

            @info "A*x totals: serial_nonzero=$n_nonzero_ser, serial_dofs=$(length(ser_mat_dict)), parallel_dofs=$(length(par_mat_dict))"
            @info "A*x: missing_from_parallel=$mat_nmiss_par, extra_in_parallel=$mat_nextra_par, mismatches=$mat_ndiff"
            @info "A*x mismatches: hanging=$mat_ndiff_hang, par<ser=$mat_ndiff_par_lt, par>ser=$mat_ndiff_par_gt"
            @info "A*x: max_abs_err=$mat_max_abs, max_rel_err=$mat_max_rel"
            mat_ok = (mat_ndiff == 0 && mat_nmiss_par == 0 && mat_nextra_par == 0)
            mat_ok ? @info("✓ MATRIX MATVEC MATCHES SERIAL") : @warn("✗ MATRIX MATVEC DOES NOT MATCH SERIAL")

            ref_mat_coords = nothing; ref_y_matvec = nothing
            ser_mat_dict = nothing; par_mat_dict = nothing; par_mat_count = nothing
            all_mat_y = nothing; hang_coords = nothing; GC.gc()
        end

        return rhs_ok && mat_ok
    end

    return nothing
end
