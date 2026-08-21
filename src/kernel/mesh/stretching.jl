function stretch_mesh!(mesh,inputs,npoin)

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)
    
    stretch_factor = inputs[:stretch_factor]
    ztop = mesh.ymax

    zsurf = zeros(Float64,npoin)
    
    for ip = 1:npoin
        stretching_factor = inputs[:stretch_factor]
        sigma = mesh.y[ip]

        # Normalize the sigma coordinate to the range [0, 1]
        sigma_normalized = sigma / ztop

        # Apply the non-linear power-law stretching function.
        # This maps the uniform [0, 1] grid to a stretched [0, 1] grid where points
        # are concentrated towards 0.
        z_normalized = sigma_normalized ^ stretching_factor

        # Map the stretched, normalized coordinate back to the physical z-domain,
        # which spans from the bottom surface zsurf[ip] to the domain top ztop.
        z = zsurf[ip] + z_normalized * (ztop - zsurf[ip])

        # Update the grid point's vertical position with the new stretched value.
        mesh.y[ip] = z
    end
    
end

function stretch_mesh_3D!(mesh,inputs, npoin)  

    comm = get_mpi_comm()
    rank = MPI.Comm_rank(comm)
    mpi_size = MPI.Comm_size(comm)
    
    stretch_factor = inputs[:stretch_factor]
    ztop = mesh.zmax
    
    zsurf = zeros(Float64,npoin)

    if inputs[:stretch_type] == "powerlaw"
        for ip = 1:npoin
            stretching_factor = inputs[:stretch_factor]
            sigma = mesh.coords[3,ip]
            
            # Normalize the sigma coordinate to the range [0, 1]
            sigma_normalized = sigma / ztop
            
            # Apply the non-linear power-law stretching function.
            # This maps the uniform [0, 1] grid to a stretched [0, 1] grid where points
            # are concentrated towards 0.
            z_normalized = sigma_normalized ^ stretching_factor
            
            # Map the stretched, normalized coordinate back to the physical z-domain,
            # which spans from the bottom surface zsurf[ip] to the domain top ztop.
            z = zsurf[ip] + z_normalized * (ztop - zsurf[ip])
            
            # Update the grid point's vertical position with the new stretched value.
            mesh.coords[3,ip] = z
        end
        
    elseif inputs[:stretch_type] == "fixed_first"
        
        # NEW INPUT: Define the desired thickness of the first grid cell from the surface.
        first_cell_size = inputs[:first_zelement_size] # For example, 1.0 meter

        #Scale it up by the order because then LGL points will scale it down again:
        first_cell_size *= (mesh.ngl-1)
        
        # --- Pre-computation Step ---

        # Store a copy of the original, uniform z-coordinates to use as the sigma field.
        # This is crucial because we will be overwriting mesh.z in the loop.
        sigma_coords = copy(mesh.coords[3,:])

        # Find the maximum height of the domain.
        ztop = mesh.zmax #maximum(sigma_coords)

        # Find the z-coordinate of the first layer in the original uniform grid (sigma_1).
        # This is the smallest positive value in the original grid coordinates.
        # We filter for > 0 to exclude any points at the z=0 boundary.
        sigma_1 = minimum(filter(x -> x > 0, sigma_coords))

        # Calculate the required stretching factor 'n' to satisfy the first_cell_size condition.
        # This is derived from the power-law equation: h = (σ/z_top)^n * z_top
        # Solving for n gives: n = log(h/z_top) / log(σ/z_top)
        computed_stretch_factor = log(first_cell_size / ztop) / log(sigma_1 / ztop)

        println("Desired resolution by the surface: ", first_cell_size/(mesh.ngl-1))
        println("Computed stretching factor: ", computed_stretch_factor)

        # --- Grid Transformation Step ---

        # Initialize the surface elevation array (assuming it's flat at z=0 for this problem).
        zsurf = zeros(Float64, npoin)

        for ip = 1:npoin
            # Get the original uniform coordinate for this point.
            sigma = sigma_coords[ip]

            if sigma == 0.0
                # A point at the bottom surface remains at the bottom surface.
                mesh.coords[3,ip] = zsurf[ip]
                continue
            end

            # Normalize the sigma coordinate to the range [0, 1].
            sigma_normalized = sigma / ztop

            # Apply the non-linear power-law stretching using the new COMPUTED factor.
            z_normalized = sigma_normalized ^ computed_stretch_factor

            # Map the stretched, normalized coordinate back to the physical z-domain.
            z = zsurf[ip] + z_normalized * (ztop - zsurf[ip])

            # Update the grid point's vertical position.
            mesh.coords[3,ip] = z
        end

    elseif inputs[:stretch_type] == "fixed_first_twoblocks_weak"
        
        # --- User Inputs ---
        # Desired thickness of the first grid cell from the surface.
        first_cell_size = inputs[:first_zelement_size] # e.g., 1.0 meter

        #Scale it up by the order because then LGL points will scale it down again:
        first_cell_size *= (mesh.ngl-1)
        
        # The z-coordinate where the stretching stops and the uniform grid begins.
        zlevel_transition = inputs[:zlevel_transition]  # e.g., 200.0 meters

        # --- Pre-computation Step ---

        # Store a copy of the original, uniform z-coordinates to use as the sigma field.
        sigma_coords = copy(mesh.coords[3,:])

        
        # 1. Initialize local minimum to positive infinity.
        # This is the identity element for the 'min' operation and elegantly
        # handles processes with no positive values.
        local_min_sigma = Inf

        # 2. Filter the local array for positive values.
        local_positive_coords = filter(x -> x > 0, sigma_coords)

        # 3. If any positive values were found locally, compute the local minimum.
        # Otherwise, local_min_sigma remains Inf.
        if !isempty(local_positive_coords)
            local_min_sigma = minimum(local_positive_coords)
        end
        
        #sigma_1 = minimum(filter(x -> x > 0, sigma_coords)) # Represents the uniform dσ
        # 4. Perform the global reduction to find the true minimum across all processes.
        # MPI.Allreduce sends the result to all processes.
        sigma_1 = MPI.Allreduce(local_min_sigma, MPI.MIN, comm)
        
        # Find the maximum height and the first layer's height in the original grid.
        ztop = mesh.zmax
      

        # 1. CALCULATE STRETCHING FACTOR for the bottom region.
        # This is the same calculation as before, ensuring the bottom layer meets the size criteria.
        computed_stretch_factor = log(first_cell_size / ztop) / log(sigma_1 / ztop)
        n = computed_stretch_factor # Use 'n' for brevity in formulas

        println("Desired resolution by the surface: ", first_cell_size/(mesh.ngl-1))
        println("Transition z-level: ", zlevel_transition)
        println("Computed stretching factor for bottom region: ", n)

        # 2. CALCULATE PROPERTIES AT THE TRANSITION POINT.
        # We need to find the point in the original grid (sigma_transition) that maps to zlevel_transition
        # in the new stretched grid, and the cell size (dz_transition) at that point.

        # Invert the stretching formula to find sigma_transition: z = ztop * (σ/ztop)^n -> σ = ztop * (z/ztop)^(1/n)
        sigma_transition = ztop * (zlevel_transition / ztop)^(1/n)

        # Calculate the cell size at the transition point using the derivative of the stretch function.
        # dz/dσ = n * (σ/ztop)^(n-1). So, dz = dσ * (dz/dσ). We use sigma_1 as dσ.
        dz_transition = (n * (sigma_transition / ztop)^(n-1)) * sigma_1

        println("Original sigma coordinate at transition: ", sigma_transition)
        println("Cell size at transition (for uniform top grid): ", dz_transition)

        # --- Grid Transformation Step ---

        # Initialize the surface elevation array.
        zsurf = zeros(Float64, npoin)

        for ip = 1:npoin
            # Get the original uniform coordinate for this point.
            sigma = sigma_coords[ip]

            # Handle the point at the surface.
            if sigma == 0.0
                mesh.coords[3,ip] = zsurf[ip]
                continue
            end

            local z::Float64 # Declare z locally for the if/else block

            # Apply transformation based on whether the point is in the bottom or top region.
            if sigma <= sigma_transition
                # --- REGION 1: STRETCHED GRID ---
                # Apply the standard power-law stretching.
                sigma_normalized = sigma / ztop
                z_normalized = sigma_normalized^n
                z = zsurf[ip] + z_normalized * (ztop - zsurf[ip])
            else
                # --- REGION 2: UNIFORM GRID ---
                # The position is zlevel_transition plus a number of uniform cells of size dz_transition.
                # Number of cells above transition = (distance in sigma) / (sigma cell size)
                num_cells_above = (sigma - sigma_transition) / sigma_1
                height_above_transition = num_cells_above * dz_transition
                z = zlevel_transition + height_above_transition
            end

            # Update the grid point's vertical position.
            mesh.coords[3,ip] = z
        end
        
    elseif inputs[:stretch_type] == "fixed_first_twoblocks_strong"

        # --- User Inputs ---
        # Desired thickness of the first grid cell from the surface.
        first_cell_size = inputs[:first_zelement_size] # e.g., 1.0 meter

        #Scale it up by the order because then LGL points will scale it down again:
        first_cell_size *= (mesh.ngl-1)
        
        # The z-coordinate where the stretching stops and the uniform grid begins.
        zlevel_transition = inputs[:zlevel_transition]  # e.g., 200.0 meters

        # --- Pre-computation Step ---

        # Store a copy of the original, uniform z-coordinates to use as the sigma field.
        sigma_coords = copy(mesh.coords[3,:])

        # Find the maximum height and the first layer's height in the original grid.
        ztop = mesh.zmax

        
        # 1. Initialize local minimum to positive infinity.
        # This is the identity element for the 'min' operation and elegantly
        # handles processes with no positive values.
        local_min_sigma = Inf

        # 2. Filter the local array for positive values.
        local_positive_coords = filter(x -> x > 0, sigma_coords)

        # 3. If any positive values were found locally, compute the local minimum.
        # Otherwise, local_min_sigma remains Inf.
        if !isempty(local_positive_coords)
            local_min_sigma = minimum(local_positive_coords)
        end
        
        #sigma_1 = minimum(filter(x -> x > 0, sigma_coords)) # Represents the uniform dσ
        # 4. Perform the global reduction to find the true minimum across all processes.
        # MPI.Allreduce sends the result to all processes.
        sigma_1 = MPI.Allreduce(local_min_sigma, MPI.MIN, comm)

        
        # 1. CALCULATE STRETCHING FACTOR for the bottom region.
        n = log(first_cell_size / ztop) / log(sigma_1 / ztop)
        println("Desired resolution by the surface: ", first_cell_size/(mesh.ngl-1))
        println("Transition z-level: ", zlevel_transition)
        println("Computed stretching factor for bottom region: ", n)

        # 2. CALCULATE THE TRANSITION POINT IN THE ORIGINAL (SIGMA) GRID.
        # This finds which original coordinate maps to our desired physical transition level.
        sigma_transition = ztop * (zlevel_transition / ztop)^(1/n)
        println("Original sigma coordinate at transition: ", sigma_transition)

        # --- Grid Transformation Step ---

        # Initialize the surface elevation array.
        zsurf = zeros(Float64, npoin)

        for ip = 1:npoin
            # Get the original uniform coordinate for this point.
            sigma = sigma_coords[ip]

            # Handle the point at the surface.
            if sigma == 0.0
                mesh.coords[3,ip] = zsurf[ip]
                continue
            end

            local z::Float64

            if sigma <= sigma_transition
                # --- REGION 1: STRETCHED GRID (Bottom) ---
                # This part is unchanged. It applies the power-law stretch.
                sigma_normalized = sigma / ztop
                z_normalized = sigma_normalized^n
                z = zsurf[ip] + z_normalized * (ztop - zsurf[ip])
            else
                # --- REGION 2: LINEARLY MAPPED GRID (Top) ---
                # This new logic ensures the grid perfectly hits ztop.
                # It performs a linear mapping from the remaining sigma-space to the remaining physical space.

                # Find the relative position of the point within the top section of the *original* grid.
                # This value (alpha) goes from 0 to 1 as sigma goes from sigma_transition to ztop.
                relative_pos = (sigma - sigma_transition) / (ztop - sigma_transition)

                # Apply that same relative position to the top section of the *new* physical grid.
                # The remaining physical space is from zlevel_transition to ztop.
                z = zlevel_transition + relative_pos * (ztop - zlevel_transition)
            end

            # Update the grid point's vertical position.
            mesh.coords[3,ip] = z
        end

    elseif inputs[:stretch_type] == "two_block uniformish"
        # --- User Inputs ---
        # Desired thickness of the first grid cell from the surface.
        first_cell_size = inputs[:first_zelement_size] # e.g., 1.0 meter

        #Scale it up by the order because then LGL points will scale it down again:
        first_cell_size *= (mesh.ngl-1)

        # The z-coordinate where the stretching stops and the uniform grid begins.
        zlevel_transition = inputs[:zlevel_transition]  # e.g., 200.0 meters
        
        # Store a copy of the original, uniform z-coordinates to use as the sigma field.
        sigma = copy(mesh.coords[3,:])

        # Find the maximum height and the first layer's height in the original grid.
        ztop = mesh.zmax

        #compute existing general cell size
        ip1 = mesh.connijk[1,1,1,1]
        ip2 = mesh.connijk[1,1,1,mesh.ngl]
        cell_size = abs(mesh.coords[3,ip1]-mesh.coords[3,ip2])
        
        #find number of vertical elements to be worked with
        total_z_elems = round(ztop/cell_size)
        
        #Determine if enough elements are present to construct uniform high resolution area
        nelem_squashed = round(zlevel_transition/first_cell_size)
    
        nelem_stretched = total_z_elems-nelem_squashed
        #compute squashing ratio:
        squash_ratio = first_cell_size/cell_size

        #Find stretched element cell sizes
        stretch_cell_size = (ztop-zlevel_transition)/nelem_stretched

        if (nelem_squashed > (9/10) * total_z_elems) || nelem_stretched < 5
            error("too few elements available for desired resolution and transition height combination")
            return nothing
        end

        #find stretching ratio

        stretch_ratio = stretch_cell_size/cell_size

        #Use half the stretched layer to transition between mesh sizes using half of the coarser area. 
        
        nelem_transition = round(nelem_stretched/2)
        #find new number of stretched elements
        nelem_stretched -= nelem_transition
        #recompute new stretched cell size
        transition_height = nelem_squashed*first_cell_size
        cumul_height = 0.0
        for e_t=1:nelem_transition
            ratio = squash_ratio + ((e_t)/nelem_transition)*(stretch_ratio-squash_ratio)
            transition_height += cell_size*ratio
            cumul_height += cell_size*ratio
        end
        stretch_ratio_old = stretch_ratio
        stretch_cell_size = (ztop-transition_height)/(nelem_stretched)
        stretch_ratio = stretch_cell_size/cell_size
        
        for ie=1:mesh.nelem
            ip1 = mesh.connijk[ie,1,1,1]
            ip2 = mesh.connijk[ie,1,1,mesh.ngl]
            zelem_extrema1 = sigma[ip1]
            zelem_extrema2 = sigma[ip2]
            zelem_max = max(zelem_extrema1, zelem_extrema2)
            zelem_min = min(zelem_extrema1, zelem_extrema2)
            elem_order = round(zelem_max/cell_size)
            if (elem_order <= nelem_squashed)
                for i=1:mesh.ngl
                    for j=1:mesh.ngl
                        for k=1:mesh.ngl
                            ip = mesh.connijk[ie,i,j,k]
                            dz = abs(sigma[ip] - zelem_min)
                            z = (elem_order-1)*first_cell_size + dz * squash_ratio
                            mesh.coords[3,ip] = z
                        end
                    end
                end
            elseif (elem_order <= (nelem_squashed + nelem_transition))
                for i=1:mesh.ngl
                    for j=1:mesh.ngl
                        for k=1:mesh.ngl
                            ip = mesh.connijk[ie,i,j,k]
                            ratio = squash_ratio + ((elem_order-nelem_squashed)/nelem_transition)*(stretch_ratio_old-squash_ratio)
                            cumul_height = 0.0
                            for e_t = 1:elem_order-nelem_squashed-1
                                ratio_temp = squash_ratio + ((e_t)/nelem_transition)*(stretch_ratio_old-squash_ratio)
                                cumul_height += cell_size * ratio_temp
                            end
                            dz = abs(sigma[ip] - zelem_min)
                            z = (nelem_squashed)*first_cell_size + cumul_height + dz*ratio 
                            mesh.coords[3,ip] = z
                        end
                    end
                end
            else
                for i=1:mesh.ngl
                    for j=1:mesh.ngl
                        for k=1:mesh.ngl
                            ip = mesh.connijk[ie,i,j,k]
                            dz = abs(sigma[ip] - zelem_min)
                            z = transition_height + (elem_order-nelem_squashed-nelem_transition - 1)*stretch_cell_size  + dz * stretch_ratio
                            mesh.coords[3,ip] = z
                        end
                    end
                end
            end
        end
    end
    
    
end
