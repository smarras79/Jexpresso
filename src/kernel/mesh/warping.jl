using Roots  # (Assuming needed elsewhere; not used here)

function warp_mesh!(mesh, inputs)
    """
    Warp a 2D mesh to follow terrain
    GMSH 2D mesh has: x (horizontal), y (vertical)
    """
    
    am = inputs[:a_mount]
    hm = inputs[:h_mount]
    xc = inputs[:c_mount]
    
    # Get GLOBAL domain top (critical for parallel runs!)
    ytop_local = maximum(mesh.y)
    ytop = MPI.Allreduce(ytop_local, MPI.MAX, MPI.COMM_WORLD)  # Global max
    
    # Compute surface elevation
    ysurf = zeros(Float64, mesh.npoin)
    
    if inputs[:mount_type] == "gauss"
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            ysurf[ip] = hm * exp(-((x - xc)^2 / (am^2)))
        end
    elseif inputs[:mount_type] == "agnesi"
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            ysurf[ip] = hm * am * am / ((x - xc) * (x - xc) + am * am)
        end
    else
        ysurf .= 0.0
    end
    
    # Apply terrain-following transformation
    for ip = 1:mesh.npoin
        sigma = mesh.y[ip]
        y_new = ysurf[ip] + (sigma / ytop) * (ytop - ysurf[ip])
        mesh.y[ip] = y_new
    end
    
    # Report statistics (each processor reports its local max)
    max_surf_local = maximum(ysurf)
    @info "Mesh warping complete:"
    @info "  Mountain type: $(inputs[:mount_type])"
    @info "  Mountain height: $(hm) m"
    @info "  Mountain width: $(2*am) m"
    @info "  Max surface elevation (local): $(max_surf_local) m"
    @info "  Domain top (global): $(ytop) m"
    
    return nothing
end







function warp_mesh_3D!(mesh, inputs)
    # Hoist common params (with defaults if missing)
    am = get(inputs, :a_mount, mesh.xmax - mesh.xmin)
    hm = get(inputs, :h_mount, 0.0)
    xc = get(inputs, :c_mount, 0.0)
    ztop = mesh.zmax
    zsurf = zeros(Float64, mesh.npoin)
    sigma = zeros(Float64, mesh.npoin)  # Always init
    
    if inputs[:mount_type] == "real topography"
        # Find surface heights by reading and interpolating real data onto grid
        fname = inputs[:topo_database]
        fname2 = inputs[:topo_geoid]
        lat_min = inputs[:read_topo_latmin]
        lat_max = inputs[:read_topo_latmax]
        lon_min = inputs[:read_topo_lonmin]
        lon_max = inputs[:read_topo_lonmax]
        zone = inputs[:read_topo_zone]
        xmin = mesh.xmin  # Fixed typo
        xmax = mesh.xmax
        ymin = mesh.ymin  # Fixed typo
        ymax = mesh.ymax
        lat, lon, z_topo = extract_region_topography_from_global_data(fname, fname2, lat_max, lon_max, lat_min, lon_min)
        
        x_topo, y_topo = Map_lat_lon_onto_simulation_domain(lat, lon, xmin, xmax, ymin, ymax, zone)
        
        interpolate_topography_onto_grid!(mesh.x, mesh.y, zsurf, x_topo, y_topo, z_topo)
        ### sigma coordinate topography
        for ip = 1:mesh.npoin
            sigma[ip] = mesh.z[ip]
            z_new = zsurf[ip] + (ztop - zsurf[ip]) / ztop * sigma[ip]  # Fixed order for consistency
            mesh.z[ip] = z_new
        end
        return nothing  # Skip damping for real topo

    elseif inputs[:mount_type] == "agnesi"
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = hm * am * am / ((x - xc) * (x - xc) + am * am)
        end
        
    elseif inputs[:mount_type] == "gauss"
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = hm * exp( -((x - xc)^2 / (am^2)) )  # Pure Gaussian
        end
        
    elseif inputs[:mount_type] == "sauer"
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = 0.5 * hm * (1.0 - cospi(2.0 * (x - xc) / am))  # Added xc shift for centering
        end
        
    elseif inputs[:mount_type] == "LESICP"
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = 0.5 * hm * (1.0 - cospi(2.0 * (x - xc) / am))  # Added xc shift
        end
                
    elseif inputs[:mount_type] == "schar"
        ac = inputs[:a_mount]
        hc = inputs[:h_mount]
        lambdac = inputs[:lambda_mount]
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = hc * exp(-(x/ac)^2) * cospi(x/lambdac)^2
        end
        
    else
        zsurf .= 0.0  # Flat
    end

    # Parameters for damping control
    z_transition_start = inputs[:z_transition_start]  # Height where damping starts (30% of domain)
    z_transition_end = inputs[:z_transition_end]    # Height where grid becomes fully flat (60% of domain)
    
    # Apply terrain-following warp with damping
    for ip = 1:mesh.npoin
        sigma[ip] = mesh.z[ip]
        
        # Original warped coordinate (follows topography)
        z_warped = zsurf[ip] + sigma[ip] * (ztop - zsurf[ip]) / ztop

        damping_factor = 1.0  # Default: full warping

        if sigma[ip] >= z_transition_start && sigma[ip] <= z_transition_end
            # Smooth transition between z_transition_start and z_transition_end
            progress = (sigma[ip] - z_transition_start) / (z_transition_end - z_transition_start)
            
            # Cosine (or 'Raised Cosine') Damping Function
            damping_factor = 0.5 * (1.0 + cospi(progress))
            
        elseif sigma[ip] > z_transition_end
            # No warping above transition end (fully flat)
            damping_factor = 0.0
        end  # Below start: full damping (1.0)

        mesh.z[ip] = sigma[ip] + damping_factor * (z_warped - sigma[ip])
    end
    
    return nothing
end

function warp_phys_grid!(x,y,z,ncol,nlay)
    if (inputs[:mount_type] == "real topography")
        # find surface heights by reading and interpolating real data onto grid
        fname = inputs[:topo_database]
        fname2 = inputs[:topo_geoid]
        lat_min = inputs[:read_topo_latmin]
        lat_max = inputs[:read_topo_latmax]
        lon_min = inputs[:read_topo_lonmin]
        lon_max = inputs[:read_topo_lonmax]
        zone = inputs[:read_topo_zone]
        xmin = minimum(x)
        xmax = maximum(x)
        ymin = minimum(y)
        ymax = maximum(y)
        lat, lon, z_topo = extract_region_topography_from_global_data(fname, fname2, lat_max, lon_max, lat_min, lon_min)

        x_topo, y_topo = Map_lat_lon_onto_simulation_domain(lat,lon,xmin,xmax,ymin,ymax,zone)
        zsurf = zeros(ncol)

        interpolate_topography_onto_grid!(x, y, zsurf, x_topo, y_topo, z_topo)
        ### sigma coordinate topography
        ztop = mesh.zmax
        sigma = zeros(nlay+1,ncol)
        for icol = 1:ncol
            for ilay = 1:nlay+1
                sigma[ilay,icol] = z[ilay,icol]
                z_new = (ztop - zsurf[icol])/ztop * sigma[ilay,icol] + zsurf[icol]
                z[ilay,icol] = z_new
            end
        end

    elseif (inputs[:mount_type] == "agnesi")
        zsurf = zeros(ncol)
        sigma = zeros(nlay+1,ncol)
        ztop = maximum(z)
        am = inputs[:a_mount]
        hm = inputs[:h_mount]
        xc = inputs[:c_mount]
        for icol = 1:ncol
            xx = x[icol]
            zsurf[icol] = hm/(1+ ((xx-xc)/am)^2)
        end
    elseif (inputs[:mount_type] == "schar")
        ac = inputs[:a_mount]
        hc = inputs[:h_mount]
        lambdac = inputs[:lambda_mount]
        for icol = 1:ncol
            xx = x[icol]
            zsurf[icol] = hc * exp(-(xx/ac)^2) * cospi(xx/lambdac)^2
        end
    end

    for icol = 1:ncol
        for ilay = 1:nlay+1
            sigma[ilay,icol] = z[ilay,icol]
            z_new = (ztop - zsurf[icol])/ztop * sigma[ilay,icol] + zsurf[icol]
            z[ilay,icol] = z_new
        end
    end
end

# --- Wrap all logic in a main function to ensure correct variable scope ---
function stretching!(z, zmax, dz1, N)
    # --- 1. Define Grid Parameters ---
    L = zmax  # Domain height in meters
    #dz1 = 10.0  # Desired spacing of the first layer (m)
    #N = 16      # Total number of grid points

    println("--- Grid Setup ---")
    println("Domain Height (L): $L m")
    println("First Spacing (dz1): $dz1 m")
    println("Number of Points (N): $N")
    println("--------------------")

    # Define the function for the root-finding algorithm
    # This function "captures" L, dz1, and N from the main function's scope
    f(r) = dz1 * (r^(N - 1) - 1) / (r - 1) - L

    # --- 2. Robustly find a bracketing interval [a, b] ---
    #
    # THE FIX IS HERE: Start 'a' slightly above 1.0 to avoid 0/0 division.
    # nextfloat(1.0) gets the smallest floating point number greater than 1.0.
    #
    a, b = nextfloat(1.0), 1.1

    # This loop will now correctly modify 'b' because it's inside a function
    while f(a) * f(b) > 0
        b += 0.1
        if b > 5.0 # Failsafe to prevent an infinite loop
            error("Could not find a bracketing interval. Check parameters.")
        end
    end

    println("Found a valid search bracket: ($a, $b)")

    # Numerically solve for 'r' using the dynamically found bracket
    r = find_zero(f, (a, b))

    println("✅ Calculated Stretching Ratio (r): $r")

    # --- 3. Generate the Grid Points ---
    #z = zeros(N)
    z[1] = zmin

    # The sum of spacings is the position of the next point.
    # The spacing itself is dz_i = dz1 * r^(i-1)
    # The position is z_k = sum_{i=1}^{k-1} dz_i
    for i in 2:N
        # More direct formula for position based on geometric series sum
        z[i] = dz1 * (r^(i - 1) - 1) / (r - 1)
    end

    # --- 4. Verification and Output ---
    println("\n--- Grid Verification ---")
    println("Top grid point: ", round(z[end], digits=2), " m (should be close to $L)")
    first_spacing = z[2] - z[1]
    last_spacing = z[end] - z[end-1]
    println("First grid spacing: ", round(first_spacing, digits=2), " m")
    println("Last grid spacing: ", round(last_spacing, digits=2), " m")

    # --- 5. Visualization ---
    p = plot(z, 0:N-1, # Plot z on x-axis, index on y-axis for vertical look
             seriestype = :scatter,
             marker = :hline,
             title = "Stretched 1D Vertical Grid",
             xlabel = "Height (z) [m]",
             ylabel = "Grid Point Index (i)",
             label = "Grid Levels",
             yflip=false, # Often vertical grids are plotted 0 at bottom
             legend = :topleft)
    display(p)
    println("\n📈 Plot generated successfully.")

end
