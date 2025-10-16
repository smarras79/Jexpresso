# Import necessary packages
using Roots

function warp_mesh!(mesh,inputs)
    am = inputs[:a_mount]
    hm = inputs[:h_mount]
    xc = inputs[:c_mount]
    ztop = maximum(mesh.y)
    zsurf = zeros(Float64,mesh.npoin)
    sigma = zeros(Float64,mesh.npoin)
    for i=1:mesh.npoin
        sigma[i] = mesh.z[i]
        
        z = (ztop - zsurf[i])/ztop * sigma[i] + zsurf[i]
        
        mesh.z[i] = z
    end
    if (inputs[:mount_type] == "agnesi")
        am = inputs[:a_mount]
        hm = inputs[:h_mount]
        xc = inputs[:c_mount]
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = hm*am*am/((x-xc)*(x-xc) + am*am)
            #zsurf[ip] = hm/(1+ ((x-xc)/am)^2)
        end
    elseif (inputs[:mount_type] == "schar")
        ac = inputs[:a_mount]
        hc = inputs[:h_mount]
        lambdac = inputs[:lambda_mount]
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = hc * exp(-(x/ac)^2) * cospi(x/lambdac)^2 
        end
    end
    
    if inputs[:lstretch]
        for ip = 1:mesh.npoin
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
    else
        for ip = 1:mesh.npoin
            sigma[ip] = mesh.y[ip]
            #if (mesh.y[ip] < 10000.0)  
            z = (ztop - zsurf[ip])/ztop * sigma[ip] + zsurf[ip]
            mesh.y[ip] = z
            #=elseif (mesh.y[ip] < 15000.0)
            factor = (15000-mesh.y[ip])/5000.0
            z = (ztop - factor*zsurf[ip])/ztop * sigma[ip] + factor*zsurf[ip]
            mesh.y[ip] = z 
            end=#
        end
    end
    
    
    #=for iedge = 1:size(mesh.bdy_edge_in_elem,1)
    iel = mesh.bdy_edge_in_elem[iedge]
    comp = mesh.bdy_edge_comp[iedge]
    for k=1:mesh.ngl
    #if (mesh.bdy_edge_type[iedge] == "free_slip")
    tag = mesh.bdy_edge_type[iedge]
    ip = mesh.poin_in_bdy_edge[iedge,k]
    @info "prewarp", mesh.y[ip]
    if (mesh.y[ip] < 10.0)
    mesh.y[ip] = (hm*(am^2))/((mesh.x[ip]-xc)^2 + am^2)
    end
    @info "postwarp", mesh.y[ip]
    #end
    end
    end=#
end

function warp_mesh_3D!(mesh,inputs)
    if (inputs[:mount_type] == "real topography")
        # find surface heights by reading and interpolating real data onto grid
        fname = inputs[:topo_database]
        fname2 = inputs[:topo_geoid]
        lat_min = inputs[:read_topo_latmin]
        lat_max = inputs[:read_topo_latmax]
        lon_min = inputs[:read_topo_lonmin]
        lon_max = inputs[:read_topo_lonmax]
        zone = inputs[:read_topo_zone]
        xmin = mesh.xmax
        xmax = mesh.xmax
        ymin = mesh.ymax
        ymax = mesh.ymax
        lat, lon, z_topo = extract_region_topography_from_global_data(fname, fname2, lat_max, lon_max, lat_min, lon_min)
        
        x_topo, y_topo = Map_lat_lon_onto_simulation_domain(lat,lon,xmin,xmax,ymin,ymax,zone)
        zsurf = zeros(mesh.npoin)
        
        interpolate_topography_onto_grid!(mesh.x, mesh.y, zsurf, x_topo, y_topo, z_topo)
        ### sigma coordinate topography
        ztop = mesh.zmax
        sigma = zeros(mesh.npoin)
        for ip = 1:mesh.npoin
            sigma[ip] = mesh.z[ip]
            z_new = (ztop - zsurf[ip])/ztop * sigma[ip] + zsurf[ip]
            mesh.z[ip] = z_new
        end

    elseif (inputs[:mount_type] == "agnesi")
        zsurf = zeros(mesh.npoin)
        sigma = zeros(mesh.npoin)
        ztop = mesh.zmax
        am = inputs[:a_mount]
        hm = inputs[:h_mount]
        xc = inputs[:c_mount]
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = hm*am*am/((x-xc)*(x-xc) + am*am)
        end
        
    elseif (inputs[:mount_type] == "sauer")
        zsurf = zeros(mesh.npoin)
        sigma = zeros(mesh.npoin)
        ztop = mesh.zmax
        #am = mesh.xmax - mesh.xmin
	am = inputs[:a_mount]
      	hm = inputs[:h_mount]
        xc = inputs[:c_mount]
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            ####zsurf[ip] = hm*(sech(x - xc)/am)^2 #not working
            zsurf[ip] = 0.5*hm*(1.0 - cospi(2.0*(x)/am))
        end
        
    elseif (inputs[:mount_type] == "LESICP")
        zsurf = zeros(mesh.npoin)
        sigma = zeros(mesh.npoin)
        ztop = mesh.zmax
        am = mesh.xmax - mesh.xmin
      	hm = inputs[:h_mount]
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = 0.5*hm*(1.0 - cospi(2.0*x/am))
        end
                
    elseif (inputs[:mount_type] == "schar")
        ac = inputs[:a_mount]
        hc = inputs[:h_mount]
        lambdac = inputs[:lambda_mount]
        for ip = 1:mesh.npoin
            x = mesh.x[ip]
            zsurf[ip] = hc * exp(-(x/ac)^2) * cospi(x/lambdac)^2
        end

    #=elseif (inputs[:mount_type] == "stretching" || inputs[:mount_type] == "bl")

         Yuo need to do detect the vertically aligned nodes to move. 
        zsurf = zeros(mesh.npoin)
        sigma = zeros(mesh.npoin)
        ztop = mesh.zmax
        
	dz1 = inputs[:a_mount]

        stretching!(zsurf, ztop, dz1, mesh.npoin)
       =#
    end

    for ip = 1:mesh.npoin
        sigma[ip] = mesh.z[ip]
        z = (ztop - zsurf[ip])/ztop * sigma[ip] + zsurf[ip]
        mesh.z[ip] = z
    end
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
