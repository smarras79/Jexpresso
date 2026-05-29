using NCDatasets

"""
    read_atmospheric_data(filename::String)

Read atmospheric data from a NetCDF file and extract coordinates and variables.
Calculates cloud liquid and ice mixing ratios from water paths and moist air density.

Returns a named tuple containing:
- x, y, z: coordinate arrays
- t_lay, p_lay, vmr_h2o: atmospheric variables
- q_liq: liquid water mixing ratio (kg/kg)
- q_ice: ice water mixing ratio (kg/kg)
- rho: moist air density (kg/m³)
"""
function read_atmospheric_data(filename::String)
    # Open the NetCDF file
    @info filename
    ds = NCDataset(filename, "r")
    
    try
        # Extract coordinate arrays
        x = Array(ds["x"])
        y = Array(ds["y"])
        z = Array(ds["z"])
        
        # Extract atmospheric variables
        t_lay = Array(ds["t_lay"])
        p_lay = Array(ds["p_lay"])
        vmr_h2o = Array(ds["vmr_h2o"])
        lwp = Array(ds["lwp"])  # liquid water path (kg/m²)
        iwp = Array(ds["iwp"])  # ice water path (kg/m²)
        
        # Calculate moist air density
        rho = calculate_moist_density(p_lay, t_lay, vmr_h2o)
        
        # Calculate mixing ratios from water paths
        q_liq, q_ice = calculate_mixing_ratios(lwp, iwp, p_lay, t_lay, vmr_h2o)
        
        # Verify dimensions
        nx, ny, nz = length(x), length(y), length(z)
        
        # Return all data as a named tuple
        return (
            x = x,
            y = y,
            z = z,
            t_lay = t_lay,
            p_lay = p_lay,
            vmr_h2o = vmr_h2o,
            q_liq = q_liq,
            q_ice = q_ice,
            rho = rho
        )
        
    finally
        close(ds)
    end
end

"""
    calculate_moist_density(p_lay, t_lay, vmr_h2o)

Calculate moist air density accounting for water vapor.

# Arguments
- `p_lay`: Pressure on layers (Pa)
- `t_lay`: Temperature on layers (K)
- `vmr_h2o`: Water vapor volume mixing ratio (mol/mol)

# Returns
- `rho`: Moist air density (kg/m³)
"""
function calculate_moist_density(p_lay, t_lay, vmr_h2o)
    # Constants
    Rd = 287.05      # specific gas constant for dry air (J/kg/K)
    Rv = 461.5       # specific gas constant for water vapor (J/kg/K)
    
    # Molecular weights
    M_dry = 28.97    # molecular weight of dry air (g/mol)
    M_h2o = 18.015   # molecular weight of water (g/mol)
    
    # Convert volume mixing ratio to mass mixing ratio
    # q = (M_h2o/M_dry) * vmr_h2o
    q_h2o = (M_h2o / M_dry) * vmr_h2o
    
    # Virtual temperature accounts for moisture effect
    # Tv = T * (1 + 0.61*q) / (1 + q)
    # Or using the exact formula with gas constants:
    # ρ = p / (Rd * Tv) where Tv = T * (1 + (Rv/Rd - 1) * q) / (1 + q)
    
    # For moist air: ρ = p / (Rd * T * (1 + (Rv/Rd) * q) / (1 + q))
    # Simplified: ρ = p / (Rd * T) * (1 + q) / (1 + (Rv/Rd) * q)
    
    # Even simpler and commonly used approximation:
    # Virtual temperature: Tv = T / (1 - (1 - Rv/Rd) * e/p)
    # But for small q: Tv ≈ T * (1 + 0.61*q)
    
    epsilon = Rd / Rv  # ≈ 0.622
    
    # Virtual temperature (accounts for water vapor being lighter than dry air)
    T_virtual = t_lay .* (1 .+ q_h2o ./ epsilon) ./ (1 .+ q_h2o)
    
    # Moist air density using virtual temperature
    rho = p_lay ./ (Rd .* T_virtual)
    
    return rho
end

"""
    calculate_mixing_ratios(lwp, iwp, p_lay, t_lay, vmr_h2o)

Convert liquid and ice water paths to mixing ratios using moist air density.

# Arguments
- `lwp`: Liquid water path (kg/m²), size (nx, ny, nz)
- `iwp`: Ice water path (kg/m²), size (nx, ny, nz)
- `p_lay`: Pressure on layers (Pa), size (nx, ny, nz)
- `t_lay`: Temperature on layers (K), size (nx, ny, nz)
- `vmr_h2o`: Water vapor volume mixing ratio (mol/mol), size (nx, ny, nz)

# Returns
- `q_liq`: Liquid water mixing ratio (kg/kg)
- `q_ice`: Ice water mixing ratio (kg/kg)
"""
function calculate_mixing_ratios(lwp, iwp, p_lay, t_lay, vmr_h2o)
    # Constants
    g = 9.80665  # gravitational acceleration (m/s²)
    
    nx, ny, nz = size(p_lay)
    
    # Initialize mixing ratio arrays
    q_liq = zeros(Float64, nx, ny, nz)
    q_ice = zeros(Float64, nx, ny, nz)
    
    # Calculate moist air density
    rho_moist = calculate_moist_density(p_lay, t_lay, vmr_h2o)
    
    # Calculate layer thickness in pressure coordinates
    for i in 1:nx
        for j in 1:ny
            # Calculate pressure thickness for each layer (absolute value)
            dp = zeros(nz)
            for k in 1:nz
                if k == 1
                    dp[k] = abs(p_lay[i, j, k+1] - p_lay[i, j, k])
                elseif k == nz
                    dp[k] = abs(p_lay[i, j, k] - p_lay[i, j, k-1])
                else
                    dp[k] = abs(p_lay[i, j, k+1] - p_lay[i, j, k-1]) / 2
                end
            end
            
            # Layer mass per unit area: dm = dp/g
            dm = dp ./ g
            
            # Total column mass
            total_mass = sum(dm)
            
            # Distribute water path proportionally to layer mass
            if total_mass > 0
                for k in 1:nz
                    # Fraction of total column mass in this layer
                    mass_fraction = dm[k] / total_mass
                    
                    # Mixing ratio = (water mass in layer) / (air mass in layer)
                    # Water mass in layer = total water path × mass fraction
                    q_liq[i, j, k] = lwp[i, j, k] * mass_fraction / dm[k]
                    q_ice[i, j, k] = iwp[i, j, k] * mass_fraction / dm[k]
                end
            end
        end
    end
    
    return q_liq, q_ice
end

# Example usage:
# data = read_atmospheric_data("your_file.nc")
# println("Moist air density range: $(extrema(data.rho)) kg/m³")
# println("Liquid water mixing ratio range: $(extrema(data.q_liq))")
# println("Ice water mixing ratio range: $(extrema(data.q_ice))")

"""
    interpolate_to_mesh(data, mesh)

Interpolate atmospheric data from regular grid to unstructured mesh using trilinear interpolation.

# Arguments
- `data`: Named tuple from `read_atmospheric_data` containing grid coordinates and fields
- `mesh`: Named tuple with fields `x`, `y`, `z` (vectors of length npoin with arbitrary ordering)

# Returns
Named tuple with interpolated fields at mesh points:
- `t_lay`, `p_lay`, `vmr_h2o`, `q_liq`, `q_ice`, `rho`: vectors of length npoin
"""
function interpolate_atmosphere_to_mesh(data, mesh)
    npoin = length(mesh.x)
    
    # Initialize output arrays
    t_lay_interp = zeros(Float64, npoin)
    p_lay_interp = zeros(Float64, npoin)
    vmr_h2o_interp = zeros(Float64, npoin)
    q_liq_interp = zeros(Float64, npoin)
    q_ice_interp = zeros(Float64, npoin)
    rho_interp = zeros(Float64, npoin)
    
    # Interpolate each point independently (no assumptions about ordering)
    for i in 1:npoin
        t_lay_interp[i] = trilinear_interpolate(
            data.x, data.y, data.z, data.t_lay,
            mesh.x[i], mesh.y[i], mesh.z[i]
        )
        
        p_lay_interp[i] = trilinear_interpolate(
            data.x, data.y, data.z, data.p_lay,
            mesh.x[i], mesh.y[i], mesh.z[i]
        )
        
        vmr_h2o_interp[i] = trilinear_interpolate(
            data.x, data.y, data.z, data.vmr_h2o,
            mesh.x[i], mesh.y[i], mesh.z[i]
        )
        
        q_liq_interp[i] = trilinear_interpolate(
            data.x, data.y, data.z, data.q_liq,
            mesh.x[i], mesh.y[i], mesh.z[i]
        )
        
        q_ice_interp[i] = trilinear_interpolate(
            data.x, data.y, data.z, data.q_ice,
            mesh.x[i], mesh.y[i], mesh.z[i]
        )
        
        rho_interp[i] = trilinear_interpolate(
            data.x, data.y, data.z, data.rho,
            mesh.x[i], mesh.y[i], mesh.z[i]
        )
    end
    
    return (
        t_lay = t_lay_interp,
        p_lay = p_lay_interp,
        vmr_h2o = vmr_h2o_interp,
        q_liq = q_liq_interp,
        q_ice = q_ice_interp,
        rho = rho_interp
    )
end

"""
    trilinear_interpolate(x_grid, y_grid, z_grid, field, x_point, y_point, z_point)

Perform trilinear interpolation of a 3D field at a single point.
Assumes x_grid, y_grid, z_grid define a regular grid (sorted coordinates).

# Arguments
- `x_grid`, `y_grid`, `z_grid`: 1D coordinate arrays defining the regular grid (must be sorted)
- `field`: 3D array of values on the grid, size (nx, ny, nz)
- `x_point`, `y_point`, `z_point`: coordinates of the point to interpolate to

# Returns
- Interpolated value at (x_point, y_point, z_point)
"""
function trilinear_interpolate(x_grid, y_grid, z_grid, field, x_point, y_point, z_point)
    # Find the indices of the surrounding grid points
    i = find_interval(x_grid, x_point)
    j = find_interval(y_grid, y_point)
    k = find_interval(z_grid, z_point)
    
    # Handle boundary cases (extrapolation - clamp to grid boundaries)
    nx, ny, nz = length(x_grid), length(y_grid), length(z_grid)
    i = clamp(i, 1, nx-1)
    j = clamp(j, 1, ny-1)
    k = clamp(k, 1, nz-1)
    
    # Get the coordinates of the surrounding cell
    x0, x1 = x_grid[i], x_grid[i+1]
    y0, y1 = y_grid[j], y_grid[j+1]
    z0, z1 = z_grid[k], z_grid[k+1]
    
    # Normalized coordinates within the cell (0 to 1)
    xd = (x_point - x0) / (x1 - x0)
    yd = (y_point - y0) / (y1 - y0)
    zd = (z_point - z0) / (z1 - z0)
    
    # Clamp to [0, 1] in case of small numerical errors
    xd = clamp(xd, 0.0, 1.0)
    yd = clamp(yd, 0.0, 1.0)
    zd = clamp(zd, 0.0, 1.0)
    
    # Get the 8 corner values
    c000 = field[i,   j,   k]
    c100 = field[i+1, j,   k]
    c010 = field[i,   j+1, k]
    c110 = field[i+1, j+1, k]
    c001 = field[i,   j,   k+1]
    c101 = field[i+1, j,   k+1]
    c011 = field[i,   j+1, k+1]
    c111 = field[i+1, j+1, k+1]
    
    # Interpolate along x
    c00 = c000 * (1 - xd) + c100 * xd
    c01 = c001 * (1 - xd) + c101 * xd
    c10 = c010 * (1 - xd) + c110 * xd
    c11 = c011 * (1 - xd) + c111 * xd
    
    # Interpolate along y
    c0 = c00 * (1 - yd) + c10 * yd
    c1 = c01 * (1 - yd) + c11 * yd
    
    # Interpolate along z
    c = c0 * (1 - zd) + c1 * zd
    
    return c
end

"""
    find_interval(grid, point)

Find the index i such that grid[i] <= point < grid[i+1].
Assumes grid is sorted in ascending order.

# Arguments
- `grid`: 1D sorted array
- `point`: value to locate

# Returns
- Index i such that grid[i] <= point < grid[i+1]
"""
function find_interval(grid, point)
    n = length(grid)
    
    # Handle out-of-bounds cases
    if point <= grid[1]
        return 1
    elseif point >= grid[end]
        return n - 1
    end
    
    # Binary search for efficiency
    i_low = 1
    i_high = n
    
    while i_high - i_low > 1
        i_mid = div(i_low + i_high, 2)
        if point < grid[i_mid]
            i_high = i_mid
        else
            i_low = i_mid
        end
    end
    
    return i_low
end

# Example usage:
# # Read NetCDF data (assumes regular grid with sorted x, y, z coordinates)
# data = read_atmospheric_data("atmospheric_data.nc")
# 
# # Define your unstructured mesh with arbitrary point ordering
# mesh = (
#     x = [2.3, 5.7, 1.2, 8.4, ...],  # npoin values in any order
#     y = [1.5, 3.2, 0.8, 4.1, ...],  # npoin values in any order
#     z = [150.0, 200.0, 100.0, 250.0, ...]  # npoin values in any order
# )
# 
# # Interpolate to mesh (works regardless of mesh point ordering)
# mesh_data = interpolate_to_mesh(data, mesh)
# 
# println("Interpolated temperature at point 1: $(mesh_data.t_lay[1])")
# println("Interpolated temperature at point 3: $(mesh_data.t_lay[3])")

