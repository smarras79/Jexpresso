"""
    ERA5 Data Reader Module for Jexpresso

    This module provides functionality to read and interpolate ERA5 reanalysis data
    from ECMWF for use in Large Eddy Simulation (LES) weather forecasting.

    ERA5 is a comprehensive global atmospheric reanalysis dataset that provides
    hourly estimates of atmospheric, land and oceanic climate variables.

    Author: Jexpresso Development Team
    Date: 2025-11-30
"""

using NCDatasets
using Interpolations
using Dates

"""
    St_ERA5Data

Structure to hold ERA5 atmospheric data fields.

Fields:
- `temperature`: Temperature [K]
- `pressure`: Pressure levels [Pa]
- `u_wind`: Zonal wind component [m/s]
- `v_wind`: Meridional wind component [m/s]
- `w_wind`: Vertical wind component [m/s] (optional)
- `specific_humidity`: Specific humidity [kg/kg]
- `geopotential`: Geopotential [m²/s²]
- `surface_pressure`: Surface pressure [Pa]
- `heights`: Geometric heights [m]
- `time`: Time coordinate
- `lat`: Latitude coordinate [degrees]
- `lon`: Longitude coordinate [degrees]
"""
Base.@kwdef mutable struct St_ERA5Data{T <: AbstractFloat}
    temperature::Array{T}          # Temperature [K]
    pressure::Array{T}              # Pressure [Pa]
    u_wind::Array{T}                # Zonal wind [m/s]
    v_wind::Array{T}                # Meridional wind [m/s]
    w_wind::Array{T}                # Vertical wind [m/s]
    specific_humidity::Array{T}     # Specific humidity [kg/kg]
    geopotential::Array{T}          # Geopotential [m²/s²]
    surface_pressure::Array{T}      # Surface pressure [Pa]
    heights::Array{T}               # Geometric height [m]
    time::Array{Float64}            # Time coordinate
    lat::Array{T}                   # Latitude [degrees]
    lon::Array{T}                   # Longitude [degrees]
end

"""
    read_era5_netcdf(filename::String;
                     var_names::Dict = Dict(),
                     time_index::Int = 1,
                     lat_range::Tuple = (nothing, nothing),
                     lon_range::Tuple = (nothing, nothing),
                     pressure_levels::Union{Vector, Nothing} = nothing)

Read ERA5 data from a NetCDF file.

# Arguments
- `filename`: Path to ERA5 NetCDF file
- `var_names`: Dictionary mapping standard names to actual variable names in the file
  - Default mappings: "T", "u", "v", "q", "z", "sp", "level"
- `time_index`: Time index to read (default: 1)
- `lat_range`: Tuple of (min_lat, max_lat) to subset data
- `lon_range`: Tuple of (min_lon, max_lon) to subset data
- `pressure_levels`: Specific pressure levels to extract [Pa]

# Returns
- `St_ERA5Data`: Structure containing ERA5 data fields

# Example
```julia
era5_data = read_era5_netcdf("era5_data.nc",
                             time_index=1,
                             lat_range=(30.0, 50.0),
                             lon_range=(-10.0, 10.0))
```
"""
function read_era5_netcdf(filename::String;
                          var_names::Dict = Dict(),
                          time_index::Int = 1,
                          lat_range::Tuple = (nothing, nothing),
                          lon_range::Tuple = (nothing, nothing),
                          pressure_levels::Union{Vector, Nothing} = nothing)

    # Default variable name mappings (can be customized)
    default_var_names = Dict(
        "T" => "t",           # Temperature
        "u" => "u",           # U wind component
        "v" => "v",           # V wind component
        "w" => "w",           # W wind component (optional)
        "q" => "q",           # Specific humidity
        "z" => "z",           # Geopotential
        "sp" => "sp",         # Surface pressure
        "level" => "level",   # Pressure levels
        "lat" => "latitude",
        "lon" => "longitude",
        "time" => "time"
    )

    # Merge user-provided variable names with defaults
    vnames = merge(default_var_names, var_names)

    println("Reading ERA5 data from: $filename")

    # Open NetCDF dataset
    ds = Dataset(filename, "r")

    try
        # Read coordinate variables
        lat_all = ds[vnames["lat"]][:]
        lon_all = ds[vnames["lon"]][:]
        time_all = ds[vnames["time"]][:]

        # Determine lat/lon indices for subsetting
        if lat_range[1] !== nothing && lat_range[2] !== nothing
            lat_idx = findall(x -> lat_range[1] <= x <= lat_range[2], lat_all)
        else
            lat_idx = 1:length(lat_all)
        end

        if lon_range[1] !== nothing && lon_range[2] !== nothing
            lon_idx = findall(x -> lon_range[1] <= x <= lon_range[2], lon_all)
        else
            lon_idx = 1:length(lon_all)
        end

        lat = lat_all[lat_idx]
        lon = lon_all[lon_idx]

        # Read pressure levels if available
        if haskey(ds, vnames["level"])
            levels_hPa = ds[vnames["level"]][:]  # Usually in hPa
            levels = levels_hPa .* 100.0         # Convert to Pa

            if pressure_levels !== nothing
                # Find closest pressure levels
                level_idx = [argmin(abs.(levels .- p)) for p in pressure_levels]
            else
                level_idx = 1:length(levels)
            end
        else
            levels = [100000.0]  # Single level surface data
            level_idx = 1:1
        end

        pressure = levels[level_idx]
        nlevels = length(level_idx)
        nlat = length(lat_idx)
        nlon = length(lon_idx)

        # Initialize arrays
        T = zeros(Float64, nlon, nlat, nlevels)
        u = zeros(Float64, nlon, nlat, nlevels)
        v = zeros(Float64, nlon, nlat, nlevels)
        w = zeros(Float64, nlon, nlat, nlevels)
        q = zeros(Float64, nlon, nlat, nlevels)
        z = zeros(Float64, nlon, nlat, nlevels)

        # Read 3D fields
        if haskey(ds, vnames["T"])
            if ndims(ds[vnames["T"]]) == 4  # (lon, lat, level, time)
                T[:,:,:] = ds[vnames["T"]][lon_idx, lat_idx, level_idx, time_index]
            elseif ndims(ds[vnames["T"]]) == 3  # (lon, lat, time) - surface only
                T[:,:,1] = ds[vnames["T"]][lon_idx, lat_idx, time_index]
            end
        else
            @warn "Temperature variable '$(vnames["T"])' not found in file"
        end

        if haskey(ds, vnames["u"])
            if ndims(ds[vnames["u"]]) == 4
                u[:,:,:] = ds[vnames["u"]][lon_idx, lat_idx, level_idx, time_index]
            elseif ndims(ds[vnames["u"]]) == 3
                u[:,:,1] = ds[vnames["u"]][lon_idx, lat_idx, time_index]
            end
        else
            @warn "U-wind variable '$(vnames["u"])' not found in file"
        end

        if haskey(ds, vnames["v"])
            if ndims(ds[vnames["v"]]) == 4
                v[:,:,:] = ds[vnames["v"]][lon_idx, lat_idx, level_idx, time_index]
            elseif ndims(ds[vnames["v"]]) == 3
                v[:,:,1] = ds[vnames["v"]][lon_idx, lat_idx, time_index]
            end
        else
            @warn "V-wind variable '$(vnames["v"])' not found in file"
        end

        if haskey(ds, vnames["w"])
            if ndims(ds[vnames["w"]]) == 4
                w[:,:,:] = ds[vnames["w"]][lon_idx, lat_idx, level_idx, time_index]
            elseif ndims(ds[vnames["w"]]) == 3
                w[:,:,1] = ds[vnames["w"]][lon_idx, lat_idx, time_index]
            end
        end

        if haskey(ds, vnames["q"])
            if ndims(ds[vnames["q"]]) == 4
                q[:,:,:] = ds[vnames["q"]][lon_idx, lat_idx, level_idx, time_index]
            elseif ndims(ds[vnames["q"]]) == 3
                q[:,:,1] = ds[vnames["q"]][lon_idx, lat_idx, time_index]
            end
        else
            @warn "Specific humidity variable '$(vnames["q"])' not found in file"
        end

        if haskey(ds, vnames["z"])
            if ndims(ds[vnames["z"]]) == 4
                z[:,:,:] = ds[vnames["z"]][lon_idx, lat_idx, level_idx, time_index]
            elseif ndims(ds[vnames["z"]]) == 3
                z[:,:,1] = ds[vnames["z"]][lon_idx, lat_idx, time_index]
            end
        else
            @warn "Geopotential variable '$(vnames["z"])' not found in file"
        end

        # Read surface pressure
        sp = zeros(Float64, nlon, nlat)
        if haskey(ds, vnames["sp"])
            if ndims(ds[vnames["sp"]]) == 3  # (lon, lat, time)
                sp[:,:] = ds[vnames["sp"]][lon_idx, lat_idx, time_index]
            elseif ndims(ds[vnames["sp"]]) == 2  # (lon, lat)
                sp[:,:] = ds[vnames["sp"]][lon_idx, lat_idx]
            end
        else
            @warn "Surface pressure variable '$(vnames["sp"])' not found, using standard pressure"
            sp .= 101325.0  # Standard sea level pressure
        end

        # Calculate geometric height from geopotential
        g0 = 9.80665  # Standard gravity [m/s²]
        heights = z ./ g0

        println("Successfully read ERA5 data:")
        println("  Dimensions: $nlon × $nlat × $nlevels")
        println("  Lat range: $(minimum(lat)) to $(maximum(lat)) degrees")
        println("  Lon range: $(minimum(lon)) to $(maximum(lon)) degrees")
        println("  Pressure levels: $(length(pressure))")
        println("  Height range: $(minimum(heights)) to $(maximum(heights)) m")

        # Create and return ERA5 data structure
        era5_data = St_ERA5Data{Float64}(
            temperature = T,
            pressure = pressure,
            u_wind = u,
            v_wind = v,
            w_wind = w,
            specific_humidity = q,
            geopotential = z,
            surface_pressure = sp,
            heights = heights,
            time = [time_all[time_index]],
            lat = lat,
            lon = lon
        )

        return era5_data

    finally
        close(ds)
    end
end

"""
    interpolate_era5_to_1d_column(era5_data::St_ERA5Data,
                                   target_lat::Real,
                                   target_lon::Real,
                                   target_heights::Vector{<:Real})

Interpolate ERA5 3D data to a 1D vertical column at specified lat/lon.

# Arguments
- `era5_data`: ERA5 data structure from read_era5_netcdf
- `target_lat`: Target latitude [degrees]
- `target_lon`: Target longitude [degrees]
- `target_heights`: Target height levels for interpolation [m]

# Returns
- Dictionary with interpolated 1D profiles for all variables
"""
function interpolate_era5_to_1d_column(era5_data::St_ERA5Data,
                                       target_lat::Real,
                                       target_lon::Real,
                                       target_heights::Vector{<:Real})

    println("Interpolating ERA5 data to 1D column at lat=$target_lat°, lon=$target_lon°")

    # Create 2D interpolators for horizontal interpolation
    lat = era5_data.lat
    lon = era5_data.lon

    nlevels = length(era5_data.pressure)
    nz_target = length(target_heights)

    # Initialize output arrays
    T_col = zeros(nz_target)
    u_col = zeros(nz_target)
    v_col = zeros(nz_target)
    w_col = zeros(nz_target)
    q_col = zeros(nz_target)
    p_col = zeros(nz_target)

    # For each pressure level, interpolate horizontally to target lat/lon
    heights_at_point = zeros(nlevels)
    T_at_point = zeros(nlevels)
    u_at_point = zeros(nlevels)
    v_at_point = zeros(nlevels)
    w_at_point = zeros(nlevels)
    q_at_point = zeros(nlevels)

    for k in 1:nlevels
        # Create 2D interpolators for this level
        itp_h = LinearInterpolation((lon, lat), era5_data.heights[:,:,k], extrapolation_bc=Flat())
        itp_T = LinearInterpolation((lon, lat), era5_data.temperature[:,:,k], extrapolation_bc=Flat())
        itp_u = LinearInterpolation((lon, lat), era5_data.u_wind[:,:,k], extrapolation_bc=Flat())
        itp_v = LinearInterpolation((lon, lat), era5_data.v_wind[:,:,k], extrapolation_bc=Flat())
        itp_w = LinearInterpolation((lon, lat), era5_data.w_wind[:,:,k], extrapolation_bc=Flat())
        itp_q = LinearInterpolation((lon, lat), era5_data.specific_humidity[:,:,k], extrapolation_bc=Flat())

        heights_at_point[k] = itp_h(target_lon, target_lat)
        T_at_point[k] = itp_T(target_lon, target_lat)
        u_at_point[k] = itp_u(target_lon, target_lat)
        v_at_point[k] = itp_v(target_lon, target_lat)
        w_at_point[k] = itp_w(target_lon, target_lat)
        q_at_point[k] = itp_q(target_lon, target_lat)
    end

    # Sort by height (ERA5 might be top-down or bottom-up)
    sort_idx = sortperm(heights_at_point)
    heights_sorted = heights_at_point[sort_idx]

    # Vertical interpolation to target heights
    if nlevels > 1
        itp_T_z = LinearInterpolation(heights_sorted, T_at_point[sort_idx], extrapolation_bc=Flat())
        itp_u_z = LinearInterpolation(heights_sorted, u_at_point[sort_idx], extrapolation_bc=Flat())
        itp_v_z = LinearInterpolation(heights_sorted, v_at_point[sort_idx], extrapolation_bc=Flat())
        itp_w_z = LinearInterpolation(heights_sorted, w_at_point[sort_idx], extrapolation_bc=Flat())
        itp_q_z = LinearInterpolation(heights_sorted, q_at_point[sort_idx], extrapolation_bc=Flat())

        for (i, z) in enumerate(target_heights)
            T_col[i] = itp_T_z(z)
            u_col[i] = itp_u_z(z)
            v_col[i] = itp_v_z(z)
            w_col[i] = itp_w_z(z)
            q_col[i] = itp_q_z(z)
        end
    else
        # Single level - use constant value
        T_col .= T_at_point[1]
        u_col .= u_at_point[1]
        v_col .= v_at_point[1]
        w_col .= w_at_point[1]
        q_col .= q_at_point[1]
    end

    # Estimate pressure from hydrostatic balance
    # p(z) = p_s * exp(-z/H) where H is scale height
    sp_itp = LinearInterpolation((lon, lat), era5_data.surface_pressure, extrapolation_bc=Flat())
    p_surface = sp_itp(target_lon, target_lat)

    R_dry = 287.05  # Gas constant for dry air [J/(kg·K)]
    g = 9.80665     # Gravity [m/s²]

    for (i, z) in enumerate(target_heights)
        # Simple barometric formula
        T_mean = mean(T_col[1:i])
        H = R_dry * T_mean / g  # Scale height
        p_col[i] = p_surface * exp(-z / H)
    end

    println("  Interpolated to $(nz_target) vertical levels")
    println("  Height range: $(minimum(target_heights)) to $(maximum(target_heights)) m")

    return Dict(
        "height" => target_heights,
        "temperature" => T_col,
        "u_wind" => u_col,
        "v_wind" => v_col,
        "w_wind" => w_col,
        "specific_humidity" => q_col,
        "pressure" => p_col
    )
end

"""
    interpolate_era5_to_mesh(backend, era5_data::St_ERA5Data,
                             mesh,
                             target_lat::Real,
                             target_lon::Real)

Interpolate ERA5 data to Jexpresso mesh coordinates.

# Arguments
- `backend`: Computation backend (CPU() or CUDABackend())
- `era5_data`: ERA5 data structure
- `mesh`: Jexpresso mesh structure with z coordinates
- `target_lat`: Target latitude [degrees]
- `target_lon`: Target longitude [degrees]

# Returns
- Dictionary with interpolated fields on mesh points
"""
function interpolate_era5_to_mesh(backend, era5_data::St_ERA5Data,
                                  mesh,
                                  target_lat::Real,
                                  target_lon::Real)

    println("Interpolating ERA5 data to Jexpresso mesh...")

    # Get mesh heights
    z_mesh = mesh.z
    npoin = mesh.npoin

    # Convert to regular array for interpolation
    z_array = if backend == CPU()
        collect(z_mesh)
    else
        Array(z_mesh)  # Copy from GPU to CPU for interpolation
    end

    # Get unique heights for initial interpolation
    z_unique = sort(unique(z_array))

    # Interpolate ERA5 to 1D column
    era5_1d = interpolate_era5_to_1d_column(era5_data, target_lat, target_lon, z_unique)

    # Create interpolators for mesh points
    itp_T = LinearInterpolation(z_unique, era5_1d["temperature"], extrapolation_bc=Flat())
    itp_u = LinearInterpolation(z_unique, era5_1d["u_wind"], extrapolation_bc=Flat())
    itp_v = LinearInterpolation(z_unique, era5_1d["v_wind"], extrapolation_bc=Flat())
    itp_w = LinearInterpolation(z_unique, era5_1d["w_wind"], extrapolation_bc=Flat())
    itp_q = LinearInterpolation(z_unique, era5_1d["specific_humidity"], extrapolation_bc=Flat())
    itp_p = LinearInterpolation(z_unique, era5_1d["pressure"], extrapolation_bc=Flat())

    # Interpolate to all mesh points
    T_mesh = zeros(npoin)
    u_mesh = zeros(npoin)
    v_mesh = zeros(npoin)
    w_mesh = zeros(npoin)
    q_mesh = zeros(npoin)
    p_mesh = zeros(npoin)

    for i in 1:npoin
        T_mesh[i] = itp_T(z_array[i])
        u_mesh[i] = itp_u(z_array[i])
        v_mesh[i] = itp_v(z_array[i])
        w_mesh[i] = itp_w(z_array[i])
        q_mesh[i] = itp_q(z_array[i])
        p_mesh[i] = itp_p(z_array[i])
    end

    # Transfer to backend if needed
    if backend != CPU()
        T_mesh_backend = KernelAbstractions.allocate(backend, Float64, npoin)
        u_mesh_backend = KernelAbstractions.allocate(backend, Float64, npoin)
        v_mesh_backend = KernelAbstractions.allocate(backend, Float64, npoin)
        w_mesh_backend = KernelAbstractions.allocate(backend, Float64, npoin)
        q_mesh_backend = KernelAbstractions.allocate(backend, Float64, npoin)
        p_mesh_backend = KernelAbstractions.allocate(backend, Float64, npoin)

        KernelAbstractions.copyto!(backend, T_mesh_backend, T_mesh)
        KernelAbstractions.copyto!(backend, u_mesh_backend, u_mesh)
        KernelAbstractions.copyto!(backend, v_mesh_backend, v_mesh)
        KernelAbstractions.copyto!(backend, w_mesh_backend, w_mesh)
        KernelAbstractions.copyto!(backend, q_mesh_backend, q_mesh)
        KernelAbstractions.copyto!(backend, p_mesh_backend, p_mesh)

        return Dict(
            "temperature" => T_mesh_backend,
            "u_wind" => u_mesh_backend,
            "v_wind" => v_mesh_backend,
            "w_wind" => w_mesh_backend,
            "specific_humidity" => q_mesh_backend,
            "pressure" => p_mesh_backend
        )
    else
        return Dict(
            "temperature" => T_mesh,
            "u_wind" => u_mesh,
            "v_wind" => v_mesh,
            "w_wind" => w_mesh,
            "specific_humidity" => q_mesh,
            "pressure" => p_mesh
        )
    end
end

"""
    export_era5_to_sounding(era5_1d::Dict, output_file::String)

Export ERA5 1D profile to Jexpresso sounding format.

# Arguments
- `era5_1d`: Dictionary with 1D profiles from interpolate_era5_to_1d_column
- `output_file`: Output file path

# Format
Space-delimited ASCII file with columns:
height[m] temperature[K] u_wind[m/s] v_wind[m/s] specific_humidity[kg/kg] pressure[Pa]
"""
function export_era5_to_sounding(era5_1d::Dict, output_file::String)

    println("Exporting ERA5 data to sounding file: $output_file")

    open(output_file, "w") do io
        # Write header
        println(io, "# ERA5 vertical profile exported for Jexpresso")
        println(io, "# Columns: height[m] temperature[K] u_wind[m/s] v_wind[m/s] q[kg/kg] pressure[Pa]")

        # Write data
        nz = length(era5_1d["height"])
        for i in 1:nz
            println(io, "$(era5_1d["height"][i]) $(era5_1d["temperature"][i]) " *
                       "$(era5_1d["u_wind"][i]) $(era5_1d["v_wind"][i]) " *
                       "$(era5_1d["specific_humidity"][i]) $(era5_1d["pressure"][i])")
        end
    end

    println("  Exported $nz vertical levels")
end
