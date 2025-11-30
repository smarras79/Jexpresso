#!/usr/bin/env julia
"""
Convert ERA5 NetCDF to Jexpresso Sounding Format

This utility script reads ERA5 NetCDF data and exports it to Jexpresso's
ASCII sounding format for easy visualization and testing.

Usage:
    julia era5_to_sounding.jl <era5_file.nc> <output_sounding.dat> [lat] [lon]

Example:
    julia era5_to_sounding.jl era5_data.nc era5_sounding.dat 40.0 -105.0

Author: Jexpresso Development Team
Date: 2025-11-30
"""

using Pkg
Pkg.activate(".")

# Load Jexpresso
include("../src/Jexpresso.jl")
using .Jexpresso

function main()
    # Parse command line arguments
    if length(ARGS) < 2
        println("Usage: julia era5_to_sounding.jl <era5_file.nc> <output.dat> [lat] [lon]")
        println("")
        println("Arguments:")
        println("  era5_file.nc  - Input ERA5 NetCDF file")
        println("  output.dat    - Output sounding file")
        println("  lat           - Target latitude [degrees] (default: center of domain)")
        println("  lon           - Target longitude [degrees] (default: center of domain)")
        println("")
        println("Example:")
        println("  julia era5_to_sounding.jl data/era5.nc data/sounding.dat 40.0 -105.0")
        exit(1)
    end

    era5_file = ARGS[1]
    output_file = ARGS[2]

    # Check if file exists
    if !isfile(era5_file)
        error("ERA5 file not found: $era5_file")
    end

    println("="^70)
    println("ERA5 to Sounding Converter")
    println("="^70)
    println("Input:  $era5_file")
    println("Output: $output_file")

    # Read ERA5 data
    println("\nReading ERA5 data...")
    era5_data = read_era5_netcdf(era5_file)

    # Determine target location
    if length(ARGS) >= 4
        target_lat = parse(Float64, ARGS[3])
        target_lon = parse(Float64, ARGS[4])
    else
        # Use center of domain
        target_lat = mean(era5_data.lat)
        target_lon = mean(era5_data.lon)
        println("No lat/lon specified, using domain center:")
    end

    println("Target location: lat = $target_lat°, lon = $target_lon°")

    # Define vertical grid (0 to max height, every 50 m)
    max_height = maximum(era5_data.heights)
    min_height = minimum(era5_data.heights)
    target_heights = collect(min_height:50.0:max_height)

    println("Vertical range: $(min_height) to $(max_height) m")
    println("Vertical levels: $(length(target_heights))")

    # Interpolate to 1D column
    println("\nInterpolating to 1D column...")
    era5_1d = interpolate_era5_to_1d_column(era5_data, target_lat, target_lon, target_heights)

    # Export to sounding format
    println("\nExporting to sounding format...")
    export_era5_to_sounding(era5_1d, output_file)

    println("\n" * "="^70)
    println("✓ Conversion complete!")
    println("="^70)
    println("\nSounding file created: $output_file")
    println("Vertical levels: $(length(target_heights))")
    println("\nFormat: height[m] T[K] u[m/s] v[m/s] q[kg/kg] P[Pa]")
    println("\nYou can now use this file with Jexpresso's sounding reader.")
    println()

    # Print summary statistics
    println("Profile Summary:")
    println("  Temperature: $(minimum(era5_1d["temperature"])) - $(maximum(era5_1d["temperature"])) K")
    println("  Wind speed:  $(minimum(sqrt.(era5_1d["u_wind"].^2 .+ era5_1d["v_wind"].^2))) - " *
            "$(maximum(sqrt.(era5_1d["u_wind"].^2 .+ era5_1d["v_wind"].^2))) m/s")
    println("  Humidity:    $(minimum(era5_1d["specific_humidity"])*1000) - " *
            "$(maximum(era5_1d["specific_humidity"])*1000) g/kg")
    println()
end

# Run main function
if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
