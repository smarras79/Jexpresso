# ERA5 Data Extension for Jexpresso

**Author:** Jexpresso Development Team
**Date:** November 30, 2025
**Version:** 1.0

## Table of Contents

1. [Overview](#overview)
2. [Features](#features)
3. [Installation Requirements](#installation-requirements)
4. [Quick Start](#quick-start)
5. [Obtaining ERA5 Data](#obtaining-era5-data)
6. [Configuration Options](#configuration-options)
7. [Example Usage](#example-usage)
8. [API Reference](#api-reference)
9. [Troubleshooting](#troubleshooting)
10. [Best Practices](#best-practices)

---

## Overview

The ERA5 extension enables Jexpresso to use real-world atmospheric reanalysis data from ECMWF's ERA5 dataset for Large Eddy Simulation (LES) weather forecasting. This extension provides:

- **Automatic reading** of ERA5 NetCDF files
- **Spatial interpolation** from ERA5 grid to Jexpresso mesh
- **Vertical interpolation** from pressure levels to geometric heights
- **Initial condition generation** from ERA5 data
- **Optional nudging/relaxation** forcing for LES simulations

### What is ERA5?

ERA5 is the fifth generation ECMWF atmospheric reanalysis of the global climate covering the period from January 1940 to present. It provides hourly estimates of atmospheric, land and oceanic climate variables at:

- **Temporal resolution:** Hourly
- **Spatial resolution:** 0.25° × 0.25° (~31 km)
- **Vertical levels:** 137 model levels or 37 pressure levels
- **Variables:** Temperature, wind (u,v,w), humidity, pressure, geopotential, and 100+ more

---

## Features

### 1. **Flexible Data Reading**
- Supports standard ERA5 NetCDF format
- Customizable variable name mappings
- Spatial subsetting (lat/lon ranges)
- Pressure level selection
- Time index selection for multi-time files

### 2. **Advanced Interpolation**
- Horizontal interpolation to specific lat/lon points
- Vertical interpolation from pressure to height coordinates
- Automatic mesh interpolation for any Jexpresso grid
- GPU-compatible data structures

### 3. **LES Integration**
- Automatic initialization from ERA5 profiles
- Optional turbulence perturbations
- Nudging/relaxation forcing toward ERA5 state
- Compatible with existing Jexpresso physics (microphysics, radiation, SGS models)

### 4. **Easy Configuration**
- Simple dictionary-based configuration
- Sensible default values
- Extensive error checking and warnings

---

## Installation Requirements

### Julia Packages

The ERA5 extension requires the following Julia packages (already included in Jexpresso):

```julia
using NCDatasets        # NetCDF file reading
using Interpolations    # Spatial/vertical interpolation
using Dates             # Time handling
using KernelAbstractions # GPU support
```

These packages are automatically loaded when you `using Jexpresso`.

### External Tools (Optional)

For downloading and preprocessing ERA5 data:

- **Python 3.7+** with `cdsapi` package (for ERA5 download)
- **CDO** (Climate Data Operators) - for data preprocessing
- **ncview** or **Panoply** - for visualizing NetCDF files

---

## Quick Start

### Step 1: Download ERA5 Data

See [Obtaining ERA5 Data](#obtaining-era5-data) section for detailed instructions.

### Step 2: Create Test Case

Copy the example case:

```bash
cp -r problems/equations/CompEuler/era5_les problems/equations/CompEuler/my_era5_case
```

### Step 3: Configure Input Parameters

Edit `problems/equations/CompEuler/my_era5_case/user_inputs.jl`:

```julia
:lERA5                => true,
:era5_file            => "./data_files/my_era5_data.nc",
:era5_target_lat      => 40.0,      # Your location
:era5_target_lon      => -105.0,
:era5_time_index      => 1,
```

### Step 4: Run Simulation

```bash
julia --project=. src/run.jl CompEuler my_era5_case
```

---

## Obtaining ERA5 Data

### Method 1: Copernicus Climate Data Store (Recommended)

1. **Create Account:** Register at [https://cds.climate.copernicus.eu](https://cds.climate.copernicus.eu)

2. **Install CDS API:**
   ```bash
   pip install cdsapi
   ```

3. **Configure API Key:** Create `~/.cdsapirc`:
   ```
   url: https://cds.climate.copernicus.eu/api/v2
   key: YOUR_UID:YOUR_API_KEY
   ```

4. **Download Data:** Use the provided Python script `scripts/download_era5.py`:

   ```python
   import cdsapi

   c = cdsapi.Client()

   c.retrieve(
       'reanalysis-era5-pressure-levels',
       {
           'product_type': 'reanalysis',
           'format': 'netcdf',
           'variable': [
               'temperature', 'u_component_of_wind', 'v_component_of_wind',
               'vertical_velocity', 'specific_humidity', 'geopotential'
           ],
           'pressure_level': [
               '1000', '950', '900', '850', '800', '750', '700',
               '650', '600', '550', '500', '450', '400', '350', '300'
           ],
           'year': '2024',
           'month': '01',
           'day': '15',
           'time': '12:00',
           'area': [45, -110, 35, -100],  # North, West, South, East
       },
       'my_era5_data.nc')
   ```

5. **Run Download:**
   ```bash
   python scripts/download_era5.py
   ```

### Method 2: Manual Download via Web Interface

1. Go to [CDS Dataset Page](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-pressure-levels)
2. Select variables, date, time, area, and format
3. Submit request and download when ready

### Required Variables

For Jexpresso LES, download these variables at minimum:

| Variable | ERA5 Name | Units |
|----------|-----------|-------|
| Temperature | `temperature` | K |
| U wind | `u_component_of_wind` | m/s |
| V wind | `v_component_of_wind` | m/s |
| Specific humidity | `specific_humidity` | kg/kg |
| Geopotential | `geopotential` | m²/s² |
| Surface pressure | `surface_pressure` | Pa |

Optional but recommended:
- `vertical_velocity` (w) - Pa/s
- `cloud_liquid_water_content` - kg/kg

### Recommended Pressure Levels

```
1000, 975, 950, 925, 900, 875, 850, 825, 800, 775, 750, 700, 650, 600,
550, 500, 450, 400, 350, 300, 250, 225, 200, 175, 150, 125, 100 hPa
```

For LES (typically 0-5 km domain), levels 1000-500 hPa are most important.

---

## Configuration Options

### Basic ERA5 Parameters

Add these to your `user_inputs.jl`:

```julia
# Enable ERA5 extension
:lERA5                => true,

# Path to ERA5 NetCDF file
:era5_file            => "./data_files/era5_data.nc",

# Target location (lat/lon in degrees)
:era5_target_lat      => 40.0,     # Latitude [degrees N]
:era5_target_lon      => -105.0,   # Longitude [degrees E]

# Time index in NetCDF file (for multi-time files)
:era5_time_index      => 1,
```

### Advanced Options

#### Custom Variable Names

If your ERA5 file uses different variable names:

```julia
:era5_var_names       => Dict(
    "T" => "temperature",      # Temperature
    "u" => "u_wind",           # U component
    "v" => "v_wind",           # V component
    "w" => "w_wind",           # W component (optional)
    "q" => "specific_humidity", # Specific humidity
    "z" => "geopotential",     # Geopotential
    "sp" => "surface_pressure", # Surface pressure
    "level" => "pressure_level", # Pressure levels
    "lat" => "latitude",
    "lon" => "longitude",
    "time" => "time"
),
```

#### Spatial Subsetting

Reduce memory by loading only a region:

```julia
:era5_lat_range       => (35.0, 45.0),   # Min/max latitude
:era5_lon_range       => (-110.0, -100.0), # Min/max longitude
```

#### Pressure Level Selection

Select specific pressure levels:

```julia
:era5_pressure_levels => [100000.0, 95000.0, 90000.0, 85000.0, 80000.0],
```

Values in Pascals. If not specified, all levels in the file are used.

#### ERA5 Nudging/Forcing

Apply relaxation forcing toward ERA5 state (experimental):

```julia
:lERA5_forcing        => true,
:era5_nudging_timescale => 3600.0,  # Relaxation timescale [seconds]
```

Nudging applies: `dq/dt = -(q - q_ERA5) / τ`

Typical values:
- `τ = 3600` (1 hour) - strong nudging
- `τ = 10800` (3 hours) - moderate nudging
- `τ = 21600` (6 hours) - weak nudging

---

## Example Usage

### Example 1: Simple LES with ERA5 Initial Conditions

```julia
# user_inputs.jl
function user_inputs()
    inputs = Dict(
        # Time integration
        :Δt                   => 1.0,
        :tend                 => 3600.0,  # 1 hour

        # ERA5 configuration
        :lERA5                => true,
        :era5_file            => "./data_files/boulder_20240115.nc",
        :era5_target_lat      => 40.0,
        :era5_target_lon      => -105.0,

        # LES domain: 5km × 5km × 3km
        :nelx                 => 20,
        :xmax                 => 5000.0,
        :nely                 => 20,
        :ymax                 => 5000.0,
        :nelz                 => 15,
        :zmax                 => 3000.0,

        # Physics
        :lvisc                => true,
        :visc_model           => SMAG(),
        :lsaturation          => true,
    )
    return inputs
end
```

### Example 2: Multi-Level Comparison

Compare forecasts from different times:

```julia
# Run 3 simulations with different ERA5 times
for time_idx in 1:3
    inputs[:era5_time_index] = time_idx
    # Run simulation
end
```

### Example 3: Regional LES with Subsetting

```julia
:lERA5                => true,
:era5_file            => "./data_files/europe_large.nc",
:era5_lat_range       => (48.0, 52.0),   # Focus on Netherlands
:era5_lon_range       => (3.0, 7.0),
:era5_target_lat      => 50.0,
:era5_target_lon      => 5.0,
```

---

## API Reference

### Main Functions

#### `read_era5_netcdf`

```julia
read_era5_netcdf(filename::String;
                 var_names::Dict = Dict(),
                 time_index::Int = 1,
                 lat_range::Tuple = (nothing, nothing),
                 lon_range::Tuple = (nothing, nothing),
                 pressure_levels::Union{Vector, Nothing} = nothing)
```

**Description:** Read ERA5 data from NetCDF file.

**Returns:** `St_ERA5Data` structure with all atmospheric fields.

**Example:**
```julia
era5 = read_era5_netcdf("data.nc", time_index=1,
                        lat_range=(30.0, 50.0),
                        lon_range=(-10.0, 10.0))
```

#### `interpolate_era5_to_1d_column`

```julia
interpolate_era5_to_1d_column(era5_data::St_ERA5Data,
                              target_lat::Real,
                              target_lon::Real,
                              target_heights::Vector{<:Real})
```

**Description:** Extract 1D vertical profile at specific location.

**Returns:** Dictionary with 1D profiles for all variables.

**Example:**
```julia
heights = 0:100:3000  # 0 to 3 km, every 100m
profile = interpolate_era5_to_1d_column(era5, 40.0, -105.0, heights)
```

#### `interpolate_era5_to_mesh`

```julia
interpolate_era5_to_mesh(backend, era5_data::St_ERA5Data,
                         mesh, target_lat::Real, target_lon::Real)
```

**Description:** Interpolate ERA5 to Jexpresso mesh points.

**Returns:** Dictionary with fields on mesh.

**Example:**
```julia
era5_mesh = interpolate_era5_to_mesh(CPU(), era5, mesh, 40.0, -105.0)
```

#### `export_era5_to_sounding`

```julia
export_era5_to_sounding(era5_1d::Dict, output_file::String)
```

**Description:** Export 1D profile to Jexpresso sounding format.

**Example:**
```julia
export_era5_to_sounding(profile, "./data_files/era5_sounding.dat")
```

### Data Structures

#### `St_ERA5Data`

Main structure holding ERA5 fields:

```julia
Base.@kwdef mutable struct St_ERA5Data{T <: AbstractFloat}
    temperature::Array{T}          # [K]
    pressure::Array{T}              # [Pa]
    u_wind::Array{T}                # [m/s]
    v_wind::Array{T}                # [m/s]
    w_wind::Array{T}                # [m/s]
    specific_humidity::Array{T}     # [kg/kg]
    geopotential::Array{T}          # [m²/s²]
    surface_pressure::Array{T}      # [Pa]
    heights::Array{T}               # [m]
    time::Array{Float64}
    lat::Array{T}                   # [degrees]
    lon::Array{T}                   # [degrees]
end
```

---

## Troubleshooting

### Common Issues

#### 1. "Variable not found in NetCDF file"

**Problem:** ERA5 file uses different variable names.

**Solution:** Specify custom variable names:
```julia
:era5_var_names => Dict("T" => "t", "u" => "u", ...)
```

Check your file with:
```bash
ncdump -h era5_file.nc
```

#### 2. "Dimension mismatch"

**Problem:** ERA5 file dimensions don't match expected format.

**Solution:**
- Ensure file has dimensions: (lon, lat, level, time)
- Check with: `ncdump -h era5_file.nc`
- Reorder if needed using CDO:
  ```bash
  cdo -f nc4 copy input.nc output.nc
  ```

#### 3. "Out of memory"

**Problem:** Large ERA5 files exceed available RAM.

**Solution:**
- Use spatial subsetting:
  ```julia
  :era5_lat_range => (min_lat, max_lat),
  :era5_lon_range => (min_lon, max_lon),
  ```
- Select fewer pressure levels
- Use smaller domain for download

#### 4. "Interpolation gives NaN values"

**Problem:** Target location outside ERA5 data range.

**Solution:**
- Verify lat/lon are within your ERA5 domain
- Check that target heights are within ERA5 vertical range
- ERA5 heights are computed from geopotential

#### 5. "File download fails"

**Problem:** CDS API issues.

**Solution:**
- Check API key in `~/.cdsapirc`
- Verify CDS Terms & Conditions are accepted
- Check CDS status: [https://cds.climate.copernicus.eu](https://cds.climate.copernicus.eu)
- Try smaller request (less data)

### Debug Tips

Enable detailed output:

```julia
# In initialize.jl
@info "ERA5 temperature range: $(minimum(era5_fields.temperature)) - $(maximum(era5_fields.temperature))"
@info "ERA5 wind range: $(minimum(era5_fields.u_wind)) - $(maximum(era5_fields.u_wind))"
```

Visualize ERA5 before simulation:

```julia
# Export to sounding format
profile = interpolate_era5_to_1d_column(era5, lat, lon, 0:50:3000)
export_era5_to_sounding(profile, "debug_sounding.dat")

# Plot with your favorite tool
```

---

## Best Practices

### 1. **Data Quality**

- Download ERA5 data closest in time to your case study
- Use 3-hourly or hourly data for transient simulations
- Include surface pressure for better vertical interpolation
- Verify data with visualization tools before running

### 2. **Domain Configuration**

- LES domain should be smaller than ERA5 spatial coverage
- Vertical domain should be within ERA5 height range (typically < 10 km)
- Use at least 5-10 ERA5 pressure levels covering your domain

### 3. **Initialization**

- Add small perturbations to trigger turbulence (done automatically)
- Consider spin-up time (1-2 hours) for LES to develop
- Save initial state for reproducibility

### 4. **Resolution Matching**

- ERA5 (~31 km) provides mesoscale background
- LES domain typically 1-10 km with ~100 m resolution
- This scale separation is physically appropriate

### 5. **Time Integration**

- Start with small time steps (Δt ~ 0.5-1 s)
- Monitor CFL condition
- Use adaptive time stepping if available

### 6. **Validation**

- Compare LES output to ERA5 at larger scales
- Check energy spectra
- Verify conservation properties
- Compare with observations if available

---

## Citation

If you use the ERA5 extension in your research, please cite:

1. **Jexpresso:**
   ```
   Marras, S., et al. (2024). Jexpresso: A Julia Framework for
   Spectral Element Methods. [Add DOI when available]
   ```

2. **ERA5:**
   ```
   Hersbach, H., et al. (2020). The ERA5 global reanalysis.
   Quarterly Journal of the Royal Meteorological Society, 146(730), 1999-2049.
   https://doi.org/10.1002/qj.3803
   ```

---

## Support and Contributing

- **Issues:** Report bugs at [GitHub Issues](https://github.com/smarras79/Jexpresso/issues)
- **Questions:** Contact Jexpresso team
- **Contributing:** Pull requests welcome!

---

## Version History

- **v1.0 (2025-11-30):** Initial release
  - Basic ERA5 NetCDF reading
  - 1D and 3D interpolation
  - LES initialization
  - Example test case

---

## License

This extension is part of Jexpresso and distributed under the same license.

---

**Happy Simulating! 🌦️**
