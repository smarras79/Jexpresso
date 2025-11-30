# ERA5 LES Test Case

This test case demonstrates using ERA5 reanalysis data for Large Eddy Simulation (LES) weather forecasting with Jexpresso.

## Overview

- **Purpose:** Real-world weather forecasting using ERA5 initial conditions
- **Domain:** 5 km × 5 km × 3 km LES domain
- **Resolution:** ~250 m horizontal, ~200 m vertical
- **Duration:** 1 hour simulation (configurable)
- **Physics:** Compressible Euler with θ equation, moisture, Smagorinsky SGS

## Quick Start

### 1. Get ERA5 Data

Download ERA5 data using the CDS API (see `scripts/download_era5.py`):

```bash
cd scripts
python download_era5.py
```

Or manually from: https://cds.climate.copernicus.eu

Place the downloaded file in `data_files/era5_example.nc`

### 2. Configure Case

Edit `user_inputs.jl` to set:
- ERA5 file path
- Target latitude/longitude
- Domain size
- Simulation duration

### 3. Run Simulation

```bash
julia --project=. src/run.jl CompEuler era5_les
```

## File Description

- `user_inputs.jl` - Configuration parameters
- `initialize.jl` - Initialize fields from ERA5 data
- `user_bc.jl` - Boundary conditions (free-slip lateral, no-slip bottom)
- `user_flux.jl` - Compressible Euler fluxes
- `user_source.jl` - Source terms (gravity)
- `user_primitives.jl` - Conservative ↔ primitive variable conversion

## Configuration Options

### Basic Setup

```julia
:lERA5                => true,
:era5_file            => "./data_files/era5_example.nc",
:era5_target_lat      => 40.0,    # Your location
:era5_target_lon      => -105.0,
```

### Advanced Options

See `../../../docs/ERA5_EXTENSION_GUIDE.md` for:
- Custom variable names
- Spatial subsetting
- Pressure level selection
- Nudging/forcing options

## Expected Output

The simulation produces VTK files in `./output/CompEuler/era5_les/` containing:
- ρ (density)
- u, v, w (velocity components)
- θ (potential temperature)
- qt (total water mixing ratio)
- ql (liquid water mixing ratio)
- P (pressure)

Visualize with ParaView, VisIt, or similar tools.

## Tips

1. **First run:** Start with a short simulation (tend = 600 s) to verify setup
2. **Resolution:** Adjust `nelx`, `nely`, `nelz` based on your computational resources
3. **Time step:** Reduce `Δt` if simulation becomes unstable
4. **Perturbations:** Small random perturbations are added automatically to trigger turbulence

## For More Information

See comprehensive documentation:
- `docs/ERA5_EXTENSION_GUIDE.md` - Complete ERA5 extension guide
- `docs/getting_started.md` - Jexpresso basics
- ERA5 documentation: https://confluence.ecmwf.int/display/CKB/ERA5

## Troubleshooting

**Problem:** "ERA5 file not found"
- Ensure ERA5 data is downloaded and path is correct in `user_inputs.jl`

**Problem:** "Variable not found in NetCDF"
- Check variable names with: `ncdump -h era5_file.nc`
- Specify custom names in `:era5_var_names`

**Problem:** Simulation crashes immediately
- Check CFL condition (reduce `Δt`)
- Verify ERA5 data quality
- Start with smaller domain

## Contact

For questions or issues with the ERA5 extension, please contact the Jexpresso development team.
