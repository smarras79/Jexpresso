# Turbulence Statistics in Jexpresso

## Overview

This module provides turbulence statistics computation capabilities for Jexpresso, analogous to the dod2d `mod_aver.f90` module. It computes time-averaged quantities and second-order statistics for turbulence analysis in LES (Large Eddy Simulation) applications.

## Features

The module computes the following turbulence statistics:

### First-Order Statistics (Mean Quantities)
- **`<u>`, `<v>`, `<w>`**: Time-averaged velocity components
- **`<ρ>`**: Time-averaged density
- **`<p>`**: Time-averaged pressure
- **`<μ_eff>`**: Time-averaged effective viscosity (molecular + SGS)

### Second-Order Statistics
- **`<u²>`, `<v²>`, `<w²>`**: Velocity variance components (for computing RMS fluctuations)
- **`<p²>`**: Pressure variance
- **Reynolds Stress Components**:
  - 2D: `<u'v'>`
  - 3D: `<u'v'>`, `<u'w'>`, `<v'w'>`

### Derived Quantities
- **RMS velocity fluctuations**: `u'_rms = sqrt(<u²> - <u>²)`
- **Reynolds stresses**: `<u'v'> = <uv> - <u><v>`

## File Structure

```
src/kernel/globalStructs.jl              # Data structure: St_turbulence_stats
src/io/turbulence_statistics.jl          # Main module with all functions
```

## Usage

### 1. Allocate Statistics Structure

Add turbulence statistics to your problem setup:

```julia
# In your problem initialization
turb_stats = allocate_turbulence_stats(SD, npoin, T, backend;
                                       l_turbulence_stats=true)
```

### 2. Initialize Statistics

Start accumulating statistics at a specific time:

```julia
# Initialize at time t_start (e.g., after spin-up)
t_start = 100.0  # Start statistics after flow develops
initialize_turbulence_stats!(turb_stats, t_start)
```

### 3. Accumulate During Time Stepping

In your time integration loop, accumulate statistics at each step:

```julia
# In the time stepping callback or RHS function
function accumulate_callback!(integrator)
    # Extract current state
    uaux = integrator.p.uaux
    press = integrator.p.q.press
    mueff = integrator.p.mueff  # Or nothing if not available

    # Get current time and time step
    t = integrator.t
    dt = integrator.dt

    # Only accumulate if we're past the start time
    if t >= turb_stats.elapsed_time[]
        accumulate_statistics!(turb_stats, uaux, press, mueff,
                              dt, mesh.npoin, ndim, SD)
    end
end
```

### 4. Normalize and Output

At the end of the simulation or at regular intervals:

```julia
# Normalize accumulated statistics
normalize_statistics!(turb_stats)

# Print summary to console
print_turbulence_stats_summary(turb_stats, ndim)

# Write to ASCII file for analysis
write_turbulence_stats_ascii(turb_stats, mesh, OUTPUT_DIR,
                             "turbulence_stats", ndim, t_current)

# Compute derived quantities
vel_rms = similar(turb_stats.avvel)
compute_rms_fluctuations!(vel_rms, turb_stats, ndim)

reynolds_stress = similar(turb_stats.avvex)
compute_reynolds_stresses!(reynolds_stress, turb_stats, ndim)
```

### 5. Reset for New Averaging Window

To compute statistics over multiple time windows:

```julia
# After writing output, reset for next window
reset_turbulence_stats!(turb_stats)
```

## Integration with Time Stepping

### Method 1: Using DiscreteCallback (Recommended)

Add a callback to the time integrator:

```julia
# Define callback for statistics accumulation
function stats_condition(u, t, integrator)
    # Accumulate every time step
    return true
end

function stats_affect!(integrator)
    accumulate_statistics!(turb_stats,
                          integrator.p.uaux,
                          integrator.p.q.press,
                          nothing,  # mueff if available
                          integrator.dt,
                          mesh.npoin,
                          ndim,
                          SD)
end

stats_cb = DiscreteCallback(stats_condition, stats_affect!)

# Add to integrator callbacks
callbacks = CallbackSet(output_cb, stats_cb)
```

### Method 2: Direct Integration in RHS

Call `accumulate_statistics!` directly in your RHS function:

```julia
function rhs_with_stats!(dU, U, params, t)
    # Compute RHS as usual
    compute_rhs!(dU, U, params, t)

    # Accumulate statistics
    if params.turb_stats.l_enabled[]
        accumulate_statistics!(params.turb_stats,
                              params.uaux,
                              params.q.press,
                              params.mueff,
                              params.dt,
                              params.mesh.npoin,
                              params.ndim,
                              params.SD)
    end
end
```

## Example: 2D Turbulent Channel Flow

```julia
using Jexpresso
include("src/io/turbulence_statistics.jl")

# Problem setup
SD = NSD_2D()
ndim = 2
npoin = mesh.npoin
T = Float64
backend = CPU()

# Allocate statistics
turb_stats = allocate_turbulence_stats(SD, npoin, T, backend;
                                       l_turbulence_stats=true)

# Initialize after spin-up time
t_spinup = 50.0
initialize_turbulence_stats!(turb_stats, t_spinup)

# Run simulation with statistics accumulation
# (integrate statistics in time loop as shown above)

# After simulation
normalize_statistics!(turb_stats)

# Compute derived quantities
vel_rms = zeros(T, npoin, ndim)
compute_rms_fluctuations!(vel_rms, turb_stats, ndim)

reynolds_stress = zeros(T, npoin, 1)  # Only <u'v'> in 2D
compute_reynolds_stresses!(reynolds_stress, turb_stats, ndim)

# Output results
write_turbulence_stats_ascii(turb_stats, mesh, OUTPUT_DIR,
                             "channel_stats", ndim, t_final)

# Print summary
print_turbulence_stats_summary(turb_stats, ndim)
```

## Output Format

### ASCII Output

The ASCII file contains space-separated columns:

**1D Output:**
```
# ip | x | <u> | <ρ> | <p> | <u*u> | <p*p> | <μ_eff>
```

**2D Output:**
```
# ip | x | y | <u> | <v> | <ρ> | <p> | <u*u> | <v*v> | <u*v> | <p*p> | <μ_eff>
```

**3D Output:**
```
# ip | x | y | z | <u> | <v> | <w> | <ρ> | <p> | <u*u> | <v*v> | <w*w> | <u*v> | <u*w> | <v*w> | <p*p> | <μ_eff>
```

## Comparison with dod2d mod_aver.f90

| dod2d Feature | Jexpresso Equivalent |
|---------------|---------------------|
| `avvel(ipoin,idime)` | `turb_stats.avvel[ip, idim]` |
| `avve2(ipoin,idime)` | `turb_stats.avve2[ip, idim]` |
| `avvex(ipoin,1:3)` | `turb_stats.avvex[ip, 1:3]` |
| `avpr`, `avpr2` | `turb_stats.avpress[ip]`, `turb_stats.avpr2[ip]` |
| `avrho` | `turb_stats.avrho[ip]` |
| `avmueff` | `turb_stats.avmueff[ip]` |
| `acutim` | `turb_stats.acutim[]` |
| `favre_average` | `accumulate_statistics!` |
| `eval_average_window` | `normalize_statistics!` |

## Performance Considerations

1. **Memory**: Statistics arrays are the same size as the solution arrays, so memory usage approximately doubles when statistics are enabled.

2. **Computational Cost**: Accumulation adds minimal overhead (~5-10% depending on problem size) as it only involves simple arithmetic operations.

3. **GPU Support**: The module uses `KernelAbstractions.zeros` for array allocation, making it compatible with GPU backends. For optimal GPU performance, consider kernel-based accumulation.

4. **MPI Parallelization**: Each MPI rank accumulates statistics for its local points. For domain-averaged quantities, add MPI reduction operations.

## Advanced Usage

### Spatial Averaging

For channel flows or other homogeneous directions:

```julia
# Average statistics in x-direction (assuming structured grid)
function compute_spanwise_average(turb_stats, mesh, ndim)
    # Implementation depends on mesh structure
    # Return 1D profiles of statistics
end
```

### Multiple Time Windows

```julia
# Accumulate statistics over multiple windows
for iwindow = 1:nwindows
    # Reset for new window
    reset_turbulence_stats!(turb_stats)

    # Run for window duration
    t_window = 10.0
    # ... time integration ...

    # Normalize and write
    normalize_statistics!(turb_stats)
    write_turbulence_stats_ascii(turb_stats, mesh, OUTPUT_DIR,
                                 "stats_window_$iwindow", ndim, t)
end
```

### Conditional Statistics

```julia
# Accumulate statistics only in specific regions
function accumulate_conditional!(turb_stats, uaux, press, mueff, dt,
                                region_mask, npoin, ndim, SD)
    # Modify accumulate_statistics! to only accumulate where region_mask[ip] == true
end
```

## References

1. **dod2d**: Original Fortran implementation in `mod_aver.f90`
2. **Pope, S.B.** (2000): *Turbulent Flows*, Cambridge University Press
3. **Sagaut, P.** (2006): *Large Eddy Simulation for Incompressible Flows*, Springer

## Support

For questions or issues:
- Check the Jexpresso documentation
- See examples in `problems/equations/CompEuler/giga_les/`
- Contact the development team

## Future Enhancements

Planned features:
- [ ] Spatial (homogeneous direction) averaging
- [ ] Higher-order moments (skewness, kurtosis)
- [ ] Two-point correlations
- [ ] Spectral analysis integration
- [ ] VTK output integration
- [ ] GPU-optimized kernel accumulation
- [ ] MPI-aware global statistics
