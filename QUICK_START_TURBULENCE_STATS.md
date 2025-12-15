# Turbulence Statistics - Quick Start Guide

## ✅ Integration Complete!

The turbulence statistics module is now **fully integrated** into Jexpresso's core infrastructure. You can enable it simply by adding parameters to your `user_inputs.jl` file.

## 🚀 How to Use (3 Simple Steps)

### Step 1: Enable in `user_inputs.jl`

Add these lines to your `user_inputs()` function (e.g., in `problems/equations/CompEuler/giga_les/user_inputs.jl`):

```julia
function user_inputs()
    inputs = Dict(
        # ... your existing parameters ...

        #---------------------------------------------------------------------------
        # Turbulence statistics parameters
        #---------------------------------------------------------------------------
        :l_turbulence_stats         => true,           # Enable statistics
        :turb_stats_start_time      => 100.0,          # Start after spin-up
        :turb_stats_output_times    => (200.0, 500.0, 1000.0, 2000.0),  # Output times
        :turb_stats_reset_after_output => false,       # Keep accumulating across outputs

        # ... rest of your parameters ...
    )
    return inputs
end
```

### Step 2: Run Your Simulation

```bash
julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "giga_les"); include("src/Jexpresso.jl")'
```

### Step 3: Find Your Results

Statistics will be automatically written to your output directory:

```
output_dir/
├── turbulence_stats_t000200.000000.dat
├── turbulence_stats_t000500.000000.dat
├── turbulence_stats_t001000.000000.dat
└── turbulence_stats_t002000.000000.dat
```

That's it! No other code changes needed.

---

## 📊 What Gets Computed

### Automatically Computed Statistics

The system computes these quantities at each timestep after `turb_stats_start_time`:

**First-Order (Mean) Quantities:**
- `<u>`, `<v>`, `<w>` - Time-averaged velocity components
- `<ρ>` - Time-averaged density
- `<p>` - Time-averaged pressure
- `<μ_eff>` - Time-averaged effective viscosity (if available)

**Second-Order Quantities:**
- `<u²>`, `<v²>`, `<w²>` - Velocity variance (for RMS)
- `<p²>` - Pressure variance
- `<u*v>`, `<u*w>`, `<v*w>` - Reynolds stress components

**Output:** All statistics are written in space-separated ASCII format for easy analysis.

---

## ⚙️ Configuration Options

### Required Parameters

| Parameter | Description | Example |
|-----------|-------------|---------|
| `:l_turbulence_stats` | Enable/disable statistics | `true` |

### Optional Parameters

| Parameter | Description | Default |
|-----------|-------------|---------|
| `:turb_stats_start_time` | Time to start accumulation | `0.0` |
| `:turb_stats_output_times` | Times to write statistics | `diagnostics_at_times` |
| `:turb_stats_reset_after_output` | Reset after each output | `false` |

---

## 📁 Output Format

### 2D Example Output

```
# Turbulence Statistics
# Time: 200.000000
# Accumulated time: 100.000000
# Number of points: 12800
# Spatial dimensions: 2
#
# ip | x | y | <u> | <v> | <ρ> | <p> | <u*u> | <v*v> | <u*v> | <p*p> | <μ_eff>
1 0.000000e+00 0.000000e+00 1.234e+01 2.345e-01 1.200e+00 1.013e+05 1.567e+02 8.901e-01 3.456e+00 1.025e+10 5.000e+02
2 1.000000e+02 0.000000e+00 1.235e+01 2.347e-01 1.201e+00 1.013e+05 1.568e+02 8.905e-01 3.458e+00 1.025e+10 5.001e+02
...
```

### 3D Example Output

```
# ip | x | y | z | <u> | <v> | <w> | <ρ> | <p> | <u*u> | <v*v> | <w*w> | <u*v> | <u*w> | <v*w> | <p*p> | <μ_eff>
```

---

## 🔬 Analysis Examples

### Load Statistics in Python

```python
import numpy as np

# Load statistics file
data = np.loadtxt('turbulence_stats_t001000.000000.dat', skiprows=6)

# Extract columns (2D example)
ip, x, y, u, v, rho, p, uu, vv, uv, pp, mu = data.T

# Compute RMS fluctuations
u_rms = np.sqrt(uu - u**2)
v_rms = np.sqrt(vv - v**2)

# Compute Reynolds stress
reynolds_stress = uv - u * v

# Plot vertical profile (assuming structured grid)
import matplotlib.pyplot as plt
plt.plot(u, y, label='Mean velocity')
plt.plot(u_rms, y, label='RMS fluctuation')
plt.xlabel('Velocity')
plt.ylabel('Height')
plt.legend()
plt.show()
```

### Load Statistics in Julia

```julia
using DelimitedFiles

# Load statistics
data = readdlm("turbulence_stats_t001000.000000.dat", ' ', Float64; skipstart=6)

# Extract columns (2D)
x, y, u, v = data[:, 2], data[:, 3], data[:, 4], data[:, 5]
uu, vv, uv = data[:, 8], data[:, 9], data[:, 10]

# Compute derived quantities
u_rms = sqrt.(uu .- u.^2)
v_rms = sqrt.(vv .- v.^2)
reynolds_stress = uv .- u .* v
```

---

## 🎯 Common Use Cases

### Case 1: Basic Statistics Collection

```julia
# In user_inputs.jl
:l_turbulence_stats    => true,
:turb_stats_start_time => 100.0,  # After spin-up
```

Statistics accumulate from t=100 onwards and output at diagnostic times.

### Case 2: Windowed Averaging

```julia
# In user_inputs.jl
:l_turbulence_stats         => true,
:turb_stats_output_times    => (100.0, 200.0, 300.0, 400.0),
:turb_stats_reset_after_output => true,  # Reset after each window
```

Computes separate statistics for windows: [0,100], [100,200], [200,300], [300,400].

### Case 3: Long-Time Averaging

```julia
# In user_inputs.jl
:l_turbulence_stats    => true,
:turb_stats_start_time => 500.0,   # Start after full spin-up
:turb_stats_output_times => (5000.0,),  # Single output at end
```

Accumulates statistics over very long time (4500 time units) for well-converged averages.

---

## 🔍 Validation Checklist

### For Homogeneous Isotropic Turbulence

- [ ] Check `<u²> ≈ <v²> ≈ <w²>` (isotropy)
- [ ] Verify energy decay: `<u²> + <v²> + <w²>` decreases with time
- [ ] Reynolds stresses `<u'v'>`, `<u'w'>`, `<v'w'>` should be ~ 0

### For Channel Flow

- [ ] Mean velocity `<u>` should follow log-law in inner region
- [ ] Reynolds stress `<u'v'>` should match DNS/LES reference data
- [ ] Peak RMS occurs near wall (y+ ~ 15)

### For Any Flow

- [ ] Mass conservation: domain-integrated `<ρ>` constant
- [ ] Statistics converge with longer averaging time
- [ ] Results independent of `turb_stats_output_times` frequency

---

## 🐛 Troubleshooting

### Problem: No statistics files generated

**Solution:** Check that:
1. `:l_turbulence_stats => true` is set
2. Simulation time reaches `:turb_stats_start_time`
3. Simulation time reaches at least one output time
4. Check terminal output for turbulence stats messages

### Problem: Statistics are all zeros

**Solution:**
- Ensure `turb_stats_start_time < tend`
- Verify flow is actually turbulent (check instantaneous fields)
- Check that pressure field is being computed

### Problem: Statistics seem wrong

**Solution:**
- Increase accumulation time (statistics need sufficient samples)
- Verify spin-up time is adequate before starting statistics
- Compare with instantaneous field values for sanity check

---

## 📚 Additional Resources

- **Full Documentation**: `docs/TURBULENCE_STATISTICS.md`
- **Examples**: `examples/turbulence_statistics_integration_example.jl`
- **Implementation**: `src/io/turbulence_statistics.jl`
- **Main README**: `TURBULENCE_STATISTICS_README.md`

---

## 💡 Performance Notes

- **Memory overhead**: ~2x solution array size when enabled
- **Computational cost**: ~5-10% for accumulation
- **I/O cost**: Minimal (ASCII files are small)
- **Recommendation**: Enable only when needed for analysis

---

## ✨ Advanced Features

For advanced usage including:
- Spatial (homogeneous direction) averaging
- Multiple time windows
- Conditional statistics
- Derived quantities (TKE, production, dissipation)
- Custom output formats

See the full documentation in `docs/TURBULENCE_STATISTICS.md` and examples in `examples/turbulence_statistics_integration_example.jl`.

---

**Ready to use!** Just add the parameters to your `user_inputs.jl` and run your simulation. 🎉
