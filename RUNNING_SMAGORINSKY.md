# How to Run Jexpresso with Smagorinsky Turbulence Model

This guide provides step-by-step instructions for running Jexpresso with the newly implemented Smagorinsky turbulence model.

## Prerequisites

1. **Julia Installation** (version 1.11.2 or higher)
2. **Jexpresso Dependencies** (see main README.md)
3. **MPI** (for parallel execution)

## Quick Start (Simplest Method)

### Option 1: Using the Run Script

From the Jexpresso root directory:

```bash
julia --project=. run_smagorinsky.jl
```

This will run the pre-configured `theta_smagorinsky` test case with Smagorinsky model enabled.

### Option 2: From Julia REPL

```julia
# Start Julia with the project environment
julia --project=.

# In the Julia REPL:
julia> push!(empty!(ARGS), "CompEuler", "theta_smagorinsky")
julia> include("src/Jexpresso.jl")
```

### Option 3: Using MPI for Parallel Execution

```bash
mpirun -np 4 julia --project=. run_smagorinsky.jl
```

Replace `4` with the number of processors you want to use.

## Understanding the Test Case

The `theta_smagorinsky` case is located at:
```
problems/equations/CompEuler/theta_smagorinsky/
```

This case solves the **compressible Euler equations** with:
- **Potential temperature formulation** (ρθ)
- **Smagorinsky LES model** for subgrid-scale turbulence
- **Viscous energy diffusion**: ∇·(κ∇θ)

## Configuration Parameters

The test case is configured in `user_inputs.jl`:

```julia
:lsmagorinsky    => true,   # Enable Smagorinsky model
:C_s             => 0.18,   # Smagorinsky constant
:Pr_t            => 0.85,   # Turbulent Prandtl number
:lmolecular_visc => true,   # Include molecular viscosity
```

### Key Parameters to Adjust

**Smagorinsky Constant (C_s)**:
- **Range**: 0.1 - 0.2
- **Atmospheric flows**: 0.18 - 0.2
- **Engineering flows**: 0.1 - 0.15
- **Default**: 0.18

**Turbulent Prandtl Number (Pr_t)**:
- **Range**: 0.7 - 1.0
- **Typical**: 0.85
- **Default**: 0.85

**Time Step (Δt)**:
- **Default**: 0.4 seconds
- **Reduce** if simulation becomes unstable
- **CFL constraint**: Δt ≤ C × h² / max(κ_total)

**Polynomial Order (nop)**:
- **Default**: 4 (5th order accuracy)
- **Higher values** (5-7) for better resolution
- **Lower values** (3) for faster runs

## Output and Results

### Output Location

Results are saved to:
```
problems/equations/CompEuler/theta_smagorinsky/output/
```

Or if you set `:output_dir` in user_inputs.jl, results go there.

### Output Files

- **VTK files**: For visualization (ParaView, VisIt)
- **Diagnostic data**: Time series, energy, etc.
- **user_inputs.jl**: Copy of configuration used

### Visualizing Results

Open VTK files in **ParaView**:
```bash
paraview problems/equations/CompEuler/theta_smagorinsky/output/*.vtu
```

**Variables to plot**:
- `rho` - Density (ρ)
- `rhou` - x-momentum (ρu)
- `rhov` - y-momentum (ρv)
- `rhotheta` - Potential temperature (ρθ)

## Step-by-Step Detailed Instructions

### 1. Navigate to Jexpresso Directory

```bash
cd /path/to/Jexpresso
```

### 2. Activate Julia Environment

```bash
julia --project=.
```

This loads all required packages from Project.toml.

### 3. (Optional) Install/Update Dependencies

If first time or dependencies changed:

```julia
julia> using Pkg
julia> Pkg.instantiate()
```

### 4. Run the Simulation

**Method A - Simple**:
```julia
julia> include("run_smagorinsky.jl")
```

**Method B - Standard Jexpresso**:
```julia
julia> push!(empty!(ARGS), "CompEuler", "theta_smagorinsky")
julia> include("src/Jexpresso.jl")
```

**Method C - Command Line**:
```bash
julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "theta_smagorinsky"); include("src/Jexpresso.jl")'
```

### 5. Monitor Progress

You'll see output like:
```
==================================================================
  Jexpresso - Spectral Element with Smagorinsky Turbulence Model
==================================================================

Running Configuration:
  Equations: CompEuler
  Case:      theta_smagorinsky
  Features:  Smagorinsky LES model enabled

Initializing mesh...
Initializing solution...
Time integration starting...
  t = 0.0 / 1000.0
  t = 100.0 / 1000.0
  ...
```

### 6. Check Results

After completion:
```bash
ls problems/equations/CompEuler/theta_smagorinsky/output/
```

## Running Different Test Cases

### Modify Existing Case

Edit parameters in:
```
problems/equations/CompEuler/theta_smagorinsky/user_inputs.jl
```

Then re-run.

### Create New Case

1. Copy the template:
```bash
cp -r problems/equations/CompEuler/theta_smagorinsky problems/equations/CompEuler/my_case
```

2. Modify `my_case/user_inputs.jl`, `initialize.jl`, etc.

3. Run:
```julia
julia> push!(empty!(ARGS), "CompEuler", "my_case")
julia> include("src/Jexpresso.jl")
```

## Integrating Smagorinsky into Existing Cases

To add Smagorinsky to an existing Jexpresso case:

### 1. Add Parameters to user_inputs.jl

```julia
:lsmagorinsky    => true,
:C_s             => 0.18,
:Pr_t            => 0.85,
:lmolecular_visc => true,
```

### 2. Include Smagorinsky Modules

In your case directory or in the driver, add:

```julia
include("../../../src/kernel/Turbulence/Smagorinsky.jl")
include("../theta/user_viscous_flux.jl")
```

### 3. Modify RHS Computation

In the RHS function (or create a wrapper), add:

```julia
if inputs[:lsmagorinsky]
    # Compute Smagorinsky viscous terms
    # (See examples/smagorinsky_integration_example.jl for details)
end
```

See `examples/smagorinsky_integration_example.jl` for complete implementation.

## Common Issues and Solutions

### Issue 1: "UndefVarError: Smagorinsky not defined"

**Solution**: The Smagorinsky module isn't loaded. Add to your case:
```julia
include("src/kernel/Turbulence/Smagorinsky.jl")
```

### Issue 2: Simulation is Unstable

**Solutions**:
- Reduce time step: `:Δt => 0.1`
- Reduce Smagorinsky constant: `:C_s => 0.1`
- Increase polynomial order: `:nop => 5`

### Issue 3: "No such file or directory: meshes/..."

**Solution**: Check mesh path in user_inputs.jl:
```julia
:gmsh_filename => "./meshes/gmsh_grids/hexa_TFI_10x10.msh"
```

Make sure the mesh file exists.

### Issue 4: Very Slow Execution

**Solutions**:
- Use coarser mesh
- Reduce polynomial order
- Disable Smagorinsky for testing: `:lsmagorinsky => false`
- Use parallel execution: `mpirun -np 4 julia ...`

### Issue 5: "LoadError: ArgumentError: Package MPI not found"

**Solution**: Install dependencies:
```julia
using Pkg
Pkg.instantiate()
```

## Performance Optimization

### For Faster Runs

1. **Reduce domain size**: Smaller mesh
2. **Lower polynomial order**: `nop = 3`
3. **Larger time step**: Increase `:Δt` (check stability)
4. **Shorter simulation**: Reduce `:tend`
5. **Less frequent output**: Increase `:diagnostics_at_times` interval

### For Better Accuracy

1. **Higher polynomial order**: `nop = 5` or `6`
2. **Finer mesh**: More elements
3. **Smaller time step**: Reduce `:Δt`
4. **Tune C_s**: Try range 0.15 - 0.20

### For Parallel Execution

```bash
# Use 4 processors
mpirun -np 4 julia --project=. run_smagorinsky.jl

# Use 8 processors
mpirun -np 8 julia --project=. run_smagorinsky.jl
```

## Advanced: GPU Execution

For GPU acceleration (requires CUDA or Metal):

In `user_inputs.jl`:
```julia
:backend => CUDABackend()  # For NVIDIA GPUs
# or
:backend => MetalBackend() # For Apple Silicon
```

**Note**: Smagorinsky GPU kernels not yet implemented. This is a future enhancement.

## Testing the Implementation

### Quick Validation Test

Run a short simulation to verify everything works:

```julia
# In user_inputs.jl, set:
:tend => 10.0              # Short run
:diagnostics_at_times => (0:5:10)
```

Then run and check for errors.

### Full Validation

For scientific validation:

1. **Energy conservation**: Monitor total energy over time
2. **Turbulence decay**: For decaying turbulence, check E(t) ∝ t^(-n)
3. **Compare with DNS**: If DNS data available
4. **Grid convergence**: Run with different mesh sizes

## Example: Complete Workflow

```bash
# 1. Clone/navigate to Jexpresso
cd ~/Jexpresso

# 2. Start Julia
julia --project=.

# 3. In Julia REPL
julia> # Install dependencies (first time only)
julia> using Pkg
julia> Pkg.instantiate()

julia> # Run simulation
julia> include("run_smagorinsky.jl")

# Wait for completion...

# 4. Visualize results
julia> exit()
$ paraview problems/equations/CompEuler/theta_smagorinsky/output/*.vtu
```

## Getting Help

1. **Documentation**:
   - `docs/spectral_element_energy_smagorinsky.md` - Mathematical details
   - `docs/SMAGORINSKY_README.md` - Implementation guide

2. **Examples**:
   - `examples/smagorinsky_integration_example.jl` - Code examples

3. **Jexpresso Main Docs**:
   - See main README.md for general Jexpresso usage

4. **Issues**:
   - Check existing GitHub issues
   - Open new issue with error details

## Summary

**Quickest way to run**:
```bash
julia --project=. run_smagorinsky.jl
```

**With MPI (4 processors)**:
```bash
mpirun -np 4 julia --project=. run_smagorinsky.jl
```

**Key files**:
- **Run script**: `run_smagorinsky.jl`
- **Test case**: `problems/equations/CompEuler/theta_smagorinsky/`
- **Configuration**: `user_inputs.jl` in test case directory
- **Output**: `theta_smagorinsky/output/`

**What it does**:
- Solves compressible Euler equations with potential temperature
- Applies Smagorinsky LES model for turbulence
- Adds viscous diffusion to energy equation: ∇·(κ∇θ)
- Outputs VTK files for visualization

Enjoy running your Smagorinsky simulations! 🚀
