# Smagorinsky Model - Quick Start Guide

## TL;DR - Run It Now!

```bash
cd /path/to/Jexpresso
julia --project=. run_smagorinsky.jl
```

That's it! The simulation will run and output results to:
```
problems/equations/CompEuler/theta_smagorinsky/output/
```

## What Just Happened?

You just ran a **compressible Navier-Stokes** simulation with **Smagorinsky LES turbulence modeling** using **spectral elements**!

The simulation:
- ✅ Solves 2D compressible Euler equations with potential temperature
- ✅ Includes Smagorinsky eddy viscosity for subgrid-scale turbulence
- ✅ Uses high-order (4th order) continuous Galerkin spectral elements
- ✅ Outputs VTK files for ParaView visualization

## With Parallel Execution (Faster!)

```bash
mpirun -np 4 julia --project=. run_smagorinsky.jl
```

Replace `4` with number of CPU cores you want to use.

## Visualize Results

```bash
paraview problems/equations/CompEuler/theta_smagorinsky/output/*.vtu
```

## Key Parameters (In user_inputs.jl)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `:C_s` | 0.18 | Smagorinsky constant (0.1-0.2) |
| `:Pr_t` | 0.85 | Turbulent Prandtl number (0.7-1.0) |
| `:Δt` | 0.4 | Time step [seconds] |
| `:tend` | 1000.0 | Final time [seconds] |
| `:nop` | 4 | Polynomial order (3-6) |

## Change Settings

Edit: `problems/equations/CompEuler/theta_smagorinsky/user_inputs.jl`

Then re-run.

## Alternative Run Methods

### From Julia REPL
```julia
julia --project=.

julia> push!(empty!(ARGS), "CompEuler", "theta_smagorinsky")
julia> include("src/Jexpresso.jl")
```

### Command Line One-Liner
```bash
julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "theta_smagorinsky"); include("src/Jexpresso.jl")'
```

## Output Files

After running, check:
```bash
ls problems/equations/CompEuler/theta_smagorinsky/output/
```

You'll see:
- `*.vtu` - VTK files (open in ParaView)
- `user_inputs.jl` - Copy of settings used
- Diagnostic output files

## Troubleshooting

**Error: Package not found**
```julia
using Pkg; Pkg.instantiate()
```

**Simulation unstable**
- Reduce time step: `:Δt => 0.1`
- Reduce C_s: `:C_s => 0.1`

**Too slow**
- Reduce domain time: `:tend => 100.0`
- Use MPI: `mpirun -np 4 julia ...`

## Learn More

- **Full guide**: `RUNNING_SMAGORINSKY.md`
- **Math details**: `docs/spectral_element_energy_smagorinsky.md`
- **Implementation**: `docs/SMAGORINSKY_README.md`
- **Code examples**: `examples/smagorinsky_integration_example.jl`

## File Structure

```
Jexpresso/
├── run_smagorinsky.jl              ← Main run script
├── QUICKSTART_SMAGORINSKY.md       ← This file
├── RUNNING_SMAGORINSKY.md          ← Detailed guide
├── docs/
│   ├── spectral_element_energy_smagorinsky.md  ← Math theory
│   └── SMAGORINSKY_README.md       ← Implementation details
├── src/kernel/Turbulence/
│   └── Smagorinsky.jl              ← Turbulence model
├── problems/equations/CompEuler/theta/
│   └── user_viscous_flux.jl        ← Viscous fluxes
├── problems/equations/CompEuler/theta_smagorinsky/
│   ├── user_inputs.jl              ← Configuration (edit this!)
│   ├── initialize.jl               ← Initial conditions
│   └── ...                         ← Other case files
└── examples/
    └── smagorinsky_integration_example.jl  ← Integration template
```

## What The Code Does

### Physics
- Compressible Navier-Stokes energy equation:
  ```
  ∂(ρθ)/∂t + ∇·(ρθu) = ∇·(κ∇θ)
  ```

### Turbulence Model
- Smagorinsky eddy viscosity:
  ```
  ν_t = (C_s Δ)² |S|
  ```
  where |S| = strain rate magnitude, Δ = filter width

### Numerical Method
- Continuous Galerkin spectral elements
- Tensor-product Lagrange basis on LGL nodes
- High-order accuracy (order = nop + 1)
- Explicit time integration (SSPRK)

## Customization

### Run Different Case

Copy template:
```bash
cp -r problems/equations/CompEuler/theta_smagorinsky problems/equations/CompEuler/my_case
```

Edit `my_case/user_inputs.jl` and run:
```julia
julia> push!(empty!(ARGS), "CompEuler", "my_case")
julia> include("src/Jexpresso.jl")
```

### Disable Smagorinsky (For Testing)

In `user_inputs.jl`:
```julia
:lsmagorinsky => false,
```

### Change Domain

Edit mesh file in `user_inputs.jl`:
```julia
:gmsh_filename => "./meshes/gmsh_grids/your_mesh.msh"
```

## Performance Tips

**Faster simulations**:
- Coarser mesh (fewer elements)
- Lower polynomial order (`:nop => 3`)
- Larger time step (`:Δt => 0.5`, check stability!)
- Shorter run (`:tend => 100.0`)

**Better accuracy**:
- Finer mesh (more elements)
- Higher polynomial order (`:nop => 5` or `6`)
- Smaller time step (`:Δt => 0.2`)

**Parallel speedup**:
- Use MPI: `mpirun -np N julia ...` where N = # cores
- Optimal N ≈ number of elements / 10

## Expected Runtime

On typical laptop (4 cores):
- Small test (10x10 mesh, 100s): ~1-2 minutes
- Medium run (20x20 mesh, 1000s): ~10-30 minutes
- Large run (40x40 mesh, 1000s): ~1-2 hours

With MPI (4 processes): ~2-3x faster

## Next Steps

1. **Run the example**: `julia --project=. run_smagorinsky.jl`
2. **Check output**: Look in `theta_smagorinsky/output/`
3. **Visualize**: Open VTK files in ParaView
4. **Modify**: Edit `user_inputs.jl` and re-run
5. **Read docs**: For deeper understanding

## Getting Help

- **Detailed instructions**: See `RUNNING_SMAGORINSKY.md`
- **Mathematical theory**: See `docs/spectral_element_energy_smagorinsky.md`
- **Code integration**: See `examples/smagorinsky_integration_example.jl`

Happy simulating! 🚀
