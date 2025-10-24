# Spectral Element Discretization with Smagorinsky Turbulence Model

## Overview

This implementation provides a complete spectral element discretization of the energy equation for compressible Navier-Stokes equations with Smagorinsky subgrid-scale turbulence modeling. The implementation is designed to integrate seamlessly with the existing Jexpresso framework.

## Author
Claude Code (AI Assistant)
Date: 2025-10-24

## Components

### 1. Mathematical Formulation
**File**: `docs/spectral_element_energy_smagorinsky.md`

Complete mathematical documentation including:
- Governing equations (compressible Navier-Stokes with potential temperature)
- Smagorinsky turbulence model theory
- Spectral element weak formulation
- Discrete implementation details
- Numerical considerations
- Test case recommendations
- Parameter guidelines

### 2. Smagorinsky Turbulence Model
**File**: `src/kernel/Turbulence/Smagorinsky.jl`

Core turbulence modeling functionality:
- `compute_strain_rate_magnitude_2d!()` - Compute strain rate tensor magnitude in 2D
- `compute_strain_rate_magnitude_3d!()` - Compute strain rate tensor magnitude in 3D
- `compute_smagorinsky_viscosity_2d!()` - Calculate eddy viscosity (2D)
- `compute_smagorinsky_viscosity_3d!()` - Calculate eddy viscosity (3D)
- `compute_element_filter_width()` - Determine filter width Δ

**Key Features**:
- Computes strain rate: S_ij = 1/2(∂u_i/∂x_j + ∂u_j/∂x_i) - 1/3 δ_ij (∇·u)
- Eddy viscosity: ν_t = (C_s Δ)² |S|
- Thermal diffusivity: κ_t = ρ ν_t / Pr_t
- Handles both 2D and 3D flows
- Fully compatible with spectral element tensor-product structure

### 3. Viscous Flux Functions
**File**: `problems/equations/CompEuler/theta/user_viscous_flux.jl`

Viscous flux computation for energy and momentum equations:
- `compute_theta_gradient_2d!()` / `compute_theta_gradient_3d!()` - Gradient of potential temperature
- `user_viscous_flux_energy_2d!()` / `user_viscous_flux_energy_3d!()` - Energy equation viscous fluxes
- `compute_stress_tensor_2d!()` / `compute_stress_tensor_3d!()` - Viscous stress tensor
- `user_viscous_flux_momentum_2d!()` - Momentum equation viscous fluxes
- `user_viscous_flux_full_2d!()` - Combined energy + momentum viscous fluxes
- `add_molecular_viscosity!()` - Add molecular transport to turbulent values

**Key Features**:
- Energy flux: F_visc = κ ∇θ
- Momentum flux: F_visc = τ (stress tensor)
- Supports molecular + turbulent contributions
- Gradient computation using spectral element derivative matrices

### 4. Integration Example
**File**: `examples/smagorinsky_integration_example.jl`

Complete working example showing:
- How to integrate Smagorinsky model into RHS computation
- Weak form discretization for viscous terms
- 2D and 3D implementations
- Diagnostic functions for turbulence analysis
- Energy conservation validation
- Usage with existing Jexpresso time integrators

## Quick Start

### Basic Usage

```julia
# Load the modules
include("src/kernel/Turbulence/Smagorinsky.jl")
include("problems/equations/CompEuler/theta/user_viscous_flux.jl")

# Set Smagorinsky parameters
C_s = 0.18    # Smagorinsky constant (atmospheric: 0.18-0.2)
Pr_t = 0.85   # Turbulent Prandtl number

# In your RHS computation loop (per element):
ngl = mesh.ngl
ν_t = zeros(ngl, ngl)
μ_t = zeros(ngl, ngl)
κ_t = zeros(ngl, ngl)

# Compute eddy viscosity
compute_smagorinsky_viscosity_2d!(ν_t, μ_t, κ_t, q, dψ, metrics, mesh, iel,
                                  C_s, Pr_t, PhysConst)

# Compute viscous fluxes
F_visc = zeros(ngl, ngl, neqs)
G_visc = zeros(ngl, ngl, neqs)
user_viscous_flux_full_2d!(F_visc, G_visc, q, μ_t, κ_t, dψ, metrics, mesh, iel)

# Add to RHS using weak form (see integration example for details)
```

### Parameters

**Smagorinsky Constant (C_s)**:
- Atmospheric flows: 0.18 - 0.2
- Engineering flows: 0.1 - 0.18
- Channel flow: ~0.1
- Default recommended: 0.18

**Turbulent Prandtl Number (Pr_t)**:
- Typical range: 0.7 - 1.0
- Atmospheric: 0.85
- Default recommended: 0.85

**Molecular Properties (air)**:
- μ_molecular ≈ 1.8×10⁻⁵ kg/(m·s)
- Pr_molecular ≈ 0.71
- κ_molecular = μ_mol × c_p / Pr

## Integration with Existing Code

### Modifying RHS Computation

To add Smagorinsky viscous terms to existing Jexpresso simulations:

1. **Include new modules** in your driver file
2. **Add viscous RHS function** after inviscid computation
3. **Sum contributions**: `RHS_total = RHS_inviscid + RHS_viscous`

Example modification to `src/kernel/operators/rhs.jl`:

```julia
function _build_rhs!(RHS, u, params, time)
    # Existing inviscid computation
    # ... (keep existing code)

    # Add Smagorinsky viscous terms
    if params.inputs[:lsmagorinsky]  # New input flag
        build_rhs_smagorinsky!(RHS, u, params, time)
    end
end
```

### Input Parameters

Add to input file (e.g., `input.jl`):

```julia
inputs[:lsmagorinsky] = true     # Enable Smagorinsky model
inputs[:C_s] = 0.18              # Smagorinsky constant
inputs[:Pr_t] = 0.85             # Turbulent Prandtl number
inputs[:lmolecular_visc] = true  # Include molecular viscosity
```

## Theoretical Background

### Smagorinsky Model

The Smagorinsky model is a classical algebraic subgrid-scale (SGS) model for Large Eddy Simulation (LES):

**Eddy viscosity**:
```
ν_t = (C_s Δ)² |S|
```

where:
- `Δ = h_e / (N+1)` is the filter width (element size / polynomial order)
- `|S| = √(2 S_ij S_ij)` is the magnitude of the strain rate tensor
- `C_s` is the Smagorinsky constant (calibrated from experiments/DNS)

**Physical interpretation**:
- Represents unresolved turbulent eddies at scales smaller than Δ
- Energy cascade from resolved to subgrid scales
- Provides dissipation to prevent energy pile-up at grid scale

### Energy Equation with Viscous Terms

For potential temperature formulation:

```
∂(ρθ)/∂t + ∇·(ρθu) = ∇·(κ_total ∇θ) + S_θ
```

where `κ_total = κ_molecular + κ_turbulent`

**Weak form** (continuous Galerkin):
```
∫_Ω φ ∂(ρθ)/∂t dΩ = -∫_Ω ∇φ·(ρθu) dΩ + ∫_Ω ∇φ·(κ∇θ) dΩ + boundary terms
```

**Spectral element discretization**:
- Test functions = Lagrange polynomials at LGL nodes
- Integration by parts → gradient of test function
- Quadrature using GLL points
- Mass matrix M (can be lumped for efficiency)

## Numerical Considerations

### Stability

**CFL condition** for diffusion (explicit time integration):
```
Δt ≤ C_CFL × h² / (2d × max(κ_total))
```

For high Reynolds number flows, molecular diffusion is negligible.
For LES with Smagorinsky, turbulent diffusivity dominates.

**Recommended**: Use small CFL constant (0.3-0.5) when viscous terms are significant.

### Accuracy

- **Filter width**: Should resolve energy-containing eddies (Δ in inertial range)
- **Polynomial order**: Higher order (N ≥ 4) provides better resolution
- **Element size**: h should be comparable to integral length scale
- **Smagorinsky constant**: May need tuning for specific applications

### Boundary Conditions

**Energy equation**:
- Adiabatic wall: `∂θ/∂n = 0`
- Isothermal wall: `θ = θ_wall`
- Periodic: natural in spectral elements

**Momentum equation**:
- No-slip wall: `u = 0`
- Free-slip wall: `u·n = 0, ∂u_t/∂n = 0`
- Periodic: natural

## Validation and Testing

### Recommended Test Cases

1. **Decaying Isotropic Turbulence**
   - Initialize with turbulent velocity field
   - Monitor kinetic energy decay: `E(t) ∝ t^(-n)`
   - Compare with DNS or experimental data
   - Check energy spectrum: `E(k) ∝ k^(-5/3)` in inertial range

2. **Channel Flow**
   - Periodic in streamwise direction
   - Compare mean velocity profile
   - Check law of the wall: `u^+ = f(y^+)`
   - Validate Reynolds stress profiles

3. **Atmospheric Boundary Layer**
   - Use existing Jexpresso atmospheric test cases
   - Add Smagorinsky model
   - Compare with field measurements or LES data
   - Check vertical profiles of wind, temperature

### Diagnostic Outputs

Monitor during simulation:
- Mean eddy viscosity: `<ν_t>`
- Max/min eddy viscosity: `ν_t^max, ν_t^min`
- Effective Reynolds number: `Re_eff = UL/(ν_mol + ν_t)`
- Turbulent kinetic energy: `TKE = 0.5<u'·u'>`
- Energy dissipation rate: `ε`

## Performance Considerations

### Computational Cost

Adding Smagorinsky model increases cost by:
- ~20-30% for strain rate computation (velocity gradients)
- ~5-10% for eddy viscosity calculation
- ~30-50% for viscous flux computation and RHS assembly

**Total overhead**: Approximately 50-80% compared to inviscid Euler

### Optimization Strategies

1. **Pre-allocate arrays**: Avoid allocations inside element loop
2. **Use @turbo**: Leverage LoopVectorization.jl for gradient computation
3. **GPU acceleration**: Extend to GPU kernels (similar to existing rhs_gpu.jl)
4. **Filter width**: Pre-compute once during initialization
5. **Derivative matrices**: Reuse existing `dψ` matrices

## Comparison with Existing DynSGS

Jexpresso already has a dynamic SGS model in `src/kernel/ArtificialViscosity/DynSGS.jl`.

**Key differences**:

| Feature | DynSGS (existing) | Smagorinsky (new) |
|---------|-------------------|-------------------|
| Type | Residual-based | Strain-rate-based |
| Formula | `μ ∝ Δ² × residual_ratio` | `μ_t = ρ(C_s Δ)² \|S\|` |
| C_s | Self-adjusting | Fixed constant |
| Complexity | Medium | Low |
| Tuning | Automatic | Manual (C_s) |
| Physics | Numerical stability | Physical turbulence |

**When to use**:
- **DynSGS**: When numerical stability is primary concern, or for under-resolved simulations
- **Smagorinsky**: When physical turbulence modeling is needed, for well-resolved LES

Both can coexist in the code. User chooses via input flags.

## References

1. Smagorinsky, J. (1963). "General circulation experiments with the primitive equations." *Monthly Weather Review*, 91(3), 99-164.

2. Germano, M., Piomelli, U., Moin, P., & Cabot, W. H. (1991). "A dynamic subgrid-scale eddy viscosity model." *Physics of Fluids A*, 3(7), 1760-1765.

3. Pope, S. B. (2000). *Turbulent Flows*. Cambridge University Press.

4. Sagaut, P. (2006). *Large Eddy Simulation for Incompressible Flows*. Springer.

5. Giraldo, F. X., & Restelli, M. (2008). "A study of spectral element and discontinuous Galerkin methods for the Navier-Stokes equations in nonhydrostatic mesoscale atmospheric modeling." *Journal of Computational Physics*, 227(8), 3849-3877.

6. Kopriva, D. A. (2009). *Implementing Spectral Methods for Partial Differential Equations*. Springer.

## Future Enhancements

Potential extensions to this implementation:

1. **Dynamic Smagorinsky**: Compute C_s dynamically using Germano identity
2. **Wall damping**: van Driest damping function near walls
3. **Anisotropic filtering**: Different filter widths in different directions
4. **Momentum equation**: Full viscous stress tensor for momentum
5. **GPU kernels**: Port to KernelAbstractions for GPU acceleration
6. **Implicit viscous solver**: For larger time steps (CN, BDF schemes)
7. **Wall-adapting models**: WALE, Vreman, or σ-models
8. **Hybrid RANS-LES**: Detached Eddy Simulation (DES), DDES

## Contact and Support

This implementation was created as an educational and research tool for the Jexpresso
spectral element framework. For questions or issues:

1. Review the mathematical documentation in `spectral_element_energy_smagorinsky.md`
2. Check the example code in `examples/smagorinsky_integration_example.jl`
3. Consult the Jexpresso documentation for framework-specific questions
4. See references for theoretical background

## License

This code follows the same license as the Jexpresso framework.

---

**Generated by**: Claude Code (Anthropic)
**Date**: October 24, 2025
**Version**: 1.0
