# Spectral Element Discretization of the Energy Equation for Compressible Navier-Stokes with Smagorinsky Turbulence Model

## 1. Governing Equations

### 1.1 Compressible Navier-Stokes Energy Equation

The total energy equation for compressible flow with viscous effects:

```
∂(ρE)/∂t + ∇·(ρEu) = -∇·(pu) + ∇·(τ·u) - ∇·q + S
```

where:
- `ρ` = density
- `E` = total specific energy = `e + (u²+v²+w²)/2`
- `e` = internal energy
- `u = (u,v,w)` = velocity vector
- `p` = pressure
- `τ` = stress tensor
- `q` = heat flux vector
- `S` = source terms

### 1.2 Potential Temperature Formulation (Used in Jexpresso)

For atmospheric applications, Jexpresso uses potential temperature `θ`:

```
∂(ρθ)/∂t + ∇·(ρθu) = ∇·(κ∇θ) + S_θ
```

where:
- `θ` = potential temperature
- `κ` = thermal diffusivity (molecular + turbulent)

## 2. Smagorinsky Turbulence Model

### 2.1 Eddy Viscosity Formulation

The Smagorinsky model computes turbulent eddy viscosity:

```
ν_t = (C_s Δ)² |S|
```

where:
- `C_s` = Smagorinsky constant (typically 0.1-0.2)
- `Δ` = filter width (element size in spectral elements)
- `|S|` = magnitude of strain rate tensor

### 2.2 Strain Rate Tensor

The strain rate tensor for compressible flow:

```
S_ij = 1/2 (∂u_i/∂x_j + ∂u_j/∂x_i) - 1/3 δ_ij (∇·u)
```

Magnitude:
```
|S| = √(2 S_ij S_ij)
```

For 2D:
```
|S| = √[2((∂u/∂x)² + (∂v/∂y)² + 1/2(∂u/∂y + ∂v/∂x)²) - 2/3(∂u/∂x + ∂v/∂y)²]
```

For 3D:
```
|S| = √[2((∂u/∂x)² + (∂v/∂y)² + (∂w/∂z)² +
         1/2((∂u/∂y + ∂v/∂x)² + (∂u/∂z + ∂w/∂x)² + (∂v/∂z + ∂w/∂y)²)) -
         2/3(∂u/∂x + ∂v/∂y + ∂w/∂z)²]
```

### 2.3 Total Viscosity and Thermal Diffusivity

```
μ_total = μ_molecular + ρ ν_t
κ_total = κ_molecular + ρ ν_t / Pr_t
```

where `Pr_t` = turbulent Prandtl number (typically 0.7-1.0)

## 3. Spectral Element Weak Formulation

### 3.1 Test Function and Integration by Parts

Multiply by test function `φ` and integrate over element `Ω_e`:

```
∫_Ω_e φ ∂(ρθ)/∂t dΩ = -∫_Ω_e ∇φ·(ρθu) dΩ + ∫_∂Ω_e φ (ρθu·n) dΓ
                        +∫_Ω_e ∇φ·(κ_total ∇θ) dΩ - ∫_∂Ω_e φ (κ_total ∇θ·n) dΓ
                        +∫_Ω_e φ S_θ dΩ
```

### 3.2 Continuous Galerkin Formulation

For continuous Galerkin (CG) spectral elements used in Jexpresso:

**Inviscid terms** (already implemented):
```
R_inv = -M^(-1) [∫_Ω_e ∇φ·F dΩ - ∫_∂Ω_e φ F*·n dΓ]
```

**Viscous terms** (to be added):
```
R_visc = M^(-1) [∫_Ω_e ∇φ·F_visc dΩ - ∫_∂Ω_e φ F_visc*·n dΓ]
```

where:
- `M` = mass matrix
- `F = ρθu` = inviscid flux
- `F_visc = κ_total ∇θ` = viscous flux
- `F*` = numerical flux (Riemann solver)
- `F_visc*` = viscous numerical flux (central or BR1/BR2)

## 4. Spectral Element Discrete Formulation

### 4.1 Tensor-Product Basis

For element `e` in 2D with Lagrange polynomials `ℓ_i(ξ)`, `ℓ_j(η)`:

```
θ(x,y) ≈ θ^e(ξ,η) = ∑_{i=1}^{N+1} ∑_{j=1}^{N+1} θ_{ij}^e ℓ_i(ξ) ℓ_j(η)
```

Test functions are the same basis functions: `φ_{kl} = ℓ_k(ξ) ℓ_l(η)`

### 4.2 Gradient Computation in Reference Element

Using chain rule:
```
∂θ/∂x = ∂θ/∂ξ · ∂ξ/∂x + ∂θ/∂η · ∂η/∂x
∂θ/∂y = ∂θ/∂ξ · ∂ξ/∂y + ∂θ/∂η · ∂η/∂y
```

In matrix form:
```
[∂θ/∂ξ]   [∂ξ/∂x  ∂ξ/∂y] [∂θ/∂x]
[∂θ/∂η] = [∂η/∂x  ∂η/∂y] [∂θ/∂y]
```

Inverse:
```
[∂θ/∂x]            [∂η/∂y  -∂η/∂x] [∂θ/∂ξ]
[∂θ/∂y] = 1/J * [-∂ξ/∂y   ∂ξ/∂x] [∂θ/∂η]
```

where `J = det(∂(ξ,η)/∂(x,y))` is the Jacobian.

### 4.3 Derivative Matrices

Reference space derivatives using differentiation matrix `D`:
```
∂θ/∂ξ|_{ξ_k,η_l} = ∑_{i=1}^{N+1} D_{ki} θ_{il}
∂θ/∂η|_{ξ_k,η_l} = ∑_{j=1}^{N+1} D_{lj} θ_{kj}
```

where `D_{ij} = dℓ_i/dξ|_{ξ_j}` (already computed in Jexpresso)

### 4.4 Strain Rate Computation

**2D velocity gradients** (element-wise):
```
∂u/∂x = (D ⊗ I) u * ∂ξ/∂x + (I ⊗ D) u * ∂η/∂x
∂u/∂y = (D ⊗ I) u * ∂ξ/∂y + (I ⊗ D) u * ∂η/∂y
∂v/∂x = (D ⊗ I) v * ∂ξ/∂x + (I ⊗ D) v * ∂η/∂x
∂v/∂y = (D ⊗ I) v * ∂ξ/∂y + (I ⊗ D) v * ∂η/∂y
```

**Strain rate components**:
```
S_xx = ∂u/∂x - 1/3(∂u/∂x + ∂v/∂y)
S_yy = ∂v/∂y - 1/3(∂u/∂x + ∂v/∂y)
S_xy = 1/2(∂u/∂y + ∂v/∂x)
```

**Magnitude**:
```
|S| = √(2(S_xx² + S_yy² + 2S_xy²))
```

### 4.5 Eddy Viscosity at Quadrature Points

For each element `e` and quadrature point `(i,j)`:

```
Δ_e = h_e / (N+1)   or   Δ_e = (volume_e)^(1/d) / (N+1)
ν_t[i,j,e] = (C_s Δ_e)² |S[i,j,e]|
μ_t[i,j,e] = ρ[i,j,e] * ν_t[i,j,e]
κ_t[i,j,e] = μ_t[i,j,e] * c_p / Pr_t
```

### 4.6 Viscous Flux Computation

**Energy equation viscous flux** (2D):
```
F_visc_x = κ_total * ∂θ/∂x
F_visc_y = κ_total * ∂θ/∂y
```

**Full momentum equation viscous stress** (2D):
```
τ_xx = 2μ_total * S_xx
τ_yy = 2μ_total * S_yy
τ_xy = 2μ_total * S_xy
```

**Viscous flux for momentum**:
```
F_visc = [0,         0        ]
         [τ_xx,      τ_xy     ]
         [τ_xy,      τ_yy     ]
         [F_visc_x,  F_visc_y ]
```

### 4.7 Weak Form Discrete Equations

**Volume integral** (for each test function node `k,l` in element `e`):
```
R_kl^e = - ∑_{i,j} ω_i ω_j J_{ij}^e [∂ℓ_k/∂ξ|_{ξ_i} ℓ_l(η_j) F_ξ_{ij} +
                                       ℓ_k(ξ_i) ∂ℓ_l/∂η|_{η_j} F_η_{ij}]
```

where:
```
F_ξ = F_x * ∂ξ/∂x + F_y * ∂ξ/∂y
F_η = F_x * ∂η/∂x + F_y * ∂η/∂y
```

**Surface integral** (continuous Galerkin):
```
R_surf = ∑_{faces} ∫_face ℓ_k ℓ_l (F* - F) · n dΓ
```

### 4.8 Time Integration

Semi-discrete form:
```
M ∂q/∂t = R(q)
```

With mass lumping or explicit mass matrix inversion:
```
∂q/∂t = M^(-1) R(q)
```

Standard time integrators: RK2, RK3, RK4, LSRK (already in Jexpresso)

## 5. Implementation in Jexpresso

### 5.1 Required Modifications

1. **Compute velocity gradients** in each element
2. **Compute strain rate magnitude** `|S|`
3. **Compute Smagorinsky eddy viscosity** `ν_t = (C_s Δ)² |S|`
4. **Compute thermal eddy diffusivity** `κ_t = ρ ν_t / Pr_t`
5. **Compute temperature/theta gradients** for viscous flux
6. **Add viscous flux** to RHS computation

### 5.2 Code Structure

Following Jexpresso's architecture:

**New file**: `problems/equations/CompEuler/theta/user_viscous_flux.jl`
```julia
function compute_strain_rate_magnitude(u, v, dψ, metrics, mesh, iel)
    # Compute velocity gradients
    # Compute strain rate tensor
    # Return |S| at each quadrature point
end

function compute_smagorinsky_viscosity(ρ, u, v, dψ, metrics, mesh, iel, C_s)
    |S| = compute_strain_rate_magnitude(u, v, dψ, metrics, mesh, iel)
    Δ = element_size(mesh, iel) / mesh.ngl
    ν_t = (C_s * Δ)^2 * |S|
    μ_t = ρ .* ν_t
    return μ_t
end

function user_viscous_flux!(F_v, G_v, SD::NSD_2D, q, ∇q, μ_t, κ_t,
                            mesh::St_mesh, PhysConst)
    # Energy equation viscous flux
    ∂θ∂x = ∇q[1,4]  # gradient of θ in x
    ∂θ∂y = ∇q[2,4]  # gradient of θ in y

    F_v[4] = κ_t * ∂θ∂x
    G_v[4] = κ_t * ∂θ∂y

    # Can also add momentum equation viscous terms
    # F_v[2] = τ_xx, F_v[3] = τ_xy
    # G_v[2] = τ_xy, G_v[3] = τ_yy
end
```

**Modified file**: `src/kernel/operators/rhs.jl`
Add viscous RHS computation after inviscid terms:
```julia
function _build_rhs_viscous!(RHS, u, params, time)
    # 1. Compute gradients of velocity and theta
    # 2. Compute Smagorinsky eddy viscosity
    # 3. Compute viscous fluxes
    # 4. Add to RHS using weak form
end
```

### 5.3 Integration with Existing DynSGS

Note: Jexpresso already has a residual-based dynamic SGS model in `DynSGS.jl`.
The standard Smagorinsky is simpler and can coexist or replace it depending on application.

Key differences:
- **Dynamic SGS** (current): `μ_SGS ∝ Δ² * max(residual ratios)` with limiter
- **Smagorinsky** (new): `μ_t = (C_s Δ)² |S|`

Both can use same framework for applying viscosity to equations.

## 6. Numerical Considerations

### 6.1 Stability

CFL condition for diffusion:
```
Δt ≤ C * h² / (2d * max(κ_total))
```

For explicit time integration, viscous terms can be very restrictive.

### 6.2 Boundary Conditions

- **Adiabatic walls**: `∂θ/∂n = 0` → no viscous flux
- **Isothermal walls**: `θ = θ_wall` → Dirichlet BC
- **Periodic**: natural in spectral elements

### 6.3 Filter Width

For spectral elements:
```
Δ = h_e / (N+1)
```
where `h_e` is characteristic element size, `N` is polynomial order.

## 7. Test Cases

### 7.1 Decaying Turbulence

- Initialize with turbulent velocity field
- Monitor kinetic energy decay
- Compare with DNS if available

### 7.2 Channel Flow

- Periodic in streamwise direction
- Compare mean profiles with experimental data
- Check law of the wall: u⁺ = f(y⁺)

### 7.3 Atmospheric Boundary Layer

- Existing Jexpresso cases can be enhanced with Smagorinsky
- Compare with measurements or LES

## 8. References

1. Smagorinsky, J. (1963). "General circulation experiments with the primitive equations"
2. Germano, M. et al. (1991). "A dynamic subgrid-scale eddy viscosity model"
3. Giraldo, F.X. (2001). "Strong and weak Lagrange-Galerkin spectral element methods"
4. Kopriva, D.A. (2009). "Implementing Spectral Methods for PDEs"

## 9. Parameter Guidelines

- Smagorinsky constant: `C_s = 0.1 - 0.2` (atmospheric: 0.18-0.2)
- Turbulent Prandtl number: `Pr_t = 0.7 - 1.0` (often 0.85)
- Molecular Prandtl number: `Pr = 0.71` (air)
- Filter width ratio: `Δ / h = 1 / (N+1)`
