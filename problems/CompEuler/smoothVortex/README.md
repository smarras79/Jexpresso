# Isentropic (Shu) vortex — the hydrodynamic control for the MHD accuracy test

The classical smooth accuracy test of the 2D compressible Euler equations, set
up here as the **control** for [`MHD/smoothVortex`](../../MHD/smoothVortex/README.md):
same box, same meshes, same error norms, same figures, same sweep driver — and
no magnetic field, so no $\nabla\cdot\mathbf B$ to be numerically non-zero and
no GLM cleaning field $\psi$ to damp. Whatever error the MHD case carries
beyond this one is the price of the divergence constraint rather than of the
discretization.

```bash
tools/smooth_vortex_mesh.sh                                   # the meshes, once
SV_CASE=CompEuler/smoothVortex tools/smooth_vortex_mpi_scan.sh
julia --project=. tools/smooth_vortex_plot.jl --case=CompEuler/smoothVortex
```

## The problem

Domain $[-5,5]^2$, doubly periodic, nondimensional ($\rho_\infty = p_\infty =
T_\infty = 1$, $R = 1$), $\gamma = c_p/c_v$, $\mathbf v_0 = (1,1)$, so the
vortex returns to its starting point at $t = 10$ and the exact solution at any
time is the initial condition translated by $\mathbf v_0 t$. With
$r^2 = x^2 + y^2$:

$$
\delta u = -\frac{\beta}{2\pi}\,y\,e^{(1-r^2)/2},\qquad
\delta v =  \frac{\beta}{2\pi}\,x\,e^{(1-r^2)/2},\qquad
\delta T = -\frac{(\gamma-1)\beta^2}{8\gamma\pi^2}e^{1-r^2},
$$

$$
T = 1 + \delta T,\qquad \rho = T^{1/(\gamma-1)},\qquad p = \rho^\gamma .
$$

The state is isentropic ($p/\rho^\gamma \equiv 1$) and steady in the co-moving
frame for any $\beta$: the swirl's centrifugal force is balanced exactly by the
pressure gradient the temperature dip produces. **No source term of any kind**
— `:lsource => false`.

$\beta = 5$ is the classical strength (the density dips to 0.494 at the core).
**$\beta = 1$ is the setting for a like-for-like comparison with the MHD
vortex**, whose perturbation amplitude ($\kappa = \mu = 1$) is comparable;
`JEXPRESSO_EV_BETA=1`.

## What the run produces

Exactly what the MHD case produces, from the same code: at the final time each
run measures the **absolute** $L^1$, $L^2$ and $L^\infty$ error of the velocity
against the exact solution, with the solver's own lumped-mass quadrature
weights, and stores it in `errors/nop<N>_nelx<M>_<visc>.dat`; the figures are
drawn from every error in that store.

| file | contents |
|---|---|
| `convergence_dsgs-it<n>.png` | the paper's Fig. 1 layout: the three norms against $1/\sqrt{\#\mathrm{DOFs}}$ on log-log axes, one line per order, slope guides, the measured rate $p$ in each legend entry |
| `convergence_dsgs_L1-it<n>.png`, `_L2`, `_Linf` | the same panels one per file, at publication size |
| `convergence_galerkin-it<n>.png` and its panels | the plain Galerkin run (`JEXPRESSO_EV_VISC=none`) |

## Overrides

| variable | default | meaning |
|---|---|---|
| `JEXPRESSO_EV_NOP` | 4 | polynomial order |
| `JEXPRESSO_EV_NELX` | 16 | elements per side; picks `MHD/smoothVortex/vortex_<n>x<n>.msh` |
| `JEXPRESSO_EV_BETA` | 5 | vortex strength; 1 matches the MHD vortex's amplitude |
| `JEXPRESSO_EV_DT` | — | time step; otherwise $\Delta t \propto 1/(n_{elx}N)$ |
| `JEXPRESSO_EV_TEND` | 1.0 | final time |
| `JEXPRESSO_EV_VISC` | `dsgs` | `none` for the plain Galerkin run |
| `JEXPRESSO_EV_SOLVER` | `ck54` | `vern9`, `vern7`, `dp8`, `ssprk54`, `tsit5` |
| `JEXPRESSO_EV_CMIN` / `_CR` / `_CMAX` | 0, 1, 0.5 | the DynSGS coefficients |
| `JEXPRESSO_EV_REL` | 1 | `:dsgs_rel`, the normalization floor |
| `JEXPRESSO_EV_NORMS` | `domain` | `:dsgs_norms` |

The meshes are the MHD case's — the same $[-5,5]^2$ doubly periodic quad grids,
and two copies would only mean two things to regenerate. Nothing else is
shared: the per-case SEM cache and the error store live next to this deck.

## Why the background floor is off here

`:dsgs_Cmin` defaults to **0**, as in the MHD case. The floor is a viscosity
$C_{min}\Delta\lambda$ applied everywhere regardless of the residual; on a
smooth solution it is an $O(h)$ error that no polynomial order can beat, so a
sweep with it on measures the floor and not the scheme.

## Files

| file | role |
|---|---|
| `user_inputs.jl` | the deck and the environment overrides |
| `initialize.jl` | the vortex, its exact-solution helper `ev_state`, and the wrap |
| `user_plot.jl` | the error norms, the store and the convergence figures |
| `user_flux.jl`, `user_primitives.jl`, `user_source.jl` | the 2D total-energy Euler system, as in `CompEuler/ffs_step` |
| `user_bc.jl` | free-slip fallback; periodicity is handled at the mesh level |
