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
weights, and stores it in `errors/nop<N>_nelx<M>_b<beta>_<visc>.dat`; the figures are
drawn from every error in that store. β is part of both the record name and
the filter that builds the figure, so the classical β = 5 sweep and a β = 1
sweep (the one amplitude-matched to the MHD vortex) neither overwrite each
other nor share a curve.

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
| `JEXPRESSO_EV_BDF` | 2 | `:dsgs_bdf_order`, the order of the residual sensor's backward difference; records `dsgs_bdf<k>` ([below](#the-order-of-the-sensors-time-difference)) |

The meshes are the MHD case's — the same $[-5,5]^2$ doubly periodic quad grids,
and two copies would only mean two things to regenerate. Nothing else is
shared: the per-case SEM cache and the error store live next to this deck.

## The order of the sensor's time difference

`JEXPRESSO_EV_BDF` (`:dsgs_bdf_order`) sets the order $k$ of the backward
difference the residual sensor takes $\partial_t q$ with (default 2, BDF2); the
option, and the study `tools/smooth_vortex_bdf_study.sh` runs on both vortices,
are described in the [MHD case's README](../../MHD/smoothVortex/README.md#the-order-of-the-sensors-time-difference-bdf-k).
On this vortex ($\beta = 5$, without the cutoff, which is this deck's default)
DynSGS is far above Galerkin from $N = 6$ on, and the order of the difference
does not change that ($32^2$, $t = 1$, $\Delta t = 1/1750$, velocity $L^1$, ratio
to Galerkin in brackets):

| $N$ | Galerkin | DynSGS, BDF2 | DynSGS, BDF5 |
|---|---|---|---|
| 4 | 3.417e-3 | 3.18e-2 (9.3) | 2.20e-2 (6.4) |
| 5 | 2.297e-3 | 1.36e-2 (5.9) | 4.52e-2 (19.7) |
| 6 | 7.223e-5 | 6.18e-2 (856) | 5.57e-2 (771) |
| 7 | 2.927e-5 | 4.22e-2 (1443) | 4.86e-2 (1661) |

$\nu$ is not sitting on the time-difference floor: at $N = 6$ it starts at
2.2e-6 right after the startup hold, the same with BDF2 and BDF5 to three digits,
and grows by about 1.7 % per step to ~1e-2 at $t \approx 0.45$ before it
saturates, with either. That growth, not the order of the difference, is what
the error above measures.

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
| `user_plot.jl` | the error norms, the record tags and the store; the convergence figure is the MHD case's [`convergence_figure.jl`](../../MHD/smoothVortex/convergence_figure.jl), so the two cases are drawn the same way |
| `user_flux.jl`, `user_primitives.jl`, `user_source.jl` | the 2D total-energy Euler system, as in `CompEuler/ffs_step` |
| `user_bc.jl` | free-slip fallback; periodicity is handled at the mesh level |
