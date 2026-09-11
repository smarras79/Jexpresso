# Smooth MHD vortex — the accuracy test of Dao & Nazarov (2022), §5.1.1

The one problem of that paper with a closed-form solution, and the one its
Fig. 1 convergence history is built on: a steady isentropic MHD vortex
carried by a uniform flow across a doubly periodic box, so the exact
solution at any time is the initial condition translated by $\mathbf v_0 t$.
Its purpose is to check that **the residual viscosity does not destroy the
accuracy of a high-order solution** where the solution is smooth — the
opposite end of the model's range from the shock tube of
[`brioWu1d`](../brioWu1d/README.md).

```bash
tools/smooth_vortex_mesh.sh     # the meshes, once
tools/smooth_vortex_scan.sh     # the sweep that builds the figure
```

## The problem

Domain $[-5, 5]^2$, doubly periodic, $\gamma = 5/3$, $\mathbf v_0 = (1,1)$,
so the vortex returns to its starting point at $t = 10$. With
$r^2 = x^2 + y^2$ and $f = e^{(1-r^2)/2}$ (Balsara, *ApJS* 151:149 (2004)):

$$
\rho = 1,\qquad
\mathbf v = \mathbf v_0 + \frac{\kappa}{2\pi} f\,(-y,\ x),\qquad
\mathbf B = \frac{\mu}{2\pi} f\,(-y,\ x),
$$

$$
p = 1 + \frac{1}{8\pi^2}\big(\mu^2(1 - r^2) - \kappa^2\big)e^{1-r^2},
$$

with $\kappa = \mu = 1$; $w$, $B_z$ and $\psi$ are zero and stay zero. That
pressure is the one that makes the vortex an **exact steady solution** in the
co-moving frame, which is what an accuracy test needs. With
$v_\theta = \kappa g r$, $B_\theta = \mu g r$ and $g = f/2\pi$, radial
momentum balance reads

$$
\frac{\rho v_\theta^2}{r} - \frac{B_\theta^2}{r}
   = \frac{d}{dr}\Big(p + \frac{B_\theta^2}{2}\Big),
$$

and the $p$ above closes it identically, both sides being
$(1/4\pi^2)\,e^{1-r^2} r\,[\kappa^2 - \mu^2(2-r^2)]$ — the derivation is in
the header of `initialize.jl`.

## What the run produces

At the final time each run measures the relative $L^1$, $L^2$ and $L^\infty$
error of the **velocity** against the exact solution, using the mesh's own
nodal quadrature weights (the lumped mass: $\omega_i\omega_j|J|$ summed over
the elements at each node; they integrate the domain area to 14 digits), and
stores it in `errors/nop<N>_nelx<M>_<visc>.dat`. The figure is then drawn
from **every** error stored there:

| file | contents |
|---|---|
| `convergence_dsgs-it<n>.png` | the paper's Fig. 1 layout: error against $1/\sqrt{\#\mathrm{DOFs}}$ on log-log axes, one line per polynomial order, with slope guides, and the measured rate of each order in its legend |
| `convergence_galerkin-it<n>.png` | the same for the plain Galerkin run (`JEXPRESSO_SV_VISC=none`), the second panel of their figure |
| `<var>-it<n>.png` | the usual field panels |

So a sweep over orders and meshes builds the whole figure and each run
replaces only its own point. `rm -r problems/MHD/smoothVortex/errors` starts
a fresh comparison; the store is git-ignored. Only errors from the same final
time are drawn together.

The hook is `user_plot_2d` in `user_plot.jl`, the 2D counterpart of the
`user_plot_1d` hook the 1D cases use, added to
`src/io/plotting/jeplots.jl`. Unlike the 1D one it is additive: the generic
field panels are still rendered.

## Running it

```bash
tools/smooth_vortex_scan.sh                          # orders 2-4, meshes 4-32
SV_NOPS="3 4" SV_NELX="8 16 32" tools/smooth_vortex_scan.sh
SV_VISC="dsgs none" tools/smooth_vortex_scan.sh      # both panels of Fig. 1
```

or one run at a time:

```bash
JEXPRESSO_SV_NOP=3 JEXPRESSO_SV_NELX=16 \
    julia --project=. src/Jexpresso.jl MHD smoothVortex
```

| variable | default | meaning |
|---|---|---|
| `JEXPRESSO_SV_NOP` | 4 | polynomial order |
| `JEXPRESSO_SV_NELX` | 16 | elements per side; picks `vortex_<n>x<n>.msh` |
| `JEXPRESSO_SV_DT` | — | time step; otherwise $\Delta t \propto 1/(n_{elx}\,N)$ |
| `JEXPRESSO_SV_TEND` | 1.0 | final time |
| `JEXPRESSO_SV_CMIN` | 0 | the DynSGS background floor `:dsgs_Cmin` |
| `JEXPRESSO_SV_VISC` | `dsgs` | `none` for the plain Galerkin run |

The time step follows the resolution so that the Courant number is the same
at every point of a sweep (about 0.05 against the fastest wave of this
initial condition) and the comparison is not contaminated by a changing time
error.

The meshes come from one parameterized `vortex_periodic.geo`:

```bash
gmsh -2 -setnumber nx 16 vortex_periodic.geo -o vortex_16x16.msh
```

which `tools/smooth_vortex_mesh.sh` wraps.

## Why the background floor is off here

`:dsgs_Cmin` defaults to **0** in this case, unlike the shock cases. The floor
is a viscosity $C_{min}\Delta\lambda$ that is applied everywhere regardless of
the residual; on a smooth solution it is an $O(h)$ error that no polynomial
order can beat, so a sweep with the floor on measures the floor and not the
scheme — every order would collapse onto slope 1. Turn it on with
`JEXPRESSO_SV_CMIN` if that is what you want to see.

## Files

| file | role |
|---|---|
| `user_inputs.jl` | the deck and the environment overrides |
| `initialize.jl` | the vortex, its exact-solution helpers `sv_state`/`sv_wrap`, and the equilibrium derivation |
| `user_plot.jl` | the error norms, the store and the convergence figure |
| `user_flux.jl`, `user_source.jl`, `user_bc.jl`, `user_primitives.jl` | the GLM-MHD equations, identical to `orszagTangBormanis2024` |
| `vortex_periodic.geo` | the parameterized periodic mesh |
