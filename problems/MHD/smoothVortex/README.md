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

At the final time each run measures the **absolute** $L^1$, $L^2$ and
$L^\infty$ error of the **velocity** against the exact solution,

$$
\int_\Omega |\mathbf u_h - \mathbf u|\,d\Omega, \qquad
\Big(\int_\Omega |\mathbf u_h - \mathbf u|^2 d\Omega\Big)^{1/2}, \qquad
\max_\Omega |\mathbf u_h - \mathbf u| ,
$$

as the paper's Fig. 1 plots them, using the mesh's own nodal quadrature
weights (the lumped mass: $\omega_i\omega_j|J|$ summed over the elements at
each node; they integrate the domain area to 14 digits), and stores it in
`errors/nop<N>_nelx<M>_<visc>.dat` — with the relative norms in the header
too, since they cost nothing and say how large the error is against the
solution it is measured on. The figures are then drawn from **every** error
stored there that carries the same final time (and the `norm=abs` marker, so a
store written before the norms became absolute is never mixed in):

| file | contents |
|---|---|
| `convergence_dsgs-it<n>.png` | the paper's Fig. 1 layout: the three norms side by side against $1/\sqrt{\#\mathrm{DOFs}}$ on log-log axes, one line per polynomial order, with slope guides, and the measured rate $p$ of each order in its legend |
| `convergence_dsgs_L1-it<n>.png`, `_L2`, `_Linf` | the same panels one per file, at publication size |
| `convergence_galerkin-it<n>.png` and its three panels | the same for the plain Galerkin run (`JEXPRESSO_SV_VISC=none`), the second panel of their figure |
| `<var>-it<n>.png` | the usual field panels |

**$p$ in the legend** is the measured convergence rate: the slope of that
order's last two points on the log-log axes,

$$
p = \frac{\log(e_{i-1}/e_i)}{\log(h_{i-1}/h_i)},\qquad h = 1/\sqrt{\#\mathrm{DOFs}},
$$

so it is the order in the mesh size $h$ between the two finest resolutions that
order has run — the number the dashed slope guides are there to be compared
against. It needs at least two resolutions per order; with one, the legend
carries the order alone. The guides are the nominal rates of the lowest and
highest orders on the figure, $\min(N)+1$ and $\max(N)+1$.

So a sweep over orders and meshes builds the whole figure and each run
replaces only its own point. **The scan clears the store before it starts**
(`SV_KEEP=1` to accumulate instead, to finish a partial sweep): since the
figure is redrawn from the whole store at every run, a leftover sweep would
otherwise appear on the comparison of a new one — with orders that this
comparison has not computed yet. `rm -r problems/MHD/smoothVortex/errors`
does the same by hand; the store is git-ignored. Only errors from the same
final time are drawn together.

The hook is `user_plot_2d` in `user_plot.jl`, the 2D counterpart of the
`user_plot_1d` hook the 1D cases use, added to
`src/io/plotting/jeplots.jl`. Unlike the 1D one it is additive: the generic
field panels are still rendered.

## Running it

```bash
tools/smooth_vortex_scan.sh                          # orders 4-7, meshes 4-32,
                                                     # DynSGS and Galerkin
SV_NOPS="3 4" SV_NELX="8 16 32" tools/smooth_vortex_scan.sh
SV_VISC=dsgs tools/smooth_vortex_scan.sh             # only the DynSGS panel
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
| `JEXPRESSO_SV_CR`, `JEXPRESSO_SV_CMAX` | 1, 0.5 | `:dsgs_CR`, `:dsgs_Cmax`; both zero keeps the sensor and applies no viscosity |
| `JEXPRESSO_SV_REL` | 1 | `:dsgs_rel`, the normalization floor; `1e-3` is what the kernels used before the fix |
| `JEXPRESSO_SV_SENSOR` | `residual` | `legacy` for the assembled-rate sensor |

The time step follows the resolution so that the Courant number is the same
at every point of a sweep (about 0.05 against the fastest wave of this
initial condition) and the comparison is not contaminated by a changing time
error.

The meshes come from one parameterized `vortex_periodic.geo`:

```bash
gmsh -2 -setnumber nx 16 vortex_periodic.geo -o vortex_16x16.msh
```

which `tools/smooth_vortex_mesh.sh` wraps.

## What this test found

It is the test that caught two defects in the DynSGS sensor itself, both fixed
in `src/kernel/physics/SGS.jl` and `src/kernel/operators/rhs.jl` (DSGS.md §4.2,
§4.4 and the defect list in §6). Before, DynSGS held $\nu$ on its cap
$C_{max}\Delta\lambda$ — an $O(h)$ viscosity, the same for every polynomial
order — and the case converged at **rate 1** while the same scheme without it
converged at its design order. Velocity $L^1$ at $t = 1$, `:nop => 4`:

| elements per side | 4 | 8 | 16 | 32 | rate |
|---|---|---|---|---|---|
| DynSGS, before | 1.004e-2 | 5.240e-3 | 2.678e-3 | 1.142e-3 | 1.0 |
| **DynSGS, after** | 4.26e-3 | 5.220e-5 | 2.208e-6 | 1.639e-7 | 4.6, 3.8 |
| plain Galerkin | 1.570e-3 | 3.541e-5 | 2.039e-6 | 1.631e-7 | 4.1, 3.6 |

The residual viscosity now costs 0.5 % of the error at the finest mesh instead
of 7000×, which is the statement this problem exists to check. (The Galerkin
rate falls from 5.5 toward 4 because $\Delta t \propto 1/(n_{elx}N)$ here, so
the fourth-order time error takes over; fix $\Delta t$ with `JEXPRESSO_SV_DT`
to see the spatial rate alone.)

The runs that isolated the cause are worth repeating after any change to the
sensor:

```bash
# the sensor with no viscosity applied: the solution is the Galerkin one and
# the printed residual is the sensor's reading of a clean solution
JEXPRESSO_SV_CR=0 JEXPRESSO_SV_CMAX=0 JEXPRESSO_SV_CMIN=0 \
JEXPRESSO_DSGS_DEBUG=1 JEXPRESSO_SV_NELX=16 \
    julia --project=. src/Jexpresso.jl MHD smoothVortex
```

`JEXPRESSO_DSGS_DEBUG=1` prints, every 200 kernel calls and for the first
twelve, $\nu_{max}$, the cap, and each equation's normalized residual and
denominator, with the node carrying the largest one.

## The time error, and why a fourth-order integrator caps the rate at 4

`CarpenterKennedy2N54` is **fourth order**, and the deck's default time step
follows the resolution, $\Delta t \propto h$, so the measured error is

$$
C_s\,h^{N+1} \;+\; C_t\,\Delta t^4 \;\propto\; h^{N+1} + h^4 ,
$$

and no order above 3 can show its own rate — $p$ saturates at 4 however fine
the mesh, on DynSGS and Galerkin alike. Measured on the plain Galerkin run at
`:nop => 4` (4/8/16/32 elements per side): 5.47, then 4.12, then 3.64, as the
second term takes over. That is the integrator, not the discretization.

Two ways out, both switchable:

| | |
|---|---|
| `JEXPRESSO_SV_DT=<fixed>` | one $\Delta t$ for the whole sweep, so the time error is a constant rather than something that shrinks at fourth order and pollutes the slope. `tools/smooth_vortex_mpi_scan.sh` sets it automatically from the finest (nelx, nop) of the sweep |
| `JEXPRESSO_SV_SOLVER=vern9` | `Vern9` (9th order), `dp8` (8th), `vern7`, `ssprk54`, `tsit5`, or `ck54` for the default |

**Measured, the fourth-order integrator is not what limits these runs.** At
32×32 elements and `:nop => 4`, absolute velocity $L^1$ at $t = 1$:

| | |
|---|---|
| `ck54`, $\Delta t = 10^{-3}$ | 2.318e-5 |
| `vern9`, $\Delta t = 10^{-3}$ | 2.320e-5 |
| `vern9`, $\Delta t = 3.3\times10^{-4}$ | 2.314e-5 |

Neither the integrator nor the step moves the error: it is purely spatial
there, and Vern9's 16 stages per step cost 3× for nothing. So
`tools/smooth_vortex_mpi_scan.sh` runs `ck54` with one fixed $\Delta t$ for
the sweep, and the honest check on any sweep is to re-run its **finest** case
at half `SV_DT` and see that the error does not move.

## The final time is part of the measurement

`SV_TEND` shortens a run, which is the right way to check that a sweep's
machinery works before committing hours to it — and the wrong thing to leave
set. At $t = 0.05$ the vortex has moved a twentieth of the box, the error is
dominated by the initial transient, and it barely converges with resolution:
measured at $\Delta t = 3.3\times10^{-4}$, `:nop => 4` on 32² gives 7.539e-6
and `:nop => 6` on 32² gives 1.433e-6, against 2.3e-5 and far less at
$t = 1$ — flat lines and negative rates that look exactly like a broken
solver. The figures therefore carry the final time and the step in their
titles.

## Comparing P1 and P3 directly, as in the paper

Dao & Nazarov's Fig. 1 is a low-order comparison. Any subset of the orders in
the store can be drawn on its own axes, either while running —

```bash
JEXPRESSO_SV_PLOT_NOPS="1 3" tools/smooth_vortex_mpi_scan.sh
```

which writes `convergence_<visc>_nop1-3-it<n>.png` beside the all-orders
figure — or afterwards, from the stored errors alone, without running
anything:

```bash
julia --project=. tools/smooth_vortex_plot.jl 1,3
julia --project=. tools/smooth_vortex_plot.jl          # every order
julia --project=. tools/smooth_vortex_plot.jl 1,3 --t=1.0 --out=figs
```

`tools/smooth_vortex_plot.jl` reads `errors/*.dat` and calls the case's own
plotting code, so its figures are the same ones a run produces.

## Running it on many cores

```bash
tools/smooth_vortex_mpi_scan.sh                          # 4 ranks per case
SV_NP=8  SV_NELX="8 16 32 64" tools/smooth_vortex_mpi_scan.sh
SV_NP=16 SV_JOBS=4 SV_NOPS="1 3" tools/smooth_vortex_mpi_scan.sh
SV_NP=1  SV_JOBS=8 tools/smooth_vortex_mpi_scan.sh       # 8 serial cases at once
```

On a SLURM cluster, `tools/smooth_vortex_slurm.sh` is the same sweep as a
batch job — as written, orders 4 and 6 on 32² and 64² elements, both panels,
16 ranks per case and 4 cases at a time:

```bash
sbatch tools/smooth_vortex_slurm.sh
```

`SV_NP × SV_JOBS` must equal `--ntasks` in its header, or the job steps queue
behind each other (too few tasks) or leave cores idle (too many). Each case is
one `srun --exclusive -n $SV_NP` step — `--exclusive` at *step* level is what
keeps four concurrent steps off each other's cores — and the script does one
serial warm-up case first so that four cases do not compile at once and
contend for the depot's precompile locks. Jexpresso's MPI.jl must be built
against the MPI the modules provide; `tools/check_mpi_setup.sh` checks that.

`SV_NP` is ranks per case — what a big grid needs, since one run must fit and
finish — and `SV_JOBS` is how many (independent) cases run at the same time;
their product is what you are asking the machine for. The norms are
MPI-correct: each unknown is weighed once, by the rank that owns it
(`mesh.gip2owner`, the map the DSS assembler uses), and the degree-of-freedom
count on the abscissa is the exact $(n_{elx}N)^2$ of this periodic box, so a
point does not move when the rank count changes.

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
