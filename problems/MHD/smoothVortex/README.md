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
tools/smooth_vortex_mpi_scan.sh                  # the sweep: DynSGS and Galerkin, every order on every mesh
julia --project=. tools/smooth_vortex_plot.jl    # the error table and the figures again, from the stored errors
```

On a laptop, the sweep of 11 September 2026 with the current setup (orders 4–7
on 4² to 32² elements, one core per case, four cases at a time):

```bash
SV_NP=1 SV_JOBS=4 SV_NOPS="4 5 6 7" SV_NELX="4 8 16 32" tools/smooth_vortex_mpi_scan.sh
```

## The problem

Domain $[-L/2, L/2]^2$, doubly periodic, with $L = 20$ by default (the
$[-10, 10]^2$ of the paper; see [The box](#the-box)), $\gamma = 5/3$,
$\mathbf v_0 = (1,1)$. With $r^2 = x^2 + y^2$ and $f = e^{(1-r^2)/2}$
(Balsara, *ApJS* 151:149 (2004)):

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

## The DynSGS setup

What `user_inputs.jl` runs, and so what a sweep measures:

| | setting | key |
|---|---|---|
| sensor | the element's own strong residual against the BDF2 history | `:dsgs_sensor => "residual"` |
| when $\nu$ is computed | once per time step, at the first stage | `:dsgs_freeze_stage`, the `DSGS_MHD` default with the residual sensor |
| **smoothness cutoff** | $c = 10^{-3}$ on the normalized residual | `:dsgs_cutoff`, **the default since 20 Oct 2026** |
| normalization | over the whole domain, each variable's spread floored at its own scale ($\bar\rho$, $\bar\rho\bar c$, $\bar\rho\bar c^2$, $\sqrt{\bar\rho}\,\bar c$) | `:dsgs_norms => "domain"`, `:dsgs_rel => 1` |
| coefficient | $\nu = \max\big(0, \min(C_{max}\Delta\lambda,\ C_R\Delta^2\max(0, r - c))\big)$ with $C_R = 1$, $C_{max} = 0.5$, $\Delta = h_e/(N+1)$ | `:dsgs_CR`, `:dsgs_Cmax` |
| background floor | none, $C_{min} = 0$ ([why](#why-the-background-floor-is-off-here)) | `:dsgs_Cmin` |
| form | physical: $\nu$ on $\rho$, $\bar\rho\nu$ in the deviatoric stress, $\kappa = \bar\rho\nu/\mathrm{Pr}$ on $\nabla T$, $\eta = \nu$ on $\mathbf B$, with the $\tau\cdot\mathbf u$ and $\eta\,\mathbf B\cdot\nabla\mathbf B$ work terms | `:dsgs_conserved => false`, `:dsgs_nazarov_energy => true` |
| startup | $\nu$ held at zero for 2 steps (the BDF2 history) | `:dsgs_hold_steps` |
| time | `CarpenterKennedy2N54` to $t = 1$ | |

**Why the cutoff.** On a smooth solution the residual sensor still reads the
truncation error of its own BDF2 time difference, $(\Delta t^2/3)\,\partial_t^3 q$,
which does not vanish on a moving vortex however well it is resolved. That
floor alone gives $\nu \propto \Delta^2\Delta t^2 \propto h^4$ at a fixed Courant
number, whatever $N$, so without a cutoff DynSGS sits on top of Galerkin up to
$N = 4$ and drifts above it as the mesh is refined from $N = 5$ on, its order
capped near 4 ([`setup_notes/dsgsCutoff.pdf`](../setup_notes/dsgsCutoff.pdf)).
The cutoff $c = 10^{-3}$ is above that floor and five decades below the
normalized residual of a shock, $O(10^2)$: with it DynSGS reproduces Galerkin
bit for bit wherever the vortex is resolved (measured here at $N = 6$, $32^2$,
$t = 0.5$; on the Alfvén wave at $N = 5$–$7$ on every mesh but the coarsest),
and still acts where it is not — on the coarsest meshes the residual exceeds
$c$ and the viscosity switches on, as it should. `JEXPRESSO_SV_CUTOFF=0` runs
without it.

**What changed since the sweep of 11 September 2026:**

| | 11 Sep 2026 | now |
|---|---|---|
| box | $[-5, 5]^2$ | $[-10, 10]^2$, the paper's (12 Sep): [the small box floors every order](#the-box) |
| normalization of the residual | per MPI rank | over the domain (12 Sep): [the rank count changed the answer 18×](#what-this-test-found) |
| form | conserved, on the physical primitives (neither form) | physical (`f9cdb02`, 13 Oct) |
| $\nu$ computed | at every Runge–Kutta stage | once per step (`b956ccd`, 14 Oct) |
| smoothness cutoff | did not exist | $10^{-3}$ (20 Oct) |
| $\Delta t$ | $\propto h$, per case (`tools/smooth_vortex_scan.sh`) | one $\Delta t$ for the sweep (`tools/smooth_vortex_mpi_scan.sh`) |
| unchanged | residual sensor, `:dsgs_rel = 1`, $C_{min} = 0$, $C_R = 1$, $C_{max} = 0.5$, hold 2, $t = 1$, the absolute velocity norms | |

## What a run produces

At the final time each run measures the **absolute** $L^1$, $L^2$ and
$L^\infty$ error of the **velocity** against the exact solution,

$$
\int_\Omega |\mathbf u_h - \mathbf u|\,d\Omega, \qquad
\Big(\int_\Omega |\mathbf u_h - \mathbf u|^2 d\Omega\Big)^{1/2}, \qquad
\max_\Omega |\mathbf u_h - \mathbf u| ,
$$

as the paper's Fig. 1 plots them, using the mesh's own nodal quadrature
weights (the lumped mass: $\omega_i\omega_j|J|$ summed over the elements at
each node; they integrate the domain area to 14 digits). It stores them in
`errors/nop<N>_nelx<M>_<tag>.dat`, whose header also carries the relative
norms and every parameter that made the number (box, $C_R$, $C_{max}$,
$C_{min}$, `:dsgs_rel`, hold, cutoff, $\Delta t$). The tag says which
experiment the record belongs to, so that no experiment overwrites another's:

| tag | the run |
|---|---|
| `dsgs` | DynSGS as the deck runs it, with the cutoff |
| `dsgs_nocut` | DynSGS without the cutoff (`JEXPRESSO_SV_CUTOFF=0`) |
| `dsgs_cut` | DynSGS with another cutoff |
| `dsgs_hold` | DynSGS with a startup hold other than 2 steps |
| `galerkin` | the plain Galerkin run (`JEXPRESSO_SV_VISC=none`) |

Every run then prints the table of **every** error in the store that carries
its final time and box, and draws the figures from them:

| file | contents |
|---|---|
| `convergence-it<n>.png` | the paper's Fig. 1 layout: the three norms side by side against #DOFs on log-log axes; one colour per polynomial order, the stabilization in the line style (DynSGS solid with a filled marker, DynSGS without the cutoff dash-dot, another cutoff or hold dash-dot-dot, Galerkin dashed with a hollow marker), dotted slope guides, and the measured rate $p$ of each curve in its legend; $t$ and $\Delta t$ in the titles |
| `convergence_L1-it<n>.png`, `_L2`, `_Linf` | the same panels one per file, at publication size |
| `convergence_nop<subset>[_L1…]-it<n>.png` | a subset of the orders on its own axes (`JEXPRESSO_SV_PLOT_NOPS`, or `tools/smooth_vortex_plot.jl 1,3`) |
| `<var>-it<n>.png` | the usual field panels |

in `output/MHD/smoothVortex/output/`. Records written before the norms became
absolute (no `norm=abs`) are never drawn; records without a box (`L=`) predate
the option and are kept; a DynSGS record without `cutoff=` was stored when the
deck ran without the cutoff (before 20 Oct 2026) and is drawn as `dsgs_nocut`.

**$p$ in the legend** is the measured convergence rate: the slope of that
curve's last two points on the log-log axes,

$$
p = \frac{\log(e_{i-1}/e_i)}{\log(h_{i-1}/h_i)},\qquad h = 1/\sqrt{\#\mathrm{DOFs}},
$$

so it is the order in the mesh size $h$ between the two finest resolutions that
curve has — the number the dotted slope guides are there to be compared
against. It needs at least two resolutions; with one, the legend carries the
order alone. The guides are the nominal rates of the lowest and highest orders
on the figure, $\min(N)+1$ and $\max(N)+1$.

So a sweep over orders and meshes builds the whole figure and each run
replaces only its own point. The hook is `user_plot_2d` in `user_plot.jl`,
the 2D counterpart of the `user_plot_1d` hook the 1D cases use, added to
`src/io/plotting/jeplots.jl`; unlike the 1D one it is additive, and the
generic field panels are still rendered.

## Running the sweep

`tools/smooth_vortex_mpi_scan.sh` runs every (stabilization, order, mesh)
case, `SV_JOBS` of them at a time, each on `SV_NP` MPI ranks (`SV_NP=1` runs
plain serial cases and needs no `mpiexec`):

```bash
tools/smooth_vortex_mpi_scan.sh                                   # orders 1-4, meshes 4-32, 4 ranks per case
SV_NP=1 SV_JOBS=4 SV_NOPS="4 5 6 7" tools/smooth_vortex_mpi_scan.sh   # serial cases, 4 at once
SV_NP=8 SV_NELX="8 16 32 64" tools/smooth_vortex_mpi_scan.sh
SV_KEEP=1 SV_VISC=dsgs SV_CUTOFF=0 SV_NP=1 SV_JOBS=4 SV_NOPS="4 5 6 7" \
    tools/smooth_vortex_mpi_scan.sh                               # add DynSGS without the cutoff to that sweep
DRYRUN=1 tools/smooth_vortex_mpi_scan.sh                          # print the plan, run nothing
```

**The store is moved aside first**, to `errors.prev`, unless `SV_KEEP=1`:
every run redraws the figures from the whole store, so a leftover sweep would
otherwise appear on the comparison of a new one. With `SV_KEEP=1` a sweep is
**resumable**: a case whose record is already there is skipped — a Galerkin
one always, a DynSGS one only if it was run with the cutoff asked for now (a
record without `cutoff=` is run again). Each case writes its log to
`logs/smoothVortex_<visc>_nop<N>_nelx<M>.log`, and the script reports every
case's error as it finishes.

| variable | default | meaning |
|---|---|---|
| `SV_NOPS` | `1 2 3 4` | polynomial orders |
| `SV_NELX` | `4 8 16 32` | elements per side |
| `SV_VISC` | `dsgs none` | DynSGS and/or the plain Galerkin run |
| `SV_NP` | 4 | MPI ranks per case |
| `SV_JOBS` | 1 | cases run at the same time; `SV_NP × SV_JOBS` is what you ask the machine for |
| `SV_DT` | one $\Delta t$ for the sweep, the deck's rule at its finest (nelx, nop) | a number fixes it; `auto` gives every case the deck's own $\Delta t \propto h$ |
| `SV_SOLVER` | `ck54` | the integrator: `vern9`, `vern7`, `dp8`, `ssprk54`, `tsit5` |
| `SV_TEND` | 1.0 | final time |
| `SV_L` | 20 | box width |
| `SV_CUTOFF` | the deck's, $10^{-3}$ | 0 runs DynSGS without the cutoff (`dsgs_nocut`) |
| `SV_HOLD` | 2 | startup hold in steps (`dsgs_hold` otherwise) |
| `SV_KEEP` | 0 | 1 adds to the store and skips the cases already in it |
| `SV_REDO` | 0 | with `SV_KEEP=1`, runs those cases again anyway |
| `SV_PLOT_NOPS` | — | a subset of the orders to draw on its own axes as well |
| `SV_CASE` | `MHD/smoothVortex` | `CompEuler/smoothVortex`, the hydrodynamic control (same box and meshes, no magnetic field, no cutoff by default) |
| `JULIA`, `MPIEXEC` | the macOS app's `julia`, else the `PATH`'s; `mpiexec` | |

`tools/smooth_vortex_scan.sh` is the serial loop of the first study: orders
4–7 on 4² to 32², one case after the other, each at the deck's own
$\Delta t \propto h$, with the macOS app's `julia`.

**One run by hand:**

```bash
JEXPRESSO_SV_NOP=4 JEXPRESSO_SV_NELX=16 julia --project=. src/Jexpresso.jl MHD smoothVortex
```

| variable | default | meaning |
|---|---|---|
| `JEXPRESSO_SV_NOP` | 4 | polynomial order |
| `JEXPRESSO_SV_NELX` | 16 | elements per side; picks `vortex_L<L>_<n>x<n>.msh` |
| `JEXPRESSO_SV_L` | 20 | box width $L$ |
| `JEXPRESSO_SV_DT` | $2\cdot10^{-3}\cdot 64/(n_{elx}N)$ | time step; the rule keeps the Courant number at about 0.2 against the fastest wave |
| `JEXPRESSO_SV_SOLVER` | `ck54` | time integrator |
| `JEXPRESSO_SV_TEND` | 1.0 | final time |
| `JEXPRESSO_SV_VISC` | `dsgs` | `none` for the plain Galerkin run |
| `JEXPRESSO_SV_SENSOR` | `residual` | `legacy` for the assembled-rate sensor |
| `JEXPRESSO_SV_CUTOFF` | $10^{-3}$ | `:dsgs_cutoff`; 0 for none |
| `JEXPRESSO_SV_HOLD` | 2 | `:dsgs_hold_steps` |
| `JEXPRESSO_SV_FREEZE` | (on) | 1 or 0 forces `:dsgs_freeze_stage` |
| `JEXPRESSO_SV_CR`, `JEXPRESSO_SV_CMAX` | 1, 0.5 | `:dsgs_CR`, `:dsgs_Cmax`; both zero keeps the sensor and applies no viscosity |
| `JEXPRESSO_SV_CMIN` | 0 | the background floor `:dsgs_Cmin` |
| `JEXPRESSO_SV_REL` | 1 | `:dsgs_rel`, the normalization floor; `1e-3` is what the kernels used before the fix of 11 Sep |
| `JEXPRESSO_SV_NORMS` | `domain` | `:dsgs_norms`; `rank` makes the answer depend on the rank count, `element` normalizes per element |
| `JEXPRESSO_SV_PLOT_NOPS` | — | as `SV_PLOT_NOPS` |

**The meshes** come from one parameterized `vortex_periodic.geo`;
`vortex_L20_<n>x<n>.msh` for $n$ = 4 … 256 and `vortex_L10_<n>x<n>.msh` for
$n$ = 4 … 64 are in the repository, and the sweep makes any other with

```bash
SV_L=20 SV_NELX="16 32 64" tools/smooth_vortex_mesh.sh     # needs gmsh
```

## What this test found

**1. Two defects in the DynSGS sensor itself** (11 Sep 2026), fixed in
`src/kernel/physics/SGS.jl` and `src/kernel/operators/rhs.jl` (DSGS.md §4.2,
§4.4 and the defect list in §6). Before, DynSGS held $\nu$ on its cap
$C_{max}\Delta\lambda$ — an $O(h)$ viscosity, the same for every polynomial
order — and the case converged at **rate 1** while the same scheme without it
converged at its design order. Velocity $L^1$ at $t = 1$, `:nop => 4`, on the
$[-5, 5]^2$ box and with the deck of that day:

| elements per side | 4 | 8 | 16 | 32 | rate |
|---|---|---|---|---|---|
| DynSGS, before | 1.004e-2 | 5.240e-3 | 2.678e-3 | 1.142e-3 | 1.0 |
| **DynSGS, after** | 4.26e-3 | 5.220e-5 | 2.208e-6 | 1.639e-7 | 4.6, 3.8 |
| plain Galerkin | 1.570e-3 | 3.541e-5 | 2.039e-6 | 1.631e-7 | 4.1, 3.6 |

The runs that isolated the cause are worth repeating after any change to the
sensor:

```bash
# the sensor with no viscosity applied: the solution is the Galerkin one
JEXPRESSO_SV_CR=0 JEXPRESSO_SV_CMAX=0 JEXPRESSO_SV_NELX=16 \
    julia --project=. src/Jexpresso.jl MHD smoothVortex
```

**2. The normalization scope changes the answer.** Normalizing each rank's
residual by its own subdomain's spread (`rank`) makes $\nu$ — and therefore
the solution — depend on how the domain was cut: at `:nop => 6`, $32^2$,
$\Delta t = 3.3\times10^{-4}$, $t = 1$, the velocity $L^1$ was 3.935e-7 on 2
ranks against 7.155e-6 on 4, an 18× difference from nothing but the partition.
The domain norm (the default) costs two small `Allreduce`s per right-hand side
and is what Dao & Nazarov's eq. 4.8 means by the norm over $\Omega$.

**3. The physical form.** Until `f9cdb02` the deck ran the conserved form on
the physical primitives $(\rho, u, v, T, w, \mathbf B, \psi)$ that
`user_primitives.jl` hands the kernel — one kinematic $\nu$ on each, no stress
and no work terms, neither form. The comparison of before, after and Galerkin
is in [`setup_notes/smoothVortex.pdf`](../setup_notes/smoothVortex.pdf).

**4. The time-difference floor, and the cutoff.** Computing $\nu$ once per
step (`b956ccd`) keeps the floor of the residual at $O(\Delta t^2)$ instead of
$O(\Delta t)$, and the cutoff removes it: see
[The DynSGS setup](#the-dynsgs-setup) above.

## The box

The vortex is a Gaussian, so on $[-L/2, L/2]^2$ the exact solution is not
periodic: the velocity perturbation at the middle of an edge is
$(L/2)\,e^{(1-(L/2)^2)/2}/2\pi$ with opposite signs on opposite edges, so the
initial condition jumps across the periodic seam by twice that — 4.9e-6 on the
$L = 10$ box of Balsara and of Dao & Nazarov's text. That jump is a
discontinuity in the *data*: it floors any accuracy study near $10^{-5}$,
whatever the order, and a high order reaches the floor sooner. $L = 15$ puts it
at $10^{-12}$ and $L = 20$ at machine zero, so $L = 20$, the $[-10, 10]^2$ the
paper's figures use, is the default (`SV_L=10` for the small box).

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
| one $\Delta t$ for the sweep | the time error is a constant rather than something that shrinks at fourth order and pollutes the slope. `tools/smooth_vortex_mpi_scan.sh` does this by default, from the finest (nelx, nop) of the sweep; `JEXPRESSO_SV_DT` for a single run |
| `JEXPRESSO_SV_SOLVER=vern9` | `Vern9` (9th order), `dp8` (8th), `vern7`, `ssprk54`, `tsit5`, or `ck54` for the default |

On the $[-5, 5]^2$ box, at 32×32 elements and `:nop => 4`, the integrator did
not move the error (absolute velocity $L^1$ at $t = 1$):

| | |
|---|---|
| `ck54`, $\Delta t = 10^{-3}$ | 2.318e-5 |
| `vern9`, $\Delta t = 10^{-3}$ | 2.320e-5 |
| `vern9`, $\Delta t = 3.3\times10^{-4}$ | 2.314e-5 |

— an error near that box's seam floor. On the default box the floor is gone,
so the honest check on any sweep is to re-run its **finest** case at half
`SV_DT` and see that the error does not move; if it does, the time error is
part of what the figure shows.

## The final time is part of the measurement

`SV_TEND` shortens a run, which is the right way to check that a sweep's
machinery works before committing hours to it — and the wrong thing to leave
set. At $t = 0.05$ the vortex has moved a twentieth of the box, the error is
dominated by the initial transient, and it barely converges with resolution:
measured on the $[-5, 5]^2$ box at $\Delta t = 3.3\times10^{-4}$, `:nop => 4`
on 32² gave 7.539e-6 and `:nop => 6` on 32² gave 1.433e-6, against 2.3e-5 and
far less at $t = 1$ — flat lines and negative rates that look exactly like a
broken solver. The figures therefore carry the final time and the step in
their titles.

## Comparing P1 and P3 directly, as in the paper

Dao & Nazarov's Fig. 1 is a low-order comparison. Any subset of the orders in
the store can be drawn on its own axes, either while running —

```bash
SV_PLOT_NOPS="1 3" tools/smooth_vortex_mpi_scan.sh
```

which writes `convergence_nop1-3-it<n>.png` beside the all-orders figure — or
afterwards, from the stored errors alone, without running anything:

```bash
julia --project=. tools/smooth_vortex_plot.jl 1,3
julia --project=. tools/smooth_vortex_plot.jl          # every order
julia --project=. tools/smooth_vortex_plot.jl 1,3 --t=1.0 --out=figs
```

`tools/smooth_vortex_plot.jl` reads `errors/*.dat`, prints the error table and
calls the case's own plotting code, so its figures are the same ones a run
produces (as `-it0`).

## Running it on many cores

```bash
SV_NP=8  SV_NELX="8 16 32 64" tools/smooth_vortex_mpi_scan.sh
SV_NP=16 SV_JOBS=4 SV_NOPS="4 6" tools/smooth_vortex_mpi_scan.sh
```

On a SLURM cluster, `tools/smooth_vortex_slurm.sh` is the same sweep as a
batch job — as written, orders 4 and 6 on 32² and 64² elements, both panels,
16 ranks per case and 4 cases at a time — and `tools/smooth_vortex_wulver.sh`
the same on NJIT's Wulver:

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

The norms are MPI-correct: each unknown is weighed once, by the rank that owns
it (`mesh.gip2owner`, the map the DSS assembler uses), and the
degree-of-freedom count on the abscissa is the exact $(n_{elx}N)^2$ of this
periodic box, so a point does not move when the rank count changes.

## Why the background floor is off here

`:dsgs_Cmin` is **0** in this case, unlike the shock cases. The floor is a
viscosity $C_{min}\Delta\lambda$ that is applied everywhere regardless of the
residual; on a smooth solution it is an $O(h)$ error that no polynomial order
can beat, so a sweep with the floor on measures the floor and not the scheme —
every order would collapse onto slope 1. Turn it on with `JEXPRESSO_SV_CMIN`
if that is what you want to see.

## Files

| file | role |
|---|---|
| `user_inputs.jl` | the deck and the environment overrides |
| `initialize.jl` | the vortex, its exact-solution helpers `sv_state`/`sv_wrap`, and the equilibrium derivation |
| `user_plot.jl` | the error norms, the record tags, the store and the convergence figures |
| `user_flux.jl`, `user_source.jl`, `user_bc.jl`, `user_primitives.jl` | the GLM-MHD equations, identical to `orszagTangBormanis2024` |
| `vortex_periodic.geo` | the parameterized periodic mesh |
| `tools/smooth_vortex_mpi_scan.sh`, `tools/smooth_vortex_scan.sh` | the sweeps |
| `tools/smooth_vortex_plot.jl` | the table and the figures from the store alone |
| `tools/smooth_vortex_mesh.sh` | the meshes |
| `tools/smooth_vortex_slurm.sh`, `tools/smooth_vortex_wulver.sh` | the sweep as a batch job |
