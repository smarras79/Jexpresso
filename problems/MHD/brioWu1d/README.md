# 1D ideal MHD: the Brio–Wu shock tube (Dao & Nazarov 2022, §5.2)

The MHD counterpart of [`problems/CompEuler/sod1d`](../../CompEuler/sod1d):
the Riemann problem of Brio & Wu (*J. Comput. Phys.* 75:400, 1988) in the
setting of Dao & Nazarov, *J. Sci. Comput.* 92:77 (2022), Section 5.2,

| | $\rho$ | $u$ | $p$ | $B_x$ | $B_y$ |
|---|---|---|---|---|---|
| $x \in [0, 0.5)$ | 1 | 0 | 1 | 0.75 | 1 |
| $x \in [0.5, 1]$ | 0.125 | 0 | 0.1 | 0.75 | −1 |

$v = w = B_z = 0$, $\gamma = 2$, magnetic pressure $\tfrac12|\mathbf{B}|^2$.
**Final time.** The paper writes $\hat t = 0.2$, but its Fig. 2 shows the
waves at the positions of $t = 0.1$ on $(0, 1)$ (rarefaction head at
$x = 0.32$, compound wave at $0.47$, contact at $0.56$, slow shock at $0.64$,
right fast rarefaction at $0.84$–$0.87$): that is Brio & Wu's original
setting, $t = 0.2$ on $(-1, 1)$, which maps to $t = 0.1$ on $(0, 1)$. The case
runs to $t = 0.1$ so that its figure and the paper's coincide. (At $t = 0.2$
on $(0, 1)$ the right fast wave has already left the domain.) The solution
carries, from left to right, a fast
rarefaction, a slow compound wave (a slow shock attached to a slow
rarefaction), a contact discontinuity, a slow shock and a fast rarefaction,
which is why it is the standard test of whether an MHD solver represents
all wave types.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD brioWu1d
```

Serial; a few seconds of time stepping after compilation. Output goes to
`./output/MHD/brioWu1d/output/`. By default (`:plot_user => true`) the case
writes the figure of the paper's Fig. 2: `density-it<n>.svg` at every
output time ($t = 0, 0.025, \dots, 0.1$), density against $x$ on the paper's
axes, the numerical solution labelled `nop <N>, <n> DOFs`, and at the final
time the reference solution with four zoom boxes: the paper's three (the
foot of the fast rarefaction, the compound wave, the contact) and one more
on the slow shock at $x \approx 0.647$ — `user_plot.jl`. Every
numerical curve is drawn with a broken line and its own marker
(nop 4 red dashed circles, nop 5 blue dash-dot squares, nop 6 green dotted
diamonds, nop 7 purple dash-dot-dot triangles); the reference is the only
solid line, so the curves are told apart in print and by line type, not by
colour alone. The background floor `:dsgs_Cmin` is written in the title when
all the curves on the figure share it, and in their legend entries when they
do not.

Beside it, at every output time, **one figure per solution quantity** —
`rho-it<n>.svg`, `u-it<n>.svg`, `v-it<n>.svg`, `w-it<n>.svg`, `p-it<n>.svg`,
`Bx-it<n>.svg`, `By-it<n>.svg`, `Bz-it<n>.svg` — each on the full canvas
against the reference, with zoom boxes on $u$, $v$, $p$ and $B_y$ where the
DynSGS solution leaves it (the slow shock at $x \approx 0.647$ and the
compound wave at $x \approx 0.47$), and the paper's four on $\rho$. $B_x$,
$w$ and $B_z$ are constant by construction, so those three are centred on
their own value with the spread in the title rather than auto-scaled to an
arbitrary window; this run reports `spread 0.0e+00` for all three.
`mu_dsgs-it<n>.svg` carries the structure of the DynSGS coefficient: $\rho$
and $B_y$ on top, $\nu$ per equation slot below on the same $x$ axis.

**Every figure is vector graphics**, because the zoom boxes exist to be looked
at closely and a raster figure is pixels at that magnification. `.svg` by
default (~80–250 kB each), `.pdf` with `JEXPRESSO_BW_FIGFMT=pdf` for
`\includegraphics`, `.png` with `JEXPRESSO_BW_FIGFMT=png` (or
`:plot_format` in the deck) if raster is wanted. `:plot_user => false` falls
back to the generic panels of the 1D plotter.

### Comparing polynomial orders, and the convergence history

At the final time a run stores its density profile and its DynSGS
coefficient $\nu$ in `curves/nop<N>_dof<M>_<form>.dat` in this case directory
(`curves/<store>/` when `JEXPRESSO_BW_STORE` names one: the scan keeps one per
$C_{min}$, `curves/Cmin0.06/`, `curves/Cmin0.03/`, …), and the figures are
drawn from **every** curve stored there whose `<form>` — the coefficient form
and the length-scale convention, the two things that make one run a different
method from another — matches the run drawing them. A sweep over orders and
resolutions therefore builds the whole comparison and each run replaces only
its own (order, DOFs) point:

```bash
BW_JOBS=4 tools/brio_wu_order_scan.sh                    # orders 4-7 at 600 and 1200 DOFs, C_min = 0.06
BW_JOBS=4 BW_CMIN="0.06 0.03 0" tools/brio_wu_order_scan.sh   # one complete sweep per C_min, to compare
BW_JOBS=4 BW_DOFS=600 tools/brio_wu_order_scan.sh        # the paper's 600 DOFs alone
BW_DOFS="150 300 600 1200" tools/brio_wu_order_scan.sh   # also 150 and 300: the convergence figure
JULIA=/path/to/julia tools/brio_wu_order_scan.sh         # default: the macOS app's julia, else the PATH's
julia --project=. tools/brio_wu_plot.jl --store=Cmin0.06 # the figures again, from the store alone
```

**Every value of `BW_CMIN` is a sweep of its own**: its runs store their
curves in `curves/Cmin<C_min>/`, and when the last of them is done the scan
draws every figure of that sweep in
`output/MHD/brioWu1d/output/Cmin<C_min>/`. With `BW_HSCALE` or `BW_NODAL`
set, both names carry it too (`Cmin0.06_nop`), so different methods never
share a figure. `BW_JOBS=n` runs up to $n$ orders at once, one Julia session
per order (each compiles the case once, then runs its resolutions); the log
of order $N$ goes to `output/MHD/brioWu1d/scan-Cmin<C_min>-nop<N>.log`.

**Each sweep starts from an empty store**: the curves and the figures of the
same $C_{min}$ from an earlier sweep are deleted first (`BW_KEEP=1` keeps
them, to finish a partial sweep), so a leftover curve can never appear on
the comparison. Inside a sweep the runs only store their curves
(`JEXPRESSO_BW_DRAWSTORE=0`); a run started by hand redraws the comparison
from the whole store itself, and says so when the store holds curves of
another method, which it leaves off the figure.

| figure | contents |
|---|---|
| `rho_By_dof<M>-it<n>.svg` | **the comparison figure, one per resolution**, in the layout of the Brio–Wu figure of Dao, Nazarov & Tomas ([arXiv:2310.18467](https://arxiv.org/abs/2310.18467), Figs. 1–2): **left** $\rho$, **right** $B_y$, of every order at $M$ DOFs against the reference, $\rho$ with the four zoom boxes of the density figure, $B_y$ with their three: the top of the compound wave ($x \in [0.416, 0.476]$, upper left), the plateau up to the slow shock ($[0.56, 0.65]$, lower left) and the right fast rarefaction ($[0.763, 0.913]$, right, under the legend) |
| `rho_By-it<n>.svg` | the same for the finest resolution in the store |
| `nu_dof<M>-it<n>.svg`, `nu-it<n>.svg`, `nu_all-it<n>.svg` | **the DynSGS coefficient $\nu$ alone**, of the same curves (the $\rho u$ slot, which in this deck's conserved form is the one $\nu$ of every equation; per element, so a staircase through the nodes). Linear axis, or $\log_{10}$ when $\nu$ spans more than two decades ($C_{min}=0$); `nu_all` holds every order at every resolution (colour the order, line type the resolution) and falls into one band per resolution ($\nu \propto \Delta \propto 1/\#\mathrm{DOFs}$). $\rho$ and $B_y$ get no such figure: with every resolution in one $B_y$ panel there are too many curves to tell which one overshoots |
| `density_dof<M>-it<n>.svg` | **one file per resolution**: every order that ran at $M$ degrees of freedom, against the reference, with the four zoom boxes. A sweep over 150, 300, 600 and 1200 DOFs leaves `density_dof150-it<n>.svg`, `density_dof300-it<n>.svg`, … side by side |
| `density-it<n>.svg` | the same figure for the finest resolution in the store, under the plain name |
| `density_nu_dof<M>-it<n>.svg`, `density_nu-it<n>.svg`, `density_nu_all-it<n>.svg` | the earlier layout, $\rho$ and $\nu$ side by side, in place of the $\rho$–$B_y$ figures: `JEXPRESSO_BW_RIGHT=nu`, or `--right=nu` for `tools/brio_wu_plot.jl`; also what a store from a sweep run before $B_y$ was stored gets |
| `convergence-it<n>.svg` | the layout of the paper's Fig. 1: the relative $L^1$, $L^2$ and $L^\infty$ error of $\rho$ against $1/\#\mathrm{DOFs}$ on log-log axes, one line per order, with slope guides and the measured rate in each legend entry |
| `convergence_smooth-it<n>.svg` | the same, restricted to $x \in (0.33, 0.41)$ inside the fast rarefaction — the one smooth, non-constant part of this solution |

The abscissa is $1/\#\mathrm{DOFs}$ because the problem is 1D and the mesh
size is $h \propto 1/\#\mathrm{DOFs}$; the paper's $1/\sqrt{\#\mathrm{DOFs}}$
is the same quantity for its 2D vortex.

The error is measured against `reference_hll.dat`, which is itself a
**first-order** finite-volume solution on 10 000 cells. Its own error is
$O(\Delta x) \approx 10^{-4}$ in the smooth fan, but a few cells of smearing at
every discontinuity — far larger there than the difference between two
spectral-element orders. That is why the smooth-window figure is the one that
can say anything about the order, and the whole-tube figure mostly measures
how wide the jumps are.

Overrides, all optional, so a sweep needs no file edits:

| variable | default | meaning |
|---|---|---|
| `JEXPRESSO_BW_NOP` | 4 | polynomial order |
| `JEXPRESSO_BW_DOFS` | 600 | degrees of freedom; sets the element count as DOFs/nop |
| `JEXPRESSO_BW_NELX` | — | element count, overriding `JEXPRESSO_BW_DOFS` |
| `JEXPRESSO_BW_DT` | — | time step, overriding the rule below |
| `JEXPRESSO_BW_CMIN` | 0.06 | the DynSGS background floor `:dsgs_Cmin` |
| `JEXPRESSO_BW_HSCALE` | `ngl` | the DynSGS length scale: `ngl` is the kernels' own $\Delta_K/(k+1) = h_e/(N+1)$, eq. (15) of the DynSGS-MHD paper; `nop` is Dao & Nazarov's $\Delta_K/k$ |
| `JEXPRESSO_BW_NODAL` | 0 | 1 for the nodal coefficient (`:ldsgs_nodal`, their eq. 4.10) |

At a fixed DOF count the orders run at equal cost: 150, 120, 100 and 86
elements for $N = 4, 5, 6, 7$ at 600 DOFs, i.e. 601, 601, 601 and 603 points.
The time step follows the resolution, $\Delta t = 5\times10^{-5}\cdot 600/\#\mathrm{DOFs}$,
so the Courant and the DynSGS parabolic numbers are the same at every point
of a sweep. At 600 DOFs the smallest LGL spacing falls from
$1.15\times10^{-3}$ at $N = 4$ to $7.5\times10^{-4}$ at $N = 7$, so the Courant
number against the fast speed of the right state runs from 0.16 to 0.25 and
the parabolic number stays below 0.25.

A sweep's directory therefore holds its complete comparison: one density
figure, one $\rho$–$B_y$ figure and one $\nu$ figure per resolution, the
$\nu$ of every resolution together, and the two convergence histories.
`tools/brio_wu_plot.jl` draws the same from a store alone (`--store=Cmin0.03`
for the sweep at $C_{min} = 0.03$, `4,7` for a subset of the orders,
`--dofs=600` for one resolution, `--fmt=pdf`, `--out=dir`); the scan calls it
at the end of each sweep. Curves stored before the coefficient or $B_y$ was
are drawn on the density figures and left off the figures that need them
(rerun them); a store with no $B_y$ at all, from a sweep run before it was
stored, gets the $\rho$–$\nu$ figures in place of the $\rho$–$B_y$ ones.
`rm -r problems/MHD/brioWu1d/curves` clears every store; it is
git-ignored. Curves whose stored final time differs from the current one are
ignored, so changing `:tend` cannot silently mix solutions from different
times.

**The two regimes, and what the paper's Fig. 2 shows.** Dao & Nazarov compare
$\mathbb{P}_1$ against $\mathbb{P}_3$ — "under the same number of degrees of
freedom, the $\mathbb{P}_3$ solution captures the compound structure more
accurately than the $\mathbb{P}_1$ solution". That is the low-order regime,
where raising the order buys a large gain in the accuracy of the underlying
Galerkin scheme.

**The length scale, and why it decides the comparison.** The kernels of
`src/kernel/physics/SGS.jl` build the coefficient on $\Delta = \Delta_K/(k+1)$
(`ngl`, the number of LGL points per element); Dao & Nazarov's eq. 4.10 uses
$h_K/k$, the polynomial **degree**. The difference is a factor $(k+1)/k$ — and
at a *fixed* number of degrees of freedom, where $n_{elx} = \#\mathrm{DOFs}/k$,

$$
\frac{\Delta_K}{k+1} = \frac{L}{\#\mathrm{DOFs}}\cdot\frac{k}{k+1},
\qquad
\frac{\Delta_K}{k} = \frac{L}{\#\mathrm{DOFs}} ,
$$

so Dao & Nazarov's convention is the same for every order while the kernels' grows
with $k$ (0.8 of it at $k=4$, 0.875 at $k=7$). That factor varies by 75 % over
$k = 1\ldots7$, it is monotone in $k$, and it enters the cap and the floor
linearly and the residual viscosity $C_R\Delta^2\mathcal{R}$ quadratically. On
a solution made of jumps, where the artificial viscosity and not the
polynomial order sets the error, it is the whole comparison. Measured here on
the nodal coefficient at 601 DOFs, relative $L^1$ error of $\rho$:

| length scale | spread over orders 1-7 | order 1 $\to$ 7 |
|---|---|---|
| $\Delta_K/(k+1)$ | 64 % | +64 % (each order worse than the last) |
| $\Delta_K/k$ | 5.5 % | −2.1 % |

With $\Delta_K/k$ the monotone penalty on the higher orders disappears, and in
the smooth fan order 1 is the worst of the seven ($5.39\times10^{-3}$ against
$3.3$–$3.5\times10^{-3}$ for orders 2-7) — which is the statement of the
paper's Fig. 2, that $\mathbb{P}_3$ captures the compound structure more
accurately than $\mathbb{P}_1$ at equal degrees of freedom.

**This case runs $\Delta_K/(k+1) = h_e/(N+1)$ by default**, the kernels' own
and the $\Delta_e$ of eq. (15) of the DynSGS-MHD paper, so that $C_R = 1$,
$C_{max} = 0.5$ and $C_{min} = 0.06$ mean here what they mean in the paper and
in every other case. Until October 2026 the default was Dao & Nazarov's
$\Delta_K/k$, for the order comparison above; `JEXPRESSO_BW_HSCALE=nop`
(`BW_HSCALE=nop` for the scan) still gives it. That is a *deck-level* switch:
`user_inputs.jl` folds the factor into this case's own coefficients —
$C_{max}$ and $C_{min}$ by $(k+1)/k$, $C_R$ by $((k+1)/k)^2$, which is
algebraically identical to changing $\Delta$ — and no kernel is touched. Two
things follow from the default. At equal DOFs the higher orders get the larger
viscosity (orders 4 to 7: 9 % more floor and cap, 20 % more residual
viscosity), the penalty the table measures. And the absolute level is lower
than with $\Delta_K/k$, which is $(k+1)/k$ *larger* (at order 4,
$7.1\to5.6\times10^{-3}$ over the whole tube); that larger scale is also what
pushed the rising bubble past the explicit viscous stability limit at
start-up, which is why the kernels use $(k+1)$.

`JEXPRESSO_BW_CMIN=0` removes the part of the viscosity that is independent of
the solution altogether.

## What is implemented

- **Equations** (`user_flux.jl`): the 1D ideal-MHD system in conservative
  form for $(\rho, \rho u, \rho v, \rho E, \rho w, B_x, B_y, B_z)$ — the 2D
  cases' slot order without the GLM field, which 1D does not need
  ($\nabla\cdot\mathbf{B} = \partial_x B_x = 0$ holds exactly with $B_x$
  constant; its slot carries a zero flux).
- **Stabilization: DynSGS** (`:visc_model => DSGS_MHD()`), with the 1D MHD
  kernel added to `kernel/physics/SGS.jl` for this case. The deck runs the
  default **element form** (one $\nu$ per element); the **nodal form**
  (`:ldsgs_nodal => true`) is Dao & Nazarov's own and gives the same profile
  once the $C_{min}$ floor is on (see below): at every node
  the element-wise residual (DSGS.md §1.2) averaged at the node, normalized by their eq. 4.7 with
  $C_l = 0.4$ (`:dsgs_Cl`), the maximum over the equations (eq. 4.8), and
  $\nu_i = \min(C_{max}h_i\lambda_i, C_R h_i^2 R_i)$ with the fast
  magnetosonic speed, $h_i = \Delta x/(k+1)$ (the element form's $\Delta$), $C_{max} = 0.5$ (`:dsgs_Cmax`), $C_R = 1$ (`:dsgs_CR`)
  (eq. 4.10); $\nu$ is a continuous field the element loop interpolates,
  so the diffusive flux has no jump at element interfaces. Applied in the
  **conserved form**, $\nabla\cdot(\nu\nabla q)$ on every slot (`user_primitives.jl`),
  i.e. exactly conservative, with the magnetic and kinetic energy removed
  from $\mathbf{B}$ and $\rho\mathbf{v}$ accounted for in $E$. This is the
  form the residual method reduces to at its cap (a Lax–Friedrichs-type
  viscosity on the conserved variables). The physical-form coefficients of
  Dao & Nazarov's eq. 4.4 ($\nu$ on $\rho$, $\rho\nu$ on $\mathbf{u}$,
  $\rho\nu/\mathrm{Pr}$ on $T$, $\nu$ on $\mathbf{B}$) are available from
  the same kernel with `:dsgs_conserved => false`, `:dsgs_nazarov_energy => true`
  and a `user_primitives!` returning $(\rho, u, v, T, w, \mathbf{B})$, but the
  1D viscous loop has no $\tau\cdot u$ or $\eta\mathbf{B}\cdot\nabla\mathbf{B}$
  work terms, so that form does not conserve total energy here.
- **Boundaries** (`user_bc.jl`): both ends held at the initial states, both
  still undisturbed at $t = 0.1$ (outermost waves at $x \approx 0.32$ and
  $0.87$). The right fast wave leaves the domain at $t \approx 0.13$; past
  that the pinned right end shows a one-node glitch (a free end is not an
  outflow condition in CG and drained the domain when tried).
- **Resolution**: 150 elements at $N = 4$ = 600 LGL points, the paper's
  "600 DOFs" ($\mathbb{P}_3$) case of its Fig. 2(a); `:nelx => 300`
  (`JEXPRESSO_BW_NELX=300`) is its 1200-DOF case. `JEXPRESSO_BW_NOP` changes
  the order at fixed DOFs, see above. $\Delta t = 5\times10^{-5}$ is a Courant number of 0.19
  against the fast speed of the right state (3.75) on the smallest LGL
  spacing; 2000 steps to $t = 0.1$.
- **Reference solution** (`reference_hll.dat`, `user_analytic.jl`): the paper
  compares against Athena at 10 000 points; the file shipped here is a
  first-order finite-volume (HLL flux) solution on 10 000 cells at $t = 0.1$,
  sampled at 2000 points, generated by `tools/brio_wu_reference.jl` (20 s).
  A first-order scheme at that resolution resolves every wave of the
  problem to well within the width of a 600-point spectral-element
  solution's shocks; it is not Athena's exact-Riemann-solver reference, and
  the compound wave's peak is the one place where the two would differ
  visibly.

## The ripples behind the compound wave, and the floor

Without a background floor the density plateau between the fast
rarefaction and the compound wave ($x = 0.42$–$0.46$ at $t = 0.1$) carries
ripples of $\pm0.5\,\%$ with one wiggle per element. They are radiated by
the compound wave — a slowly moving slow shock, the classical source of
post-shock oscillations — into the plateau; they are present at
$\mathbb{P}_3$ (200 elements, the paper's 600-DOF layout) and $\mathbb{P}_4$,
and with the element and the nodal coefficient alike, so neither the order
nor the continuity of $\nu$ is their cause. The residual viscosity cannot
remove them: $\nu = C_R h^2 R$ and the residual of a ripple is its
amplitude over $h$, so $\nu \sim C_R h\,\times$ amplitude, far below the
first-order cap for a $0.5\,\%$ ripple. The case therefore runs with
`:dsgs_Cmin => 0.06`, a floor of 6 % of the first-order viscosity
$C_{max}h(|u| + c_f)$: it damps an element-scale mode at a rate
$\nu(\pi/h)^2 \approx 600$ per unit time and diffuses a resolved profile by
$\sqrt{2\nu t} \approx 0.006$ over the run, less than one element; with
0.03 the ripples remain, with 0.06 none (measured), at the price of a
contact about one node wider; 0.03 leaves a faint trace of them (±0.1 %)
with sharper shocks, 0 gives the sharpest shocks with the ripples in full.
Those were measured on Dao & Nazarov's $\Delta_K/k$. On the default
$h_e/(N+1)$ the same 0.06 is an effective 0.048–0.053 on $\Delta_K/k$
($N = 4$–7), and re-measured at 600 DOFs nothing is visible at the figure's
scale: between $x = 0.425$ and $0.455$ successive extrema of $\rho$ differ by
at most 0.4 % on either length (dip at the foot of the compound wave
included), with a few more small extrema at $N = 4, 5$ on $h_e/(N+1)$.
(With the residual and the time stencil of the runs before September 2026,
DSGS.md §1.2 and §4.4, 0.03 was enough: the old sensor over-fired on the
moving waves and supplied the rest of the damping itself;
`:dsgs_sensor => "legacy"` reproduces that sensor.) Dao & Nazarov's $\mathbb{P}_3$
solution shows no ripples with no floor; their elements integrate the
nonlinear flux exactly on uniform nodes, this code's collocated LGL flux
does not, and that is the remaining difference between the two
discretizations.

## Results

Serial, 600 points, $\Delta t = 5\times10^{-5}$, 2000 steps: a few seconds
of time stepping. At $t = 0.1$ (`density-it5.svg`; the per-variable figures with
`:plot_user => false`):

- every wave is where the reference and the paper's Fig. 2 put it: fast
  rarefaction from $x = 0.32$ to $0.42$, compound wave at $0.47$, contact at
  $0.56$, slow shock at $0.64$, right fast rarefaction at $0.84$–$0.87$;
- the compound wave's density peak is $0.84$ against the (first-order,
  smoothed) reference's $0.80$ and the paper's $\mathbb{P}_3$ value of
  $\approx 0.82$; the contact and the slow shock are 2–3 nodes wide, the
  paper's zoom boxes show the same smearing for its $\mathbb{P}_3$ solution;
  $p$, $B_y$, $u$, $v$ follow the reference to plotting accuracy;
- the DynSGS coefficient (`mu_dsgs-it5.svg`, and `nu-it5.svg`) is on its floor $C_{min}\Delta_e(|u| + c_f)$
  everywhere, at 600 DOFs and $N = 4$–7: $1.4$–$1.6\times10^{-4}$ on the left
  state, $2.9$–$3.2\times10^{-4}$ on the right one and at most
  $3.4$–$3.8\times10^{-4}$ behind the slow shock, the higher orders a little
  higher (the equal-DOF penalty of $h_e/(N+1)$ above). The residual part never
  exceeds it: with the floor off (`JEXPRESSO_BW_CMIN=0`) it peaks at the
  compound wave and at the slow shock at about the floor's size and falls to
  $10^{-5}$–$10^{-6}$ elsewhere (measured on $\Delta_K/k$);
- the relative $L^1$ error of $\rho$ over the tube is $6.1$, $6.6$, $6.7$ and
  $7.0\times10^{-3}$ for $N = 4$–7; on $\Delta_K/k$ it was $7.6$–$8.0\times10^{-3}$,
  and $2.7$–$3.7\times10^{-3}$ there with the floor off.

## Files

| file | contents |
|---|---|
| `user_inputs.jl` | solver parameters (Δt, tend, DynSGS constants, mesh, output) |
| `initialize.jl` | the left/right states; output variable list |
| `user_flux.jl` | 1D ideal-MHD fluxes, `γ_mhd = 2`, `pressure_mhd1d` |
| `user_bc.jl` | Dirichlet ends |
| `user_primitives.jl` | conserved variables as DynSGS primitives; `user_uout!` (ρ, u, v, p, By) |
| `user_analytic.jl` | interpolates `reference_hll.dat` for the overlay at t = 0.1 |
| `user_plot.jl` | the paper's Fig. 2 density figure with zoom boxes (`density-it<n>.svg`), one figure per solution quantity, and the DynSGS coefficient figure |
| `reference_hll.dat` | the reference solution |
| `EQUATIONS.md` | the equations and the wave structure |
