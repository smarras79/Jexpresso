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
writes the figure of the paper's Fig. 2: `density-it<n>.png` at every
output time ($t = 0, 0.025, \dots, 0.1$), density against $x$ on the paper's
axes, the numerical solution in red labelled with its polynomial order and
number of degrees of freedom, and at the final time the reference solution
in black with the paper's three zoom boxes (the foot of the fast
rarefaction, the compound wave, the contact) — `user_plot.jl`. With
`:plot_user => false` the output is the format of `CompEuler/sod1d`: one
figure `fields-it<n>.png` per output time with a panel per output variable
($\rho$, $u$, $v$, $p$, $B_y$; the reference dashed at the final time) and a
last panel with the DynSGS coefficient per element (`:plot_matrix => false`
writes those panels as separate files).

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
  the assembled lumped-mass BDF2 residual, normalized by their eq. 4.7 with
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
  "600 DOFs" ($\mathbb{P}_3$) case of its Fig. 2(a); `:nelx => 300` is its
  1200-DOF case. $\Delta t = 5\times10^{-5}$ is a Courant number of 0.19
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
`:dsgs_Cmin => 0.03`, a floor of 3 % of the first-order viscosity
$C_{max}h(|u| + c_f)$: it damps an element-scale mode at a rate
$\nu(\pi/h)^2 \approx 300$ per unit time and diffuses a resolved profile by
$\sqrt{2\nu t} \approx 0.004$ over the run, less than one element; with
0.01 a trace of the ripples remains, with 0.03 none (measured), at the
price of a contact about one node wider. Dao & Nazarov's $\mathbb{P}_3$
solution shows no ripples with no floor; their elements integrate the
nonlinear flux exactly on uniform nodes, this code's collocated LGL flux
does not, and that is the remaining difference between the two
discretizations.

## Results

Serial, 600 points, $\Delta t = 5\times10^{-5}$, 2000 steps: a few seconds
of time stepping. At $t = 0.1$ (`density-it5.png`; `fields-it5.png` with
`:plot_user => false`):

- every wave is where the reference and the paper's Fig. 2 put it: fast
  rarefaction from $x = 0.32$ to $0.42$, compound wave at $0.47$, contact at
  $0.56$, slow shock at $0.64$, right fast rarefaction at $0.84$–$0.87$;
- the compound wave's density peak is $0.84$ against the (first-order,
  smoothed) reference's $0.80$ and the paper's $\mathbb{P}_3$ value of
  $\approx 0.82$; the contact and the slow shock are 2–3 nodes wide, the
  paper's zoom boxes show the same smearing for its $\mathbb{P}_3$ solution;
  $p$, $B_y$, $u$, $v$ follow the reference to plotting accuracy;
- the DynSGS coefficient (last panel of `fields-it5.png`) is at its floor,
  $\approx 1\times10^{-4}$, in the smooth regions, $2$–$4\times10^{-4}$ at the
  compound wave and the contact and $1\times10^{-3}$ at the slow shock, a
  third of its cap $C_{max}h(|u| + c_f) \approx 3\times10^{-3}$.

## Files

| file | contents |
|---|---|
| `user_inputs.jl` | solver parameters (Δt, tend, DynSGS constants, mesh, output) |
| `initialize.jl` | the left/right states; output variable list |
| `user_flux.jl` | 1D ideal-MHD fluxes, `γ_mhd = 2`, `pressure_mhd1d` |
| `user_bc.jl` | Dirichlet ends |
| `user_primitives.jl` | conserved variables as DynSGS primitives; `user_uout!` (ρ, u, v, p, By) |
| `user_analytic.jl` | interpolates `reference_hll.dat` for the overlay at t = 0.1 |
| `user_plot.jl` | the paper's Fig. 2 density figure with zoom boxes (`density-it<n>.png`) |
| `reference_hll.dat` | the reference solution |
| `EQUATIONS.md` | the equations and the wave structure |
