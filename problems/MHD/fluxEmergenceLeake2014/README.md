# 3D visco-resistive GLM-MHD: emergence of a twisted flux tube into a coronal dipole (Leake, Linton & Antiochos 2014)

This case runs the three-dimensional emergence of a buoyant, twisted
convection-zone flux tube through the photosphere into a pre-existing coronal
dipole, the simulations SD, MD, WD and ND of

> J. E. Leake, M. G. Linton, S. K. Antiochos,
> *Simulations of Emerging Magnetic Flux. II. The Formation of Unstable Coronal
> Flux Ropes and the Initiation of Coronal Mass Ejections*,
> ApJ **787**:46 (2014). https://doi.org/10.1088/0004-637X/787/1/46
> (arXiv:1402.2645) — **Paper II**

whose equations, domain, grid, boundaries and initial state are "exactly the
same as used in Paper I",

> J. E. Leake, M. G. Linton, T. Török,
> *Simulations of Emerging Magnetic Flux. I. The Formation of Stable Coronal
> Flux Ropes*, ApJ **778**:99 (2013). https://doi.org/10.1088/0004-637X/778/2/99
> (arXiv:1308.6204) — **Paper I**, Section 2

except that the dipole is reversed so that it reconnects with the top of the
emerging tube. Paper II's MD erupts: its coronal flux rope rises to the
damping layer at $z = 180L_0$ between $t \approx 120$ and $160\,t_0$ (its
Fig. 11); ND (no dipole) does not erupt.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2014
mpiexec -n 128 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2014
```

(use the `mpiexec` that matches the MPI your `MPI.jl` is built against, see
`INSTALL.md`), or interactively `Jexpresso.run_case("MHD", "fluxEmergenceLeake2014")`.
Output goes to `./output/MHD/fluxEmergenceLeake2014/<run>/` every $5t_0$ up
to $t = 160t_0$.

The simulation is chosen in `user_inputs.jl`:

| key | values | meaning |
|---|---|---|
| `:lk_dipole` | `"SD"`, `"MD"` (default), `"WD"`, `"ND"` | Paper II's dipoles, $\Phi_{dip}/\Phi_{tube} = 0.13,\ 0.1,\ 0.067,\ 0$ |
| `:lk_flux_def` | `"whole"` (default) / `"top"` | $\Phi_{tube}$ of the whole tube or of its top half (Deviations, item 3) |
| `:lk_Bd` | a number | overrides `:lk_dipole` with an explicit signed $B_d$ (negative = Paper II's orientation, positive = Paper I's) |
| `:lk_tube` | `true` / `false` | `false` runs the background alone (atmosphere + dipole), which must stay at rest |
| `:lk_cz_gradient` | `1.0` (default) | $a$ of the convection-zone temperature gradient (below) |
| `:lk_mf` | `1.0` (default) / `1.25` | $p = \rho T/m_f$ (Deviations, item 4) |
| `:lk_damp_xd`, `:lk_damp_zd` | 96, 180 | start of the side and top velocity-damping layers (Paper I used 170 at the top) |
| `:lk_glm_alpha` | 1.0 | Dedner damping $\alpha_p$ of $\psi$, rate $\alpha_pc_h/\Delta h$ (Deviations, item 8) |
| `LK_REF` (environment variable) | 1, 2 (default), 3 | mesh and time step (below); also `LK_TEND`, `LK_OUT_DT` and the restart (Restart) |

## Units

Both papers normalize with $L_0 = 170$ km, $B_0$ and $g_0 = g_\odot = 274$
m s$^{-2}$, hence $P_0 = B_0^2/\mu_0$, $\rho_0 = P_0/(L_0g_0)$,
$v_0 = \sqrt{L_0g_0} = 6.82$ km s$^{-1}$, $t_0 = \sqrt{L_0/g_0} = 24.9$ s and
$T_0 = m_pL_0g_0/k_B = 5640$ K. In these units $g = 1$, the magnetic pressure
is $\tfrac12|\mathbf B|^2$ (Heaviside–Lorentz form, as in every MHD case of
this tree) and the gas law is $p = \rho T$; $\gamma = 5/3$ (Deviations, item 4).
(Paper I prints $B_0 = 0.12$ T; Paper II, and Paper I's own dimensional
$B_d$, give 0.13 T. The case is dimensionless and does not depend on it.)

## The problem

Coordinates: $x$ across the tube, $y$ along it, $z$ up.

**Atmosphere** (Paper I Eqs. 11–14), hydrostatic with $\rho(0) = 1$:

$$
T(z) = \begin{cases}
1 - a\,\dfrac{\gamma-1}{\gamma}\,z & z \le 0 \quad\text{(convection zone)}\\
1 & 0 < z < 10 \quad\text{(photosphere/chromosphere)}\\
150^{(z-10)/10} & 10 \le z < 20 \quad\text{(transition region)}\\
150 & z \ge 20 \quad\text{(corona)}
\end{cases}
$$

$dp/dz = -\rho = -p/T$ integrates in closed form in every layer
(`lk_pressure`): a polytrope $p = T^{\gamma/(a(\gamma-1))}$ below
$z = 0$ (exactly adiabatic for $a = 1$), exponentials above. The density
falls by nine and a half decades, from $\rho(-30) = 46.9$ to
$\rho(209.8) = 1.2\times10^{-8}$, and the pressure from 609 to
$6.3\times10^{-6}$ in the corona, as the papers' Fig. 1(a).

**Dipole** (Paper I Eq. 15), invariant along $y$:
$\mathbf B = \nabla\times(A_y\mathbf e_y)$,
$A_y = B_d(z - z_d)/r_1^3$, $r_1^2 = x^2 + (z - z_d)^2$, $z_d = -100$.
Paper II reverses Paper I's dipole ($B_x < 0$ above the tube, antiparallel to
the top of the tube) and specifies its strength by the azimuthal flux per unit
length above and below the separatrix on the $z$ axis (its Eqs. 9–10),

$$
\frac{\Phi_{dip}}{\Phi_{tube}} = \frac{\left|\int_{z_{sep}}^{z_{top}} B_x(0,0,z)\,dz\right|}{\int_{z_{tube}}^{z_{sep}} B_x(0,0,z)\,dz} = 0.13,\ 0.1,\ 0.067 \quad\text{(SD, MD, WD)}.
$$

`lk_solve_dipole` finds $B_d$ from that ratio by bisection (see Deviations).

**Flux tube** (Paper I Eqs. 16–18), along $y$ at $x = 0$, $z_t = -12$:
$B_y = B_te^{-r^2/R^2}$, $B_\theta = qrB_y$, $r^2 = x^2 + (z - z_t)^2$,
$B_t = 5$, $R = 2.5$, $q = 1/R$ (right-handed), so
$B_x = q(z - z_t)B_y$, $B_z = -qxB_y$. Its pressure deficit balances its own
Lorentz force,
$p_1 = \tfrac12B_t^2e^{-2r^2/R^2}\left[q^2(R^2/2 - r^2) - 1\right]$,
and the buoyancy is set by (Paper I Eqs. 19–20)

$$
p = p_0(z) + p_1, \qquad \rho = \rho_0(z)\left[1 + \frac{p_1}{p_0(z)}\,e^{-(y/\lambda)^2}\right],\quad \lambda = 10,
$$

so the middle of the tube ($y = 0$) is isothermal with its surroundings and
rises, and its ends are neutrally buoyant.

**Dissipation** (Paper I Sec. 2.2, Eq. 9): Re = Rm = 100, i.e.
$\nu = \eta = 0.01$ in units of $L_0v_0$, with the resistivity tapered off
towards the side boundaries,
$\eta = 0.01\left[\tanh\!\left(-(r_\eta - 100)/5\right) + 1\right]/2$,
$r_\eta = \sqrt{x^2 + y^2}$.

**Velocity damping** (Paper I Eq. 10): $D\mathbf v/Dt$ gets $-N\mathbf v$
where $|\kappa| > \kappa_d$, for $\kappa = x, y$ ($\kappa_d = 96$) and $z$
($\kappa_d = 180$ in Paper II), with
$N = (|\kappa| - \kappa_d)/(\max|\kappa| - \kappa_d)$ rising linearly to 1 at
the boundary; the three directions add up.

**Boundaries** (Paper I Sec. 2.3): on all six faces all components of the
velocity are zero and the normal gradients of $\mathbf B$, $\rho$ and the
specific internal energy are zero.

**Domain and grid** (Paper I Eqs. 7–8): $x, y \in [-127.6, 127.6]$,
$z \in [-30, 209.8]$, 304 cells per direction, stretched by

$$
x = \pm\left[(1 + f_h)\chi + f_hw\ln\frac{\cosh((\chi - L_h)/w)}{\cosh(-L_h/w)}\right],\ \chi\in[0, 100],
\qquad
z = -30 + (1 + f_v)\chi + f_vw\ln\frac{\cosh((\chi - L_v)/w)}{\cosh(-L_v/w)},\ \chi\in[0, 130],
$$

$f_h = 2.1$, $f_v = 1.83$, $L_h = 95$, $L_v = 100$, $w = 10$: uniform
spacing 0.658 horizontally inside $|x| < 95$ and 0.428 vertically below
$z \approx 70$, three times coarser near the boundaries.

## Meshes and time step

`LK.geo` places the **element** boundaries on the same map, so at `:nop => 4`
each element spans four of the paper's cells and the nodal spacing follows
the paper's stretching:

| `LK_REF` | elements | file | nodes ($N = 4$) | `:Δt` | |
|---|---|---|---|---|---|
| 1 | $16^3$ | `LK_16x16x16.msh` (shipped) | $65^3$ | $3\times10^{-3}$ | smoke test |
| 2 | $38^3$ | `LK_38x38x38.msh` (shipped) | $153^3$ | $1.2\times10^{-3}$ | half the paper's resolution (default) |
| 3 | $76^3$ | `LK_76x76x76.msh` (generate) | $305^3$ | $6\times10^{-4}$ | the paper's $304^3$ |

The $76^3$ mesh is not shipped (about 50 MB); generate it with

```bash
cd problems/MHD/fluxEmergenceLeake2014
gmsh -3 LK.geo -setnumber ref 3 -o LK_76x76x76.msh
```

The time step is set by the fastest wave on the smallest LGL gap. The fast
speed of the initial state peaks just above the transition region, where the
dipole is still strong and the density has dropped by $10^{7.4}$:
$c_h = \max\sqrt{(\gamma p + |\mathbf B|^2)/\rho}$ (printed at start-up,
held constant as the GLM cleaning speed). It grows with the dipole:

| | $B_d$ (solved) | $z_{sep}$ | $c_h$ |
|---|---|---|---|
| SD | $-16006$ | $-6.01$ | 92 |
| MD | $-12351$ | $-5.86$ | 72 |
| WD | $-8314$ | $-5.64$ | 50 |
| ND | 0 | – | 16 |

The smallest LGL gap at $N = 4$ is $0.173$ of the element height at the
bottom, where the elements are $130/n_v$ tall: $\Delta h = 1.40,\ 0.591,\ 0.295$
for `LK_REF` = 1, 2, 3. The step is fixed while the speeds grow (Paper II's
reconnection outflows reach $30v_0$, its Fig. 8), so the `:Δt` of the mesh
table give $(|\mathbf v| + c_f)\Delta t/\Delta h \approx 0.25$ at
$|\mathbf v| + c_f = 120$, the initial fast speed of SD plus the outflows; the
CFL monitor (Notes) prints the actual value as it grows.

## Deviations from the papers, and why

1. **Kinematic viscosity.** Paper I's $\nu = 3.35\times10^3$ kg m$^{-1}$
   s$^{-1}$ is dynamic ($0.01\nu_0$, $\nu_0 = \rho_0L_0v_0$): in the corona,
   where $\rho \sim 10^{-8}$, it would be a kinematic viscosity of $10^6$
   (Re $\sim 10^{-6}$, and an explicit viscous time step $10^7$ times smaller
   than the wave one). The papers' "Re = 100" is evaluated with the
   normalizing density. This case uses the kinematic $\nu = 0.01$
   everywhere, i.e. Re = 100 at every height.
2. **Twist sign.** Paper I's tube is right-handed with $q = 1/R$; Paper II
   prints $q = -1/a$ but states that its configuration is Paper I's with only
   the dipole reversed, and that the dipole's $B_x < 0$ is antiparallel to
   the top of the tube. With $B_y > 0$ only the right-handed tube has
   $B_x > 0$ at its top, so $q = -1/a$ is read as a different sign
   convention for $\theta$ and the case keeps $q = 1/R$.
3. **Dipole strength from Paper II's flux ratios.** Paper II defines its
   simulations by $\Phi_{dip}/\Phi_{tube} = 0.13, 0.1, 0.067$, so $B_d$ is solved
   from them (Paper I's $B_d = [10, 7.5, 5]\times10^3$ belong to its parallel
   dipole). Its Eq. (10) integrates $\Phi_{tube}$ over the top half of the tube,
   which gives $B_d = 7887, 6115, 4138$; but the $\beta(0, 0, z)$ curves of its
   Fig. 1(a), read off the figure at $z = 20$–$64$, need $B_d$ larger by
   2.08, 2.01 and 1.96 at every height (the profiles' shape is matched
   exactly: the atmosphere and the field form are right, only the amplitude
   differs). With $\Phi_{tube}$ the azimuthal flux of the whole tube,
   $\int|B_{x,\rm tube}|\,dz = B_tR$, the same ratios give $B_d = 16006, 12351,
   8314$, within 3 % of the figure. The case uses that (`:lk_flux_def =>
   "whole"`; `"top"` gives the printed one) and prints both ratios at
   start-up. `:lk_Bd` sets $B_d$ directly.
4. **Gas law in normalized units.** Both papers give $\mu_m = 1.25m_p$ and
   $T_0 = m_pL_0g_0/k_B$, i.e. $p = \rho T/1.25$. Their own Fig. 1(a) does not
   follow that: the pressure is $\approx 5\times10^2$ at $z = -30$ and
   $\approx 5\times10^{-6}$ in the corona, the values of $p = \rho T$ (609 and
   $6.3\times10^{-6}$; $p = \rho T/1.25$ gives 2422 and $2.5\times10^{-7}$), and
   Paper I's statement that $T = 1 - (\gamma - 1)/\gamma\cdot z$ is exactly
   adiabatic holds only for $p = \rho T$. Leake et al. (2022)'s Fig. 1 shows
   the same atmosphere in kg m$^{-3}$. The reduced mass sits in the code's
   temperature unit, so `:lk_mf => 1.0` (`initialize.jl`); 1.25 gives the
   literal-text atmosphere, whose corona is 20 times less dense, with Alfvén
   speeds 4.5 times higher and time steps 4.5 times smaller.
5. **Numerics.** The papers use Lare3d (staggered Lagrangian-remap, shock
   viscosity, constrained transport). Here: continuous spectral elements at
   $N = 4$ with GLM divergence cleaning ($c_h$ above, Dedner damping
   $\alpha_p = 1$, item 8), DynSGS shock capturing in conserved form on the
   departure from the reference state (`DSGS_MHD`, legacy sensor, fast-speed
   floors, thermal-energy residual), the explicit $\nu$ and $\eta$ added to
   the DynSGS coefficient element by element (`:mhd_dissipation`), and
   positivity floors $\rho \ge 0.2\rho_e(z)$, $p \ge 0.2p_e(z)$ as the
   Runge–Kutta stage limiter (every hit is counted and printed).
6. **Discrete equilibrium (`:lwell_balanced`).** Neither the stratification
   nor the dipole is an exact equilibrium of the discrete operator, and in the
   corona ($\rho \sim 10^{-8}$, $\beta \sim 10^{-1}$) their interpolation errors
   are large accelerations. Without correction the background alone (no tube,
   $16^3$, measured with the earlier $m_f = 1.25$ atmosphere) reaches
   $|\mathbf v| = 54v_0$ at the transition region by $t = 0.05t_0$ and
   $130v_0$ by $0.1t_0$, and the positivity floors fire from the first step. The deck sets `:lwell_balanced => true`: the kernel
   evaluates the element RHS of the reference state (atmosphere + dipole,
   fluxes and sources) once and subtracts it at every stage, so the background
   is an exact discrete steady state (unchanged to the last bit at
   $t = 0.1t_0$) and only the departure from it (the tube) drives the flow.
7. **Boundaries.** CG imposes the three velocity components and $\psi = 0$
   strongly on all faces; $\rho$, $E$ and $\mathbf B$ on the boundary nodes
   evolve with their own element derivatives, and the natural condition of the
   dissipative terms is a zero normal gradient. This is the CG counterpart of
   the papers' zero-gradient conditions, not a ghost-cell copy.
8. **GLM transient and $\alpha_p$.** Constrained transport keeps Lare3d's
   $\nabla\cdot\mathbf B$ at round-off. Here the tube's nodal field has a
   discrete divergence that GLM carries away at $c_h$ (72 for MD), i.e. from
   the tube to the transition region in $0.4t_0$. There $\beta$ is low
   (0.06 for MD), and in total-energy form the gas
   pressure is the small remainder
   $E - \tfrac12\rho|\mathbf v|^2 - \tfrac12|\mathbf B|^2 - \tfrac12\psi^2$:
   a field error $\delta B$ changes it by about $B\,\delta B$, a sizeable
   fraction of the local $p$ for $\delta B/B \sim 10^{-1}$, and drives
   Alfvénic flows $\sim v_A\,\delta B/B$. The damping rate $\alpha_pc_h/\Delta h$ attenuates a
   $\psi$ front by $e^{-\alpha_p L/\Delta h}$ over a distance $L$; this case
   uses $\alpha_p = 1$, the top of the range of Dedner et al. (2002) and
   Mignone & Tzeferacos (2010), rather than the 0.2 of
   `fluxEmergenceSon2025` (Verification). The discrete divergence of the tube
   at $t = 0$, estimated elementwise on the three meshes, is 0.82, 0.082 and
   0.018 for `LK_REF` = 1, 2, 3 (the tube has $R = 2.5$, and the $16^3$
   elements are 12.5 wide), so the transient falls fast with resolution.
9. **The papers' dipole is not current-free.** Paper II's Eq. 8 (Paper I's
   Eq. 15) prints the dipole as "translationally invariant along y",
   $A_y = B_d(z - z_d)/r_1^3$ with $r_1^2 = x^2 + (z - z_d)^2$, and this case
   uses it as printed. A field invariant along $y$ is potential only if $A_y$
   is harmonic in $(x, z)$, and $s/r_1^3$ ($s = z - z_d$) is not:
   $\nabla\times\mathbf B = -\nabla^2A_y\,\hat{\mathbf y} = -3B_ds/r_1^5\,\hat{\mathbf y}$,
   i.e. $|\mathbf j| = 1.5|\mathbf B|/r_1$ on the axis (checked by finite
   differences of the coded field). For MD at the base of the corona
   $(0, 0, 20)$ the Lorentz force $|\mathbf j\times\mathbf B| = 2.6\times10^{-6}$
   is 60 times the gravity $\rho g = 4.2\times10^{-8}$: the background is not
   an equilibrium of the continuous equations either, and only
   `:lwell_balanced` (item 6) holds it at rest, by subtracting that force with
   the rest of the background's discrete RHS at every stage; a departure from
   the background (the emerging tube, the GLM transient) does not feel the
   subtraction and moves in a corona that, physically, would not be at rest.
   The potential field invariant along $y$ is the line dipole
   $A_y \propto s/r_1^2$ (used by `fluxEmergenceLeake2022_2D`); it falls off
   as $1/s^2$ rather than $1/s^3$, so it does not reproduce the shape of the
   $\beta(0, 0, z)$ curves of Paper II's Fig. 1, which this field matches.

## Output

`ρ, u, v, w, p, Bx, By, Bz, ψ` and the derived `T` ($= m_fp/\rho$, in $T_0$),
`β` ($= 2p/|\mathbf B|^2$), `Bmag` and `log10ρ`, as `.pvtu` per snapshot.
Mass and total energy are logged every 10 steps (`:conservation_slots`).

**Disk.** A snapshot holds about 240 bytes per node: 65 MB at `LK_REF = 1`,
0.85 GB at 2 and 6.7 GB at 3, so the default 33 snapshots (every $5t_0$)
take 28 GB at `LK_REF = 2`. Paper II's rise speed (Fig. 12) and Fig. 9
($t = 122.5t_0$) need a denser cadence late in the run, for example

```julia
:diagnostics_at_times => vcat(0.0:5.0:100.0, 101.0:1.0:165.0, [122.5]),
```

(any list of times works; 87 snapshots, 74 GB at `LK_REF = 2`).

## Restart

Any snapshot of this case restarts the run: it holds the nine primitives in
Float64 at every node, and `initialize.jl` rebuilds the conserved state from
them exactly ($E = p/(\gamma-1) + \tfrac12\rho|\mathbf v|^2 +
\tfrac12|\mathbf B|^2 + \tfrac12\psi^2$), keeps the analytic reference state
and the $c_h$ of $t = 0$, and starts the clock at the snapshot's time. Each
rank reads its own piece, so **use the rank count of the run that wrote it**.
The DynSGS residual history starts afresh, so the continuation is close to,
not bitwise equal to, the uninterrupted run. Environment variables, so the
deck need not be edited:

| variable | meaning |
|---|---|
| `LK_RESTART_DIR` | output directory of the earlier run (`output/MHD/fluxEmergenceLeake2014/output-<date>`) |
| `LK_RESTART_TIME` | snapshot to start from (nearest in its `simulation.pvd`; default the last) |
| `LK_TEND` | `:tend` |
| `LK_OUT_DT` | a snapshot every `LK_OUT_DT` (default 5) |
| `LK_REF` | mesh (must be the one of the snapshot) |
| `LK_DIAG` | the stage diagnostics (Notes) every `LK_DIAG` stages |

For example, to rerun the last stretch before a blow-up at $t = 5.87$ with a
snapshot every $0.05t_0$ and the diagnostics every 5 steps:

```bash
LK_REF=1 LK_RESTART_DIR=output/MHD/fluxEmergenceLeake2014/output-<date> LK_RESTART_TIME=5 \
LK_TEND=6 LK_OUT_DT=0.05 LK_DIAG=25 mpiexec -n 14 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2014
```

The new run writes to a new output directory, starting with the restart
state itself.

## Figures

`tools/plot_leake2014.py` draws Paper II's figures from the output
(`numpy` and `matplotlib`; the 3D field-line figures also need
`pip install pyvista`, which renders off-screen, so it also runs on a cluster
node). `RUN` is an output directory `output/MHD/fluxEmergenceLeake2014/output-<date>`;
a restarted run is `DIR1,DIR2`. Every command takes `--outdir` (default
`figs_leake2014`), `--times` (nearest snapshots) and `-h`.

```bash
T=tools/plot_leake2014.py
python3 $T ic MD=$MD                                              # Fig. 1(a), with the nodal values of MD at t = 0
python3 $T render $MD --preset ic        --times 0                # Fig. 1(b)
python3 $T render $MD --preset emergence --times 5,40,60          # Fig. 2(a-c)
python3 $T render $MD --preset emergence --times 5,40,60 --view top   # Fig. 2(d-f)
python3 $T render $MD --preset arcade    --times 100,110,120      # Fig. 3
python3 $T render $SD --preset dipoles   --times 100              # Fig. 4 (and $MD, $WD)
python3 $T render $MD --preset rotation  --times 70,85,100        # Fig. 5(a-c)
python3 $T slice  $MD --plane z=0 --var Bz --clim -0.2,0.2 --cmap jet --vectors \
                      --xlim 0,30 --ylim -20,20 --times 70,85,100 # Fig. 5(d-f): zoom on the positive polarity
python3 $T render $MD --preset fluxrope  --times 100,110,120      # Fig. 6
python3 $T slice  $MD --plane y=0 --var jy --fieldlines 1.5 --nulls --line-color k \
                      --xlim -60,60 --zlim -10,80 --times 100,110,120   # the (Bx, Bz) lines of Figs. 6-7, O/X-points
python3 $T render $MD --preset current   --times 100,110,120      # Fig. 7(a-c)
python3 $T render $MD --preset current   --times 100,110,120 --view top --bz-cmap gray  # Fig. 7(d-f)
python3 $T slice  $MD --plane y=0 --var j --cmap gist_heat --clim 0,0.01 \
                      --xlim -40,40 --zlim 15,85 --times 140,150,160    # Fig. 8(a-c)
python3 $T line   $MD --var w --zlim 15,45 --times 140,150,160    # Fig. 8(d-f)
python3 $T render $MD --preset eruption  --times 110,120,130,140,150,160  # Fig. 10
python3 $T height SD=$SD MD=$MD WD=$WD ND=$ND                     # Figs. 11 and 12 (+ one CSV per run)
python3 $T height MD=$MD MD2=$MD2 MD3=$MD3 MD4=$MD4 --mono        # Fig. 13
python3 $T axflux SD=$SD MD=$MD WD=$WD                            # Fig. 14(b) (+ CSVs)
```

How the quantities are computed:

- **Current density.** $\mathbf j = \nabla\times\mathbf B$ in $j_0 = B_0/(\mu_0L_0)$,
  second-order differences on the node lattice (the output is in the papers'
  units; `LK.geo` makes the nodes a tensor-product lattice, which the script
  rebuilds from the rank pieces).
- **O-point** (Figs. 6, 10, 11). The nulls of $(B_x, B_z)$ in $y = 0$ are found
  cell by cell (bilinear), classified by the sign of the Jacobian determinant
  (O-point $> 0$, X-point $< 0$); O–X pairs closer than `--pair-dist`
  ($3L_0$) are grid noise and dropped; the flux-rope axis is the highest
  O-point with $|x| \le$ `--xwin` (30) and $z \le 180$, as in Paper II
  ("the higher one is used"). The position is accurate to the local node
  spacing (at $t = 0$ it sits at $z = -11.6$ on the $16^3$ mesh, $-12$
  exactly). `height` writes every time, height and the O/X counts to a CSV:
  check the counts, and the `slice --nulls` plot, where the curve jumps.
- **Rise speed** (Fig. 12): $\partial z/\partial t$ of the height by centred
  differences, in `--speed-range` (120–165); `--smooth N` averages over N
  snapshots first.
- **Axial flux** (Fig. 14b): $\int\!\!\int_{z > z_h} B_y\,dx\,dz \,/
  \int\!\!\int B_y\,dx\,dz$ in $y = 0$ (trapezoids), for `--heights` 0 and 10.
- **3D figures.** Field lines are traced (VTK) from the papers' seeds: the
  bottom boundary for the gray dipole lines, circles around the tube axis on
  the side boundaries $y = \pm127.6$ (radii 0.5, 1, 2 in Fig. 4), the axis
  itself (black/purple), the O-point in $y = 0$ (orange) and a disk of radius
  10 around it (blue, Fig. 10). The $B_z$ plane at $z = 0$ uses the papers'
  rainbow and range (±0.5 Fig. 2, ±0.2 Figs. 5–7 and 10, ±7.5 G in Fig. 1b);
  `--box`, `--azimuth`, `--elevation`, `--zoom`, `--dipole-seeds`,
  `--tube-radii`, `--bz-clim` override a preset. Paper II's seed points and
  camera are not published, so the presets match its figures by eye. Fig. 9
  (hand-picked seeds) has no preset.

## Verification here

Short runs on 4 threads; the full runs are for a cluster.

- **Dipole strength.** MD solves to $B_d = -12351$, $z_{sep} = -5.86$ ($\Phi_{dip}/\Phi_{tube} =
  0.1000$ with the whole tube, 0.207 with its top half, both printed at
  start-up); SD and WD as in the table above. With them the $\beta$ curves
  of `tools/plot_leake2014.py ic` have the minima 0.032, 0.054, 0.12 above
  $z = 0$, those of the paper's Fig. 1(a).
- **Atmosphere.** The nodal $\rho$, $p$ take the closed-form values
  (checked with $m_f = 1.25$: 233 to $4.28\times10^{-10}$ and 2422 to
  $5.14\times10^{-8}$; with $m_f = 1$: 46.9 to $1.18\times10^{-8}$ and 609 to $1.76\times10^{-6}$).
- **Background at rest** (`:lk_tube => false`, `LK_REF = 1`): with
  `:lwell_balanced` every field at $t = 0.1t_0$ equals the initial one to the
  last bit, and no floor fires. Without it, see Deviations, item 6.
- **Full problem** (MD, `LK_REF = 1`): runs to $t = 0.4t_0$ with no NaN. The
  $16^3$ mesh does not resolve the tube (12.5-wide elements for $R = 2.5$),
  so this is a test of the pipeline, not of the physics. The GLM transient of
  item 8 at $t = 0.2t_0$:

  | setup | $\alpha_p$ | max $\lvert\mathbf v\rvert$ | max $p/p(t=0)$ | max $\lvert\psi\rvert$ at $z = 12$–$19$ |
  |---|---|---|---|---|
  | $m_f = 1.25$, top-half $B_d$, $2\times10^{-3}$ | 0.2 | 18.3 | 66 | $1.8\times10^{-3}$ |
  | $m_f = 1.25$, top-half $B_d$, $2\times10^{-3}$ | 1 | 6.1 | 6.2 | $6.1\times10^{-5}$ |
  | $m_f = 1$, whole-tube $B_d$, $3\times10^{-3}$ (default; $t = 0.21$) | 1 | 0.83 | 1.10 | $1.5\times10^{-4}$ |

  With $m_f = 1.25$ the positivity floors hit $\sim10^2$ nodes per stage
  from $t \approx 0.02t_0$ ($10^5$ hits by $0.1t_0$); with the defaults they
  start at $t = 0.115t_0$ ($10^4$ hits by 0.15), and the CFL monitor stays at
  its initial 0.153.

- **Restart** (`LK_REF = 1`, Restart). A run with snapshots at $t = 0.03$
  and $0.06$ restarted from its $t = 0.03$ snapshot, on 1 rank and on 2 MPI
  ranks: the restart state equals the snapshot to $5\times10^{-16}$ of each
  field's maximum ($T$ $7\times10^{-15}$), piece by piece. Ten steps later
  (1 rank) $\rho$, $p$ and $T$ agree with the uninterrupted run to
  $2\times10^{-4}$ and $\mathbf v$, $\mathbf B$ to 0.7–5 % of their maxima
  ($\psi$ 9 %): the DynSGS residual history restarts from zero.

Nothing was run here at `LK_REF` = 2 or 3.

**Cost.** `LK_REF = 1` takes about 6 s per time step on 4 threads
($65^3$ nodes). Scaling by nodes and steps, `LK_REF = 2` ($153^3$ nodes,
$1.3\times10^5$ steps to $160t_0$) is $\approx 1.2\times10^4$ core-hours and
`LK_REF = 3` ($305^3$, $2.7\times10^5$ steps) $\approx 1.8\times10^5$, before
parallel losses. In Paper II the O-point of MD jumps at $t \approx 115t_0$
and reaches the damping layer before $160t_0$; a shorter `:tend` covers the
formation of the flux rope.

## Notes

- **Partition.** The deck sets `:lxy_partition => false`: a balanced 3D
  partition (METIS, through GridapGmsh's distributed reader, which reads the
  mesh file on every rank). Jexpresso's default (`true`) bins the $x$–$y$
  bounding box into a uniform $n_x\times n_y$ grid of columns, one per rank;
  on these horizontally stretched meshes the outer elements are several bins
  wide, so from a few hundred ranks on some ranks get no element and the run
  stops in `mod_mesh_read_gmsh!` with `reducing over an empty collection`
  (at $38^3$: 32 of 448 and 116 of 896 ranks empty), and at any count the
  central ranks carry 1.3–1.9 times the average load.
- **Progress.** Snapshots are $5t_0$ apart (1667 steps at `LK_REF = 1`), and
  the floor and CFL lines below are printed only when their totals pass a
  power of ten or grow by 20 %, so a healthy run can print nothing for tens of
  minutes. The deck sets `:lstep_heartbeat => true` (`#   step N   t = …`
  every 100 steps), and `conservation.dat` in the run's output directory gets
  $t$, the step, the total mass and energy every 10 steps
  (`tail -f output/MHD/fluxEmergenceLeake2014/<run>/conservation.dat`).
  `LK_REF = 1` takes about 9 s per step on 4 threads.
- **The positivity-floor line.** `# positivity floors: N node hits so far on
  all ranks (n in this stage) at t = T` is printed by rank 0 each time the
  total N passes a power of ten. T is the Runge–Kutta stage time
  ($t_n + c_i\Delta t$, in $t_0$), not an output time: snapshots are written at
  `:diagnostics_at_times` only. A node on a partition interface counts once
  per rank that holds it.
- **Stage diagnostics** (`LK_DIAG=n` in the environment, off by default):
  every `n`-th Runge–Kutta stage rank 0 prints `#DIAG t=…` with the largest
  $|\mathbf v|$ and $|\mathbf v| + c_f$, the smallest $\rho/\rho_e$, the
  largest and smallest $p/p_e$, $|\psi|$, $T$ and $|\delta\mathbf B|$, the
  number of nodes at the density and pressure floors and the heights they span,
  the largest DynSGS + explicit coefficient of the $\rho$, momentum and field
  slots, and the largest explicit-diffusion number $\mu\Delta t/\Delta h^2$
  (at $N = 4$ the explicit Laplacian goes unstable from $\approx 0.3$–$0.5$ at
  this integrator), then a line per extreme saying where it is. `LK_DIAG=25`
  is every 5 steps.
- **CFL monitor.** `# CFL: max (|v| + c_f) Δt/Δh = …` is printed by rank 0
  every time the largest value over the nodes ($c_f = \sqrt{(\gamma p +
  |\mathbf B|^2)/\rho}$, $\Delta h$ the smallest LGL gap of the elements around
  the node, checked every 100 Runge–Kutta stages) has grown by 20 %, as a
  `WARNING` above 0.5: halve `:Δt` and restart from the last snapshot if it
  keeps growing.
- **The shared CFL printout.** Its advective and acoustic lines
  (`src/kernel/physics/soundSpeed.jl`) assume the CompEuler layout and report
  nonsense for this case (max|u| $\sim 10^4$, $c = 0$). The viscous line is
  the DynSGS coefficient actually applied.
- **Constants shared by name** (`γ_mhd`, `c_h_mhd`, `pressure_mhd`) with the
  other MHD cases: `user_flux.jl` refuses to run if a case with another
  $\gamma$ was loaded earlier in the same Julia session. Restart Julia.

## Files

| file | content |
|---|---|
| `user_inputs.jl` | `LK_REF`, Δt, the problem keys, DynSGS, positivity floors and CFL monitor, $\nu$/$\eta$ |
| `initialize.jl` | atmosphere, dipole ($B_d$ from the flux ratio), tube, $c_h$, Δh |
| `user_flux.jl` | 3D ideal GLM-MHD fluxes, $\gamma = 5/3$ |
| `user_source.jl` | gravity, GLM damping, velocity damping |
| `user_bc.jl` | $\rho\mathbf v = 0$, $\psi = 0$ on all faces |
| `user_primitives.jl` | DynSGS acts on $q - q_e$; derived output fields |
| `LK.geo` | the stretched hexahedral mesh (`ref` = 1, 2, 3) |
