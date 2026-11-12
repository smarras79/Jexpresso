# 2.5D visco-resistive GLM-MHD: flux emergence into a line dipole (the cross-section of Leake et al. 2022)

This case is the plane perpendicular to the flux rope of
[`fluxEmergenceLeake2022`](../fluxEmergenceLeake2022/README.md),

> J. E. Leake, M. G. Linton, S. K. Antiochos, *The role of reconnection in the
> onset of solar eruptions*, ApJ **934**:10 (2022),

with everything invariant along the rope: the same atmosphere, the same rope
cross-section, the same explicit viscosity and resistivity, the same numerics
(DynSGS, positivity floors, GLM cleaning, well balance, velocity damping), on
the paper's stretched grid in the plane. It costs a few percent of the 3D run
at the same resolution, so it is a testbed for what the 3D case does in the
transition region, at the floors and in the damping layers, and a 2.5D study
of a twisted rope emerging into an overlying field. It is not a reproduction
of the paper, whose results are 3D (Deviations).

## Run

```bash
julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022_2D
mpiexec -n 14 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022_2D
```

Output goes to `./output/MHD/fluxEmergenceLeake2022_2D/<run>/`, a snapshot
every $5t_0$ to $t = 400t_0$ (plus 142 and 184, the times of the 3D paper's
Figs. 3–4). The keys of `user_inputs.jl`:

| key | values | meaning |
|---|---|---|
| `:lk2d_orientation` | `"antiparallel"` (default), `"parallel"` | the rope's twist above its axis against the dipole field: the 3D paper's θ = π (simulation 4, which erupts in 3D) or θ = 0 (stable in 3D) |
| `:lk2d_D`, `:lk2d_hd` | 150, −100 | line-dipole strength and depth (below) |
| `:lk2d_buoyancy` | 1.0 | $b$ of the rope's density deficit (below) |
| `:lk2d_dipole`, `:lk2d_rope` | `true` | `false` removes the coronal field or the rope (the background alone must stay at rest) |
| `:mhd_dissipation` | `lk2d_dissipation` | the paper's $\nu = \eta = 0.01$; delete the line for ideal MHD + DynSGS |
| `:lk2d_damp_x_frac`, `:lk2d_damp_z_frac` | 0.75, 0.875 | start of the side and top velocity damping (the paper's Eqs. 21–22) |
| `:lk2d_mf`, `:lk2d_cz_gradient`, `:lk2d_glm_alpha` | 1, 1, 1 | as the 3D case |

and the environment variables `LK2D_REF` (mesh, below), `LK2D_TEND`,
`LK2D_OUT_DT` (snapshot interval), `LK2D_DIAG` (Notes) and the restart's
(Restart).

## Coordinates

The code's 2D plane is $(x, y)$ with $y$ **up**: it is the paper's $(x, z)$.
The rope axis, the paper's $y$, is the invariant direction, the code's third
component. In the output, `u`, `Bx` are horizontal, `v`, `By` vertical and
`w`, `Bz` along the rope, so `Bz` is the rope's axial field.

## The problem

Units of the 3D case: $L_0 = 170$ km, $B_0 = 0.13$ T, $g_0 = 274$ m s$^{-2}$,
$v_0 = 6.82$ km s$^{-1}$, $t_0 = 24.9$ s; $g = 1$, magnetic pressure
$\tfrac12|\mathbf B|^2$, $p = \rho T$ ($m_f = 1$, the 3D case's README),
$\gamma = 5/3$. With $h$ the height (the code's $y$):

**Atmosphere** (the paper's Eqs. 23–25), hydrostatic with $\rho(0) = 1$:
$T = 1 - \tfrac{\gamma-1}{\gamma}h$ below $h = 0$, 1 to $h = 10$,
$150^{(h-10)/10}$ to $h = 20$, 150 above; closed-form $p$ in every layer.

**Coronal field.** A line dipole at $(0, h_d)$ below the domain, the 2D
potential field

$$
A = -D\,\frac{s}{r^2},\qquad \mathbf B = \left(\frac{\partial A}{\partial y}, -\frac{\partial A}{\partial x}\right) = D\left(\frac{s^2 - x^2}{r^4},\ -\frac{2sx}{r^4}\right),\qquad s = h - h_d,\ r^2 = x^2 + s^2,
$$

pointing along $+x$ above the dipole, as the paper's point dipole (moment
$-\hat{\mathbf x}$) does above the rope. $h_d = -100$ and $D = 150$ give
$|\mathbf B| = D/h_d^2 = 0.015B_0 = 19.5$ G at $(0, 0)$, the 3D dipole's value
there.

**Flux rope** (the paper's Eqs. 28–34) along the invariant direction at
$(0, -12)$, width $w = 2.5$, strength $B_t = 5$, right-handed twist:

$$
B_z = \sigma B_te^{-r^2/w^2},\qquad (B_x, B_y) = \frac{\sigma B_te^{-r^2/w^2}}{w}\,(\zeta, -x),\qquad r^2 = x^2 + \zeta^2,\ \zeta = h + 12,
$$

with $\sigma = +1$ for `"parallel"` and $-1$ for `"antiparallel"`. The
pressure deficit $p_1 = \tfrac12B_t^2e^{-2r^2/w^2}\left[(w^2/2 - r^2)/w^2 - 1\right]$
balances the rope's Lorentz force, and $p = p_0(h) + p_1$,
$\rho = \rho_0(h)\left[1 + b\,p_1/p_0(h)\right]$: $b = 1$ is the middle of the
3D rope (isothermal with its surroundings, $y_r = 0$), $b = 0$ its neutrally
buoyant ends.

**Dissipation.** $\nu = \eta = 0.01$ (Re = Rm = 100), kinematic, added to the
DynSGS coefficient of the momentum and field (`:mhd_dissipation`). This is
the paper's physics, not a stabilization: the eruption is reconnection between
the rope and the overlying field, whose rate depends on $\eta$. DynSGS keeps
the discretization stable, as in every MHD case; without
`:mhd_dissipation` the case is ideal MHD + DynSGS.

**Velocity damping** (the paper's Eqs. 21–22): $-N\mathbf v$ for
$|x| > 3L_x/8$ and $h > \tfrac78h_{\rm max}$, $N$ rising linearly to 1 at the
boundary.

**Boundaries:** all three velocity components and $\psi$ are zero on the four
sides; $\rho$, $E$, $\mathbf B$ are free (the natural condition of the
dissipative terms is a zero normal gradient).

## Deviations from the 3D case, and why

1. **Line dipole.** The paper's point dipole falls off as $1/r^3$ in every
   direction and has no counterpart invariant along the rope. The line dipole
   is the potential field that has; it falls off as $1/s^2$, so above the
   photosphere it is stronger than the 3D field with the same value at
   $(0, 0)$ (on the axis $x = 0$): 0.0104 against 0.0087 at $h = 20$, 0.0067
   against 0.0044 at $h = 50$, 0.0038 against 0.0019 at $h = 100$. (The 2014 papers' 2D formula $A \propto s/r^3$ is not potential
   in 2D: it carries a current $|\mathbf j| \approx 1.5|\mathbf B|/r$, see the
   2014 case's README.)
2. **Two orientations.** Only the twist above the axis parallel or
   antiparallel to the dipole field, the paper's θ = 0 and π, exist in 2.5D;
   its intermediate angles need the third dimension.
3. **Buoyant along the whole rope.** In 3D only the middle of the rope
   ($|y_r| \lesssim \lambda = 10$) is buoyant, its ends anchor it, and it
   emerges as an Ω loop that drains mass down its legs. The 2.5D rope is the
   3D one's middle cross-section everywhere, with no legs to drain into; $b$
   lowers its buoyancy uniformly.
4. **Reconnection in the plane.** The rope's axial field is a guide field:
   2.5D reconnection only involves the twist field and the dipole.
5. **Reference-weighted DynSGS** (`:dsgs_ref_weight`, default on;
   `LK2D_REF_WEIGHT=0` for the 3D decks' default). The conserved-form operator
   diffuses $\rho$, $\rho\mathbf v$, $E$ as the relative departure,
   $\nabla\cdot(\mu\rho_e\nabla((q - q_e)/\rho_e))$, instead of
   $\nabla\cdot(\mu\nabla(q - q_e))$ (the fields keep the absolute form), as
   `fluxEmergenceSon2025DSGS`. Both vanish at rest, but across the
   $10^3$–$10^7$ density drop of the chromosphere and transition region only
   the relative form obeys a maximum principle: with the absolute one a
   departure of the dense side, diffused upward, exceeds the whole state of
   the light side, and the transition region above the rope is evacuated and
   heated to $10^3$–$10^4$ within $10t_0$ at `LK2D_REF = 1`. With the
   elements half as tall (`LK2D_REF = 4`) the absolute form lasts to
   $t \approx 13$, then evacuates the chromosphere above the rope; the
   relative form, from the same state, touches the floors at up to 20 nodes
   and recovers (Verification).

## Meshes and time step

`LK22_2D.geo` is the $x$–$z$ plane of the 3D `LK22.geo` (the paper's
stretched cells, $x \in [-277.2, 277.2]$, $h \in [-30, 530.5]$), with the
element boundaries every $k$-th cell ($k_x$ in $x$, $k_z$ in $z$ for 4):

| `LK2D_REF` | elements | $k$ | file | nodes ($N = 4$) | `:Δt` |
|---|---|---|---|---|---|
| 1 (default) | $32\times48$ | 16 | `LK22_2D_32x48.msh` | 24,897 | $2\times10^{-3}$ |
| 2 | $64\times96$ | 8 | `LK22_2D_64x96.msh` | 98,945 | $10^{-3}$ |
| 3 | $128\times192$ | 4 | `LK22_2D_128x192.msh` | 394,497 | $5\times10^{-4}$ |
| 4 | $32\times96$ | 16, 8 | `LK22_2D_32x96.msh` | 49,665 | $10^{-3}$ |

`LK2D_REF = 1` is the cross-section of the 3D `LK22_REF = 2` mesh, 3 the
paper's resolution, 4 the `LK2D_REF = 1` columns with the elements of 2 in
the vertical (Notes, "Vertical resolution"). All four ship (0.1, 0.4, 1.8, 0.2 MB); regenerate with
`gmsh -2 LK22_2D.geo -setnumber ref N -o <file>`. The smallest LGL gap is
0.173 of the bottom elements' height (4, 2, 1, 2 $L_0$), and `:Δt` gives
$(|\mathbf v| + c_f)\Delta t/\Delta h \approx 0.25$ at $|\mathbf v| + c_f = 80$,
as the 3D deck.

## Restart

As the 3D Leake cases: any snapshot restarts the run (same mesh and rank
count; `initialize.jl` rebuilds the conserved state exactly from the nine
primitives and keeps the $c_h$ of $t = 0$; the DynSGS history restarts, so the
continuation is close to, not bitwise equal to, an uninterrupted run).

```bash
LK2D_RESTART_DIR=output/MHD/fluxEmergenceLeake2022_2D/output-<date> LK2D_RESTART_TIME=100 LK2D_TEND=200 \
    mpiexec -n 14 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022_2D
```

`LK2D_RESTART_TIME` picks the nearest snapshot (default the last); the new run
writes to a new output directory, starting with that snapshot.

## Output and figures

`ρ, u, v, w, p, Bx, By, Bz, ψ, T, β, Bmag` and `log10_ρ`, `log10_p`,
`log10_β` per snapshot (`.pvtu`), and the DynSGS coefficients actually applied
(`mu_dsgs_*`); mass and total energy in `conservation.dat` every 10 steps.
ParaView reads `simulation.pvd`. `tools/plot_fluxemergence_vtu.py` draws a
field with the field lines in the plane (isolines of the flux function of
$(B_x, B_y)$) and, with `--vectors`, the velocity; `--xlim`/`--ylim` zoom on
the emergence region:

```bash
RUN=output/MHD/fluxEmergenceLeake2022_2D/output-<date>
P="tools/plot_fluxemergence_vtu.py --xlim -60 60 --ylim -30 90 --xlabel x/L0 --ylabel z/L0 --time-unit t0 --outdir figs2d"
python3 $P $RUN                                   # log10 rho and field lines, every snapshot
python3 $P $RUN --var Bz --symmetric --steps 1,11,21    # the axial field (iter_N numbers)
python3 $P $RUN --var T --cmap inferno --vectors --vec-ref 10
```

## Notes

- **Progress.** The deck sets `:lstep_heartbeat => true` (`#   step N   t = …`
  every 100 steps); `conservation.dat` in the output directory gets $t$, the
  step, mass and energy every 10 steps.
- **Positivity floors and CFL monitor**: as the 3D case, `# positivity
  floors: N node hits …` each time the total passes a power of ten, `# CFL: max
  (|v| + c_f) Δt/Δh = …` each time the largest value has grown by 20 % (a
  `WARNING` above 0.5).
- **Stage diagnostics** (`LK2D_DIAG=n`, off by default): every `n`-th
  Runge–Kutta stage, `#DIAG t=…` with the largest $|\mathbf v|$ and
  $|\mathbf v| + c_f$, the smallest $\rho/\rho_e$, the largest and smallest
  $p/p_e$, $T$, the nodes at the floors and the heights they span, the largest
  DynSGS + explicit coefficient of the $\rho$, momentum and field slots and the
  largest explicit-diffusion number $\mu\Delta t/\Delta h^2$, then where each
  extreme is. `LK2D_DIAG=25` is every 5 steps.
- **Vertical resolution.** At $N = 4$ an element couples all its nodes: a
  disturbance at its bottom node also forces its top node, which, with $\rho$
  falling by a factor $C$ across the element, gets of order $C/10$ times the
  bottom node's velocity, and the next element does the same (the 2014 case's
  README, "Vertical resolution", where it is measured in 3D). Through the
  chromosphere ($T = 1$, scale height 1) the 4 $L_0$ elements of
  `LK2D_REF = 1` have $C = 55$ (42 at the base of the transition region), the
  2 $L_0$ elements of `LK2D_REF = 4` $C = 7.4$ (9.6): about the 3D 2014
  case's `LK_REF = 2` (3.4 $L_0$, $C$ up to 41) and `LK_REF = 5` (2 $L_0$,
  9.5) meshes. The two runs share the horizontal mesh (Verification).
- **Partition.** `:lxy_partition => false` (balanced METIS partition), as the
  3D Leake cases.
- **The shared CFL printout** (`Advective CFL`, `max|u|`) assumes the
  CompEuler layout and is wrong for MHD; the CFL monitor above is the one to
  read.

## Verification here

On 4 threads at `LK2D_REF = 1`, and 4 where stated; nothing was run here at
2 or 3.

- **Background at rest** (`:lk2d_rope => false`): with `:lwell_balanced`
  every field at $t = 0.1t_0$ (50 steps) is bitwise equal to the initial
  one, and no floor fires.
- **Explicit dissipation.** The coefficient of the momentum and field slots
  is $0.01$ from the first step, the $\rho$ slot's 0 (the DynSGS residual is
  zero there): the 2D `:mhd_dissipation` hook adds $\nu$, $\eta$ as in 3D.
- **Startup.** $c_h = 53.4$, CFL monitor 0.155; mass constant to the printed
  digits over the first 10 steps.
- **Restart** from the $t = 5$ snapshot: the restart state equals the
  snapshot to $4\times10^{-16}$ of each field's maximum ($T$
  $4\times10^{-15}$).
- **The full problem (antiparallel), with the absolute DynSGS form
  (`LK2D_REF_WEIGHT=0`), to $t = 20t_0$** shows the same transition-region
  failure as the 3D runs, at 1/30 of their cost:
  - from $t = 2.3$ the pressure floor holds one to three nodes of the
    photosphere above the rope, $(0, 2)$, where the 1 $L_0$ scale height
    has node gaps of 0.69 and 1.31 at this resolution;
  - from $t \approx 5.2$ the chromosphere and transition region above the
    rope are evacuated: 50–75 nodes at the density floor at $h = 8$–20,
    $T$ up to $7\times10^3$ (150 in the corona), $|\mathbf v|$ up to 38,
    the DynSGS coefficient up to 69 and the explicit-diffusion number to
    0.29, the CFL monitor to 0.69 at $t = 8.75$ (`WARNING`);
  - by $t = 10.4$ it relaxes (CFL 0.4, diffusion number 0.04), with the
    density floor still holding 74 nodes at $h = 15$–20;
  - a second burst at $t = 13$–17 ($|\mathbf v|$ to 290, $T$ to $10^4$,
    CFL to 1.28) leaves a hot dome ($T = 10^3$–$2\times10^3$ for a coronal
    150) over $|x| < 70$, $h = 25$–100 at $t = 20$, with 205 and 193 nodes at
    the density and pressure floors.

  The run is about 0.3 s per step ($10^4$ steps to $t = 20$ in 50 minutes).
- **The same with the reference-weighted operator** (`:dsgs_ref_weight`, the
  default) **to $t = 10t_0$**: no transition-region evacuation. The one
  photospheric node above the rope touches the pressure floor at
  $t = 2.4$–2.8 and 5.2–5.4 only; at $t = 6$ (8, 10) $\rho/\rho_e \ge 0.81$
  (0.66, 0.78), $T \le 166$ (166, 167), $|\mathbf v| \le 1.7$ (1.5, 1.1), the
  DynSGS coefficient $\le 1$ (1.3, 2.1), the explicit-diffusion number
  0.004, and the CFL monitor never passes its initial 0.155 (the absolute
  form at $t = 6$: 43 nodes at the density floor, $T$ 665, $|\mathbf v|$ 6.3,
  coefficient 46).
- **The absolute form at `LK2D_REF = 4`** (the same columns, elements 2 tall
  instead of 4; Notes, "Vertical resolution") **to $t = 10t_0$**: no floor
  fires. Over $0 \le t \le 10$, $\rho/\rho_e \ge 0.58$, $p/p_e \ge 0.74$,
  $T \le 183$, $|\mathbf v| \le 1.4$, the DynSGS coefficient $\le 1.6$, the
  explicit-diffusion number $\le 0.004$, and the CFL monitor never passes its
  initial 0.155. At $t = 10$ the coarse mesh had $T$ to $6.6\times10^3$ (434
  nodes above 300), 9 nodes at the density floor and the coefficient at 5.9.
  The largest $p/p_e$, 2.8 at $(0, 9)$ at $t = 10$, is the chromosphere above
  the rope lifted by about one scale height, as the rope's top
  ($|B_z| = 0.1$ on $x = 0$) rose from $h = -7.0$ to $-6.0$. It rises at the
  same rate in all three runs ($-6.1$ at $t = 10$ on the coarse mesh, with
  either operator); the refined one keeps more of its axial field (largest
  $|B_z|$ 3.8 at $t = 10$, against 3.2; 5 at $t = 0$).

  None of the three runs has emerged any flux by $t = 10$ (no $B_z$ above
  $h = 0$; on the coarse mesh the rope's top is at $h = -2.9$ at $t = 20$), and
  sound from the rope reaches the photosphere at $t \approx 5$: the coarse
  run's floors from $t = 2.3$ are numerical, and these runs test the start-up
  transient and the early rise, not the emergence. 0.54 s per step
  ($10^4$ steps in 90 minutes).
- **From $t = 10$ to 20 at `LK2D_REF = 4`, both forms**, each restarted
  from the $t = 10$ snapshot of the run above (Restart):
  - **absolute form:** the pressure floor fires at $(0, 7)$, in the
    chromosphere above the rope, from $t = 12.74$; the nodes at the floors
    grow to 28 at $t = 14.2$, at $h = 7.7$–9.7, and 93 at $t = 15.1$, at
    $h = 7.7$–16.3, over $|x| \lesssim 8$: the chromosphere and the base of
    the transition region above the rope, with the corona falling into them at
    up to 13.5. The DynSGS coefficient reaches 8.9, the explicit-diffusion
    number 0.074 and $T$ 225; the run was stopped at $t = 15.2$;
  - **reference-weighted operator** (`LK2D_REF_WEIGHT=1`): the pressure
    floor fires at $(0, 7.7)$ from $t = 14.06$. At most 20 nodes, all at
    $h = 7.7$–10.3, are at the floors until $t = 17.3$, and none after.
    Over $10 \le t \le 20$, $|\mathbf v| \le 3.4$, the coefficient $\le 3.1$,
    the explicit-diffusion number $\le 0.008$, $T \le 182$, and the CFL
    monitor never passes its initial 0.159. At $t = 20$, $\rho/\rho_e \ge 0.76$
    and $p/p_e \ge 0.70$.

  Both start the same way. At $t = 10$ the rope has lifted the chromosphere
  over $|x| \lesssim 10$ into a dome: $\rho/\rho_e = 2.5$ at $(0, 9)$, rising
  at 0.75. From $t \approx 11.5$ the dome's material drains sideways into
  ridges at $|x| \approx 15$–18, where $\rho/\rho_e$ at $h = 11$ is 4.9
  (absolute) or 3.2 (weighted) at $t = 15$. Its middle thins and the corona
  falls into it, at up to 13.5 with the absolute form and 3.4 with the
  weighted one. With the weighted operator the infall refills the middle: at $t = 20$
  the column over $x = 0$ falls at up to 0.75, with $\rho/\rho_e$ up to 4.7
  at $T \approx 1.7$ around $h = 12$, and $T = 12.25$ (at $h = 15$ in the
  reference) is at $h = 16.8$. The drainage is smooth over several elements
  and appears with both operators, so it may be the chromosphere's own
  response to the rope's rise. How deep the cavity gets depends on the
  operator, and where the floors hold nodes they, not the equations, set the
  state. The rope's top is at $h = -4.8$ at $t = 15$ in both runs and $-3.9$
  at $t = 20$ (the coarse absolute run's: $-4.7$, $-2.9$). No flux has
  emerged by $t = 20$; at the rope's 0.2 per $t_0$ it would reach $h = 0$ at
  $t \approx 40$.
- **The reference-weighted operator at `LK2D_REF = 1` from $t = 10$ to 20**,
  restarted from the $t = 10$ snapshot of its run above: no evacuation above
  the rope and no heating. $T \le 178$, $|\mathbf v| \le 1.2$, the
  coefficient $\le 5$, the explicit-diffusion number $\le 0.011$, and the CFL
  monitor never passes its initial 0.161. Two nodes touch the pressure floor
  at $t = 13.1$–13.3, at $(\pm171, 6)$, and at $t = 19.6$–19.7, at
  $(\pm144, 6)$, far from the rope. From $t \approx 12$ the chromosphere at
  $|x| > 100$ has node-to-node cold spots on the element face $h = 6$:
  $p/p_e$ down to 0.29 and $T/T_e$ to 0.35 at $t = 20$, 72 nodes off $p_e$
  by more than 10 %. At `LK2D_REF = 4` there are 6, with $p/p_e \ge 0.86$;
  the absolute form on this mesh, with the bursts above the rope, is at the
  floors there from $t \approx 15$, $T/T_e$ to 0.05. The dome above the rope
  is weaker than at
  `LK2D_REF = 4` (the largest $p/p_e$ is 1.7, at the photosphere above the
  rope, against 3–4 in the chromosphere), and the rope is more diffused, as
  in the absolute run on this mesh: at $t = 20$ its largest $|B_z|$ is 2.1
  against 2.5, and its top is at $h = -2.9$ against $-3.9$.

## Files

| file | content |
|---|---|
| `user_inputs.jl` | mesh and Δt, the problem keys, DynSGS, positivity floors, CFL monitor, stage diagnostics, $\nu$/$\eta$ |
| `initialize.jl` | atmosphere, line dipole, rope, $c_h$, Δh, the restart |
| `user_flux.jl` | 2D ideal GLM-MHD fluxes (three components of $\mathbf v$, $\mathbf B$), $\gamma = 5/3$ |
| `user_source.jl` | gravity, GLM damping, velocity damping |
| `user_bc.jl` | $\rho\mathbf v = 0$, $\psi = 0$ on all sides |
| `user_primitives.jl` | DynSGS acts on $q - q_e$; output fields |
| `LK22_2D.geo` | the stretched quadrilateral mesh (`ref` = 1–4) |
