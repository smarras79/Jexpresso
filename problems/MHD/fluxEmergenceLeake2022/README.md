# 3D visco-resistive GLM-MHD: flux emergence into a coronal dipole, the eruption parameter study (Leake, Linton & Antiochos 2022)

This case runs the emergence of a buoyant, twisted convection-zone flux rope
into the field of a coronal point dipole, the eight simulations and the
"No dipole" run of

> J. E. Leake, M. G. Linton, S. K. Antiochos,
> *The Role of Reconnection in the Onset of Solar Eruptions*,
> ApJ **934**:10 (2022). https://doi.org/10.3847/1538-4357/ac74b7

The rope is the one of Leake et al. (2013, 2014)
(`problems/MHD/fluxEmergenceLeake2014`), rotated about the vertical by
$\theta = n\pi/4$, $n = 0, \dots, 7$ (the paper's simulations 0–7). The
angle sets how the rope's twist field above its axis lies against the
dipole field above it: parallel for $\theta = 0$, antiparallel for
$\theta = \pi$. The paper finds:

| simulations | $\theta$ | outcome |
|---|---|---|
| 7, 0, 1 | twist mostly parallel to the dipole | stable: a coronal flux rope forms and rises to ≈ 18 Mm by $600t_0$ |
| 3, 4, 5 | twist mostly antiparallel | breakout reconnection, then flare reconnection below the rope and an eruption at $150$–$200t_0$, confined at 50–70 Mm |
| 6 | twist perpendicular | early eruption, as 3–5 |
| 2 | twist perpendicular | slow rise, a later, weaker confined eruption ($t \approx 350t_0$) |
| No dipole | (rope of simulation 0) | no eruption |

Simulation 4 ($\theta = \pi$, the orientation of the 2014 case) is the
eruption the paper describes in detail (its Figs. 3, 4, 12) and is the default.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022
mpiexec -n 128 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022
```

(use the `mpiexec` that matches the MPI your `MPI.jl` is built against, see
`INSTALL.md`), or `Jexpresso.run_case("MHD", "fluxEmergenceLeake2022")`.
Output goes to `./output/MHD/fluxEmergenceLeake2022/<run>/` every $10t_0$
(and at the paper's $142$ and $184t_0$) up to $t = 400t_0$.

The simulation is chosen in `user_inputs.jl`:

| key | values | meaning |
|---|---|---|
| `:lk22_sim` | 0 … 7 (default 4) | $\theta = n\pi/4$, the paper's simulation number |
| `:lk22_dipole` | `true` / `false` | `false` with `:lk22_sim => 0` is the paper's "No dipole" run |
| `:lk22_rope` | `true` / `false` | `false` runs the background alone (atmosphere + dipole), which must stay at rest |
| `:lk22_d`, `:lk22_zd` | 15000, −100 | dipole strength and depth (Deviations, items 2–3) |
| `:lk22_mf` | 1.0 (default) / 1.25 | $p = \rho T/m_f$ (Deviations, item 5) |
| `:lk22_cz_gradient` | 1.0 | convection-zone temperature gradient factor (Eq. 23 as printed) |
| `:lk22_glm_alpha` | 1.0 | Dedner damping $\alpha_p$ of $\psi$ |
| `:lk22_damp_x_frac`, `:lk22_damp_z_frac` | 0.75, 0.875 | the damping layers start at $|x|, |y| = 3L_x/8$ and $z = 7\max(z)/8$ |
| `LK22_REF` (environment variable) | 1, 2 (default), 3, 4 | mesh and time step (below); also `LK22_TEND`, `LK22_OUT_DT` and the restart (Restart) |
| `LK22_REF_WEIGHT`, `LK22_ENERGY_SPLIT` (environment variables) | 0 (default), 1 | the reference-weighted DynSGS operator of `fluxEmergenceLeake2022_2D` (its Deviations, item 5) and the split energy flux (`src/kernel/physics/SGS.jl`); `LK22_REF_WEIGHT=1` on any mesh (Meshes and time step, "Vertical resolution") |

## Units

$L_0 = 170$ km, $B_0 = 0.13$ T, $g_0 = 274$ m s$^{-2}$, hence
$P_0 = B_0^2/\mu_0 = 1.345\times10^4$ Pa, $\rho_0 = 2.887\times10^{-4}$ kg m$^{-3}$,
$v_0 = 6.825$ km s$^{-1}$, $t_0 = 24.909$ s, $T_0 = 5637$ K,
$J_0 = B_0/(\mu_0L_0) = 0.609$ A m$^{-2}$, $E_0 = \rho_0v_0^2L_0^3 = 6.607\times10^{19}$ J.
The output is in these units: lengths × 0.17 give Mm, $B$ × 1300 gives G.
In them $g = 1$, the magnetic pressure is $\tfrac12|\mathbf B|^2$, $\gamma = 5/3$
and $p = \rho T$ (item 5).

## The problem

**Atmosphere** (Eqs. 23–25, those of the 2014 case): $T = 1 - (\gamma-1)/\gamma\cdot z$
below $z = 0$, $T = 1$ for $0 < z < 10$, $T = 150^{(z-10)/10}$ for
$10 \le z < 20$, $T = 150$ above; hydrostatic with $\rho(0) = 1$.
From 0.0135 kg m$^{-3}$ at $-5.1$ Mm to $1.2\times10^{-11}$ at the base of the
corona (3.4 Mm) and $4\times10^{-13}$ at the top (90.2 Mm), the paper's Fig. 1.

**Dipole** (Eqs. 26–27): a point dipole at $(0, 0, z_d)$ below the domain,
moment along $x$:
$\mathbf B_{\rm dip} = d\,[3(\hat{\mathbf m}\cdot\hat{\mathbf r})\hat{\mathbf r} - \hat{\mathbf m}]/r^3$,
$\mathbf r = (x, y, z - z_d)$, $d = 15000$, $z_d = -100$, $\hat{\mathbf m} = -\hat{\mathbf x}$:
$|\mathbf B| = 0.015B_0 = 19.5$ G at the photosphere above the dipole.

**Flux rope** (Eqs. 28–34): in the frame rotated by $\theta$,
$x_r = x\cos\theta + y\sin\theta$, $y_r = -x\sin\theta + y\cos\theta$,
the rope lies along $y_r$ at $x_r = 0$, $z = -12$, with
$B_{\rm axial} = B_te^{-r^2/w^2}$ along $\hat{\mathbf y}_r$ and a right-handed twist
$B_{\rm twist} = (r/w)B_te^{-r^2/w^2}$, $r^2 = x_r^2 + (z + 12)^2$,
$B_t = 5$ (6500 G), $w = 2.5$. Its pressure deficit
$p_1 = \tfrac12B_t^2e^{-2r^2/w^2}[(w^2/2 - r^2)/w^2 - 1]$ balances its Lorentz
force, and $\rho = \rho_0(z)[1 + (p_1/p_0(z))e^{-(y_r/\lambda)^2}]$,
$\lambda = 10$, makes its middle buoyant.

**Dissipation** (Sec. 2.1): $\eta = \nu = 0.01$ in normalized units
(Re = Rm = 100), uniform.

**Boundaries** (Sec. 2.3): all velocity components and the normal
gradients of the other variables zero on all six faces; velocity damping
in the outer $L_x/8$ of each side ($|x|, |y| > 207.9 = 35.3$ Mm) and above
$z = 7\max(z)/8 = 464.2$ (78.9 Mm): every step $\mathbf v \to \mathbf v/F$,
$F = 1 + \Delta t\,N$, $N$ rising linearly from 0 to 1 at the boundary
(Eqs. 21–22); here the equivalent source $-N\mathbf v$.

**Domain and grid** (Sec. 2.3): $x, y \in [-277.2, 277.2]$,
$z \in [-30, 530.5]$ ($[-47.1, 47.1]^2\times[-5.1, 90.2]$ Mm) on
$512\times512\times768$ cells stretched from a uniform $0.5\times0.5\times0.25$
grid by Eqs. (19)–(20) (item 1): 0.085 Mm horizontally in the centre to
0.714 Mm at the sides, 0.0425 Mm vertically up to $z \approx 40$ to 0.2125 Mm at
the top.

## Meshes and time step

`LK22.geo` places the element boundaries on every $k$-th face of the paper's
grid, so with `:nop => 4` the nodes follow its stretching:

| `LK22_REF` | elements | $k$ | file | nodes | `:Δt` |
|---|---|---|---|---|---|
| 1 | $16\times16\times24$ | 32 | `LK22_16x16x24.msh` (shipped) | 0.41 M | $4\times10^{-3}$ |
| 2 | $32\times32\times48$ | 16 | `LK22_32x32x48.msh` (shipped) | 3.2 M | $2\times10^{-3}$ |
| 3 | $64\times64\times96$ | 8 | generate | 25 M | $10^{-3}$ |
| 4 | $128\times128\times192$ | 4 | generate | 202 M | $5\times10^{-4}$ |

```bash
cd problems/MHD/fluxEmergenceLeake2022
gmsh -3 LK22.geo -setnumber ref 3 -o LK22_64x64x96.msh
gmsh -3 LK22.geo -setnumber ref 4 -o LK22_128x128x192.msh
```

`LK22_REF = 4` is the paper's resolution. The smallest LGL gap is 0.173 of
the bottom elements' height ($8, 4, 2, 1\,L_0$). The fast speed of the initial
state peaks at 45 just above the transition region (the GLM speed $c_h$,
printed at start-up), and the paper's reconnection outflows reach tens of
$v_0$; `:Δt` gives $(|\mathbf v| + c_f)\Delta t/\Delta h \approx 0.25$ at
$|\mathbf v| + c_f = 80$. The step is fixed (Lare3D's is not): the CFL
monitor (Notes) prints the actual value as it grows.

**Vertical resolution.** The bottom elements of `LK22_REF = 2` are 4 $L_0$
tall. In the 2.5D cross-section of this case (`fluxEmergenceLeake2022_2D`,
its README's Verification) that mesh, with the absolute DynSGS form of this
deck, hits the pressure floors at $t = 2.3$ and evacuates and heats the
transition region above the rope from $t \approx 5$, before any flux has
emerged; the same columns with elements 2 tall run to $t = 10$ without a
floor, the rope rising at the same rate, and evacuate the chromosphere above
the rope from $t = 12.7$. With the reference-weighted operator
(`LK22_REF_WEIGHT=1`) neither mesh is evacuated by $t = 20$. The elements 2
tall touch the floors at up to 20 nodes over $t = 14$–17. Those 4 tall touch
them at 2 nodes, but grow node-to-node cold spots in the chromosphere far from
the rope ($T$ down to $0.35T_e$ at $t = 20$) and diffuse the rope more. The
2014 case's README ("Vertical resolution") has the mechanism. `LK22_REF = 3`
has elements 2 tall, and eight times the elements of 2.

**Cost.** `LK22_REF = 1` takes about 10 s per time step on 4 threads
(0.41 M nodes; 13 s measured with other jobs on the machine). Scaling by
nodes and steps to $t = 400t_0$: $\approx 1.1\times10^3$ core-hours at
`LK22_REF = 1` ($10^5$ steps), $1.7\times10^4$ at 2, $2.7\times10^5$ at 3 and
$4\times10^6$ at 4 (the paper's resolution), before parallel losses.
Simulation 4 erupts by $200t_0$: `:tend => 250.0` covers it for 60 % of
that.

## Deviations from the paper, and why

1. **Grid stretching.** Eqs. (19)–(20) print $d = d_0 + f\tanh((|s_0| - l)/w)$,
   which is negative near the centre. The form
   $d = d_0[1 + f(1 + \tanh((|s_0| - l)/w))]$, evaluated at each cell's outer
   face, with $l_x = L_{x,0}/2.1$, $w_x = L_{x,0}/10$ for the full width
   $L_{x,0} = 256$ and $l_z = 70$, $w_z = L_{z,0}/20$ for $L_{z,0} = 192$,
   reproduces every number the paper quotes: 0.085 and 0.714 Mm, 0.0425 and
   0.2125 Mm, $\pm47.12$ Mm and 90.19 Mm.
2. **Dipole depth and normalization.** The paper does not print $z_d$. With
   $z_d = -100$ (Leake et al. 2013, 2014) and the normalized
   $\mathbf A = d\,\hat{\mathbf m}\times\mathbf r/r^3$, $d = 15000$ gives
   19.5 G at the photosphere, the paper's "representative 20 G". The literal
   $(\mu_0/4\pi)\,\mathbf m\times\mathbf r/r^3$ with
   $m = 15000B_0L_0^3/\mu_0 = 7.6\times10^{24}$ A m$^2$ would give 1.55 G there
   (or need $z_d = -42.7$). `:lk22_d`, `:lk22_zd` change them.
3. **Dipole sign.** The printed $\mathbf m = d\hat{\mathbf x}$ puts the dipole
   field above the rope along $-x$, while the printed rope ($q = -1/w$ in its
   left-handed $(r_1, \phi_1, y_1)$ frame, a right-handed twist, "for this
   right-handed twist" in Sec. 3.2) has its twist above the axis along
   $+x_r$: $\theta = 0$ would be antiparallel. The paper defines $\theta = 0$ as
   parallel and stable (Figs. 2, 10), so the moment is $-d\hat{\mathbf x}$
   here.
4. **Density perturbation sign.** Eq. (34) prints $\rho_1 = -\rho_0p_1/p_0\,e^{-(y_1/\lambda)^2}$,
   which with $p_1 < 0$ is positive; the text says "$\rho_1$ is maximally
   negative at the middle of the rope", as in Leake et al. (2013, 2014). The
   case uses $\rho_1 = \rho_0p_1/p_0\,e^{-(y_r/\lambda)^2}$.
5. **Gas law in normalized units.** Sec. 2.1–2.2 give $\mu_m = 1.25m_p$ and
   $T_0 = m_pL_0g/k_B$, i.e. $p = \rho T/1.25$. The paper's own Fig. 1 does
   not follow that: its densities (0.012 kg m$^{-3}$ at $-5.1$ Mm,
   $10^{-11}$ at the base of the corona, a 25 Mm coronal scale height) are those
   of $p = \rho T$, and so are Leake et al. (2014)'s pressures (Fig. 1a); with
   $p = \rho T$ the printed convection-zone gradient is also exactly the
   adiabatic one, as the 2013 paper says it is. The reduced mass sits in
   the code's temperature unit, so `:lk22_mf => 1.0`; 1.25 gives the
   literal-text atmosphere, whose corona is 20 times less dense.
6. **Kinematic viscosity**, as in the 2014 case: the dynamic
   $\nu = 3.35\times10^3$ kg m$^{-1}$ s$^{-1}$ would be a kinematic viscosity of
   $7\times10^6$ at the top of the corona; $\nu = 0.01$ is kinematic here.
7. **Numerics.** The paper uses Lare3D (staggered Lagrangian-remap,
   constrained transport, shock viscosity with $c_1, c_2 = 0.5, 0.8$). Here:
   continuous spectral elements at $N = 4$, GLM cleaning ($c_h$ fixed at the
   initial fast-speed maximum, Dedner damping $\alpha_p = 1$), DynSGS in
   conserved form on the departure from the reference state (legacy sensor,
   fast-speed floors, thermal residual) with the explicit $\nu$, $\eta$ added,
   positivity floors $\rho \ge 0.2\rho_e$, $p \ge 0.2p_e$, and the discrete
   reference RHS subtracted (`:lwell_balanced`), as in the 2014 case (its
   README, items 6 and 8).
8. **Boundaries.** CG imposes the velocity and $\psi = 0$ strongly; $\rho$,
   $E$, $\mathbf B$ on the boundary evolve with their element derivatives (the
   natural condition of the dissipative terms is a zero normal gradient).

## Output

`ρ, u, v, w, p, Bx, By, Bz, ψ, T, β, Bmag, log10ρ` per snapshot (`.pvtu`);
mass and total energy every 10 steps.

**Disk.** About 240 bytes per node: 0.1 GB per snapshot at `LK22_REF = 1`,
0.77 GB at 2, 6 GB at 3, 48 GB at 4. The default 43 snapshots take 33 GB at
`LK22_REF = 2`. The paper's figures use $t = 70, 100, 120, 130, 142, 160, 170,
184, 190, 240, 260, 380, 400t_0$; all but 142 and 184 are on the default
10 $t_0$ cadence.

## Restart

Any snapshot of this case restarts the run, so a long run can be split into
jobs that each fit a cluster's wall-time limit. A snapshot holds the nine
primitives in Float64 at every node, and `initialize.jl` rebuilds the
conserved state from them exactly ($E = p/(\gamma-1) + \tfrac12\rho|\mathbf
v|^2 + \tfrac12|\mathbf B|^2 + \tfrac12\psi^2$), keeps the analytic reference
state and the $c_h$ of $t = 0$, and starts the clock at the snapshot's time.
Each rank reads its own piece, so **every job must use the rank count and the
mesh of the run that wrote the snapshot**. The DynSGS residual history starts
afresh, so the continuation is close to, not bitwise equal to, an
uninterrupted run (checked on the 2014 case, which has the same restart:
its README, Verification). Environment variables, so the deck need not be
edited:

| variable | meaning |
|---|---|
| `LK22_RESTART_DIR` | output directory of the earlier job (`output/MHD/fluxEmergenceLeake2022/output-<date>`) |
| `LK22_RESTART_TIME` | snapshot to start from (nearest in its `simulation.pvd`; default the last) |
| `LK22_TEND` | `:tend` (default 400) |
| `LK22_OUT_DT` | a snapshot every `LK22_OUT_DT` (default 10), plus 142 and 184 |
| `LK22_REF` | mesh (default 2; must be the one of the snapshot) |

For example, simulation 4 to $t = 250t_0$ in three jobs of the same size:

```bash
LK22_TEND=90  mpiexec -n 512 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022
LK22_TEND=170 LK22_RESTART_DIR=output/MHD/fluxEmergenceLeake2022/output-<job 1> mpiexec -n 512 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022
LK22_TEND=250 LK22_RESTART_DIR=output/MHD/fluxEmergenceLeake2022/output-<job 2> mpiexec -n 512 julia --project=. src/Jexpresso.jl MHD fluxEmergenceLeake2022
```

End each job on a snapshot time (a multiple of `LK22_OUT_DT`), or the next
one restarts from the last snapshot before it and repeats the stretch in
between. Each job writes to a new output directory, starting with the
snapshot it restarted from; `tools/plot_leake2022.py` takes the pieces as
`DIR1,DIR2,DIR3`.

## Figures

`tools/plot_leake2022.py` draws the paper's figures from the output (numpy,
matplotlib; pyvista for the 3D figures). `RUN` is an output directory; its
copy of `user_inputs.jl` tells the script $\theta$ (`:lk22_sim`) and whether
there is a dipole, and the curves get the paper's labels and colors.

```bash
T=tools/plot_leake2022.py
python3 $T ic                                                   # Fig. 1 (right)
python3 $T ic $RUN4                                             # ... with the nodal values at t = 0
python3 $T render $RUN4 --times 0                               # Fig. 1 (left)
python3 $T render $RUN4 --times 100,142,184 --jb                # Fig. 3 (top)
python3 $T render $RUN4 --times 100,142,184 --view top          # Fig. 3 (bottom)
python3 $T xsec $RUN4 --var JB --style fig4 --photosphere --times 100,142,184 --xr-max 40 --zlim -5,60   # Fig. 4 (bottom)
python3 $T xsec $RUN4 --var Baxial --photosphere --times 100,142,184 --xr-max 40 --zlim -5,60          # Fig. 4 (top)
python3 $T render $RUN3 --times 200                             # Fig. 5 (and $RUN4, $RUN5, $RUN2, $RUN6)
python3 $T xsec $RUN0 --var Baxial --times 100,170,400 --title "Sim 0: "   # Fig. 6, one row per simulation:
python3 $T xsec $RUN2 --var Baxial --times 100,240,380 --title "Sim 2: "   #   0, 2, 4, 6
python3 $T xsec $RUN4 --var Baxial --times 130,160,190 --title "Sim 4: "
python3 $T xsec $RUN6 --var Baxial --times 120,190,260 --title "Sim 6: "
python3 $T xsec $RUN4 --var JB --times 130,160,190 --title "Sim 4: "      # Fig. 7 (and the other rows)
python3 $T height $RUN7 $RUN0 $RUN1 $RUN2 $RUN3 $RUN4 $RUN5 $RUN6 $RUNND  # Fig. 8 (+ CSVs)
python3 $T energy $RUN7 $RUN0 $RUN1 $RUN2 $RUN3 $RUN4 $RUN5 $RUN6 $RUNND  # Fig. 9 (+ CSVs)
python3 $T xsec $RUN4 --var JB --contours 12 --xr-max 10 --zlim 0,20 --times 130,160      # Fig. 12 (left)
python3 $T lorentz $RUN4 --times 130,160                        # Fig. 12 (right)
```

How the quantities are computed:

- **Cross-section** (Figs. 4, 6, 7, 12): the vertical plane $y_r = 0$
  halfway along the rope, sampled bilinearly from the node lattice;
  $B_{\rm axial} = \mathbf B\cdot\hat{\mathbf y}_r$ in G, $|\mathbf J|/|\mathbf B|$
  with $\mathbf J = \nabla\times\mathbf B$ by second-order differences.
- **Rope height** (Fig. 8): on the line $x = y = 0$, (1) the local maximum of
  $B_{\rm axial}$ nearest to (2) the highest height where the horizontal field
  orthogonal to the axis, $\mathbf B\cdot\hat{\mathbf x}_r$, turns from $-$ to $+$
  going up (the right-handed rope's axis), with the sign holding on two
  nodes either side (`--persist`; skips the node-to-node zigzag of an
  under-resolved field); the band spans the two, the line is their mean, the
  dashed line is the lowest such height when there are several (the original
  rope after flare reconnection), as the paper.
- **Energies** (Fig. 9) in $z \ge 10$: kinetic $\int\tfrac12\rho|\mathbf v|^2dV$;
  free magnetic $\int\tfrac12|\mathbf B|^2dV$ minus that of the potential field
  with the same normal component on the four sides, the top and the plane
  $z = 10$, solved by a vertex-centred finite-volume Laplacian (conjugate
  gradients, Jacobi preconditioned); the net flux the discrete
  $\nabla\cdot\mathbf B$ leaves on the boundary (printed) is spread over it
  first. `--stride 2` uses every second node (8 times fewer unknowns) on the
  big meshes.
- **Lorentz force** (Fig. 12): $(\mathbf J\times\mathbf B)_z$ split into the
  tension $(\mathbf B\cdot\nabla)B_z$ and the pressure $-\partial_z|\mathbf B|^2/2$
  on $x = y = 0$; the dotted line marks the axial-field maximum.
- **3D field lines** (Figs. 1, 3, 5): seeded every 1.5 Mm on $x = y = 0$
  (`--seeds z0,z1,dz`), colored by where their two ends are: purple both on
  the dipole (bottom, top or a side above the surface), yellow both on the
  rope (a side boundary below `--rope-zmax`, 0 Mm), cyan one of each; gray
  $B_z$ at $z = 0$ in $\pm100$ G; `--jb` adds the $|\mathbf J|/|\mathbf B|$
  cross-section. Seeds and camera are set by eye from the paper's figures.

Figs. 2 (schematic), 10 and 11 (hand-picked field lines and arrows) are not
reproduced.

## Verification here

Short runs on 4 threads at `LK22_REF = 1`; the full runs are for a cluster.

- **Grid.** `LK22.geo` puts the faces at $\pm47.12$ Mm and $z = -5.10$ to
  90.19 Mm, the paper's extents.
- **Atmosphere** (`tools/plot_leake2022.py ic`): the nodal $T$ and $\rho$ on
  $x = y = 0$ fall on the closed-form curves, $\rho$ from 0.0135 kg m$^{-3}$
  at $-5.1$ Mm to $1.2\times10^{-11}$ at 3.4 Mm and $4.0\times10^{-13}$ at the
  top, the paper's Fig. 1.
- **Dipole.** 19.5 G along $+x$ at $(0, 0, 0)$; divergence- and curl-free to
  round-off (checked offline).
- **Rope orientation.** Offline, the twist above the axis points along $+x$,
  $+y$, $-x$ for $\theta = 0, \pi/2, \pi$ (parallel, perpendicular,
  antiparallel to the dipole). In the code: simulation 1 ($\theta = \pi/4$)
  puts the rope along $y = -x$ (top-view render) with its axis found at
  $-1.99$ Mm on $x = y = 0$ (`height`), simulation 4 at $-2.0$ to $-2.15$ Mm.
- **Background at rest** (`:lk22_rope => false`): with `:lwell_balanced` every
  field at $t = 0.05t_0$ equals the initial one to the last bit (checked with
  $m_f = 1.25$; the subtraction does not depend on it).
- **Simulation 4 to $t = 0.4t_0$**: no NaN; the CFL monitor stays at its
  initial 0.126; the positivity floors start at $t = 0.14t_0$ ($10^4$ hits by
  0.19); the GLM transient (2014 case, item 8) reaches $|\mathbf v| = 1.5v_0$
  at $z \approx 18$ and $p/p(0) = 2.2$ in the transition region. The coronal
  free energy is $-0.12E_0$ at $t = 0$ (the discretization error of the
  potential-field solve, out of 60.8 for the dipole) and 0.63 at $0.4t_0$, the
  kinetic energy 0.06. The $16\times16\times24$ elements are 16 $L_0$ wide in the
  centre for a rope of width 2.5, so this tests the pipeline, not the physics.
- **Figures.** Every command above ran on these outputs.

- **Restart** (Restart). A run with snapshots at $t = 0.02$ and $0.04$
  restarted from its $t = 0.02$ snapshot: the restart state equals the
  snapshot to $4\times10^{-16}$ of each field's maximum ($T$
  $2\times10^{-15}$); five steps later $\rho$, $p$ and $T$ agree with the
  uninterrupted run to $2\times10^{-4}$, $\mathbf v$ to 2–3.5 % and
  $\mathbf B$, $\psi$ to 0.8–1.8 % of their maxima (the DynSGS residual
  history restarts from zero). The same code was checked on 2 MPI ranks in the
  2014 case.

Nothing was run here at `LK22_REF` ≥ 2.

## Notes

- **Partition.** The deck sets `:lxy_partition => false`: a balanced 3D
  partition (METIS, through GridapGmsh's distributed reader, which reads the
  mesh file on every rank). Jexpresso's default (`true`) bins the $x$–$y$
  bounding box into a uniform $n_x\times n_y$ grid of columns, one per rank;
  on this horizontally stretched mesh the outer elements are several bins
  wide, so from a few hundred ranks on some ranks get no element and the run
  stops in `mod_mesh_read_gmsh!` with `reducing over an empty collection`
  (at `LK22_REF = 2`: 96 of 448 and 412 of 896 ranks empty), and at any
  count the central ranks carry more than the average (1.75 times at 14
  ranks, 5.2 at 896). With METIS the useful rank count is limited by the
  elements per rank (49 152 elements at `LK22_REF = 2`).
- **Progress.** Snapshots are $10t_0$ apart (5000 steps at `LK22_REF = 2`,
  3.2 M nodes; the 2014 case's similar $38^3$ mesh takes about 8 s per step
  on 14 cores), and the floor and CFL lines below are
  printed only when their totals pass a power of ten or grow by 20 %, so a
  healthy run can print nothing for hours. The deck sets
  `:lstep_heartbeat => true` (`#   step N   t = …` every 100 steps), and
  `conservation.dat` in the run's output directory gets $t$, the step, the
  total mass and energy every 10 steps
  (`tail -f output/MHD/fluxEmergenceLeake2022/<run>/conservation.dat`).
- **CFL monitor.** `# CFL: max (|v| + c_f) Δt/Δh = …` is printed by rank 0
  every time the largest value over the nodes has grown by 20 % (checked
  every 100 Runge–Kutta stages), as a `WARNING` above 0.5: halve `:Δt` and
  restart from the last snapshot if it keeps growing. The advective and
  acoustic lines of the shared CFL printout assume the CompEuler layout and
  do not apply.
- **The positivity-floor line**, `# positivity floors: N node hits so far on
  all ranks …`, is printed each time the total passes a power of ten; its
  `t` is the Runge–Kutta stage time.
- **Constants shared by name** (`γ_mhd`, `c_h_mhd`, `pressure_mhd`) with the
  other MHD cases: restart Julia between cases with a different $\gamma$.

## Files

| file | content |
|---|---|
| `user_inputs.jl` | `LK22_REF`, Δt, the problem keys, DynSGS, positivity floors and CFL monitor, $\nu$/$\eta$ |
| `initialize.jl` | atmosphere, point dipole, rotated rope, $c_h$, Δh, damping layers |
| `user_flux.jl` | 3D ideal GLM-MHD fluxes, $\gamma = 5/3$ |
| `user_source.jl` | gravity, GLM damping, velocity damping |
| `user_bc.jl` | $\rho\mathbf v = 0$, $\psi = 0$ on all faces |
| `user_primitives.jl` | DynSGS acts on $q - q_e$; derived output fields |
| `LK22.geo` | the stretched hexahedral mesh (`ref` = 1–4) |
