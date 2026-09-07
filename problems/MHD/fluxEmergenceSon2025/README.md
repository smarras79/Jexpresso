# 2D ideal GLM-MHD: flux emergence in the solar atmosphere (Son, Jang & Magara 2025)

This case runs the two-dimensional emergence of a horizontal magnetic flux
sheet through a two-temperature (chromosphere + corona) stratified
atmosphere — the nonlinear Parker instability of Shibata et al. (1989a) — with
the setup of

> D. Son, Y. Jang, T. Magara,
> *A Comparative Analysis of High-resolution Shock-capturing Schemes for
> Two-dimensional Magnetohydrodynamic Simulation of Flux Emergence in the
> Solar Atmosphere*, ApJS **277**:46 (2025) — **Sections 2.1–2.2** (model)
> and **4.1** (results). https://doi.org/10.3847/1538-4365/adb617
> (open access; `problems/MHD/Son_2025_ApJS_277_46.pdf`)

on $[0, 80H_0]\times[0, 35H_0]$, $\gamma = 1.05$, over $t\in[0, 54\tau_0]$.

A full write-up of the governing equations, the units, the initial
magnetostatic atmosphere, the boundary conditions and what to expect from the
run is in [EQUATIONS.md](EQUATIONS.md); this file covers the implementation
choices and how to run.

## Run

Ten MPI ranks (the configuration this case was sized for):

```bash
mpiexec -n 10 julia --project=. src/Jexpresso.jl MHD fluxEmergenceSon2025
```

(use the `mpiexec` that matches the MPI your `MPI.jl` is built against, see
`INSTALL.md`; `julia --project=. -e 'using MPI; println(mpiexec())'` prints
the bundled one). Serially:

```bash
julia --project=. src/Jexpresso.jl MHD fluxEmergenceSon2025
```

or, interactively, `Jexpresso.run_case("MHD", "fluxEmergenceSon2025")`.

Output goes to `./output/MHD/fluxEmergenceSon2025/<run>/`.

## Sizing: "coarsest admissible" at N = 4

| | value |
|---|---|
| elements | $80\times35$, each $1H_0\times1H_0$ (`FE_80x35.geo/.msh`, shipped; `FE_80x70` also ships) |
| polynomial order | 4 (LGL) → $320\times140$ unique points |
| nodal spacing | 0.25 $H_0$ mean, 0.173 $H_0$ smallest |
| time step | $2.5\times10^{-3}\tau_0$, Courant ≈ 0.07 initially, ≈ 0.15 at late times |
| steps to $54\tau_0$ | 21,600 |
| output | every $\tau_0$ (55 snapshots) |

The tanh transitions of the problem are $0.5H_0$ (flux-sheet edges) and
$0.6H_0$ (chromosphere–corona), so $1H_0$ elements at $N = 4$ are the
coarsest grid that keeps three nodes across them. It is comparable to the
paper's coarsest $300^2$ grid horizontally ($\Delta x = 0.27H_0$) and two
times coarser vertically ($\Delta z = 0.12H_0$). To match the paper's
vertical spacing point `:gmsh_filename` at the shipped `FE_80x70.msh`
($0.5H_0$ tall elements, $320\times280$ points) and halve `:Δt`; the cost
quadruples (twice the elements, twice the steps).

## What is (and is not) implemented

- **Conservative GLM-MHD fluxes**, 9 unknowns `(ρ, ρu, ρv, ρE, ρw, Bx, By,
  Bz, ψ)`: the same `user_flux.jl` as `orszagTangBormanis2024` and
  `kelvinHelmholtzChan2022`, with **γ = 1.05** (the paper's value) and a
  pressure floor inside the flux (see EQUATIONS.md §2). The GLM cleaning
  speed `c_h` is the maximum wave speed of the initial condition
  (= the coronal sound speed, $5\,C_s$), constant in time.
- **Sources** (`user_source.jl`): gravity $g_0 = 1/\gamma$ in the vertical
  momentum and the energy; Dedner's parabolic ψ damping with the paper's
  Mignone parametrization $\alpha_p = 0.2$ (`glm_alpha_p_mhd`); the
  absorbing layer of the top boundary, a Rayleigh damping toward the
  initial state above $z_s = 30H_0$ with $\sigma_{max} = 2/\tau_0$
  (`sponge_zs_mhd`, `sponge_sigma_mhd`). All three are `Ref`s tunable from
  the REPL.
- **Boundaries** (`user_bc.jl`): periodic in $x$ (mesh tags), symmetric
  (free-slip, $B_z = 0$) bottom wall, free-slip top wall behind the
  absorbing layer.
- **Initial condition** (`initialize.jl`): the paper's Eqs. (1)–(5). The
  magnetostatic equilibrium is integrated numerically on a $10^{-4}H_0$
  table and interpolated at the nodes, and the residual of the *discrete*
  vertical balance of that interpolant (at most $4\times10^{-4}$ of $\rho g$
  on the shipped mesh) is subtracted as a static source (`fe_well_balanced`),
  so the initial state is an exact equilibrium of the discrete operator. **The paper never states the
  plasma beta $\beta_*$ of the sheet**; $\beta_* = 1$ was inferred from
  its Fig. 1(b) (table in EQUATIONS.md §3) and is a `Ref`
  (`fe_beta_star`).
- **Stabilization: DynSGS** (`:visc_model => DSGS_MHD()`, see
  [DSGS.md](../../../DSGS.md) §4), with `:dsgs_gamma => 1.05` matching the
  flux and the stratified-atmosphere variants of DSGS.md §4.5 switched on:
  `:dsgs_local_norms` with `:dsgs_local_rel = 1` (the residual is measured
  against the element's own $\rho c/\tau$; with the domain norm the
  $10^8$-times denser photosphere hid the corona from the sensor and a
  grid-scale sawtooth grew across the transition region, with a $10^{-3}$
  floor the sensor saturated over the whole quiet chromosphere and eroded
  the sheet), and `:dsgs_conserved` — one kinematic coefficient and a
  Laplacian on the conserved variables, applied to the departure from the
  magnetostatic reference state $q - q_e$ (`user_primitives.jl`). That
  last choice is what keeps the density positive at the transition region
  without breaking the thermodynamics: the physical form of the
  Orszag–Tang case (κ∇T on the energy, no mass diffusion) smeared the
  temperature jump and let the contact undershoot below its $10^{-8}$ light
  side; mass diffusion under the $T$ closure drove $p$ negative; the
  conserved Laplacian on $q$ itself is a steady mass source in an
  exponential atmosphere. On $q - q_e$ it vanishes at rest: measured to
  $t = 8\tau_0$, horizontal-mean vertical velocity $\lesssim 10^{-4}$, sheet
  field unchanged to three digits, $\mu \approx 10^{-4}$–$10^{-3}$ in the
  sheet and $10^{-2}$ only at the transition region, the Parker crest rising
  at $0.1\,C_s$. No filter, no entropy-stable/KEP fluxes. All multipliers
  `:μ` are 1.
- **Positivity floors** (`fe_positivity_limiter!`, top of `user_inputs.jl`,
  passed to `CarpenterKennedy2N54` as its stage limiter): $\rho \ge 0.2\rho_e(z)$
  and $p \ge 0.2p_e(z)$ at every Runge–Kutta stage, $p$ raised through the
  total energy at fixed velocity and field. Why: the sheet's initial
  adjustment launches an acoustic wave (period ≈ 6 τ₀, above the isothermal
  cutoff) whose amplitude grows as $\rho^{-1/2}$ up the 8-decade
  stratification — ≈ 0.02 C_s at the sheet, 0.1 at $z = 13H_0$ ($t = 8$),
  0.4 at the transition region ($t = 12$) — lifts the corona, which falls
  back under gravity at $2\,C_s$ ($t = 15$) and corrugates the 25× density
  contact; in the rarefied troughs $\rho$ fell below the coronal value and
  $T$ past 300 with the DynSGS coefficient already near its cap there, and the
  run broke (measured). Dissipation cannot prevent an evacuation; the
  floors are the density/pressure floor of the finite-volume solar codes
  (the paper itself refers to its schemes' behaviour "in near-vacuum
  regions with high-Mach-number flows"). The run prints the first floor
  hit and every decade of hits with time and place (`# positivity
  floors:` lines), so the log says whether they were ever needed.
- **A shared-code fix this case needed**: the boundary routine only imposes
  a value that differs from the current one, and that test used an
  absolute tolerance of $2\times10^{-6}$ — with $\rho = 7\times10^{-9}$ at the
  top of the corona the wall-normal momentum never reached it, the
  free-slip top was silently never applied and the atmosphere drained
  through it with an exponentially growing downflow. The test is now
  relative (`bc_value_changed` in `src/kernel/boundaryconditions/BCs.jl`).
- **NOT implemented**: the paper's per-step update of $c_h$ (the initial
  value already bounds the coronal sound speed and is of the order of the
  late-time Alfvén speed); the Powell/Galilean-GLM non-conservative term
  (only needed by split-form schemes).

## Figures

The solver writes PNGs directly (`:outformat => "png"`); nothing has to be
post-processed. At every output time, in the run directory:

| file | contents | paper |
|---|---|---|
| `ρ-it<n>.png` | $\log_{10}(\rho/\rho_0)$, `jet` colormap fixed to $[-8.5, 0]$, black magnetic field lines (isocontours of $A_y$), white velocity vectors with the "= 5.0" reference arrow | Fig. 2 (and the Fig. 14 sequence) |
| `v-it<n>.png`, `vA-it<n>.png`, `Bx-it<n>.png` | $V_z/C_s$, $V_A/C_s$, $B_x/B_0$ maps | — |
| `p-it<n>.png`, `T-it<n>.png`, `β-it<n>.png` | $\log_{10}p$, $T/T_0$, $\log_{10}\beta$ maps | Fig. 6(e) for $\beta$ |
| `profile-it<n>.png` | $V_z/C_s$, $V_A/C_s$, $\log_{10}(B_x/B_0)$, $\log_{10}(\rho/\rho_0)$ along $x = 40H_0$, axes fixed to $[0,1.7]$, $[0,4.5]$, $[-2,1]$, $[-8.5,1]$, $z_{cor} = 18H_0$ marked | Fig. 5 (paper times 33, 40, 47, 51–54 $\tau_0$) |

`iter` counter `n` = output slot: `it1` is $t = 0$, `it52` is $t = 51\tau_0$
(paper Fig. 2), `it55` is $t = 54\tau_0$. Under MPI the nodal data is gathered
on rank 0, which renders the whole domain — one file per variable, not one
per rank.

These figures use the generic options of the 2D PNG writer added with this
case (`:plot_vars`, `:plot_log10`, `:plot_clims`, `:plot_fieldlines`,
`:plot_vectors`, `:plot_overlay_on`, `:plot_profile_*`, …; defaults and
meaning in `src/io/mod_inputs.jl`), so they can be re-styled from
`user_inputs.jl` alone. Set `:plot_dsgs => true` to also get one
`μ_dsgs_<var>-it<n>.png` panel per equation with the DynSGS viscosity that
was applied.

`:outformat => "vtk"` writes the same twelve output variables
(`ρ u v w p Bx By Bz ψ T vA β`) plus the `mu_dsgs_<var>` fields to
`iter_<n>.pvtu` for ParaView.

## Results (shipped configuration, run to $54\tau_0$)

Measured on the run that produced this commit: 4 MPI ranks, `FE_80x35`,
$N = 4$, $\Delta t = 2.5\times10^{-3}\tau_0$, 21,600 steps, PNG output every
$\tau_0$. It ran to $t = 54\tau_0$ without an abort in 5194 s of wall
clock on four cores of the development container, of which ≈ 5 min were
start-up and compilation (≈ 0.23 s per step, ≈ 1.6 min per $\tau_0$), so
ten cores should need roughly 35–40 min.

**Against the paper** (IMWENO-P at $300^2$ unless noted; the PNGs named
below are in the run directory):

| quantity | paper | this run |
|---|---|---|
| crest height at $t = 33\tau_0$ (Fig. 14b) | $\approx 10H_0$ | $11.5H_0$ (`ρ-it34.png`) |
| crest at $t = 43$ / $47\tau_0$ (Fig. 14c,d) | $\approx 17$ / $21H_0$ | $18.5H_0$ at $t = 44$ / $21H_0$ at $47$ (`ρ-it45/48.png`) |
| loop at $t = 51\tau_0$ (Fig. 2) | top $\approx 23$–$25H_0$, corona lifted to $27H_0$, legs at $x \approx 25$, $55H_0$, dense pockets ($\log_{10}\rho \approx -3$) at $z \approx 7H_0$ | the same: top $25H_0$, lifted corona to $27H_0$, legs at $25$, $55H_0$, pockets at $z \approx 7H_0$ (`ρ-it52.png`) |
| centerline $V_z$, $t = 51\tau_0$ (Fig. 5c, 6a) | peak $\approx 1.2\,C_s$ at $z \approx 26H_0$ | $1.15\,C_s$ at $z \approx 27H_0$ (`profile-it52.png`) |
| centerline $V_A$, $t = 51\tau_0$ (Fig. 5g, 6b) | peak $\approx 3.9\,C_s$ at $z \approx 21H_0$ | $3.2\,C_s$ at $z \approx 20H_0$ |
| centerline $\log_{10}\rho$, $t = 51\tau_0$ (Fig. 5o) | $-4$ at $z = 15$, $-5$ at $20$, $-5.3$ at $25H_0$, drop to $-8$ at $27$–$28H_0$ | $-4.2$, $-5$, $-5.4$; drop to $-8$ over $z = 28$–$30H_0$ |
| centerline $\log_{10}B_x$, $t = 51\tau_0$ (Fig. 5k) | $-0.65$ at $z = 5$, $-1.5$ at $18H_0$ (Gaussian units, see EQUATIONS.md §6) | $-1.2$ at $z = 5$, $-2.0$ at $18H_0$: the paper's curves minus $\log_{10}\sqrt{4\pi} = 0.55$ |
| downflows at $t = 51\tau_0$ (paper text) | $4$–$5\,C_s$ along the loop sides, $2$–$3\,C_s$ near the footpoints | $5\,C_s$ at $(x, z) \approx (5, 16)$ and $(75, 16)$, $3\,C_s$ at $(25, 13)$ and $(55, 13)$ (`v-it52.png`) |
| centerline $\beta$ inside the loop (Fig. 6e, TENO-LAD $2400^2$) | $\approx 0.15$–$0.2$ | $\approx 0.03$–$0.1$ (`β-it52.png`) |
| $t = 54\tau_0$ (Fig. 6a,b,d, TENO-LAD $2400^2$) | $V_z$ peak $\approx 1.25\,C_s$ at $z \approx 28H_0$, $V_A$ peak $\approx 4.1\,C_s$ at $z \approx 21H_0$, $\rho$ drop at $z \approx 29$–$30H_0$ | $1.28\,C_s$ at $z \approx 32H_0$, $3.65\,C_s$ at $z \approx 21H_0$, drop at $z \approx 31H_0$; loop top at $29H_0$, inside the absorbing layer (`profile-it55.png`, `ρ-it55.png`) |

The emergence itself — timing, height, loop shape, rise speed, density
inside and above the loop, the downflows — is reproduced at this
"coarsest admissible" resolution. The departures are: the Alfvén-speed
peak 15–20 % low and the loop $\beta$ 2–5 times low (both point to the loop
interior being slightly denser in gas pressure terms / weaker in field
than the paper's; the eddy viscosity acts on $\mathbf{B}$ too); the contact
at the loop top smeared over $\approx 2H_0$ against $\approx 1H_0$ in the
paper; and, before the emergence, a transition-region disturbance the
paper's snapshots do not show — the coronal fall-back described under
"Positivity floors" leaves $\pm1\,C_s$ vertical oscillations above
$z = 18H_0$ at $t = 25$–$33\tau_0$ (`v-it26/34.png`) and hot ($T \approx 90$)
floored pockets at the transition region (`T-it26.png`). They are gone by
the time the loop passes through ($t \gtrsim 40\tau_0$) and do not visibly
alter it.

**Floors.** First hit at $t = 13.85\tau_0$ at $(x, z) = (27.5, 17.8)$, i.e.
the transition region under the perturbation; $10^4$ stage-node hits by
$t = 14.5$, $10^6$ per central rank by $t = 32\tau_0$, all at
$z \approx 17$–$18.5H_0$ (the corrugated contact). From $t \approx 50\tau_0$
they also fire, at ≈ 200 nodes per stage, in the evacuated downflow
regions at the loop sides ($(4, 16)$, $(62, 8)$) — the paper's
"near-vacuum regions with high-Mach-number flows". They never fire inside
the loop.

## Notes

- **Coordinates**: Jexpresso's `y` is the paper's `z`; `v` is $V_z$ and
  `By` is $B_z$. `w` and `Bz` must stay at machine zero for this problem —
  nothing in the 2D system couples them once they vanish; they are a free
  consistency check.
- **Units**: $H_0 = C_s = \rho_0 = 1$, $p = \rho T/\gamma$, $g_0 = 1/\gamma$,
  Heaviside–Lorentz $\mathbf{B}$ ($\tfrac12B^2$ magnetic pressure). The
  output `T` is $T/T_0 = \gamma p/\rho$, i.e. exactly 1 in the chromosphere
  and 25 in the corona at $t = 0$.
- **The CFL printout.** The advective and acoustic lines of the shared
  diagnostic in `src/kernel/physics/soundSpeed.jl` are written for
  CompEuler (γ = 1.4, pressure without the magnetic energy) and are not
  meaningful for this case — same as for the other two MHD cases (they
  report a "sound speed" of 16–35 where the true maximum is 5–10). The
  viscous line is meaningful: it uses the DynSGS coefficient actually
  applied, kinematic on every slot in the conserved form of this case
  (0.007 at $t = 0$, at most ≈ 0.09 during the run).
- **Low β, γ close to 1.** With $\gamma - 1 = 0.05$ the gas pressure is
  $(\gamma - 1)$ times the internal energy: at $\beta \approx 0.05$ inside
  the emerged loop $p$ is a 2 % residue of the total energy.
  `initialize.jl` prints the smallest initial pressure; if a run reports
  `p` at the floor of `user_flux.jl` ($10^{-9}$) in the loop, the
  regularization is too weak there — raise `:dsgs_C1`, or `:μ[4]`.
- **Constants shared by name** (`γ_mhd`, `c_h_mhd`, …) with the other MHD
  cases: `user_flux.jl` refuses to run if a γ = 5/3 case was loaded earlier
  in the same Julia session. Restart Julia.

## Files

| file | contents |
|---|---|
| `user_inputs.jl` | solver parameters (Δt, tend, mesh, DynSGS, PNG styling) |
| `initialize.jl` | atmosphere, flux sheet, perturbation, `c_h`, Δh, β_* |
| `user_flux.jl` | GLM-MHD fluxes `F`, `G`, `γ_mhd = 1.05`, `g_mhd`, `pressure_mhd` |
| `user_source.jl` | gravity, GLM ψ damping (α_p), absorbing layer |
| `user_bc.jl` | symmetric bottom, free-slip top; periodic x at mesh level |
| `user_primitives.jl` | conserved → primitive/output mapping (adds `vA`, `β`) |
| `FE_80x35.geo`, `.msh` | gmsh geometry of the $[0,80]\times[0,35]$ grid and the generated mesh (default) |
| `FE_80x70.geo`, `.msh` | the same with $0.5H_0$ tall elements (the paper's vertical spacing) |
| `EQUATIONS.md` | equations, units, initial condition, expected results |
