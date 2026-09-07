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
| elements | $80\times35$, each $1H_0\times1H_0$ (`FE_80x35.geo/.msh`, shipped) |
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
vertical spacing set `ny = 70` in `FE_80x35.geo`, regenerate the mesh
(`gmsh -2 FE_80x35.geo -o FE_70.msh`, point `:gmsh_filename` at it) and halve
`:Δt`; cost doubles.

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
  table and interpolated at the nodes. **The paper never states the
  plasma beta $\beta_*$ of the sheet**; $\beta_* = 1$ was inferred from
  its Fig. 1(b) (table in EQUATIONS.md §3) and is a `Ref`
  (`fe_beta_star`).
- **Stabilization: DynSGS** (`:visc_model => DSGS_MHD()`, see
  [DSGS.md](../../../DSGS.md) §4), with `:dsgs_gamma => 1.05` matching the
  flux, and the two stratification variants of DSGS.md §4.5 switched on:
  `:dsgs_local_norms` (per-element normalization of the residual — with
  the domain norm the $10^8$-times denser photosphere hides the corona from
  the sensor, and a grid-scale sawtooth grew across the transition region
  at $t \approx 1.5\tau_0$) and `:dsgs_nodal_rho` (dynamic coefficient with
  the nodal density — the element mean over-diffuses the light side of a
  transition-region element by a factor 25 and breaks the viscous CFL).
  The paper's shocks (fast and intermediate along the loop sides, slow
  near the footpoints) are captured by its HLLD/WENO machinery; here the
  residual-based eddy viscosity regularizes them. No filter, no
  entropy-stable/KEP fluxes.
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

## Notes

- **Coordinates**: Jexpresso's `y` is the paper's `z`; `v` is $V_z$ and
  `By` is $B_z$. `w` and `Bz` must stay at machine zero for this problem —
  nothing in the 2D system couples them once they vanish; they are a free
  consistency check.
- **Units**: $H_0 = C_s = \rho_0 = 1$, $p = \rho T/\gamma$, $g_0 = 1/\gamma$,
  Heaviside–Lorentz $\mathbf{B}$ ($\tfrac12B^2$ magnetic pressure). The
  output `T` is $T/T_0 = \gamma p/\rho$, i.e. exactly 1 in the chromosphere
  and 25 in the corona at $t = 0$.
- **Ignore the CFL printout.** The shared diagnostic in
  `src/kernel/physics/soundSpeed.jl` is written for CompEuler (γ = 1.4, no
  magnetic pressure, multiplier instead of the actual eddy viscosity) and is
  not meaningful for this case — same as for the other two MHD cases.
- **Low β.** With $\gamma - 1 = 0.05$ and $\beta \sim 10^{-4}$ inside the
  emerged loop, the pressure is a $10^{-5}$ residue of the total energy.
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
| `FE_80x35.geo` | gmsh geometry of the $[0,80]\times[0,35]$ grid |
| `FE_80x35.msh` | the generated mesh |
| `EQUATIONS.md` | equations, units, initial condition, expected results |
