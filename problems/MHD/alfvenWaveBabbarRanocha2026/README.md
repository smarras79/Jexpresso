# Smooth Alfvén wave: the convergence test of Babbar & Ranocha (2026), §5.5.1

A circularly polarized Alfvén wave (Tóth, *J. Comput. Phys.* 161 (2000) 605)
crosses a doubly periodic box at 45°. The paper uses it to show that its scheme
for GLM-MHD converges at the optimal order N + 1 (its Fig. 13). The wave is an
exact nonlinear solution of ideal MHD, so the error can be measured against it.

> A. Babbar, H. Ranocha, *Compact Runge-Kutta flux reconstruction methods for
> non-conservative hyperbolic equations*, J. Comput. Phys. 563 (2026) 115060.

```bash
tools/alfven_wave_scan.sh                      # N = 3 ... 7, DSGS and Galerkin, then the figure
AW_JOBS=5 tools/alfven_wave_scan.sh            # the same, the five orders at once
julia --project=. tools/alfven_wave_plot.jl    # the figure and the table again, from the stored errors
JEXPRESSO_AW_NOP=3 JEXPRESSO_AW_NELX=32 \
    julia --project=. src/Jexpresso.jl MHD alfvenWaveBabbarRanocha2026    # one run
```

The sweep ends with `output/MHD/alfvenWaveBabbarRanocha2026/output/convergence-it0.png`
and the table of errors and orders; see [The sweep](#the-sweep).

## The problem

With α = π/4, x⊥ = x cos α + y sin α and B⊥ = 0.1 sin(2π(x⊥ + t)):

$$
\rho = 1,\qquad p = 0.1,\qquad \gamma = 5/3,\qquad \psi = 0,
$$

$$
\mathbf v = \big(-B_\perp \sin\alpha,\ B_\perp\cos\alpha,\ 0.1\cos(2\pi(x_\perp + t))\big),\qquad
\mathbf B = (\cos\alpha + v_1,\ \sin\alpha + v_2,\ v_3).
$$

The box is [0, 1/cos α] × [0, 1/sin α] = [0, √2]², periodic, and the run ends at t = 2.

- |B|² = 1.01 and p are uniform, so the total pressure is too, and the wave keeps its
  shape.
- The wave moves along −(cos α, sin α) at the Alfvén speed B∥/√ρ = 1.
- The box holds one wavelength along each side, so the solution is back to the initial
  condition at every integer time.
- Half of the wave is out of the plane (w = B_z ≠ 0), so slots 5 and 8 of the state
  (ρw, B_z) carry it as much as the in-plane ones do.

## Where each setting comes from

Each setting was checked against the code that produced the paper's figure, not only
against the paper's text: the reproducibility repository
<https://github.com/Arpit-Babbar/paper_crk_nonconservative>, commit `19a6900`,
`run/convergence_alfven_wave/`.

| | Paper's code | Here |
|---|---|---|
| exact solution | Tenkai.jl `ExactSolutionAlfvenWave` → Trixi.jl `initial_condition_convergence_test(::IdealGlmMhdEquations2D)` | `initialize.jl`, the same formula |
| γ, box, final time | 5/3, [0, √2]², 2.0 | the same |
| orders, meshes | N = 1, 2, 3 on 16², 32², 64², 128² elements | N = 3 to 7 (Jexpresso is run at N ≥ 3), each on its own meshes between 4² and 64² ([The sweep](#the-sweep)) |
| error | Tenkai `compute_error`: sqrt((1/\|Ω\|)∫(ρ_h − ρ)²dΩ), N + 4 Gauss points per direction, **density only** | the same norm, for every conserved variable |

Two things are easy to misread in the paper:

- **The third velocity component.** The printed formula reads v = v⊥(−sin α, cos α,
  0.1 cos(2π(x⊥ + t))). In the code, the third component is 0.1 cos(2π(x⊥ + t)) itself,
  not v⊥ times it.
- **The error of Fig. 13 is the error in ρ.** Tenkai's `compute_error` loops over the
  first conserved variable only. The exact ρ is 1, so that number measures how far the
  scheme moves ρ off a constant (through the nonlinear coupling), not the error of the wave
  itself. `user_plot.jl` stores both: the ρ error, which is the paper's measure, and the
  error of the whole conserved state, sqrt(Σ_k e(q_k)²), which measures the wave.

## The equations

The flux (`user_flux.jl`) and the ψ damping (`user_source.jl`) are those of
`orszagTangBormanis2024` and `smoothVortex`: conservative GLM-MHD with c_h = 1.185 (the
largest |v| + c_f of the initial condition) and S_ψ = −(c_h/0.18)ψ.

The paper's system (its eq. 81) also has the non-conservative Godunov-Powell term and the
GLM transport v·∇ψ. Both vanish on this solution (∇·B = 0, ψ = 0), so it is an exact
solution of the conservative system too. The two systems differ only in how they carry
the discretization's own ∇·B and ψ errors.

## Settings

- **Stabilization.** DynSGS is set up as for `smoothVortex`, the other accuracy test of
  this directory:
  - residual sensor, `:dsgs_Cmin => 0`, `:dsgs_CR => 1`, `:dsgs_Cmax => 0.5`;
  - domain-wide normalization;
  - physical form (`:dsgs_conserved => false`) with the split energy flux;
  - the coefficient computed once per step (`:dsgs_freeze_stage`, the default for
    `DSGS_MHD` with the residual sensor; `JEXPRESSO_AW_FREEZE=0` computes it at every
    stage).

  `JEXPRESSO_AW_VISC=none` gives the plain Galerkin run.
- **Time step.** Δt gives a Courant number of 0.1 on the mean node spacing √2/(nelx N),
  against |v| + c_f = 1.185, and is rounded down so that a whole number of steps reaches
  t = 2. On the smallest LGL spacing that is a Courant number of 0.12 for N = 3, 0.14 for
  N = 4 and 0.22 for N = 7.
- **Integrator.** CarpenterKennedy2N54 up to N = 5, Vern7 from N = 6 on. CK54 is 4th
  order: its stability polynomial ends in z⁵/200, so it shifts the phase of the wave
  (ω = 2π) by (ωΔt)⁵/300 per step, which by t = 2 is an error of about 9Δt⁴ in the state.
  - Up to N = 5 that stays below a percent of the spatial error.
  - At N = 7 on 16² it is twice the spatial error: 1.32e-11 with CK54 against 5.53e-12 with
    Vern7 or Vern9, which agree to four digits.
  - The density error, the paper's measure, does not see it at all, since ρ = 1 whatever
    the phase of the wave.
  - Vern7 costs 10 right-hand sides per step, against 5 for CK54. With DynSGS the
    integrator does not change the error: 1.0438e-6 (Vern7) against 1.0428e-6 (CK54) at
    N = 6 on 8². The frozen coefficient works with Vern7's and Vern9's stages at t + Δt:
    6.760e-8 with Vern9 against 6.750e-8 with CK54 at N = 4 on 32².

| Variable | Meaning | Default |
|---|---|---|
| `JEXPRESSO_AW_NOP` | polynomial order N | 3 |
| `JEXPRESSO_AW_NELX` | elements per side (picks `alfven_<n>x<n>.msh`) | 16 |
| `JEXPRESSO_AW_DT` | time step, instead of the rule above | |
| `JEXPRESSO_AW_TEND` | final time | 2.0 |
| `JEXPRESSO_AW_VISC` | `dsgs`, or `none` for Galerkin | `dsgs` |
| `JEXPRESSO_AW_SENSOR` | `:dsgs_sensor`, `residual` or `legacy` | `residual` |
| `JEXPRESSO_AW_CMIN` | `:dsgs_Cmin` | 0 |
| `JEXPRESSO_AW_FREEZE` | `1` or `0` forces `:dsgs_freeze_stage` (the coefficient once per step, or at every stage) | unset: `mod_inputs.jl`'s default, on |
| `JEXPRESSO_AW_SOLVER` | `ck54`, `vern7`, `vern9`, `dp8`, `ssprk54` | `ck54` up to N = 5, `vern7` above |

The paper's N = 1 and 2 are not part of the sweep, since Jexpresso is run at N ≥ 3.
(N = 1 would not build a periodic 2D mesh in any case: `add_high_order_nodes_edges!`
returns before it fills `poin_in_edge` when nop < 2.)

## What a run produces

At t = 2, a run stores the error of every conserved variable in
`errors/nop<N>_nelx<M>_<visc>.dat`, in this directory. Each record holds the L² error
above and the largest pointwise error at the Gauss points.

`<visc>` is `dsgs` for the deck as committed and `galerkin` for `JEXPRESSO_AW_VISC=none`.
Anything else a run changes is added to it, so that the record does not replace the
deck's: `perstage` (`JEXPRESSO_AW_FREEZE=0`), `legacy`, `cmin`, an integrator other than
the deck's for that N (`vern9`, ...), or a time step of its own (`dt0.00124`).

From every record stored for the same final time, the run then:

- prints a table with the observed orders;
- draws `output/convergence-it<n>.png` in the layout of the paper's Fig. 13 (L² error
  against #DOFs), with panel (a) for ρ and panel (b) for the whole state. Each order has
  its own colour and marker (filled for DynSGS, hollow and dashed for Galerkin), and a
  dotted O(h^(N+1)) line, labelled at its end.

`tools/alfven_wave_plot.jl` draws the same figure and table from the store without
running anything, as `convergence-it0.png` (`3,5,7` for a subset of the orders,
`--visc=dsgs,galerkin` for a subset of the records, `--out=DIR` elsewhere).

#DOFs is the number of distinct unknowns per variable of this continuous Galerkin
discretization on the periodic box, (nelx N)². The paper's flux reconstruction has
(nx (N + 1))² on nx² elements. On the #DOFs axis the two methods are therefore compared
at equal cost, not on equal meshes.

To start a fresh comparison, run `rm -r problems/MHD/alfvenWaveBabbarRanocha2026/errors`.

## The sweep

`tools/alfven_wave_scan.sh` runs every order on its meshes, with DynSGS and with Galerkin,
and ends with `tools/alfven_wave_plot.jl`: the figure of the whole sweep,
`output/MHD/alfvenWaveBabbarRanocha2026/output/convergence-it0.png`, and its table.

Each order gets the meshes that take its Galerkin error from the coarse-mesh range down to
about 1e-12, short of round-off, so the five curves cover about the same span of #DOFs:

| N | meshes | #DOFs | Galerkin error on the finest |
|---|---|---|---|
| 3, 4 | 8², 16², 32², 64² | 24² ... 256² | 4e-7, 2e-9 |
| 5, 6 | 4², 8², 16², 32² | 20² ... 192² | 5e-10, 4e-12 |
| 7 | 4², 8², 16² | 28² ... 112² | 6e-12 (32² would be at round-off) |

| Setting | Effect |
|---|---|
| `AW_NOPS="3 4"` | the orders (default `3 4 5 6 7`) |
| `AW_NELX="16 32 64 128"` | the meshes of every order (the paper's, here) |
| `AW_NELX_7="4 8 16 32"` | the meshes of one order |
| `AW_VISC=dsgs` | only the DynSGS runs (default `dsgs none`) |
| `AW_JOBS=5` | up to 5 orders at once, one Julia session each (default 1) |
| `AW_KEEP=1` | keep the stored errors (they are cleared first by default) |
| `JULIA=/path/to/julia` | default: the macOS app's julia if it starts, else the one on the PATH (a juliaup launcher is bypassed: the runs call the binary itself) |

The sweep as committed, measured here with `AW_JOBS=4` on four cores, took 83 min. That
is the time of the N = 4 session, whose two 64² runs take 67 of its 82 min. The five
sessions add up to 3.75 h. `AW_JOBS=1` takes a little less, about 3.5 h, because running
side by side slowed them down: alone, N = 4 on 64² with DynSGS takes 38 min instead of 48.
Three quarters of the total is the finest mesh of each order: 64² for N = 3 and 4, 32² for
N = 5 and 6.

With `AW_JOBS`, each session writes its log to
`output/MHD/alfvenWaveBabbarRanocha2026/scan-nop<N>.log`. The order is the unit of
parallelism for a reason. Two sessions on the same mesh and order would write the same
`.jexpresso_cache` file at the same time, and JLD2, which maps that file into memory, then
dies with a bus error (it happened while this sweep was being set up).

## What this test found

Measured here with `tools/alfven_wave_scan.sh` as committed, t = 2: DynSGS (the deck) and
Galerkin (`JEXPRESSO_AW_VISC=none`), CarpenterKennedy2N54 up to N = 5 and Vern7 above.
Item 5 compares the coefficient computed at every stage (`JEXPRESSO_AW_FREEZE=0`).

- "ρ" is the paper's measure.
- "state" is sqrt(Σ_k e(q_k)²) over the nine conserved variables.
- Each cell is the error, with the rate between consecutive meshes in brackets.
- The paper's values are those of its Fig. 13 for N = 3 on the same number of elements. They
  were read from the PDF's vector paths (marker centres, calibrated on the gridlines); the
  calibration puts the points on exactly 16, 32, 64, 128 elements.

| N | elements | #DOFs | ρ, DSGS | ρ, Galerkin | state, DSGS | state, Galerkin | DSGS/Galerkin (state) | paper ρ, (a) GLL g₂ | paper ρ, (b) GL Radau |
|---|---|---|---|---|---|---|---|---|---|
| 3 | 8² | 24² | 9.67e-04 | 1.22e-03 | 2.00e-03 | 1.23e-03 | 1.6 |  |  |
| 3 | 16² | 48² | 1.17e-04 (3.0) | 1.19e-04 (3.4) | 1.20e-04 (4.1) | 1.20e-04 (3.4) | 1.0 | 7.60e-05 | 2.91e-07 |
| 3 | 32² | 96² | 1.02e-05 (3.5) | 1.02e-05 (3.5) | 1.03e-05 (3.6) | 1.03e-05 (3.5) | 1.0 | 7.18e-06 | 1.06e-08 |
| 3 | 64² | 192² | 3.63e-07 (4.8) | 3.63e-07 (4.8) | 3.66e-07 (4.8) | 3.65e-07 (4.8) | 1.0 | 2.44e-07 | 5.89e-10 |
| 4 | 8² | 32² | 1.24e-05 | 6.90e-06 | 6.57e-05 | 7.98e-06 | 8.2 |  |  |
| 4 | 16² | 64² | 4.02e-07 (4.9) | 3.93e-07 (4.1) | 1.93e-06 (5.1) | 4.28e-07 (4.2) | 4.5 |  |  |
| 4 | 32² | 128² | 1.88e-08 (4.4) | 2.33e-08 (4.1) | 6.75e-08 (4.8) | 2.60e-08 (4.0) | 2.6 |  |  |
| 4 | 64² | 256² | 1.33e-09 (3.8) | 1.45e-09 (4.0) | 2.43e-09 (4.8) | 1.63e-09 (4.0) | 1.5 |  |  |
| 5 | 4² | 20² | 1.34e-04 | 1.28e-04 | 4.96e-04 | 1.29e-04 | 3.8 |  |  |
| 5 | 8² | 40² | 2.33e-06 (5.8) | 2.38e-06 (5.7) | 3.42e-06 (7.2) | 2.40e-06 (5.8) | 1.4 |  |  |
| 5 | 16² | 80² | 4.52e-08 (5.7) | 4.40e-08 (5.8) | 1.10e-07 (5.0) | 4.43e-08 (5.8) | 2.5 |  |  |
| 5 | 32² | 160² | 6.35e-10 (6.2) | 4.74e-10 (6.5) | 6.39e-09 (4.1) | 4.76e-10 (6.5) | 13.4 |  |  |
| 6 | 4² | 24² | 3.68e-06 | 1.04e-06 | 2.28e-05 | 1.08e-06 | 21.0 |  |  |
| 6 | 8² | 48² | 2.64e-07 (3.8) | 1.83e-08 (5.8) | 1.04e-06 (4.4) | 2.04e-08 (5.7) | 51.3 |  |  |
| 6 | 16² | 96² | 1.08e-08 (4.6) | 2.09e-10 (6.5) | 5.42e-08 (4.3) | 2.36e-10 (6.4) | 229.9 |  |  |
| 6 | 32² | 192² | 8.37e-10 (3.7) | 3.30e-12 (6.0) | 3.20e-09 (4.1) | 3.75e-12 (6.0) | 854.4 |  |  |
| 7 | 4² | 28² | 5.07e-07 | 4.51e-07 | 7.94e-06 | 4.54e-07 | 17.5 |  |  |
| 7 | 8² | 56² | 8.53e-08 (2.6) | 2.10e-09 (7.7) | 6.45e-07 (3.6) | 2.12e-09 (7.7) | 305.0 |  |  |
| 7 | 16² | 112² | 6.17e-09 (3.8) | 5.42e-12 (8.6) | 4.22e-08 (3.9) | 5.53e-12 (8.6) | 7626.2 |  |  |

1. **Galerkin converges at N + 1 or better for odd N, and at N for even N.**
   - Odd N, on the finest pair: 4.8 for N = 3, 6.5 for N = 5, 8.6 for N = 7.
   - Even N: 4.0 for N = 4 on every pair, and 6.0 for N = 6 on the finest pair.
   - Order N is what the standard error estimate guarantees for unstabilized Galerkin on a
     hyperbolic problem; the odd orders do better on these uniform meshes.
   - N = 7 reaches 5.5e-12 on 16², and N = 6 reaches 3.75e-12 on 32².
2. **DynSGS costs nothing at N = 3, little at N = 4, and converges at about 4 from N = 5 on.**
   - N = 3: on 32² and 64² the ρ error is the same to three digits with and without DynSGS,
     and it converges at 4.8 on the finest pair, as the paper's Fig. 13(a) does (4.9).
   - N = 4: on 64² the state error is within 1.5× of Galerkin's, at order 4.8 (item 5).
   - N = 5, 6, 7: the state error converges at 4.1, 4.1 and 3.9 on the finest pairs. There
     it is 13×, 850× and 7600× the Galerkin error (32², 32², 16²). At equal #DOFs the three
     DynSGS curves lie within 1.5× of each other from about 48² DOFs on (panel (b) of the
     figure), so with DynSGS a higher N buys no accuracy on this wave.
   - The paper's measure shows it less at N = 5 (ρ: 6.3e-10 against 4.7e-10 on 32²), but at
     N = 6 and 7 the ρ error is 250× and 1100× Galerkin's.
3. **Why order 4: the time derivative in the frozen residual.** At the stage where the
   coefficient is computed, the residual's ∂ₜq is the BDF2 of three step-level states, whose
   truncation error is (Δt²/3)∂ₜ³q. A travelling wave has ∂ₜ³q ≠ 0 however well it is
   resolved, so that is the floor of the residual, and ν ≈ C_R Δ² O(Δt²) ∝ h⁴ at a fixed
   Courant number, whatever N.
   - The floor predicts the error with no fitted constant. The residual is normalized by
     ρ̄c̄ = 0.42 (the floor of the momentum and B norms, above the wave's spread), so
     ν ≈ C_R Δ² (Δt²/3) ω³A/(ρ̄c̄), with ω = 2π and A = 0.1. The wave is damped by ν k² t
     (k = 2π), so the state error is 0.14 ν k² t. Predicted against measured:
     N = 5 on 32², 6.6e-9 against 6.4e-9; N = 6 on 16² and 32², 5.4e-8 and 3.4e-9 against
     5.4e-8 and 3.2e-9; N = 7 on 8² and 16², 4.9e-7 and 3.0e-8 against 6.5e-7 and 4.2e-8.
   - Δ = Δ_el/(N + 1) and Δt ∝ Δ_el/N, so at equal #DOFs ν hardly depends on N: that is why
     the DynSGS curves of N = 5, 6, 7 coincide. At N = 3 the same floor is a hundredth of the
     Galerkin error, so DynSGS does not show.
   - Halving Δt alone (`JEXPRESSO_AW_DT`) divides the DynSGS error by 3.9 at N = 6 and 4.0 at
     N = 7, on 16²: 5.42e-8 → 1.38e-8 and 4.22e-8 → 1.05e-8. The Galerkin error does not move
     (N = 6 on 8²: 2.036e-8 at Δt and at Δt/2).
   - On 8² the dependence is not clean (N = 6: 1.04e-6, 3.96e-7, 4.61e-7 at Δt, Δt/2, Δt/4):
     on that mesh the spatial part of the residual still counts.
   - A backward difference of higher order for ∂ₜq would lower that floor to O(Δt^k), at one
     more stored state per order; the deck does not do it.
4. **Against the paper**, at equal #DOFs rather than on equal meshes, N = 3 is about 2.7×
   below the GLL/g₂ variant of Fig. 13(a) and about 190× above the GL/Radau variant of (b).
   At N = 3, from 16² on, the ρ error is the largest error. The errors of the wave variables
   (ρu, ρv, ρw, B) are 8–16× smaller, so the state error follows ρ. On 8² with DynSGS, they
   are as large as the ρ error. On 64², DynSGS's B error is twice the Galerkin one (7.4e-9
   against 3.8e-9).
5. **At N = 4, DynSGS converges at 4.8 because its coefficient is frozen per step.**
   - With the frozen coefficient, on 64² the state error is within 1.5× of Galerkin's
     (2.4e-9 against 1.6e-9), and the ρ error is below Galerkin's (1.3e-9 against 1.4e-9).
   - Computed at every stage (`JEXPRESSO_AW_FREEZE=0`, the default before b956ccd), the
     coefficient caps the rate at about 3. `src/io/mod_inputs.jl` (`:dsgs_freeze_stage`)
     explains why: at an intermediate RK stage, the three-point time difference of the
     residual sensor uses an internal stage whose error is O(Δt). So ν ~ C_R h² O(Δt) ∝ h³
     when Δt ∝ h.
   - State error and order, N = 4:

     | N = 4 | 16² | 32² | 64² | orders |
     |---|---|---|---|---|
     | DSGS, frozen (default) | 1.94e-06 | 6.75e-08 | 2.43e-09 | 4.84, 4.80 |
     | DSGS, per stage (`JEXPRESSO_AW_FREEZE=0`) | 3.89e-06 | 4.29e-07 | 5.24e-08 | 3.18, 3.03 |
     | DSGS, per stage, `JEXPRESSO_AW_SOLVER=vern9` | 2.02e-06 | 6.88e-08 | | 4.87 |
     | Galerkin | 4.28e-07 | 2.60e-08 | 1.63e-09 | 4.04, 4.00 |

     A higher-order integrator lifts the per-stage rate too, because its internal stages are
     more accurate. The freeze removes the internal stage from the stencil.
   - The cap also showed at N = 3 per stage, in the wave variables: the B error converged at
     3.2 on the finest pair (4.3e-8 on 64²). Frozen, it converges at 4.8 (7.4e-9).
   - The legacy sensor is not frozen (it would read −½∂ₜq at the frozen stage). On this wave
     it is far more dissipative: N = 4 on 16² gives a state error of 3.5e-3, against 1.9e-6
     with the residual sensor.

## Files

| File | Contents |
|---|---|
| `initialize.jl` | the exact solution `aw_state(x, y, t, γ)`, the initial state, c_h |
| `user_flux.jl`, `user_source.jl` | conservative GLM-MHD flux, ψ damping (as `smoothVortex`) |
| `user_primitives.jl`, `user_bc.jl` | DynSGS primitives (physical form), periodic BCs (as `smoothVortex`) |
| `user_plot.jl` | the error at the final time, its store, the convergence figure and table |
| `user_inputs.jl` | the deck and its `JEXPRESSO_AW_*` overrides |
| `alfven_periodic.geo`, `alfven_<n>x<n>.msh` | the mesh definition and the meshes, n = 4 ... 128 |
| `../../../tools/alfven_wave_mesh.sh` | meshes from the `.geo` (gmsh) |
| `../../../tools/alfven_wave_scan.sh` | the sweep, N = 3 ... 7, in one Julia session or one per order |
| `../../../tools/alfven_wave_plot.jl` | the figure and the table from the stored errors, without running |
