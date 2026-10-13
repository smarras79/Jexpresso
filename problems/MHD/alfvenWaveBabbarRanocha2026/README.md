# Smooth Alfvén wave: the convergence test of Babbar & Ranocha (2026), §5.5.1

A circularly polarized Alfvén wave (Tóth, *J. Comput. Phys.* 161 (2000) 605)
crosses a doubly periodic box at 45°. The paper uses it to show that its scheme
for GLM-MHD converges at the optimal order N + 1 (its Fig. 13). The wave is an
exact nonlinear solution of ideal MHD, so the error can be measured against it.

> A. Babbar, H. Ranocha, *Compact Runge-Kutta flux reconstruction methods for
> non-conservative hyperbolic equations*, J. Comput. Phys. 563 (2026) 115060.

```bash
tools/alfven_wave_scan.sh                      # N = 3, 4 on 16² ... 128² elements, DSGS and Galerkin
JEXPRESSO_AW_NOP=3 JEXPRESSO_AW_NELX=32 \
    julia --project=. src/Jexpresso.jl MHD alfvenWaveBabbarRanocha2026    # one run
```

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
| orders, meshes | N = 1, 2, 3 on 16², 32², 64², 128² elements | N = 3 and 4 (Jexpresso is run at N ≥ 3), the same meshes plus 8² |
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
  - physical form (`:dsgs_conserved => false`) with the split energy flux.

  `JEXPRESSO_AW_VISC=none` gives the plain Galerkin run.
- **Time.** CarpenterKennedy2N54. Δt gives a Courant number of 0.1 on the mean node
  spacing √2/(nelx N), against |v| + c_f = 1.185, and is rounded down so that a whole number
  of steps reaches t = 2. On the smallest LGL spacing that is a Courant number of 0.12 for
  N = 3 and 0.14 for N = 4. The integrator's own error is far below the spatial one here
  (ωΔt ≈ 6·10⁻³ at N = 4 on 32²). The integrator matters through the DynSGS coefficient
  instead; see the results below.

| Variable | Meaning | Default |
|---|---|---|
| `JEXPRESSO_AW_NOP` | polynomial order N | 3 |
| `JEXPRESSO_AW_NELX` | elements per side (picks `alfven_<n>x<n>.msh`) | 16 |
| `JEXPRESSO_AW_DT` | time step, instead of the rule above | |
| `JEXPRESSO_AW_TEND` | final time | 2.0 |
| `JEXPRESSO_AW_VISC` | `dsgs`, or `none` for Galerkin | `dsgs` |
| `JEXPRESSO_AW_SENSOR` | `:dsgs_sensor`, `residual` or `legacy` | `residual` |
| `JEXPRESSO_AW_CMIN` | `:dsgs_Cmin` | 0 |
| `JEXPRESSO_AW_FREEZE` | `1`: `:dsgs_freeze_stage`, the coefficient once per step (see below) | 0 |
| `JEXPRESSO_AW_SOLVER` | `ck54`, `ssprk54`, `vern9`, `dp8` | `ck54` |

The paper's N = 1 and 2 are not part of the sweep, since Jexpresso is run at N ≥ 3.
(N = 1 would not build a periodic 2D mesh in any case: `add_high_order_nodes_edges!`
returns before it fills `poin_in_edge` when nop < 2.)

## What a run produces

At t = 2, a run stores the error of every conserved variable in
`errors/nop<N>_nelx<M>_<visc>.dat`, in this directory. Each record holds the L² error
above and the largest pointwise error at the Gauss points.

From every record stored for the same final time, the run then:

- prints a table with the observed orders;
- draws `output/convergence-it<n>.png` in the layout of the paper's Fig. 13 (L² error
  against #DOFs), with panel (a) for ρ and panel (b) for the whole state.

#DOFs is the number of distinct unknowns per variable of this continuous Galerkin
discretization on the periodic box, (nelx N)². The paper's flux reconstruction has
(nx (N + 1))² on nx² elements. On the #DOFs axis the two methods are therefore compared
at equal cost, not on equal meshes.

To start a fresh comparison, run `rm -r problems/MHD/alfvenWaveBabbarRanocha2026/errors`.

## What this test found

Measured here, t = 2, with the deck as committed (DynSGS per stage, CarpenterKennedy2N54)
and with `JEXPRESSO_AW_VISC=none` (Galerkin):

- "ρ" is the paper's measure.
- "state" is sqrt(Σ_k e(q_k)²) over the nine conserved variables.
- "order" is the rate between consecutive meshes.
- The paper's values are those of its Fig. 13 for N = 3 on the same number of elements. They
  were read from the PDF's vector paths (marker centres, calibrated on the gridlines); the
  calibration puts the points on exactly 16, 32, 64, 128 elements.

| N | elements | #DOFs | ρ, DSGS | order | ρ, Galerkin | order | state, DSGS | order | state, Galerkin | order | paper (a) GLL, g₂ | paper (b) GL, Radau |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 3 | 8² | 24² | 9.67e-04 |  | 1.22e-03 |  | 1.99e-03 |  | 1.23e-03 |  |  |  |
| 3 | 16² | 48² | 1.17e-04 | 3.04 | 1.19e-04 | 3.36 | 1.21e-04 | 4.04 | 1.20e-04 | 3.36 | 7.60e-05 | 2.91e-07 |
| 3 | 32² | 96² | 1.02e-05 | 3.53 | 1.02e-05 | 3.55 | 1.03e-05 | 3.55 | 1.03e-05 | 3.55 | 7.18e-06 | 1.06e-08 |
| 3 | 64² | 192² | 3.63e-07 | 4.81 | 3.63e-07 | 4.81 | 3.82e-07 | 4.75 | 3.65e-07 | 4.81 | 2.44e-07 | 5.89e-10 |
| 4 | 8² | 32² | 1.08e-05 |  | 6.90e-06 |  | 6.82e-05 |  | 7.98e-06 |  |  |  |
| 4 | 16² | 64² | 4.20e-07 | 4.69 | 3.93e-07 | 4.13 | 3.89e-06 | 4.13 | 4.28e-07 | 4.22 |  |  |
| 4 | 32² | 128² | 4.06e-08 | 3.37 | 2.33e-08 | 4.07 | 4.29e-07 | 3.18 | 2.60e-08 | 4.04 |  |  |
| 4 | 64² | 256² | 4.12e-09 | 3.30 | 1.45e-09 | 4.01 | 5.24e-08 | 3.03 | 1.63e-09 | 4.00 |  |  |

1. **At N = 3, DynSGS costs nothing in the paper's measure.**
   - On 32² and 64², the ρ error is the same to three digits with and without DynSGS.
   - On the finest pair it converges at order 4.8, as the paper's Fig. 13(a) does (4.9).
   - The ρ error is the largest error. The errors of the wave variables (ρu, ρv, ρw, B) are
     at least 6× smaller (7–8× with DynSGS on 64²), so the state error follows ρ.
   - DynSGS raises the wave-variable errors 2.4–11.5× above Galerkin on 64², which is still
     below the ρ error.
2. **Against the paper**, at equal #DOFs rather than on equal meshes, N = 3 is about 2.7×
   below the GLL/g₂ variant of Fig. 13(a) and about 190× above the GL/Radau variant of (b).
3. **Galerkin converges at order N = 4 for N = 4** (4.0 on every pair). **For N = 3 it
   reaches 4.8 on the finest pair**, after 3.4 and 3.6. Order N is what the standard error
   estimate guarantees for unstabilized Galerkin on a hyperbolic problem. N = 3 does better
   on these uniform meshes, and N = 4 does not.
4. **At N = 4, the per-stage DynSGS coefficient caps the rate at about 3.**
   - The wave-variable errors are 19–148× the Galerkin ones, and the gap widens with
     refinement. The state error converges at 3.2, then 3.0.
   - `src/io/mod_inputs.jl` (`:dsgs_freeze_stage`) explains why. At an intermediate RK stage,
     the three-point time difference of the residual sensor uses an internal stage whose
     error is O(Δt). So ν ~ C_R h² O(Δt) ∝ h³ when Δt ∝ h.
   - Measured on 16² → 32², state error and order:

     | N = 4, DynSGS | 16² | 32² | order |
     |---|---|---|---|
     | per stage, CarpenterKennedy2N54 (default) | 3.89e-06 | 4.29e-07 | 3.18 |
     | `JEXPRESSO_AW_FREEZE=1`, CarpenterKennedy2N54 | 1.94e-06 | 6.75e-08 | 4.84 |
     | per stage, `JEXPRESSO_AW_SOLVER=vern9` | 2.02e-06 | 6.88e-08 | 4.87 |
     | Galerkin | 4.28e-07 | 2.60e-08 | 4.04 |

   Freezing the coefficient at the tⁿ stage restores the rate, as does a higher-order
   integrator. With the freeze, the 32² error is 6× smaller, and DynSGS's ρ error falls
   below the Galerkin one (1.9e-8 against 2.3e-8).
   - At N = 3 the cap shows only in the wave variables: the B error converges at 3.2 on the
     finest pair, against 3.8 for Galerkin. ρ, which DynSGS does not change, dominates the
     state error, so neither ρ nor the state shows the cap.

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
| `../../../tools/alfven_wave_scan.sh` | the sweep, in one Julia session |
