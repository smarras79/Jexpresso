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
- **Time.** CarpenterKennedy2N54 is fourth order, enough for N ≤ 3. Δt gives a Courant
  number of 0.1 on the mean node spacing √2/(nelx N), against |v| + c_f = 1.185, and is
  rounded down so that a whole number of steps reaches t = 2. On the smallest LGL spacing
  that is a Courant number of 0.12 for N = 3 and 0.14 for N = 4.

| Variable | Meaning | Default |
|---|---|---|
| `JEXPRESSO_AW_NOP` | polynomial order N | 3 |
| `JEXPRESSO_AW_NELX` | elements per side (picks `alfven_<n>x<n>.msh`) | 16 |
| `JEXPRESSO_AW_DT` | time step, instead of the rule above | |
| `JEXPRESSO_AW_TEND` | final time | 2.0 |
| `JEXPRESSO_AW_VISC` | `dsgs`, or `none` for Galerkin | `dsgs` |
| `JEXPRESSO_AW_SENSOR` | `:dsgs_sensor`, `residual` or `legacy` | `residual` |
| `JEXPRESSO_AW_CMIN` | `:dsgs_Cmin` | 0 |
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
