# 2D ideal GLM-MHD: magnetized Kelvin-Helmholtz instability

This case implements the two-dimensional ideal GLM-MHD equations of

> J. Chan, H. Ranocha, A. M. Rueda-Ramírez, G. Gassner, T. Warburton,
> *On the entropy projection and the robustness of high order entropy stable
> discontinuous Galerkin schemes for under-resolved flows*,
> Frontiers in Physics 10:898028 (2022) — **Section 3.2, Eq. (6)**,

with the magnetized Kelvin-Helmholtz initial condition of **Section 3.2.1,
Eq. (10)** on the doubly periodic domain [-1, 1]², T_final = 15.

A full write-up of the governing equations, the GLM divergence cleaning,
the Smagorinsky stabilization and the test setup is in
[EQUATIONS.md](EQUATIONS.md); this file covers the implementation choices.

## What is (and is not) implemented

- **Conservative GLM-MHD fluxes** (Eq. 6): 9 unknowns
  `(ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)`. The 2D system keeps the third
  components of velocity and magnetic field, and drops the z-flux.
  `ψ` is the GLM divergence-cleaning field; the cleaning speed `c_h` is set
  to the maximum wave speed of the initial condition and kept constant
  (computed in `initialize.jl`, stored in the `c_h_mhd` Ref).
- **Mixed (damped) GLM cleaning**: `user_source.jl` adds Dedner's ψ-damping
  source `S_ψ = -(c_h/c_r) ψ` (`c_r = 0.18`, tunable via the `glm_cr_mhd`
  Ref). Without it, the undamped ψ waves have no dissipation mechanism in a
  continuous Galerkin discretization (no interface Riemann fluxes) and
  accumulate grid-scale noise on a periodic domain.
- **Stabilization: Smagorinsky SGS** (`:visc_model => SMAG()`), not the
  entropy-stable / KEP two-point fluxes of the paper: `:lkep => false`,
  `:entropy_variables => false`. Momentum gets the full deviatoric stress,
  the energy gets `κ∇T` plus the `τ·u` viscous work, and `ρw`, `B`, `ψ` are
  diffused as scalars (a turbulent resistivity for `B`).
- **NOT implemented (yet)**: the non-conservative term Υ of Eq. (7) — the
  Powell term `(∇·B)(0, B, v·B, v, 0)` plus the Galilean GLM term
  `(0, 0, ψ v·∇ψ, 0, v·∇ψ)`. It is only required for entropy stability of
  the split-form discretizations, which are not used here.

## Notes

- **Variable ordering**: the paper writes the state as `(ρ, ρv, E, B, ψ)`;
  here the total energy sits in **slot 4** and `ρw` in slot 5 because
  Jexpresso's shared 2D kernels (viscous τ·u augmentation, sound-speed/CFL
  diagnostic) assume slot 4 of a 2D system is the energy.
- **γ**: the paper does not state a value; `γ = 5/3` (standard ideal-MHD
  choice) is set in `user_flux.jl` (`γ_mhd`). Use 1.4 to match the companion
  Euler case `problems/CompEuler/kelvinHelmholtzChan2022`.
- The sound-speed/CFL *printout* uses the Euler assumptions of the shared
  kernel (γ = 1.4, no magnetic pressure); it is diagnostic-only.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD kelvinHelmholtzChan2022
```

The mesh is the same doubly periodic 32×32 [-1,1]² grid used by the
companion Euler KHI case (`hexa_TFI_32x32_unitsquare.msh` from
JexpressoMeshes).
