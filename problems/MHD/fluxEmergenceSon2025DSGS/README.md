# 2D ideal GLM-MHD: flux emergence (Son, Jang & Magara 2025) — DynSGS only, no limiter

The same problem as [`fluxEmergenceSon2025`](../fluxEmergenceSon2025/README.md)
— identical physics, mesh, initial and boundary conditions, time step, output
and DynSGS constants — with one difference: **nothing but the Marras–Nazarov
DynSGS dissipation keeps the solution admissible.** The sibling case keeps
its density and pressure positive with a nodal floor applied at every
Runge–Kutta stage (`fe_positivity_limiter!`), which is not conservative and
is not part of the discretization. Here that limiter is gone: the time
integrator is the plain `CarpenterKennedy2N54()`, no variable is clipped,
no mass or energy is added anywhere, and every dissipative term is in
divergence form.

Read the sibling's README and EQUATIONS.md for the problem itself; this
file only covers what differs.

## Run

```bash
mpiexec -n 10 julia --project=. src/Jexpresso.jl MHD fluxEmergenceSon2025DSGS
```

Output goes to `./output/MHD/fluxEmergenceSon2025DSGS/<run>/`, the same PNG
set as the sibling (`ρ-it<n>.png` with field lines and vectors on the paper's
colour scale, `profile-it<n>.png` on the axes of the paper's Fig. 5, …).

## What changes, and why it is enough

The sibling's DynSGS operator is a Laplacian on the departure of the
conserved variables from the magnetostatic reference state, $\nabla\cdot(\mu\nabla(q - q_e))$.
It vanishes at rest and is conservative, but across the 25× reference jump
of the chromosphere–corona transition it does not obey a maximum
principle: a departure of $-10\%$ on the dense side ($\rho_e = 2\times10^{-8}$)
is $-2\times10^{-9}$ in absolute terms, more than the entire density of the
light side ($8\times10^{-10}$), and diffusing it across the contact drives
the light side negative. That is the mechanism by which the sibling, with
$\mu$ already at its wave-speed cap, still evacuated the coronal foot of the
transition region at $t \approx 14\tau_0$ and came to need floors. A
first-order (Lax–Friedrichs-strength) viscosity is positivity-preserving
only when it acts on a variable whose admissible set is preserved by convex
combinations, and $\rho - \rho_e$ across a jump of $\rho_e$ is not such a
variable.

This case therefore runs the operator on the **relative** departure with
the reference density as weight (`:dsgs_ref_weight => true`,
[DSGS.md §4.5](../../../DSGS.md)):

$$
\partial_t q_i = \dots + \nabla\cdot\!\left(\mu\,\rho_e\,\nabla\frac{q_i - q_{e,i}}{\rho_e}\right),
\qquad i = \rho,\ \rho u,\ \rho v,\ E,\ \rho w,
$$

and the unchanged $\nabla\cdot(\mu\nabla(q_i - q_{e,i}))$ on $B_x, B_y, B_z, \psi$.
Properties:

- **zero at rest** ($q = q_e$), exactly as the sibling's, so the
  magnetostatic atmosphere and the sheet are not eroded;
- **conservative**: divergence form, and the well-balanced correction and
  boundary treatment of the sibling apply unchanged;
- **maximum principle across the reference jump**: $r = \rho/\rho_e$ is
  diffused with a positive weight, so it is pulled toward its neighbours'
  values and no further; $\rho = r\rho_e$ stays positive wherever the
  neighbours are. Momentum and energy take the same weight, so the fluid
  slots are diffused as one scaled state and the pressure of a diffused
  node is that of a convex combination of admissible states. (The usual
  SEM caveat applies: the high-order stiffness matrix is not an M-matrix,
  so this is a maximum principle up to its positive off-diagonal entries,
  not a theorem.)
- **the same shock capturing in the loop**: where $q_e$ is negligible
  ($\rho_e \approx 10^{-8}$ in the corona against $10^{-5}$ inside the
  emerged loop) the operator is the conserved-variable Laplacian
  $\nabla\cdot(\mu\nabla q)$ with the same residual-based $\mu$, same
  $C_1$, $C_2$, $C_0$ and local norms as the sibling.

The coefficient $\mu$ itself is the sibling's — residual sensor, wave-speed
cap $C_2\Delta(\lVert\mathbf{v}\rVert + c_f)$, background floor $C_0$ —
untouched. No state-dependent trigger, no positivity sensor, no extra term
of any kind was added to it.

## Implementation

| where | what |
|---|---|
| `user_primitives.jl` | slots 1–5: $(q - q_e)/\rho_e$; slots 6–9: $q - q_e$; slot 10 (the spare `neqs+1` slot of `uprimitive`): $\rho_e$ |
| `src/kernel/operators/rhs.jl`, `_expansion_visc!` (2D CG) | multiplies the coefficient of slots 1–5 by slot 10 at the quadrature point when `dsgs_ref_weight[]` is set |
| `src/kernel/physics/SGS.jl` | the `dsgs_ref_weight` Ref (set by rhs.jl from `:dsgs_ref_weight`) |
| `user_inputs.jl` | the sibling's, minus the limiter, plus `:dsgs_ref_weight => true` |
| `initialize.jl`, `user_flux.jl`, `user_source.jl`, `user_bc.jl` | one-line files that `include` the sibling's |

The mesh is the sibling's `FE_80x35.msh` (shared, not copied).

## Results

(filled in from the validation run below)
