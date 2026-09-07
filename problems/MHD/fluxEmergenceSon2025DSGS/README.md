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
colour scale, `profile-it<n>.png` on the axes of the paper's Fig. 5, …) plus
`log10_μ_dsgs_ρ-it<n>.png`, the DynSGS coefficient actually applied
($\log_{10}$ of the kinematic $\mu$, floored at $10^{-6}$; all nine slots
carry the same $\mu$ in the conserved form). This is the only dissipation in
the run, so that panel is the whole story of where and how much the scheme
regularizes.

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

## Results (4 ranks, run to $54\tau_0$)

The run completed without an abort, with nothing clipped anywhere: the
density and pressure stayed positive through the transition-region
disturbance of $t \approx 13$–$18\tau_0$ that made the sibling fire its floors
$10^6$ times, and the temperature never left its initial range $[1, 25]$
during that phase (the sibling's floored pockets reached $T \approx 90$).
The operation count per step is the sibling's (one extra multiply per
quadrature point); the measured wall clock, 2012 s for the 21,600 steps on
four cores of the development container, is not comparable with the
sibling's 5194 s, which was taken in an earlier, slower session of the same
container.

**Against the sibling and the paper** (same comparison points as the
sibling's README; the paper's numbers are IMWENO-P at $300^2$ unless noted):

| quantity | paper | sibling (floors) | this case (DynSGS only) |
|---|---|---|---|
| max $\lvert\mathbf{v}\rvert$ at $t = 15$ / $18\tau_0$ | — | $2.06$ / $1.61\,C_s$ | $0.97$ / $0.82\,C_s$ |
| crest at $t = 30$–$33\tau_0$ (Fig. 14b) | $\approx 10H_0$ | $11.5H_0$ ($t = 33$) | $10.5H_0$ ($t = 30$) |
| loop top at $t = 47\tau_0$ (Fig. 14d) | $\approx 21H_0$ | $21H_0$ | $21H_0$ |
| loop at $t = 51\tau_0$ (Fig. 2) | top $\approx 25H_0$, legs at $x \approx 25$, $55H_0$, pockets at $z \approx 7H_0$ | the same | the same |
| centerline $V_z$ peak, $t = 51\tau_0$ (Fig. 5c) | $\approx 1.2\,C_s$ at $z \approx 26H_0$ | $1.15\,C_s$ at $27H_0$ | $1.23\,C_s$ at $27H_0$ |
| centerline $V_A$ peak, $t = 51\tau_0$ (Fig. 5g) | $\approx 3.9\,C_s$ at $z \approx 21H_0$ | $3.2\,C_s$ | $3.1\,C_s$ at $20H_0$ |
| centerline $\log_{10}\rho$, $t = 51\tau_0$ (Fig. 5o) | $-5$ at $z = 20$, drop to $-8$ at $27$–$28H_0$ | $-5$, drop over $28$–$30H_0$ | $-5$, drop over $28$–$30H_0$ |
| downflows at $t = 51\tau_0$ | $4$–$5\,C_s$ at the loop sides, $2$–$3\,C_s$ near the footpoints | $5$ / $3\,C_s$ | $3$ / $3\,C_s$ |
| centerline $\beta$ in the loop (Fig. 6e) | $\approx 0.15$–$0.2$ | $0.03$–$0.1$ | $0.03$–$0.1$ |
| $t = 54\tau_0$ | $V_z \approx 1.25\,C_s$ at $z \approx 28H_0$ | $1.28\,C_s$ at $32H_0$ | $1.57\,C_s$ at $32H_0$, $V_A$ $3.45\,C_s$ at $21H_0$ |

The emergence is the same as the sibling's and the paper's — timing,
height, loop shape, rise speed, density inside and above the loop. The
pre-emergence phase is cleaner than the sibling's: with the relative
operator the transition-region corrugation of $t = 13$–$18\tau_0$ diffuses
without evacuating, so there are no floored pockets, no $T \approx 90$ spots,
and the coronal fall-back peaks at $1\,C_s$ instead of $2\,C_s$. The
lateral downflows at $t = 51$ are weaker than the sibling's ($3$ against
$5\,C_s$; the paper has $4$–$5$).

**The coefficient.** `log10_μ_dsgs_ρ-it<n>.png` shows where the run's only
dissipation acts. At rest it is the $C_0$ floor, $0.03\Delta c_f$:
$10^{-1.4} = 0.037$ in the corona, $10^{-2.1} = 0.008$ in the chromosphere,
slightly more in the sheet. The residual sensor lifts it above the floor
in three places only: the transition-region band $z = 17$–$19H_0$ under
the perturbation ($0.05$–$0.12$ at $t = 10$, up to $0.3$ at $t = 33$ where
the contact is displaced), the coronal fall-back region above it at
$t = 15$ ($0.1$, $z = 18$–$27H_0$), and the absorbing layer. Inside the
rising sheet and loop it stays at $0.01$–$0.03$, and the wave-speed cap
($C_2\Delta c_f \approx 0.6$ in the corona) is never reached before the
emergence. Two consequences worth knowing: the coronal smoothing seen in
the $\beta$ maps (blurred edges of the low-$\beta$ region above $z = 18H_0$)
is the floor, which alone spreads a structure by $\sqrt{2\mu t} \approx 1.6H_0$
by $t = 33\tau_0$, in this case and the sibling alike; and the relative
operator, written out, is $\nabla\cdot(\mu\nabla\delta q) - \nabla\cdot(\mu\,\delta q\,\nabla\ln\rho_e)$,
i.e. a Laplacian plus a drift of the departure toward higher $\rho_e$ at
speed $\mu\lvert\mathrm{d}\ln\rho_e/\mathrm{d}z\rvert$ — $5.4\mu$ at the
transition region, so $0.2\,C_s$ at the floor and $1.6\,C_s$ where the
sensor fires there. That drift is the positivity mechanism; it also acts
on loop mass crossing the fixed height $z_{cor}$, a bias of $0.1$–$0.2\,C_s$
against a rise speed of $0.5$–$0.8\,C_s$ (the measured crest heights agree
with the sibling's within $1H_0$, so it is modest, but it is there).
During the emergence ($t = 51\tau_0$, `log10_μ_dsgs_ρ-it52.png`) the
sensor fires where the paper places its shocks: along the oblique loop-side
downflow fronts from $(x, z) \approx (5, 20)$ up to $(30, 28)H_0$ and their
mirror ($0.1$–$0.3$), at the loop top ($z \approx 30H_0$, $0.2$), and at the
transition-region band $z = 17$–$19H_0$ outside the loop where the lateral
downflows land, the only place the cap ($0.5$–$0.8$) is reached; the loop
interior stays at $0.03$–$0.05$ and the chromosphere at the floor. The
rerun that produced these panels reproduced the first run (same crest
heights and centerline values to the plotted precision) in 1842 s.
