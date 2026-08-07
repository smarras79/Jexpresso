# DynSGS — residual-based Dynamic Sub-Grid Scale stabilization in Jexpresso

This document describes the Marras–Nazarov Dynamic SGS model (`DSGS`), how it is
formulated for a general conservation law, and how it is implemented for each of
the three equation sets that use it in Jexpresso:

1. [The general conservation law](#1-the-general-conservation-law)
2. [1D CompEuler — `sod1d`, `case1`](#2-1d-compeuler--sod1d-case1)
3. [2D CompEuler θ-form — `theta_dsgs`](#3-2d-compeuler-θ-form--theta_dsgs)
4. [2D ideal GLM-MHD — `orszagTangBormanis2024`](#4-2d-ideal-glm-mhd--orszagtangbormanis2024)
5. [Code map, inputs and output](#5-code-map-inputs-and-output)
6. [Known issues](#6-known-issues)

**References**

- S. Marras, M. Nazarov, F. X. Giraldo, *Stabilized high-order Galerkin methods
  based on a parameter-free dynamic SGS model for LES*,
  J. Comput. Phys. **301** (2015) 77–101.
- M. Nazarov, J. Hoffman, *Residual-based artificial viscosity for simulation of
  turbulent compressible flow using adaptive finite element methods*,
  Int. J. Numer. Meth. Fluids **71** (2013) 339–357.
- S. Marras, M. A. Kopera, F. X. Giraldo, *Simulation of shallow-water jets with
  a unified element-based continuous/discontinuous Galerkin model with grid
  flexibility on the sphere*, Q. J. R. Meteorol. Soc. **141** (2015) 1727–1739.

---

## 1. The general conservation law

Take a system of conservation laws

$$
\frac{\partial \mathbf{q}}{\partial t} + \nabla\cdot\mathbf{F}(\mathbf{q}) = \mathbf{s}(\mathbf{q}),
\qquad \mathbf{q} = (q_1,\dots,q_{neqs})^T .
$$

A continuous-Galerkin spectral element discretization of this system carries no
numerical dissipation: there are no interface Riemann fluxes, and the weak form
is energy-neutral. Under-resolved features — shocks, sharp fronts, the small
scales of a turbulent cascade — therefore produce Gibbs oscillations that grow
until the solution leaves the realizable state space (negative density or
pressure) and the run dies.

The classical fix is an eddy-viscosity closure such as Smagorinsky,
$\mu_t = \rho\,C_s^2\Delta^2|S|$. Its defect for this purpose is that it depends
only on the *strain rate*: it cannot distinguish a well-resolved shear layer
from an under-resolved shock, so it damps both. Making it strong enough to
survive the shocks over-damps everything else.

**DynSGS instead measures how badly the discrete solution fails to satisfy the
PDE, and puts viscosity exactly there.** Define the residual of equation $i$,

$$
R_i = \frac{\partial q_i}{\partial t} + \nabla\cdot\mathbf{F}_i - s_i .
$$

For an exact solution $R_i \equiv 0$. For a discrete solution $R_i$ is small
wherever the solution is resolved and spikes by orders of magnitude at shocks
and under-resolved features. That makes $|R_i|$ a natural, *parameter-free*
sensor — no tuned constant decides where the model is active.

### 1.1 The model

Per element $e$:

$$
\boxed{\;
\mu_{res}\big|_e = C_1\,\Delta_e^2\,
\max_i \frac{\lVert R_i\rVert_{\infty,e}}{\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}},
\qquad
\mu_{max}\big|_e = C_2\,\Delta_e\,\big(\lVert\mathbf{v}\rVert + c\big)_{\infty,e},
\qquad
\mu\big|_e = \max\!\big(0,\ \min(\mu_{max},\ \mu_{res})\big)\;}
$$

with $C_1 \approx 1$, $C_2 \approx 0.5$, and $\Delta_e$ the element length scale
divided by the polynomial count, $\Delta_e = \Delta_{elem}/(N+1)$.

The three ingredients:

- **Normalization.** Each residual is divided by the *global* spread of its own
  variable, $\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$. This makes
  every ratio a frequency $[1/\mathrm{time}]$ regardless of the physical units
  or magnitude of $q_i$, so the $\max_i$ over equations compares like with like,
  and the whole model is dimensionally consistent without any tuned scale.

- **Units.** $R_i/\lVert\cdot\rVert$ has units $1/T$ and $\Delta^2$ has units
  $L^2$, so $\mu$ comes out as $L^2/T$ — a **kinematic** viscosity. The
  wave-speed cap $C_2\Delta(\lVert v\rVert+c)$ is $L \cdot L/T$, the same.
  Whether the applied coefficient must be multiplied by $\rho$ depends on which
  primitive variable the diffusion operator acts on — see §4.3.

- **The cap.** $\mu_{max}$ is the first-order upwind viscosity, the most any
  sane scheme should ever add at this resolution. Taking the min means DynSGS
  degrades gracefully to first-order-upwind-like dissipation at the strongest
  discontinuities and switches itself off in smooth regions, where $\mu_{res}$
  is tiny.

### 1.2 The residual in practice

$\partial q_i/\partial t$ is not directly available inside a Runge–Kutta stage,
so it is approximated by the second-order backward difference

$$
\frac{\partial q_i}{\partial t}\Big|^n \approx
\frac{3q_i^n - 4q_i^{n-1} + q_i^{n-2}}{2\Delta t},
$$

and the spatial part $\nabla\cdot\mathbf{F}_i - s_i$ is read off the RHS the
solver has just assembled. Since Jexpresso's `params.RHS` holds
$-(\nabla\cdot\mathbf{F} - \mathbf{s})$ in **weak form** at the point where the
viscous term is built — the mass-matrix division happens later in `rhs!` — the
strong-form residual is

$$
R_i = \frac{3q_i^n - 4q_i^{n-1} + q_i^{n-2}}{2\Delta t} - M^{-1}\,\mathrm{RHS}_i .
$$

The $M^{-1}$ is **required** for the units above to work out. See
[§6](#6-known-issues) for where the code does and does not do this.

### 1.3 Per-equation split

Marras et al. eq. (10) applies the single element coefficient $\mu|_e$ to each
equation with an equation-dependent factor: mass is normally left untouched so
the scheme stays strictly conservative, momentum takes $\mu$, and the
thermodynamic equation takes $\mu$ scaled by a turbulent Prandtl number. In
Jexpresso each slot is additionally multiplied by the user's `inputs[:μ][ieq]`
vector, so any equation's contribution can be scaled or switched off from the
case file without touching the kernel.

---

## 2. 1D CompEuler — `sod1d`, `case1`

`compute_dsgs_viscosity!(::DSGS, ::NSD_1D)` in `src/kernel/physics/SGS.jl`.

**State.** $\mathbf{q} = (\rho,\ \rho u,\ \rho E)$, `neqs = 3`, total-energy form.

**Coefficients.** $C_1 = 1$, $C_2 = 0.5$, $\gamma = 1.4$ (hardcoded in the
function).

**Element scale.** $\Delta = \Delta x_e/n_{gl}$, from `mesh.Δx`.

**Residual.** All three equations enter the max, and the strong form is used:

```julia
R1 = abs((3*q[ip,1] - 4*q1[ip,1] + q2[ip,1])*inv2Δt - Mi*rhs[ip,1])
```

**Denominators.** $\lVert q_i - \langle q_i\rangle\rVert_\infty$ plus machine
`eps` only — no physical floor. This is adequate for the shock-tube cases, whose
initial condition already has an O(1) jump in every variable, so no denominator
starts at zero.

**Wave speed.** `uTmx = max(|u| + sqrt(γ*Tl))` with
`Tl = max(e - ½u², 0)`, i.e. the specific internal energy.

**Split.** One coefficient shared by all three equations, scaled per equation:

$$
\mu_{dsgs}[e,i] = \texttt{inputs[:μ][i]}\cdot\mu\big|_e ,\qquad i=1,2,3
$$

Note `sod1d` runs with `:μ => [1.0, 1.0, 1.0]` — mass diffusion is **on** in 1D,
unlike the 2D paths.

**Output.** `μ_dsgs_pnode` feeds the 1D PNG writer, which renders the viscosity
as an extra panel alongside the solution.

---

## 3. 2D CompEuler θ-form — `theta_dsgs`

`compute_dsgs_viscosity!(::DSGS, ::NSD_2D)` in `src/kernel/physics/SGS.jl`.

**State.** $\mathbf{q} = (\rho,\ \rho u,\ \rho v,\ \rho\theta)$, `neqs = 4`.

**Coefficients.** $C_1 = 1$, $C_2 = 0.5$; $\gamma$ and $C_0$ from `PhysConst`.

**Element scale.** $\Delta = \Delta_{elem,e}/n_{gl}$, from `mesh.Δelem` (the
smallest corner-to-corner distance in the element).

**Equation of state.** The θ-form closure $p = C_0(\rho\theta)^\gamma$, hence
$c = \sqrt{\gamma p/\rho}$.

**Denominators.** Machine `eps` on all four, plus a **momentum floor**
$10^{-3}\,\rho_{avg}c_{avg}$ on the two momentum slots. This exists because the
rising-bubble initial condition is globally at rest, so
$\lVert\rho u - \langle\rho u\rangle\rVert_\infty$ starts at exactly zero; with
only `eps` to absorb it the ratio runs away and pins $\mu$ at the wave-speed cap
before any flow exists, which drives $\rho\theta$ negative in the first RK
substage.

**Residual.** Uses the **weak-form** `rhs[ip,i]` directly, *without* $M^{-1}$ —
see [§6](#6-known-issues).

**Split.** As committed:

$$
\mu_{dsgs}[e,1] = 0,\quad
\mu_{dsgs}[e,2] = 0,\quad
\mu_{dsgs}[e,3] = 0,\quad
\mu_{dsgs}[e,4] = \texttt{inputs[:μ][4]}\cdot\frac{Pr}{\gamma-1}\,\mu\big|_e
$$

The two momentum slots are zeroed by an explicit diagnostic block left in the
source (`# DIAG: momentum DSGS forced to zero to isolate whether …`), so **only
the $\rho\theta$ equation is currently stabilized** on this path. That is a
debugging state, not the intended model — Marras eq. (10) puts $\mu$ on both
momentum components.

---

## 4. 2D ideal GLM-MHD — `orszagTangBormanis2024`

`compute_dsgs_viscosity!(::DSGS_MHD, ::NSD_2D)` in `src/kernel/physics/SGS.jl`,
selected with `:visc_model => DSGS_MHD()`.

This is a separate `AbstractVT` tag rather than an extension of `DSGS()` because
the residual set, the equation of state and the wave speed all differ from the
θ-system, and because the θ-path carries the debugging state described above.

**State.** $\mathbf{q} = (\rho,\ \rho u,\ \rho v,\ E,\ \rho w,\ B_x,\ B_y,\ B_z,\ \psi)$,
`neqs = 9`. See
[`problems/MHD/orszagTangBormanis2024/EQUATIONS.md`](problems/MHD/orszagTangBormanis2024/EQUATIONS.md)
for the equation set.

**Coefficients.** From `inputs`, defaulted in `mod_inputs.jl`:
`:dsgs_C1` (1.0), `:dsgs_C2` (0.5), `:dsgs_gamma` (5/3), `:dsgs_Prt` (0.7).
$\gamma$ is a case input rather than `PhysConst.γ` because the latter is air's
1.4, not the monatomic plasma's 5/3.

**Equation of state and wave speed.**

$$
p = (\gamma-1)\Big(E - \tfrac12\rho\lVert\mathbf{v}\rVert^2 - \tfrac12\lVert\mathbf{B}\rVert^2 - \tfrac12\psi^2\Big),
\qquad
c_f \le \sqrt{\frac{\gamma p}{\rho} + \frac{\lVert\mathbf{B}\rVert^2}{\rho}}
$$

the fast magnetosonic speed bounded over all propagation directions.

### 4.1 Residual set

The max runs over the **eight genuine conservation laws** and excludes $\psi$.
The GLM field is a numerical constraint carrier, not a conserved quantity: its
residual is dominated by the Dedner damping source rather than by any
under-resolution of the flow, so including it would inject viscosity for a
reason unrelated to the sensor's purpose. $\psi$ still *receives* viscosity.

### 4.2 Denominator floors

Every field is uniform at $t=0$, and $\rho w$ and $B_z$ are identically zero for
all time in this problem, so $\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$
is exactly zero for them. Each denominator is floored at $10^{-3}$ of that
field's natural scale, built from the domain-mean state:

| slot | floor |
|---|---|
| $\rho$ | $10^{-3}\bar\rho$ |
| $\rho u,\ \rho v,\ \rho w$ | $10^{-3}\bar\rho\,\bar c$ |
| $E$ | $10^{-3}\bar\rho\,\bar c^2$ |
| $B_x, B_y, B_z$ | $10^{-3}\sqrt{\bar\rho}\,\bar c$ |

A degenerate field then contributes $0/\text{floor} = 0$ rather than
$0/\texttt{eps} = $ garbage.

### 4.3 Per-equation split and units

$\mu$ from the model is **kinematic**. `_expansion_visc!` applies
`visc_coeff·∇²(primitive)` to each equation, and `user_primitives!` hands it
$u,v,w$ and $T$ for the momentum and energy slots but the raw $B$ components for
the magnetic slots. So:

| slot | coefficient | units |
|---|---|---|
| $\rho$ | $0$ | — (mass stays conservative) |
| $\rho u,\rho v,\rho w$ | $\texttt{:μ}[i]\cdot\bar\rho\,\mu$ | dynamic |
| $E$ | $\texttt{:μ}[4]\cdot\bar\rho\,\mu\cdot\dfrac{\gamma}{(\gamma-1)Pr_t}$ | dynamic |
| $B_x,B_y,B_z$ | $\texttt{:μ}[i]\cdot\mu$ | kinematic (turbulent resistivity) |
| $\psi$ | $\texttt{:μ}[9]\cdot\mu$ | kinematic |

$\bar\rho$ is the element-mean density. The energy factor follows from slot 4's
primitive being $T = p/\rho\ (= R\,T_{phys})$: the physical flux is
$\nabla\cdot(k\nabla T_{phys})$ with $k = \mu_{dyn}c_p/Pr_t$, so rewriting in
terms of $T$ gives $k/R = \mu_{dyn}\gamma/((\gamma-1)Pr_t)$ since
$c_p = \gamma R/(\gamma-1)$.

### 4.4 Step-cadenced history

`params.qp.qnm1/qnm2` are advanced on every RK **stage** (`rhs.jl`, in
`_build_rhs!`), so they hold consecutive *stage* snapshots. A BDF2 stencil built
on them differences intermediate stage states over a full $\Delta t$ and is not
an approximation of $\partial q/\partial t$ at all.

DynSGS-MHD therefore carries its own pair, `params.dsgs_qnm1/dsgs_qnm2`, rolled
exactly once per time step by a gate in `rhs!`:

```julia
if time - params.dsgs_thist[] >= 0.999*params.Δt
    params.dsgs_qnm1 .= params.dsgs_qnm2
    params.dsgs_qnm2 .= params.uaux
    params.dsgs_thist[] = time
end
```

`time` sweeps $t + c_i\Delta t$ within a step, so the gate fires once per step
regardless of the stage layout. Both buffers are initialized to the initial
state in `params_setup.jl`, so the first residual is identically zero rather
than $3q/(2\Delta t)$. They are shaped from `size(qp.qn)`, not `(npoin, neqs)` —
`uaux` carries one extra trailing column beyond the `neqs` solution slots.

### 4.5 MPI

$\langle q_i\rangle$ and $\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$
are **domain** norms by definition, so both reductions are `MPI.Allreduce`d. A
rank-local version would make the eddy viscosity depend on the partitioning. The
cost is two small collectives per RHS call.

### 4.6 Measured effect

On the Orszag–Tang vortex at $128^2$, run to $t = 1$ (see
`problems/MHD/orszagTangBormanis2024/README.md` for the full table):

| | Smagorinsky ×8 | **DynSGS** |
|---|---|---|
| $\max\lvert B_x\rvert$ @ $t=0.55$ | 0.577 | **0.607** |
| $\max\lvert B_x\rvert$ @ $t=1$ | 0.448 | **0.476** |
| $\max\lvert B_y\rvert$ @ $t=1$ | 0.504 | **0.538** |
| $\rho$ range @ $t=1$ | [0.084, 0.364] | **[0.060, 0.360]** |
| $\min p$ @ $t=0.55$ | 4.3e-2 | 3.1e-2 |
| $\max\lvert\psi\rvert$ | 3.2e-3 | 3.4e-3 |

For reference the *under-dissipated* Smagorinsky run (`:μ = 1`), which aborts at
$t\approx0.55$, reaches $\max|B_x| = 0.648$ just before dying. DynSGS recovers
most of the magnetic field strength that the 8× multiplier was destroying while
staying stable, at the cost of a thinner (but still comfortable) pressure
margin — which is exactly what putting less viscosity in the smooth regions
should look like.

---

## 5. Code map, inputs and output

| file | contents |
|---|---|
| `src/kernel/abstractTypes.jl` | `struct DSGS`, `struct DSGS_MHD` |
| `src/kernel/physics/SGS.jl` | `compute_dsgs_viscosity!` (1D, 2D-θ, 2D-MHD), `broadcast_dsgs_to_nodes!`, the `SGS_diffusion` accessors |
| `src/kernel/operators/rhs.jl` | dispatch in `viscous_rhs_el!`, `_viscous_rhs_el_2d_dsgs!`, the step-cadenced history gate in `_build_rhs!` |
| `src/kernel/infrastructure/params_setup.jl` | `μ_dsgs`, `μ_dsgs_pnode`, `visc_coeff_dsgs`, `dsgs_qnm1/2`, `dsgs_avg/denom`, `dsgs_thist` |
| `src/io/mod_inputs.jl` | `:dsgs_C1`, `:dsgs_C2`, `:dsgs_gamma`, `:dsgs_Prt` defaults |
| `src/io/write_output.jl` | the `mu_dsgs_*` VTK fields |

**Case inputs**

```julia
:lvisc      => true,
:visc_model => DSGS(),        # 1D CompEuler / 2D CompEuler θ
:visc_model => DSGS_MHD(),    # 2D ideal GLM-MHD
:μ          => [0.0, 1.0, …], # per-equation multipliers, length neqs
:dsgs_C1    => 1.0,           # DSGS_MHD only
:dsgs_C2    => 0.5,
:dsgs_gamma => 5.0/3.0,
:dsgs_Prt   => 0.7,
```

**Output.** The per-element coefficients are broadcast to nodes by
`broadcast_dsgs_to_nodes!` and written to VTK as one point-data field per
equation, named after the solution variable each one damps:
`mu_dsgs_ρ`, `mu_dsgs_ρu`, …, `mu_dsgs_ψ`. They are piecewise constant per
element by construction, and shared (DSS) nodes take the value of the last
element that writes them — fine for visualization, not a nodal field.

⚠ The slots are **not all in the same units** (see §4.3): momentum and energy
carry the dynamic $\bar\rho\mu$, the magnetic and $\psi$ slots the kinematic
$\mu$. Compare a slot against itself over time, not against a different slot.

In 1D the same data also drives an extra panel in the PNG writer.

---

## 6. Known issues

These are pre-existing defects in the `DSGS()` (Euler) paths, recorded here
because they matter to anyone reading or extending the model. The `DSGS_MHD()`
path does not share them.

1. **The 2D θ-path drops $M^{-1}$ from the residual.** It uses the weak-form
   `rhs[ip,i]` directly, which is dimensionally inconsistent with the
   $\partial q/\partial t$ term it is subtracted from — the two differ by a
   factor of the mass matrix. A source comment records that the
   dimensionally-correct $M^{-1}\cdot$RHS "shrinks the residual by ~10³ on 2D
   atmospheric meshes and effectively turns DSGS off", and the inconsistent
   form was kept for that reason. The 1D path uses $M^{-1}$ correctly. If the
   correct residual really does switch the model off, the coefficient $C_1$ is
   the thing to raise, not the units.

2. **The 2D θ-path has its momentum slots zeroed by a leftover diagnostic**
   (`# DIAG: momentum DSGS forced to zero …`), so only $\rho\theta$ is
   stabilized there. The source comment asks for
   `visc_coeff[2:3] * μ` to be restored "once the bad actor is identified".

3. **`SGS_diffusion(::DSGS, ::NSD_2D)` defines only the `inputs` signature.**
   The generic 2D `_expansion_visc!` calls `SGS_diffusion` two different ways —
   `(…, PhysConst, Δ2, VT, SD)` from the momentum and scalar branches, and
   `(…, PhysConst, Δ2, inputs, VT, SD)` from the τ·u viscous-work term — so the
   shorter call has no matching method.

   Verified against a live run of `problems/CompEuler/theta_dsgs`: it raises

   ```
   MethodError: no method matching SGS_diffusion(::Vector{Float64}, ::Int64,
     ::Float64, ::Float64, ::Float64, ::Float64, ::Float64,
     ::PhysicalConst{Float64}, ::Float64, ::DSGS, ::NSD_2D)
   ```

   at `rhs.jl:2207` — the *"other scalars"* branch, reached by `ieq = 1`
   (density) on the very first RHS call. So the 2D `DSGS()` path does not fail
   partway through; **it cannot start at all**, and `theta_dsgs` has evidently
   not been run since the calling convention diverged.

   (`DSGS_MHD` defines both signatures; `SMAG`/`VREM` in 2D are unaffected
   because they go through the separate cache-reading `_expansion_visc!` at
   `rhs.jl:2249`, which calls the `(…, ρ, ip, sgs, ltheta_eqn, SD)` accessor
   instead.)

4. **The BDF2 history is stage-cadenced for the `DSGS()` paths.** As described
   in §4.4, `qp.qnm1/qnm2` advance every RK stage, so the 1D and 2D-θ residuals
   difference stage snapshots over a full $\Delta t$. `DSGS_MHD` avoids this
   with its own buffers; the Euler paths would need the same treatment.

5. **1D wave speed.** The cap uses `sqrt(γ·Tl)` with `Tl` the specific internal
   energy $e - \tfrac12u^2$. For a perfect gas $a^2 = \gamma(\gamma-1)e_{int}$,
   so unless the case is nondimensionalized with $c_v = 1$ this overestimates
   the sound speed by $1/\sqrt{\gamma-1}$ (≈ 1.58 at $\gamma = 1.4$). The effect
   is a looser cap, so $\mu_{res}$ governs more often than it should.
