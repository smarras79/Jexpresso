# DynSGS — residual-based Dynamic Sub-Grid Scale stabilization in Jexpresso

This document describes the Marras–Nazarov Dynamic SGS model (`DSGS`), how it is
formulated for a general conservation law, and how it is implemented for each of
the four equation sets that use it in Jexpresso:

1. [The general conservation law](#1-the-general-conservation-law)
2. [1D CompEuler — `sod1d`, `case1`](#2-1d-compeuler--sod1d-case1)
3. [2D CompEuler θ-form — `theta_dsgs`](#3-2d-compeuler-θ-form--theta_dsgs)
4. [2D ideal GLM-MHD — `orszagTangBormanis2024`](#4-2d-ideal-glm-mhd--orszagtangbormanis2024)
5. [Shallow water on the spherical shell — `SWsphereDSGS`](#5-shallow-water-on-the-spherical-shell--swspheredsgs)
6. [Code map, inputs and output](#6-code-map-inputs-and-output)
7. [Defects found and fixed](#7-defects-found-and-fixed)

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

The $M^{-1}$ is **required** for the units above to work out; all three paths
apply it (see [§7](#7-defects-found-and-fixed) — the 2D θ-path did not until
recently).

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

**Wave speed.** `uTmx = max(|u| + sqrt(γ(γ-1)·e_int))` with
`e_int = max(e - ½u², 0)` the specific internal energy. For a perfect gas
$p = (\gamma-1)\rho e_{int}$, hence $a^2 = \gamma p/\rho = \gamma(\gamma-1)e_{int}$.

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

**Residual.** $R_i = (3q^n-4q^{n-1}+q^{n-2})/(2\Delta t) - M^{-1}\mathrm{RHS}_i$, as
in §1.2. All four equations enter the max.

**Split.** `user_primitives!` hands `_expansion_visc!` the set
$(\rho, u, v, \theta)$, so the momentum and $\theta$ slots both take the
**dynamic** coefficient $\bar\rho\,\mu$ ($\bar\rho$ = element-mean density):

$$
\mu_{dsgs}[e,1] = 0,\quad
\mu_{dsgs}[e,2] = \texttt{:μ}[2]\,\bar\rho\mu,\quad
\mu_{dsgs}[e,3] = \texttt{:μ}[3]\,\bar\rho\mu,\quad
\mu_{dsgs}[e,4] = \texttt{:μ}[4]\,\frac{Pr}{\gamma-1}\,\bar\rho\mu
$$

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

## 5. Shallow water on the spherical shell — `SWsphereDSGS`

`compute_dsgs_viscosity!(::DSGS_SW, ::NSD_2D)` in `src/kernel/physics/SGS.jl`,
selected with `:visc_model => DSGS_SW()`. The case is
[`problems/ShallowWater/SWsphereDSGS`](problems/ShallowWater/SWsphereDSGS/README.md),
whose README carries the measurements and the deck; this section is the model.

**State.** $\mathbf{q} = (\varphi,\ \varphi u,\ \varphi v,\ \varphi w)$,
`neqs = 4`, with $\varphi = gh$ the geopotential and $(u,v,w)$ the full
**Cartesian** velocity of the Marras/Kopera/Giraldo Eq. (8) formulation.

**Coefficients.** `:dsgs_C1` (1.0), `:dsgs_C2` (0.5). `:dsgs_gamma` and
`:dsgs_Prt` are ignored — this system has neither an equation of state nor an
energy equation for them to act on.

**Element scale.** $\Delta = \Delta_{elem,e}/n_{gl}$ with $\Delta_{elem}$ the
shorter of the two element chords, built by `build_sphere_element_size` in
`sphere_rhs.jl`. `mesh.Δelem` is not reused: it is the flat routine's
`min(dx, dy)` over Cartesian axis extents, which on a manifold measures the
projection of an element onto a coordinate plane rather than its size.

**Wave speed.** The **gravity** wave, $c = \sqrt{\varphi} = \sqrt{gh}$ — the
same speed `sphere_cfl_dt` uses — so $\mu_{max} = C_2\Delta(\lVert u\rVert +
\sqrt{\varphi})$. There is no sound speed here.

**Residual.** $R_i = (3q^n-4q^{n-1}+q^{n-2})/(2\Delta t) - M^{-1}\mathrm{RHS}_i$,
as in §1.2, with all four equations in the max. Unlike $\psi$ in the GLM-MHD
system, $\varphi w$ is not a constraint carrier: the three Cartesian momentum
components are one vector equation and none of them is distinguished. What *is*
special about the normal direction is removed by the Lagrange projection, not by
this model.

### 5.1 One shared denominator for the three momentum slots

This is the one structural difference from the other three paths, and it is not
cosmetic.

On a sphere the split of one tangent momentum vector into three Cartesian
components is a property of the **frame**. A purely zonal jet — the Galewsky
test, and most classic shallow-water test cases — has
$\mathbf{u} = u\,\mathbf{e}_\lambda$ with
$\mathbf{e}_\lambda = (-\sin\lambda, \cos\lambda, 0)$, so $\varphi w \equiv 0$
everywhere and for all time. Its component-wise
$\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$ is exactly zero.

Its residual is not: the $\varphi w$ equation carries a large Coriolis source,
$-f(\mathbf{x}\times\varphi\mathbf{u})|_z \approx -2\Omega\sin\phi\,r\cos\phi\,\varphi u$,
that the flux divergence cancels only to discretization error. Normalizing that
by a floor inflates the ratio by the ratio of the floor to the real momentum
scale — measured on this case, **~260×** — so $\max_i$ is won by $\varphi w$ on
every element, $\mu_{res}$ saturates the cap everywhere, and the model
degenerates into uniform first-order-upwind dissipation. On the shipped grid
that also sits at the explicit diffusive stability limit, and the run dies at
$t \approx 1100$ s.

The three momentum slots therefore share

$$
\lVert \varphi\mathbf{u} - \langle\varphi\mathbf{u}\rangle\rVert_{\infty,\Omega}
= \max_x \big|(\varphi\mathbf{u} - \langle\varphi\mathbf{u}\rangle)(x)\big|
$$

the Euclidean norm of the deviation of the one vector they are components of,
which is frame-independent. With it, $\nu$ sits at ~1.5% of the cap and the
Galewsky run completes six days.

### 5.2 Floors

| slot | floor |
|---|---|
| $\varphi$ | $10^{-3}\bar\varphi$ |
| $\varphi u,\varphi v,\varphi w$ | $10^{-3}\bar\varphi\,\bar c$, $\bar c = \sqrt{\bar\varphi}$ |

Both are non-binding on the Galewsky jet; they exist so a shell case started
from rest (a resting-balance test, a bump on a flat ocean) does not divide a
zero residual by machine `eps` in its first stage.

### 5.3 Per-equation split, and no ρ factor

$\mu$ from the model is **kinematic**, and on the shell it is applied by
`_sphere_visc_el!` directly to the **conservative** variables:
$\nu\nabla_s^2(\varphi\mathbf{u})$, which is literally the
$\delta\nu\nabla^2(\varphi\mathbf{u})$ of Eq. (8b). So, unlike §3 and §4, there
is **no** $\bar\rho$ factor — the Euler-θ and MHD paths need one only because
their `user_primitives!` hands the viscous operator $(u,v,\theta)$ instead.

$$
\mu_{dsgs}[e,i] = \texttt{:μ}[i]\cdot\mu\big|_e,\qquad i = 1\dots4
$$

with `:μ` here a **dimensionless** per-equation multiplier (see §6, "Case inputs"), built from
the deck by `build_sphere_viscosity` out of `:μ` and `:ivisc_equations`.
`[2,3,4]` reproduces the paper's placement, momentum only; `[1,2,3,4]` also
damps $\varphi$, which is what a shell run with **no modal filter** needs, since
continuous Galerkin gives $\varphi$ no upwinding and hence no dissipation of its
own.

### 5.4 Assembly order

`sphere_rhs!` assembles the inviscid RHS over **all** elements, sizes $\nu$,
then runs a **second** element loop for the viscous term. The order is forced:
the residual is $\partial q/\partial t - M^{-1}\mathrm{RHS}_{inviscid}$, and RHS
is the complete inviscid right-hand side only after every element has
contributed through direct stiffness summation. Folding the viscous term into
the first loop — which the constant-$\nu$ path does, correctly, because its
$\nu$ does not depend on RHS — would size $\nu$ from a partially assembled RHS
whose value at a shared node depends on the element numbering.

### 5.5 Step-cadenced history and $\Delta t$

The buffers `sp.qnm1/qnm2` are rolled once per time step by a gate in
`_sphere_ode_rhs!` that fires on the first stage of each step (`t` sweeps
$t + c_i\Delta t$ within a step), for the reason in §4.4. Both start at the
initial state.

$\Delta t$ is published to `sp.Δt[]` **after** the "land exactly on `tend`"
adjustment, since a residual built on the requested step rather than the taken
one is wrong by that rounding on every step. `sphere_rhs!` errors rather than
proceeding if it is unset: an unset $\Delta t$ does not degrade the model, it
inverts it — $R\to\infty$ everywhere and $\mu$ pins at the cap.

### 5.6 MPI

The shell path is single-rank by construction (`drivers.jl` refuses
`MPI.Comm_size > 1`), so $\langle q_i\rangle$ and
$\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$ are computed without
collectives. They are **domain** norms by definition, so the day the shell path
is parallelised those two passes need `MPI.Allreduce` exactly as §4.5 already
does.

### 5.7 Measured effect

Galewsky jet, nop = 5, 10 elements per cubed-sphere panel edge, `:cfl => 0.35`,
at $t = 5.90$ d (full table in the case README):

| | filter + $\nu = 2\times10^5$ (`SWsphere`) | **DynSGS** (`SWsphereDSGS`) |
|---|---|---|
| $\max\lvert\zeta\rvert$ | 6.79e-5 | **9.25e-5** |
| $\max\lvert\zeta-\zeta_0\rvert$ | 8.88e-5 | **1.50e-4** |
| $\delta E/E$ | −8.84e-4 | −8.14e-4 |
| $\nu$ (momentum) | 2e5 everywhere | mean 8.6e4, max 4.5e5 (1.6% of cap) |

The initial jet carries $\max\lvert\zeta_0\rvert = 1.12\times10^{-4}$, so the
filtered constant-$\nu$ deck has destroyed 39% of the vorticity by day 6 against
DynSGS's 17%, while carrying a *larger* mean viscosity and decaying slightly
more energy. The difference is not how much dissipation each adds but where.

$\mu_{res}$ governs throughout — $\nu$ never approaches the
$C_2\Delta(\lVert u\rVert+c)$ cap — which is the same regime §7 reports for the
other two 2D paths.

---

## 6. Code map, inputs and output

| file | contents |
|---|---|
| `src/kernel/abstractTypes.jl` | `struct DSGS`, `struct DSGS_MHD`, `struct DSGS_SW` |
| `src/kernel/physics/SGS.jl` | `compute_dsgs_viscosity!` (1D, 2D-θ, 2D-MHD, 2D-SW-sphere), `broadcast_dsgs_to_nodes!`, the `SGS_diffusion` accessors |
| `src/kernel/operators/rhs.jl` | dispatch in `viscous_rhs_el!`, `_viscous_rhs_el_2d_dsgs!`, the step-cadenced history gate in `_build_rhs!` |
| `src/kernel/operators/sphere_rhs.jl` | the shell path: `sphere_dsgs_requested`, `build_sphere_element_size`, the DynSGS pass in `_sphere_rhs_kernel!`, `sphere_dsgs_nodal!`/`sphere_dsgs_extra` |
| `src/kernel/solvers/sphere_time_loop.jl` | `sphere_dsgs_mu_bound` (the CFL bound), `sphere_dsgs_init_history!`, `sphere_dsgs_roll_history!`, the banner and diagnostics lines |
| `src/kernel/infrastructure/params_setup.jl` | `μ_dsgs`, `μ_dsgs_pnode`, `visc_coeff_dsgs`, `dsgs_qnm1/2`, `dsgs_avg/denom`, `dsgs_thist` (flat paths) |
| `src/io/mod_inputs.jl` | `:dsgs_C1`, `:dsgs_C2`, `:dsgs_gamma`, `:dsgs_Prt` defaults |
| `src/io/write_output.jl` | the `mu_dsgs_*` VTK fields |
| `test/test_sphere_dsgs.jl` | the shell model: cap, formula, locality, selectivity, the zonal regression, deck switches |
| `tools/plot_orszag_tang.jl` | off-line figures from a finished MHD run, including the viscosity map |
| `tools/vtu_reader.jl` | the minimal `.pvtu`/`.vtu` reader that script uses |

**Case inputs**

```julia
:lvisc      => true,
:visc_model => DSGS(),        # 1D CompEuler / 2D CompEuler θ
:visc_model => DSGS_MHD(),    # 2D ideal GLM-MHD
:visc_model => DSGS_SW(),     # shallow water on the spherical shell
:μ          => [0.0, 1.0, …], # per-equation multipliers, length neqs
:dsgs_C1    => 1.0,           # DSGS_MHD and DSGS_SW
:dsgs_C2    => 0.5,
:dsgs_gamma => 5.0/3.0,       # DSGS_MHD only
:dsgs_Prt   => 0.7,           # DSGS_MHD only
```

⚠ On the shell, `:μ` is the **dimensionless** multiplier of §5.3, not a
viscosity in m²/s — the model supplies $\nu$ itself. It defaults to 1.0 there
(rather than 0.0), and `build_sphere_viscosity` rejects values above 100 so that
a constant-$\nu$ value left over from a `SWsphere` deck is an error rather than
a $10^5\times$ overdose. Which equations it applies to is chosen with
`:ivisc_equations`, as on the constant-$\nu$ path.

**Output.** The per-element coefficients are broadcast to nodes by
`broadcast_dsgs_to_nodes!` and written to VTK as one point-data field per
equation, named after the solution variable each one damps:
`mu_dsgs_ρ`, `mu_dsgs_ρu`, …, `mu_dsgs_ψ`. They are piecewise constant per
element by construction, and shared (DSS) nodes take the value of the last
element that writes them — fine for visualization, not a nodal field.

⚠ The slots are **not all in the same units** (see §4.3): momentum and energy
carry the dynamic $\bar\rho\mu$, the magnetic and $\psi$ slots the kinematic
$\mu$. Compare a slot against itself over time, not against a different slot.

In 1D the same data also drives an extra panel in the PNG writer. In 2D, setting
`:outformat => "png"` renders one `μ_dsgs_<var>` filled-contour panel per output
time alongside the solution fields; with the default `:outformat => "vtk"` the
fields go to the `.pvtu` only, and

```bash
julia --project=. tools/plot_orszag_tang.jl
```

renders the multi-time composite used in the documentation
(`assets/MHD_OT_mu_dsgs_By.png`, the viscosity of the $B_y$ equation with $B_y$
isolines on top) from a finished Orszag–Tang run.

---

## 7. Defects found and fixed

Writing this document meant reading the three paths side by side, which turned
up five defects. All are now fixed; they are recorded here because the fixes
change the behaviour of `sod1d`, `case1` and `theta_dsgs`, and because the
reasoning matters if the model is revisited.

1. **The 2D θ-path dropped $M^{-1}$ from the residual.** It differenced the
   BDF2 $\partial q/\partial t$ against the raw **weak-form** RHS, two
   quantities separated by a factor of the mass matrix, so the "residual" was
   not one. A source comment defended this: the dimensionally-correct
   $M^{-1}\cdot$RHS "shrinks the residual by ~10³ on 2D atmospheric meshes and
   effectively turns DSGS off".

   *Fixed*, and the concern turns out not to apply to the corrected code —
   see the measurements below. The most likely explanation is that the old
   comment was written while defect 4 was also present: with a stage-cadenced
   history the BDF2 term is not a time derivative at all, so removing $M^{-1}$
   was compensating for a broken numerator. Had it still under-stabilized, the
   fix would have been to raise $C_1$, not to restore wrong units.

2. **The 2D θ-path had its momentum slots zeroed** by a leftover diagnostic
   (`# DIAG: momentum DSGS forced to zero …`), so only $\rho\theta$ was
   stabilized. *Fixed* — slots 2 and 3 now carry Marras eq. (10a).

   Doing so exposed a units bug that would otherwise have been silent: the
   θ-path had the same kinematic-vs-dynamic mismatch as the MHD one. $\mu$ from
   the model is kinematic, but `user_primitives!` supplies $(\rho,u,v,\theta)$,
   so the momentum and θ slots need $\bar\rho\mu$. Restoring the momentum
   slots without that factor would have applied a coefficient wrong by $\rho$.

3. **`SGS_diffusion(::DSGS, ::NSD_2D)` defined only the `inputs` signature.**
   The generic 2D `_expansion_visc!` calls it both ways —
   `(…, PhysConst, Δ2, VT, SD)` from the momentum and scalar branches,
   `(…, PhysConst, Δ2, inputs, VT, SD)` from the τ·u viscous-work term. Verified
   against a live run: `theta_dsgs` raised

   ```
   MethodError: no method matching SGS_diffusion(…, ::DSGS, ::NSD_2D)
   ```

   at `rhs.jl:2207` — the *"other scalars"* branch, reached by `ieq = 1`
   (density) on the very first RHS call. The path did not fail partway through;
   **it could not start at all**, so `theta_dsgs` had evidently not been run
   since the calling convention diverged. *Fixed* by adding the missing method;
   the case now runs to its final time.

   (`SMAG`/`VREM` in 2D were unaffected — they go through the separate
   cache-reading `_expansion_visc!` at `rhs.jl:2249`, which calls the
   `(…, ρ, ip, sgs, ltheta_eqn, SD)` accessor instead.)

4. **The BDF2 history was stage-cadenced.** `qp.qnm1/qnm2` advance on every RK
   *stage*, so the 1D and 2D-θ residuals differenced consecutive stage snapshots
   over a full $\Delta t$ — not an approximation of $\partial q/\partial t$.
   *Fixed*: the `dsgs_qnm1/qnm2` buffers built for `DSGS_MHD` (§4.4) are now
   allocated and used for `DSGS()` as well.

5. **The 1D wave-speed cap used `sqrt(γ·e_int)`.** For a perfect gas
   $p = (\gamma-1)\rho e_{int}$, so $a^2 = \gamma(\gamma-1)e_{int}$; the cap was
   inflated by $1/\sqrt{\gamma-1} \approx 1.58$ at $\gamma = 1.4$, letting
   $\mu_{res}$ govern more often than the Marras bound intends. *Fixed.*

### Verification

All three cases run to completion with the fixes in place:

| case | result |
|---|---|
| `CompEuler/sod1d` (1D) | completes, $t = 0.2$ |
| `CompEuler/theta_dsgs` (2D θ) | completes, $t = 1000$ — previously could not start |
| `MHD/orszagTangBormanis2024` (2D GLM-MHD) | completes, $t = 1$ |

The corrected model is **active and residual-governed**, not switched off — the
worry attached to defect 1. Measured DynSGS coefficients:

| case | $\mu$ range | vs. the $C_2\Delta(\lVert v\rVert+c)$ cap |
|---|---|---|
| MHD, $t=0.55$ | mean 1.6e-4, max 9.9e-4 (kinematic) | ≈ 10% of cap |
| `theta_dsgs`, $t=1000$ | mean 74, max 321 (dynamic $\bar\rho\mu$) | ≈ 10% of cap |

In both cases $\mu$ sits an order of magnitude below the first-order-upwind
bound, i.e. $\mu_{res}$ is doing the work and $\min(\mu_{max},\mu_{res})$ is not
saturating — exactly the regime the model is designed for. `theta_dsgs` shows
the expected cold-start spike (max $\mu \approx 9.7\times10^3$ at $t=0$, when
every denominator is at its floor) decaying to $O(10^2)$ once the flow develops,
and the rising thermal bubble itself behaves physically ($w$ up to 10 m/s, the
2 K perturbation intact at 1.45 K after 1000 s).
