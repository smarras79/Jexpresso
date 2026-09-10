# DynSGS — residual-based Dynamic Sub-Grid Scale stabilization in Jexpresso

> The equations exactly as implemented, with a term-by-term comparison against Dao & Nazarov (2022), are in [`docs/DSGS.tex`](docs/DSGS.tex) (compile with `pdflatex`).

This document describes the Marras–Nazarov Dynamic SGS model (`DSGS`), how it is
formulated for a general conservation law, and how it is implemented for each of
the three equation sets that use it in Jexpresso:

1. [The general conservation law](#1-the-general-conservation-law)
2. [1D CompEuler — `sod1d`, `case1`](#2-1d-compeuler--sod1d-case1)
3. [2D CompEuler θ-form — `theta_dsgs`](#3-2d-compeuler-θ-form--theta_dsgs)
4. [2D ideal GLM-MHD — `orszagTangBormanis2024`](#4-2d-ideal-glm-mhd--orszagtangbormanis2024)
5. [Code map, inputs and output](#5-code-map-inputs-and-output)
6. [Defects found and fixed](#6-defects-found-and-fixed)

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
\mu_{res}\big|_e = C_R\,\Delta_e^2\,
\max_i \frac{\lVert R_i\rVert_{\infty,e}}{\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}},
\qquad
\mu_{max}\big|_e = C_{max}\,\Delta_e\,\big(\lVert\mathbf{v}\rVert + c\big)_{\infty,e},
\qquad
\mu\big|_e = \max\!\big(0,\ \min(\mu_{max},\ \mu_{res})\big)\;}
$$

with $C_R \approx 1$, $C_{max} \approx 0.5$ (Dao & Nazarov's names, kept throughout), and $\Delta_e$ the element length scale
divided by the polynomial count, $\Delta_e = \Delta_{elem}/(N+1)$.

The three ingredients:

- **Normalization.** Each residual is divided by the *global* spread of its own
  variable, $\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$. This makes
  every ratio a frequency $[1/\mathrm{time}]$ regardless of the physical units
  or magnitude of $q_i$, so the $\max_i$ over equations compares like with like,
  and the whole model is dimensionally consistent without any tuned scale.

- **Units.** $R_i/\lVert\cdot\rVert$ has units $1/T$ and $\Delta^2$ has units
  $L^2$, so $\mu$ comes out as $L^2/T$ — a **kinematic** viscosity. The
  wave-speed cap $C_{max}\Delta(\lVert v\rVert+c)$ is $L \cdot L/T$, the same.
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
apply it (see [§6](#6-defects-found-and-fixed) — the 2D θ-path did not until
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

**Coefficients.** $C_R = 1$, $C_{max} = 0.5$, $\gamma = 1.4$ (hardcoded in the
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
selected with `:visc_model => DSGS_MHD()`. Also used, with `:dsgs_gamma => 1.05`,
by the stratified flux-emergence case
[`problems/MHD/fluxEmergenceSon2025`](problems/MHD/fluxEmergenceSon2025/README.md),
where the same model regularizes the fast/slow/intermediate shocks of an
emerging loop across eight decades of density.

This is a separate `AbstractVT` tag rather than an extension of `DSGS()` because
the residual set, the equation of state and the wave speed all differ from the
θ-system, and because the θ-path carries the debugging state described above.

**State.** $\mathbf{q} = (\rho,\ \rho u,\ \rho v,\ E,\ \rho w,\ B_x,\ B_y,\ B_z,\ \psi)$,
`neqs = 9`. See
[`problems/MHD/orszagTangBormanis2024/EQUATIONS.md`](problems/MHD/orszagTangBormanis2024/EQUATIONS.md)
for the equation set.

**Coefficients.** From `inputs`, defaulted in `mod_inputs.jl`:
`:dsgs_CR` (1.0), `:dsgs_Cmax` (0.5), `:dsgs_gamma` (5/3), `:dsgs_Prt` (0.7).
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
is exactly zero for them. Each denominator is bounded from below ("floored":
`denom = max(denom, a)`, so it can never be smaller than `a`) at $10^{-3}$ of that
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
| $\rho$ | $\texttt{:μ}[1]\cdot\mu$ | kinematic (a conservative mass diffusion $\nabla\cdot(\mu\nabla\rho)$, the $\nu\nabla\rho$ of Dao & Nazarov's eq. 4.4 that keeps $\rho$ positive; on in the Orszag–Tang and flux-emergence cases, 0 in KH) |
| $\rho u,\rho v,\rho w$ | $\texttt{:μ}[i]\cdot\bar\rho\,\mu$ | dynamic |
| $E$ | $\texttt{:μ}[4]\cdot\bar\rho\,\mu\cdot\dfrac{\gamma}{(\gamma-1)Pr_t}$, or $\texttt{:μ}[4]\cdot\bar\rho\,\mu/Pr_t$ with `:dsgs_nazarov_energy` | dynamic |
| $B_x,B_y,B_z$ | $\texttt{:μ}[i]\cdot\mu$ | kinematic (turbulent resistivity) |
| $\psi$ | $\texttt{:μ}[9]\cdot\mu$ | kinematic |

$\bar\rho$ is the element-mean density. The energy factor follows from slot 4's
primitive being $T = p/\rho\ (= R\,T_{phys})$: the physical flux is
$\nabla\cdot(k\nabla T_{phys})$ with $k = \mu_{dyn}c_p/Pr_t$, so rewriting in
terms of $T$ gives $k/R = \mu_{dyn}\gamma/((\gamma-1)Pr_t)$ since
$c_p = \gamma R/(\gamma-1)$. Dao & Nazarov (2022, §4.4) instead take
$\kappa = \mu_{dyn}/Pr$ directly on $T = p/\rho$, i.e. $c_p = 1$ in these
units, a conduction $\gamma/(\gamma-1)$ times weaker (2.5× at $\gamma = 5/3$);
`:dsgs_nazarov_energy => true` selects that, and the MHD cases use it with
$Pr = 1$ as in their runs. With `:μ[1] = 1` the slot coefficients are then
exactly their set by equation: $\nu$ on $\rho$, $\rho\nu$ on $\mathbf{u}$,
$\rho\nu/Pr$ on $T$, $\nu$ on $\mathbf{B}$. The energy equation of the
physical form also carries, next to the $\tau\cdot\mathbf{v}$ viscous work, the
resistive work $\nabla\cdot(\eta\,\mathbf{B}\cdot\nabla\mathbf{B})$ that matches the
component-Laplacian induction term, so that the magnetic energy removed from
$\mathbf{B}$ reappears as heat and the total energy is conserved (their
eq. 4.4 has the same term for their curl-curl form).

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

### 4.5 Stratified atmospheres: `:dsgs_norms => "element"`, `:dsgs_conserved`, `:dsgs_ref_weight`, `:dsgs_nazarov_energy`, `:dsgs_nodal_rho`

Two opt-in variants of the MHD kernel, both `false` by default (the
Orszag–Tang results below are unchanged), added for
[`fluxEmergenceSon2025`](problems/MHD/fluxEmergenceSon2025/README.md), whose
density spans eight decades between the photosphere and the corona:

- **`:dsgs_norms => "element"`** normalizes the residual of equation $i$ in
  element $e$ by the spread of $q_i$ over that element,
  $\lVert q_i - \langle q_i\rangle_e\rVert_{\infty,e}$, bounded from below at
  `:dsgs_local_rel` (default 1) times the *element-mean* scales of §4.2
  ($\rho_e$, $\rho_e c_e$, $\rho_e c_e^2$, $\sqrt{\rho_e}c_e$), instead of the
  domain spread. Unlike the domain norms, whose $10^{-3}$ floors only guard
  against a degenerate spread, here the floors *are* the normalization of a
  quiescent element (the spread of $\rho\mathbf{v}$ in an atmosphere at rest
  is zero): with a $10^{-3}$ floor a residual of $2.5\times10^{-3}\rho c$ per
  unit time already drove $\mu$ to the wave-speed cap over the whole quiet
  chromosphere of the flux-emergence case, whose sheet then eroded by
  resistive diffusion (12% of its peak field in $2\tau_0$) and sank at
  $0.05\,C_s$; at the default value the same run keeps $\mu \approx 4\times10^{-4}$
  in the sheet, its peak field to four digits, and a horizontal-mean vertical
  velocity below $10^{-4}$, while a grid-scale sawtooth (ratio $\sim v_{saw}/\Delta$)
  or a shock ($\sim c/\Delta$) still saturates the cap. With the domain norm the dense bottom of the atmosphere sets
  the scale of $\rho$, $\rho\mathbf{v}$ and $E$, and a residual in the corona
  — where those fields are $10^{-8}$ of it — is invisible: in the flux-emergence
  run a grid-scale sawtooth grew across the chromosphere–corona transition with
  $\mu = 10^{-11}$ there. The element spread of a smooth stratified field is
  $O(q_i)$ ($\rho$ changes by $e^{-1}$ across a $1H_0$ element), so the ratio
  keeps the meaning of a relative under-resolution rate.
- **`:dsgs_conserved => true`** hands every slot the same kinematic $\mu$ and
  drops the $\tau\cdot u$ term; with a `user_primitives!` that returns the
  conserved variables the operator becomes a Laplacian on
  $(\rho, \rho\mathbf{v}, E, \mathbf{B}, \psi)$, the form in which an isobaric
  contact (the 25× density drop of the solar transition region) diffuses
  consistently: $\rho$ spreads, $E$ (constant across it) does not, $p$ stays
  what it was. Diffusing $\rho$ alone under the $T$-based energy closure drove
  $p$ negative within a few $\tau_0$. The flux-emergence case applies it to the
  *departure from its magnetostatic reference state*, $q - q_e$ (see its
  `user_primitives.jl`): on $q$ itself the Laplacian of the exponential
  stratification is a steady mass source that the sensor feeds on (sheet
  sinking at $0.3\,C_s$ by $t = 8\tau_0$, measured), on $q - q_e$ the operator
  vanishes at rest and reduces to the plain conserved-variable Laplacian
  wherever the state has left the reference.
- **`:dsgs_Cmin`** (default 0) floors the coefficient at $C_{min}\Delta(\lVert\mathbf{v}\rVert + c_f)$,
  a fraction of the $C_{max}$ cap (not in Dao & Nazarov). The residual cannot see a node-to-node mode
  (the discrete operator returns nearly nothing on it), and the CG
  discretization leaves such a mode undamped; the flux-emergence case uses
  $C_{min} = 0.03$, which damps it at a rate $C_{min}\Delta c(\pi/\Delta)^2 \approx 7/\tau_0$ in
  its corona while spreading a resolved structure by $\sqrt{C_{min}\Delta c\,t} \approx 1H_0$
  over the whole run.
- **`:dsgs_ref_weight => true`** (with `:dsgs_conserved`) diffuses the fluid
  slots as the *relative* departure from the reference state with the
  reference density as weight, $\nabla\cdot(\mu\rho_e\nabla((q - q_e)/\rho_e))$
  on $(\rho, \rho\mathbf{v}, E)$, the magnetic slots as before. The case's
  `user_primitives!` returns $(q - q_e)/\rho_e$ in slots 1–5 and stores
  $\rho_e$ in the spare slot `neqs+1` of `uprimitive`, which
  `_expansion_visc!` multiplies into the coefficient at the quadrature
  point. Like the $q - q_e$ form it vanishes at rest and is conservative;
  unlike it, it obeys a maximum principle across a reference jump: a
  $-10\%$ departure on the dense side of the solar transition region is,
  in absolute terms, more than the whole density of the light side, and
  the $q - q_e$ Laplacian carries it across — that is how the
  flux-emergence case lost positivity at $t \approx 14\tau_0$ with $\mu$
  already at its cap and came to need its floors. The weighted form pulls
  $\rho/\rho_e$ toward its neighbours and no further. Where $q_e$ is
  negligible (the emerged loop) it reduces to the conserved-variable
  Laplacian. Used by
  [`fluxEmergenceSon2025DSGS`](problems/MHD/fluxEmergenceSon2025DSGS/README.md),
  the limiter-free version of the flux-emergence case.
- **`:dsgs_nazarov_energy => true`** (with `:dsgs_conserved`) gives the
  slots the coefficients of Dao & Nazarov (2022, *J. Sci. Comput.* 92:77,
  §4.4): one kinematic $\nu$ from the maximum of the normalized residuals
  (their eq. 4.8, the $\mu$ of §4.1), then $\nu$ on $\nabla\rho$, $\rho\nu$ in
  the momentum stress, $\kappa = \rho\nu/\mathrm{Pr}$ on the temperature and
  $\eta = \nu$ on $\mathbf{B}$. In conserved variables the $\rho$, $\rho\mathbf{v}$
  and $\mathbf{B}$ slots already are that ($\nu\nabla(\rho\mathbf{v}) \approx \rho\nu\nabla\mathbf{v}$),
  the energy slot is not: $\nu\nabla E$ conducts the internal energy
  $\rho T/(\gamma(\gamma-1))$ at $\nu$, a heat conduction
  $\rho\nu/(\gamma(\gamma-1))$ that is $19\times$ Nazarov's at $\gamma = 1.05$,
  $\mathrm{Pr} = 1$. The option therefore **splits the energy flux**: the
  case's `user_primitives!` hands slot 4 the non-thermal part
  $\delta(\tfrac12\rho|\mathbf{v}|^2 + \tfrac12|\mathbf{B}|^2)$, which the
  kernel diffuses with the $\rho\mathbf{v}$ slot's $\nu$ so that its
  magnetic and kinetic fluxes keep matching the $\mathbf{B}$ and
  $\rho\mathbf{v}$ Laplacians, and the spare slot `neqs+2` the thermal part
  $\delta(p/(\gamma-1))$, diffused with $\max(\gamma(\gamma-1)/\mathrm{Pr}_t\cdot\nu_{res},\ \nu_{floor})$
  — Nazarov's $\kappa = \rho\nu/\mathrm{Pr}$ on the residual viscosity with
  the $C_{min}$ floor of §4.5 kept in full (`dsgs_split_energy` in SGS.jl,
  `_expansion_visc!` in rhs.jl). Scaling the whole energy slot instead —
  the first implementation — let $\nu\nabla\mathbf{B}$ spread the flux
  sheet's field while 95 % of its magnetic energy stayed put, and cut the
  floor that damps the node-to-node mode: on the flux-emergence case the
  sheet core overheated, a temperature sawtooth grew across the corona by
  $t = 10\tau_0$ and the emergence stalled (measured; with the slot back
  at $\nu$ the validated result returned). Their $\kappa\nabla T$ itself cannot be used in a
  two-temperature atmosphere: at rest it conducts across the reference
  temperature jump of the transition region, and on $T - T_e$ it heats loop
  gas crossing the fixed height of that jump; $E$ has no such jump (the
  contact is isobaric), which is why the conserved form is kept and only
  the coefficient is Nazarov's.
- **`:dsgs_nodal_rho => true`** forms the dynamic coefficient of the momentum
  and energy slots with the density of the quadrature point (in
  `SGS_diffusion(::DSGS_MHD)`) instead of the element mean $\bar\rho$, i.e.
  the viscous flux is $\nabla\cdot(\rho\mu\nabla u)$. With $\bar\rho$ the
  effective diffusivity at the light side of an element is $(\bar\rho/\rho)\mu$
  — up to $25\mu$ across the transition region — and breaks the explicit
  viscous stability limit as soon as the model switches on there. With this
  option the `mu_dsgs_ρu`, `mu_dsgs_ρv`, `mu_dsgs_ρw` and `mu_dsgs_ρE` output
  fields hold the kinematic coefficient, like the magnetic slots.

### 4.6 The 1D kernel: `brioWu1d`

`compute_dsgs_viscosity!(::DSGS_MHD, ::NSD_1D)` is the same model for the
8-variable 1D system $(\rho, \rho u, \rho v, \rho E, \rho w, B_x, B_y, B_z)$
(no GLM field in 1D): the same residual, normalization, cap with the
1D fast speed, floor and slot assignment as §4.1–4.5, called from the 1D
viscous path (`viscous_rhs_el!(…, ::NSD_1D)`), which then applies one scalar
Laplacian per slot exactly as the Euler `DSGS()` path does. Used by
[`problems/MHD/brioWu1d`](problems/MHD/brioWu1d/README.md) in the conserved
form; the physical-form coefficients are available but the 1D loop carries no
$\tau\cdot u$ or $\eta\mathbf{B}\cdot\nabla\mathbf{B}$ work terms, so
only the conserved form conserves total energy in 1D.

### 4.7 Element or nodal coefficient: `:ldsgs_nodal`

Every kernel above computes **one coefficient per element** (Marras's form):
the maximum of the normalized residual over the element's nodes, applied
as a constant over the element, so that $\nu$ is a staircase with a jump
at every element interface. `:ldsgs_nodal => true` (default `false`, i.e.
element form) selects instead the **nodal form**, which is Dao & Nazarov's
(2022) formulation itself, for the 1D and 2D kernels, `DSGS` (θ and
total-energy forms) and `DSGS_MHD` alike (`compute_dsgs_viscosity_nodal!`):

- the residual is the assembled lumped-mass nodal residual
  $R_i = |\mathrm{BDF2}(q)_i - M_i^{-1}\,\mathrm{rhs}_i|$;
- it is normalized by $n(w)_i = \bar S(w)\,(1 - C_l\,(\max_{I(i)} w - \min_{I(i)} w)/(\max w - \min w))$,
  their eq. 4.7: $\bar S$ the global spread of §4.2 (with its floors), $I(i)$
  the support of node $i$ (the elements containing it), `:dsgs_Cl` their
  $C_l$ (0 = the classical $\bar S$, 0.4 in their runs), with the
  $n^2/(n^2+\epsilon)$ guard of eq. 4.8;
- $\nu_i = \min(C_{max}h_i\lambda_i,\ C_R h_i^2 R_i)$ at every node (eq. 4.10),
  $h_i = \max \Delta_K/(k+1)$ over the support (the element form's $\Delta$;
  the paper's $h_K/k$ with the circumradius is $\Delta_K/(\sqrt2 k)$ on a
  square, within 12 % of it at $k = 4$, while $\Delta_K/k$ put the
  rising-bubble case past the explicit viscous limit at start-up),
  $C_{max} =$ `:dsgs_Cmax`, $C_R =$ `:dsgs_CR`,
  bounded from below by $C_{min} h_i\lambda_i$ (a `max`, the counterpart of the `min` cap);
- the slot coefficients from $\nu_i$ exactly as in the element kernels, with
  the **nodal** density in the dynamic coefficients;
- $\nu$ is a continuous ($C^0$, DSS'd) field: the element loop gathers the
  element's nodal values (`params.dsgs_μloc`, preallocated) and the viscous
  expansion uses $\nu$ at each quadrature point, so the diffusive flux has
  no jump at element interfaces. `μ_dsgs_pnode` holds the nodal field
  itself and `μ_dsgs[ie,:]` its element means (for the staircase output).

The whole path is allocation-free (`params.dsgs_qmin/qmax/nmin/nmax/hnod`
are its scratch). The element form remains the default of every case,
`brioWu1d` included (its deck carries the nodal switch commented out; with
the $C_{min}$ floor both forms give the same profile). There is no 3D DynSGS kernel (the 3D viscous
path dispatches the Smagorinsky/Vreman caches only), so the switch has no
3D counterpart yet. Measured on the Brio–Wu tube: both forms give the same
solution, and the element-scale ripples the compound wave radiates into
the plateau behind it are damped by neither — the residual viscosity
scales with their amplitude — and need the $C_{min}$ floor (3 % there, see the
case README).

### 4.8 MPI

$\langle q_i\rangle$ and $\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$
are **domain** norms by definition, so both reductions are `MPI.Allreduce`d. A
rank-local version would make the eddy viscosity depend on the partitioning. The
cost is two small collectives per RHS call.

### 4.9 Measured effect

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
| `src/io/mod_inputs.jl` | `:dsgs_CR`, `:dsgs_Cmax`, `:dsgs_gamma`, `:dsgs_Prt` defaults |
| `src/io/write_output.jl` | the `mu_dsgs_*` VTK fields |
| `tools/plot_orszag_tang.jl` | off-line figures from a finished MHD run, including the viscosity map |
| `tools/vtu_reader.jl` | the minimal `.pvtu`/`.vtu` reader that script uses |

**Case inputs**

```julia
:lvisc      => true,
:visc_model => DSGS(),        # 1D CompEuler / 2D CompEuler θ
:visc_model => DSGS_MHD(),    # 2D ideal GLM-MHD
:μ          => [0.0, 1.0, …], # per-equation multipliers, length neqs
:dsgs_CR    => 1.0,           # DSGS_MHD only
:dsgs_Cmax    => 0.5,
:dsgs_gamma => 5.0/3.0,
:dsgs_Prt   => 0.7,
```

**Output.** The per-element coefficients are broadcast to nodes by
`broadcast_dsgs_to_nodes!` and written to VTK as one point-data field per
equation, named after the solution variable each one damps:
`mu_dsgs_ρ`, `mu_dsgs_ρu`, …, `mu_dsgs_ψ`. Slots that carry the same
coefficient are written once, under a name listing them
(`mu_dsgs_ρ_ρu_ρv_ρw_Bx_By_Bz_ψ` and `mu_dsgs_ρE` in the conserved form of
§4.5, or a single `mu_dsgs` when all nine agree). They are piecewise constant per
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

## 6. Defects found and fixed

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
   fix would have been to raise $C_R$, not to restore wrong units.

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

| case | $\mu$ range | vs. the $C_{max}\Delta(\lVert v\rVert+c)$ cap |
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
