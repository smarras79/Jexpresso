# DynSGS — residual-based Dynamic Sub-Grid Scale stabilization in Jexpresso

This document describes the Marras–Nazarov Dynamic SGS model (`DSGS`), how it is
formulated for a general conservation law, and how it is implemented for each of
the four equation sets that use it in Jexpresso:

1. [The general conservation law](#1-the-general-conservation-law)
2. [1D CompEuler — `sod1d`, `case1`](#2-1d-compeuler--sod1d-case1)
3. [2D CompEuler θ-form — `theta_dsgs`](#3-2d-compeuler-θ-form--theta_dsgs)
4. [2D ideal GLM-MHD — `orszagTangBormanis2024`](#4-2d-ideal-glm-mhd--orszagtangbormanis2024)
5. [3D CompEuler θ-form — `rtb3d_dsgs`, `LESICP2-64x64x60-dynsgs`](#5-3d-compeuler-θ-form--rtb3d_dsgs-lesicp2-64x64x60-dynsgs)
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
divided by the polynomial count, $\Delta_e = \Delta_{elem}/(N+1)$. Which element
length that is differs between the 1D/2D paths and the 3D one — see
[§5.3](#53-filter-width).

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

The $M^{-1}$ is **required** for the units above to work out; every path
applies it (see [§7](#7-defects-found-and-fixed) — the 2D θ-path did not until
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

## 5. 3D CompEuler θ-form — `rtb3d_dsgs`, `LESICP2-64x64x60-dynsgs`

`compute_dsgs_viscosity!(::SGS_DSGS, ::NSD_3D)` and
`compute_sgs_cache!(::SGS_DSGS, …, ::NSD_3D)` in `src/kernel/physics/SGS.jl`.

**State.** $\mathbf{q} = (\rho,\ \rho u,\ \rho v,\ \rho w,\ \rho\theta)$,
`neqs = 5`, plus any moisture slots.

**Coefficients.** $C_1$, $C_2$ from `:dsgs_C1` / `:dsgs_C2` (1.0 and 0.5).

### 5.1 Why this path goes through a closure struct

The 1D and 2D paths bypass `params.sgs`: they compute one coefficient per
element, pack it into `visc_coeff_dsgs` and hand it to a scalar
$\nabla\cdot(\mu\nabla q)$ kernel. That is not what the 3D viscous kernel is.
The 3D one assembles the full stress tensor

$$
\tau_{ij} = \mu\Big(\frac{\partial u_i}{\partial x_j} + \frac{\partial u_j}{\partial x_i}\Big)
          - \tfrac23\mu\,(\nabla\cdot\mathbf{u})\,\delta_{ij}
$$

and it is reachable only through the `sgs::AbstractSGSModel` method of
`_expansion_visc!`, which reads its coefficient from `SGS_diffusion`. So
DynSGS-3D is an `AbstractSGSModel`, `SGS_DSGS`, whose `μ_turb` is filled from
the residual instead of from the strain rate. Four things come with that and
none of them is available on the 1D/2D path:

| consumer | file | what it reads |
|---|---|---|
| the viscous RHS | `operators/rhs.jl` | the full $\tau_{ij}$, not a Laplacian |
| the CFL table's viscous row | `hevi/cfl_diagnostics.jl` | `sgs.μ_turb` |
| implicit vertical diffusion | `hevi/vdiffusion.jl` | `SGS_diffusion` on this struct |
| LES subfilter stresses | `io/les_statistics.jl` | `sgs.μ_turb`, `sgs.S11…S23` |

**Per-equation split.** DynSGS supplies $\mu_t$; the split is
`SGS_diffusion`'s, i.e. the standard LES one and identical to Smagorinsky's:

$$
\text{momentum } (\mu_{mol}+\mu_t)\cdot\texttt{:μ}[i],\qquad
\theta:\ \frac{\mu_t}{Pr_t}\cdot\texttt{:μ}[5],\qquad
\rho:\ \texttt{:μ}[1]\cdot\frac{\mu_t}{Sc_t}\ (=0)
$$

This is deliberately **not** the $Pr/(\gamma-1)$ factor the 2D θ-path applies
to slot 4. That factor is Nazarov & Hoffman's artificial heat conduction on
the internal energy of a shock-capturing scheme; $\theta$ is not that
variable, and an atmospheric LES wants a turbulent Prandtl number.

### 5.2 The normalisation is eq. (9) to the letter

$\lVert q_i - \langle q_i\rangle\rVert_{\infty,\Omega}$ is taken on the
**stored prognostic variable**, with no `qe` arithmetic — the same as every
other path in `SGS.jl`.

This is worth stating explicitly because the opposite was tried, shipped, and
reverted. The argument for taking the norms on $q - q_e$ was that over a 5 km
column the total field's spread is the hydrostatic background rather than the
turbulence, so a total-field denominator divides the residual by something
100–300× too large. That argument is wrong, and a rising bubble shows why in
one number: the bubble is constructed at nearly constant **pressure**, and
$p = C_0(\rho\theta)^\gamma$ pins $\rho\theta$ to the pressure — so a 2 K
perturbation in $\theta$ is carried almost entirely by $\rho$ going down, and
the **conserved** variable barely moves:

| | field norm | perturbation norm |
|---|---|---|
| $\lVert\rho\theta - \langle\cdot\rangle\rVert_\infty$ | ≈ 110 | ≈ 0.42 |

a factor of 260 on the coefficient. Measured on `CompEuler/rtb3d_dsgs`, the
perturbation form put $\nu$ at $1.0\times10^4\ \mathrm{m^2/s}$ within two
steps and at the $C_2\Delta(\lVert v\rVert+c)$ cap by the third, and drove
$\rho\theta$ negative. The field form gives $\nu \approx 90$–$120$, in the
same range as the 2D case on the same bubble. `test/sgs/test_dsgs3d.jl` pins
it.

The general lesson: the denominator is what makes $C_1 = 1$ a universal
constant. Redefining it silently redefines $C_1$. If a case is
under-stabilized, `:dsgs_C1` is the lever — not the norm.

Note that under `:SOL_VARS_TYPE => PERT()` the stored prognostic variable
**is** the perturbation, and the field norm is then automatically taken on it.
That is not an inconsistency: it is the norm of whatever the scheme is
actually advancing, which is what eq. (9) asks for.

### 5.3 Filter width

$\Delta = \texttt{Δelem\_filter}[e]/\texttt{nop}$ — the same width SMAG and VREM
get from `compute_element_size_driver` (`:les_filter_width`, default `:max`)
and the same one `les_statistics.jl` reports against. Not `Δelem/ngl` as in
2D: on a 160 × 160 × 40 m element that is 8 m against 40 m, a factor 25 in
$\mu$, and a closure and its own diagnostic disagreeing about the filter width
is how a run ends up unable to explain its output.

### 5.4 Denominator floors

Every field of a PBL case starts horizontally uniform, so every
$\lVert q'_i-\langle q'_i\rangle\rVert_\infty$ is exactly zero at $t=0$. Each
denominator is floored at $10^{-3}$ of that field's natural scale, built from
the domain-mean **total** state:

| slot | floor |
|---|---|
| $\rho$ | $10^{-3}\bar\rho$ |
| $\rho u,\rho v,\rho w$ | $10^{-3}\bar\rho\,\bar c$ |
| $\rho\theta$ | $10^{-3}\bar\rho\,\bar\theta$ |
| other scalars | $10^{-3}\langle\lvert q_i\rvert\rangle$ |

### 5.5 The BDF2 stencil is evaluated as $3(q-q^{n-1}) - (q^{n-1}-q^{n-2})$

Algebraically identical to $3q^n-4q^{n-1}+q^{n-2}$, and not the same rounding.
The direct form builds $3q$ and $4q^{n-1}$ and then cancels them, so its error
is $\varepsilon\lvert q\rvert$ — taken on the total state, which here is
$\rho\theta \approx 360$, not the perturbation. Divided by $2\Delta t$ and by a
**floored** denominator and multiplied by $\Delta^2 = 1600\,\mathrm{m}^2$, that
noise comes out as an eddy viscosity on a fluid at rest:

```
3*1.2 - 4*1.2 + 1.2 = -4.4e-16   ->   nu ~ 6e-10 m^2/s
```

Negligible in magnitude and still wrong in kind: the whole claim of a
residual-based model is that an exact solution gets **zero** viscosity.
Differencing consecutive states instead gives an exact zero on a steady state
and noise that scales with the change rather than with the background. The 1D,
2D and MHD paths still use the direct form; changing them would perturb their
committed reference solutions at round-off for no benefit their denominators
(which are not floored against a stratified background) need.

### 5.6 The cap is inert in a low-Mach flow

$\mu_{max} = C_2\Delta(\lVert v\rVert + c)$ with $c \approx 340$ m/s and
$\Delta = 40$ m is $\approx 7\times10^{3}\ \mathrm{m^2/s}$ — two to three orders
above anything a PBL closure produces, so $\min(\mu_{max},\mu_{res})$ never
binds and $\mu_{res}$ governs alone. That is the regime the model is designed
for, but it does mean `:dsgs_C2` is not a working knob here unless it is
dropped by orders of magnitude.

### 5.7 `:dsgs_add_smagorinsky`

Off by default. DynSGS is a **stabilization**: it puts viscosity where the
discrete solution fails to satisfy the PDE, which is nearly nothing in a
locally smooth region. In the surface layer of a wall-modelled PBL the
subfilter stress carries most of the momentum flux and has to be there whether
or not the solution is smooth — that is what makes the log law, and no
residual sensor produces it. With the switch on, $\mu_t = \mu_{smag}+\mu_{dsgs}$:
the closure keeps its wall behaviour and DynSGS adds dissipation only where the
residual asks for it. The Richardson stability function and the near-wall
mixing-length limit multiply the **Smagorinsky part only** — a residual is not
a strain rate.

The diagnostic that decides whether to switch it on is the near-surface
`*_sfs` columns of the LES statistics: if $\langle u'w'\rangle_{sfs}$ collapses
in the first few hundred metres and $\langle u'w'\rangle_{res}$ does not rise
to compensate, the surface layer is under-stressed.

### 5.8 Stability

DynSGS does not exempt a case from the parabolic step limit: it produces an
eddy viscosity like any other closure and $2\nu_{eff}/h_z^2$ applies to it the
same way. What changes is only *where* the viscosity sits. On the LESICP2
grid ($h_z = 6.91$ m) that rate is what `:implicit_vdiff` exists to remove, and
the DynSGS deck carries it on by default for the same reason the Smagorinsky
one does.

### 5.9 Cases

| case | what it is for |
|---|---|
| `CompEuler/rtb3d_dsgs` | the small one. A rising thermal bubble that is a **cylinder along y**, on a mesh one element thick in y: 10 × 1 × 10 elements over 10000 × 1000 × 10000 m, 8405 gridpoints, one core, minutes. |
| `CompEuler/LESICP2-64x64x60-dynsgs` | the production one. 15.9 M nodes, the turbulent PBL. |

`rtb3d_dsgs` is the one to run first, and it is a *test* rather than a small
run for three reasons:

1. **The solution must be y-invariant.** The bubble radius is measured in the
   x–z plane only and the y faces are free-slip, so $v \equiv 0$ exactly. `v`
   in the VTU is then a free correctness check on the whole 3D path — the
   stress tensor including its $-\tfrac23\mu\nabla\cdot\mathbf{u}$ term, the
   metrics, the DSS, the coefficient itself — and far easier to read than a
   slightly-too-diffuse bubble.
2. **It should reproduce the 2D case.** Same geometry, bubble and resolution
   as `CompEuler/theta` and `theta_dsgs`. Not bit-for-bit: the 2D path gives
   slot 4 the $Pr/(\gamma-1)$ artificial conduction ($0.25\mu$ at
   `:Pr => 0.1`) where the 3D path gives $\mu/Pr_t = 1.43\mu$, a factor 5.7 on
   the θ diffusivity. `:μ[5] => 0.175` reproduces the 2D coefficient exactly
   if you want them on equal footing.
3. **The §5.2 perturbation normalisation is visible in miniature.** Over the
   10 km column the reference state runs $\rho$ 1.17 → 0.4 and $\rho\theta$
   350 → 120 while the bubble's departure from it is $O(2)$: normalising on the
   total field would divide the residual by ~230 instead of ~2 and switch the
   model off.

**Mesh note.** `rtb3d_dsgs_10x1x10.msh` is isotropic — $dx = dy = dz$ =
1000 m — *deliberately*, so `:les_filter_width` cannot be set wrongly. A dummy
direction that is one element deep but many elements *wide* is harmless for an
inviscid run and not harmless for an LES closure: the default `:max` would take
the filter width from it and $\nu \propto \Delta^2$. The `.geo` header carries
the arithmetic.

### 5.10 Tests

`test/sgs/test_dsgs3d.jl` (standalone, no Jexpresso load, same discipline as
`test/sgs/test_closures.jl`): the coefficient on a manufactured residual is
the value eq. (9) predicts; an exact solution gets identically zero; the cap
binds only from above; the perturbation normalisation is invariant under
adding a sounding; $\mu \propto \Delta^2$ below the cap and $\propto\Delta$ at
it; the floors keep a uniform field at zero; every equation can drive the max.

`test/sgs/test_dsgs3d_wiring.jl` loads the **real** `sgsStructs.jl` and
`SGS.jl` against stubs for KernelAbstractions/MPI and runs them, so a
signature that does not match its call site, a field that is not on the
struct, or an allocator that hands 2D a closure it must not have fails there
rather than three hours into a 15.9 M-node run.

---

## 6. Code map, inputs and output

| file | contents |
|---|---|
| `src/kernel/abstractTypes.jl` | `struct DSGS`, `struct DSGS_MHD` |
| `src/kernel/physics/sgsStructs.jl` | `struct SGS_DSGS` and its `allocate_SGS` (3D only) |
| `src/kernel/physics/SGS.jl` | `compute_dsgs_viscosity!` (1D, 2D-θ, 2D-energy, 2D-MHD, **3D**), `compute_sgs_cache!(::SGS_DSGS)`, `broadcast_dsgs_to_nodes!`, the `SGS_diffusion` accessors |
| `src/kernel/operators/rhs.jl` | dispatch in `viscous_rhs_el!`, `_viscous_rhs_el_2d_dsgs!`, the step-cadenced history gate in `_build_rhs!` |
| `src/kernel/infrastructure/params_setup.jl` | `μ_dsgs`, `μ_dsgs_pnode`, `visc_coeff_dsgs`, `dsgs_qnm1/2`, `dsgs_avg/denom`, `dsgs_thist` |
| `src/io/mod_inputs.jl` | `:dsgs_C1`, `:dsgs_C2`, `:dsgs_gamma`, `:dsgs_Prt`, `:dsgs_add_smagorinsky` defaults |
| `src/kernel/physics/SGS.jl` | `dsgs_monitor` — `JEXPRESSO_DSGS_MONITOR=1` |
| `src/io/write_output.jl` | the `mu_dsgs_*` VTK fields |
| `tools/plot_orszag_tang.jl` | off-line figures from a finished MHD run, including the viscosity map |
| `tools/vtu_reader.jl` | the minimal `.pvtu`/`.vtu` reader that script uses |

**Case inputs**

```julia
:lvisc      => true,
:visc_model => DSGS(),        # 1D CompEuler / 2D CompEuler θ or E / 3D CompEuler θ
:visc_model => DSGS_MHD(),    # 2D ideal GLM-MHD
:μ          => [0.0, 1.0, …], # per-equation multipliers, length neqs
:dsgs_C1    => 1.0,           # DSGS_MHD and 3D DSGS
:dsgs_C2    => 0.5,           # DSGS_MHD and 3D DSGS
:dsgs_gamma => 5.0/3.0,       # DSGS_MHD only
:dsgs_Prt   => 0.7,           # DSGS_MHD only
:dsgs_add_smagorinsky => false, # 3D DSGS only — see §5.7
:dsgs_residual        => :tendency,  # | :strict — see defect 6
:ldsgs_global_norms   => false, # domain vs rank-local ⟨q⟩ and ‖q−⟨q⟩‖
```

**`JEXPRESSO_DSGS_MONITOR=1`** prints one line per step with $\nu$ max/mean,
the normalising denominators, and the coordinates of the node that produced
the worst ratio — flagged if that node is on a boundary. It is what found
defect 9, and it is the first thing to turn on when a DynSGS run misbehaves.

**On `:dsgs_C1`.** Because the default sensor is a *tendency* sensor and not
the paper's residual, the paper's $C_1 = 1$ carries no authority. It over-fires
during a spin-up from rest by roughly 4× — measured on both bubble cases —
and `CompEuler/rtb3d_dsgs` therefore runs at `:dsgs_C1 => 0.25`, with the
sweep that established it recorded in the deck.

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

6. **The sensor is a TENDENCY sensor, not the strict residual — and that is
   now said out loud rather than assumed.** `_build_rhs!` advances the
   step-cadenced pair *above* the viscous block and writes `uaux` (the current
   state) into `dsgs_qnm2`, while every call site passes
   `(q, q1, q2) = (uaux, dsgs_qnm2, dsgs_qnm1)`. So `q1` holds the same
   contents as `q` and the stencil collapses:

   $$3q^n - 4q^1 + q^2 \;\to\; 3q^n - 4q^n + q^{n-1} = -(q^n - q^{n-1})$$

   giving $R \approx 1.5\lvert\partial q/\partial t\rvert$.

   That was found and "fixed" — and the fix was reverted, because it makes the
   model worse and the measurement says so. Rolling *after* the read gives the
   literal BDF2, so $R = \lvert\partial q/\partial t - M^{-1}\mathrm{RHS}\rvert$;
   but `params.RHS` at that point holds the **inviscid** flux divergence only
   (the viscous term is added a few lines later) while $\partial q/\partial t$
   is the derivative of the solution the viscous term helped produce. So the
   strict residual equals **the artificial viscosity's own contribution**, and
   the model reads itself. On a weak-form Galerkin RHS that is a contraction.
   Measured on `CompEuler/sod1d` at $t = 0.2$ (element 85 is the shock,
   element 1 is $x=0$):

   | | $\mu_{max}$ | where |
   |---|---|---|
   | tendency (default) | $1.83\times10^{-3}$ | element 85 — **the shock** |
   | strict BDF2 | $1.19\times10^{-4}$ | element 1 — the boundary |

   with visible oscillations in $\rho$ and $\rho E$ in the second. So the
   default is `:dsgs_residual => :tendency`, named for what it is, and
   `:strict` is kept as an option because the question of how to get a true
   residual out of a weak-form RHS is worth returning to — the answer is
   probably to compare against `RHS_inviscid + RHS_visc` from the *previous*
   step, which is the consistency error of the scheme rather than either of
   these.

7. **There was no warm-up guard.** The stencil needs three step-cadenced
   levels and a run starts with one, so the first two steps produce nothing
   meaningful. `params.dsgs_hist` holds the model at $\mu = 0$ until
   $t \ge t_{init} + 2\Delta t$. Written in `time` rather than as a roll
   counter because `time` is already snapshotted and restored by the
   integrator warm-up, the precompile pass and the restart path.

8. **The denominator floors were used as values where they should be gates
   (3D only).** The floors exist so a field with $\lVert q-\langle q\rangle\rVert = 0$
   gives $0/\mathrm{floor} = 0$ rather than $0/\varepsilon$. But for the
   momenta the floor is an **acoustic** scale, $10^{-3}\bar\rho\bar c = 0.25$
   kg/(m²s), i.e. $\lvert u\rvert = 0.33$ m/s in an atmosphere whose flow
   speeds are 10 m/s. A bubble starting from rest therefore divides a real
   vertical-momentum tendency by a floor 40× below the scale the flow will
   reach. *Fixed* in the 3D workers by setting $\mathrm{denom}_i = \infty$
   while the measured norm is still below the floor, so that equation is
   excluded from the max until its own field develops structure. **Not**
   applied to the 1D or 2D paths: 1D carries no floors at all, and the 2D ones
   are what `sod1d`, `theta_dsgs` and `ffs_step` were tuned against.

9. **The residual was taken at boundary nodes, where it is the boundary
   condition.** The largest single defect in the *picture* the model produces,
   and it affects every path.

   `apply_boundary_conditions_dirichlet!` runs at the **top** of every RHS
   call and projects the wall-normal momentum out of `uaux` at every free-slip
   node. The inviscid RHS then puts it straight back, and the next call
   projects it out again. So at a wall node

   | | |
   |---|---|
   | $\mathrm{BDF2}(q^n, q^{n-1}, q^{n-2})$ | sees the **projected** states |
   | $M^{-1}\mathrm{RHS}$ | does **not** contain the projection |

   and the two differ by the whole projected flux, every step, for ever.

   Measured on `CompEuler/rtb2d_dsgs`: the worst residual in the domain was on
   a boundary node at **every step of the run**, always on the wall-normal
   momentum equation, driving the element coefficient to ≈ 700 m²/s against
   ≈ 14 in the interior. Since the coefficient is constant per element, every
   element touching a wall lit up — a picture of red squares along the walls,
   where there is no gradient at all, with the bubble the model is supposed to
   be tracking two colour-bar decades down and invisible.

   *Fixed*: `params.dsgs_wres` (built in `params_setup.jl` from the mesh's own
   boundary-node lists) masks those nodes out of the per-element $L^\infty$.
   The **denominators are not masked** — they are domain norms of the
   *solution*, which is perfectly well defined on a boundary. An element keeps
   its interior nodes (16 of 25 in 2D at `nop = 4`, 75 of 125 in 3D), and only
   the max over the element is being taken.

   After the fix, on the same case, $\nu_{max}$ drops from ≈ 700–1300 to
   ≈ 130–210 and the worst residual moves to the bubble's edge — which is
   where the initial condition's linear taper $\theta_c(1-r/r_0)$ has its
   slope discontinuity, i.e. the actual under-resolved feature.

> ⚠ Defect 9 (the boundary mask) changes the coefficient on **every** DynSGS
> path, including `CompEuler/sod1d`, which has a committed reference in
> `test/CI-ref/`. That reference has to be regenerated
> (`.github/workflows/generate-ci-ref.yml`, or `test/generate_ci_ref.jl`)
> before the CI suite will pass. Defect 8 is 3D-only. Defects 6 and 7 leave
> the 1D and 2D coefficient unchanged in the default configuration — verified
> by running `sod1d` and reproducing its published $\mu = 1.8\times10^{-3}$ at
> the shock. The measurements in the tables below predate all of this.

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
