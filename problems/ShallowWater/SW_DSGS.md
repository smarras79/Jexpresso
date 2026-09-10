# DynSGS for the shallow-water equations (`:visc_model => DSGS_SW()`)

What the residual-based artificial viscosity does in the shallow-water
cases (`SoliWaveIslandDSGS`), written from the implementation in
`src/kernel/physics/SGS.jl` (`compute_dsgs_viscosity!(::DSGS_SW, ::NSD_2D)`
and its nodal form) and `src/kernel/operators/rhs.jl`. The general theory,
the other systems and the comparison with Dao & Nazarov (2022) are in
[`DSGS.md`](../../DSGS.md) and `docs/DSGS.tex`; this note keeps to what the
shallow-water kernel computes and where it is applied.

## 1. The equations and the terms that receive the viscosity

The solver advances the conserved variables

$$
\mathbf q = (H,\ Hu,\ Hv),
$$

$H$ the water depth above the bathymetry $H_b(x,y)$, $\mathbf v = (u, v)$
the depth-averaged velocity, with the well-balanced split of the case
(`user_flux.jl`, `user_source.jl`): pressure flux $g(H^2 - H_e^2)/2$ and
bathymetry source $-g(H - H_e)\nabla H_b$, where $\mathbf q_e = (H_e, 0, 0)$
is the lake at rest stored as the reference state. The DynSGS term is added
to the right-hand side of every equation as the divergence of a viscous flux
$\mathbf F_V$ (weak form, no boundary integral):

| equation | viscous flux added | primitive it acts on (`user_primitives.jl`) |
|---|---|---|
| $\partial_t H$ | $\nu_1\,\nabla (H - H_e)$ | the depth departure from the lake at rest, so that the rest state stays an exact equilibrium (the depth itself is cone-shaped over the island) |
| $\partial_t (Hu)$ | $\big(\tau_{xx},\ \tau_{xy}\big)$ | $Hu$ |
| $\partial_t (Hv)$ | $\big(\tau_{xy},\ \tau_{yy}\big)$ | $Hv$ |

with the stress built on the conserved momenta $(m_1, m_2) = (Hu, Hv)$ in
the same form the code uses for every 2D system (`_expansion_visc!`):

$$
\tau_{xx} = 2\nu_2\,\partial_x m_1 - \tfrac23\nu_2\,(\partial_x m_1 + \partial_y m_2),\qquad
\tau_{yy} = 2\nu_3\,\partial_y m_2 - \tfrac23\nu_3\,(\partial_x m_1 + \partial_y m_2),\qquad
\tau_{xy} = \nu_{2,3}\,(\partial_y m_1 + \partial_x m_2).
$$

The three coefficients are the **same kinematic viscosity** $\nu$ (m²/s) of
the element (or node), times the deck's per-equation multipliers `:μ`:

$$
\nu_q = \mu_q\,\nu,\qquad q = 1,2,3,\qquad \texttt{:μ => [1.0, 1.0, 1.0]}.
$$

There is no energy equation and no density factor: the shallow-water
momenta are already "$H$-weighted", so a kinematic $\nu$ on $\nabla(H\mathbf v)$
plays the role $\rho\nu\nabla\mathbf v$ plays in the Euler kernels.

## 2. How $\nu$ is built

One value per element (the default, `:ldsgs_nodal => false`). At every
evaluation of the right-hand side (every Runge–Kutta stage):

### 2.1 Residual of each equation

At each node $i$ of element $K$ and for each of the three equations,

$$
R^{q,K}_i = \Big|\ \underbrace{w_1 q_i + w_2 q_{A,i} + w_3 q_{B,i}}_{\partial_t q_i}
\;-\; \frac{\mathrm{rhs}^{q}_{K,i}}{m^K_i}\ \Big| ,
$$

where

- $\partial_t q_i$ is estimated from the current stage state and the stored
  step states $q^n, q^{n-1}, q^{n-2}$ with the stage-consistent stencil
  (BDF2 at the first stage of a step, the second-order three-point formula
  at $\tau = t - t^n$ afterwards; `DSGS.md` §4.4);
- $\mathrm{rhs}^{q}_{K,i}$ is **element $K$'s own** weak inviscid right-hand
  side at node $i$ (its fluxes and bathymetry source, before the assembly
  over elements) and $m^K_i = \omega_i\omega_j J_{K,i}$ its lumped mass
  entry. At interior nodes this is the strong nodal divergence; at element
  interfaces it differs from the assembled rate by the jump of the flux
  divergence, which is the under-resolution the sensor is meant to see.
- At the Dirichlet boundary nodes (the free-slip walls of the basin) the
  residual is set to zero: the wall constrains the assembled rate but not
  the element flux, so the difference there is the constraint force, not a
  resolution error.

### 2.2 Normalization

Each residual is divided by a scale of its variable. The scales are taken on
the **departure from the lake at rest**, $\delta\mathbf q = \mathbf q - \mathbf q_e$,
over the rank's part of the basin (`:dsgs_norms => "rank"`, the default: no
MPI reduction; `"domain"` reduces them over the ranks, the same numbers on
one rank):

$$
n^{q} = \max\Big(\ \max_i\big|\delta q_i - \langle\delta q\rangle\big|,\ 10^{-3}\,s^{q}\Big) + \epsilon,
\qquad
s^{H} = \bar H,\quad s^{Hu} = s^{Hv} = \bar H\sqrt{g\bar H},
$$

$\bar H$ the mean depth, $\epsilon = 10^{-16}$. The departure is used
because the spread of $H$ itself would be the island, not the wave; the
lower bound $10^{-3}s^q$ keeps the ratio finite where a field is uniform or
at rest (the momenta before the wave arrives).

### 2.3 Element residual ratio, wave speed, coefficient

$$
\mathcal R_K = \max_{i\in K}\ \max_{q\in\{H,Hu,Hv\}} \frac{R^{q,K}_i}{n^{q}},
\qquad
\lambda_K = \max_{i\in K}\Big(\sqrt{u_i^2 + v_i^2} + \sqrt{g\max(H_i,0)}\Big),
\qquad
(u_i, v_i) = \frac{(Hu, Hv)_i}{\max(H_i, h_{min})},
$$

with $h_{min}$ = `:dsgs_swe_hmin` (the case's wet/dry threshold, 1 mm) so a
thin film does not produce a velocity $|Hu|/\varepsilon$, and $g$ =
`:dsgs_swe_g`. Then, with $\Delta = \Delta_K/(k+1)$ ($\Delta_K$ the shortest
side of the element, $k$ the polynomial degree),

$$
\boxed{\ \nu_K = \max\Big(\ C_{min}\,\Delta\,\lambda_K,\ \ \max\big(0,\ \min(\,C_{max}\,\Delta\,\lambda_K,\ \ C_R\,\Delta^2\,\mathcal R_K\,)\big)\Big)\ }
$$

| constant | deck key | `SoliWaveIslandDSGS` | role |
|---|---|---|---|
| $C_R$ | `:dsgs_CR` | 1.0 | residual viscosity $C_R\Delta^2\mathcal R$ |
| $C_{max}$ | `:dsgs_Cmax` | 0.5 | first-order (Lax–Friedrichs-like) cap $C_{max}\Delta\lambda$ |
| $C_{min}$ | `:dsgs_Cmin` | 0.05 | background floor, a fraction of the cap: $\nu \ge 0.05\,\Delta\lambda \approx 0.02$ m²/s in the still basin (the sibling `SoliWaveIsland` uses a constant 0.05); it keeps an element-scale transverse mode off the wave front, which the residual only catches once it has grown |

Where the solution is resolved $\mathcal R_K$ is small and $\nu_K$ sits on
the floor; at the wave front, the run-up and the wet/dry ring it rises
towards the cap (0.14 m²/s measured at the run-up on the island).

### 2.4 Nodal form (`:ldsgs_nodal => true`, optional)

The same quantities per node instead of per element: the residual of a node
is the mass-weighted average of its elements' residuals, the normalization
is Dao & Nazarov's eq. 4.7 with `:dsgs_Cl`, $h_i = \Delta_K/(k+1)$ of the
elements containing the node, $\lambda_i$ the nodal wave speed, and
$\nu_i$ is a continuous field interpolated inside each element. The default
of the case is the element form.

## 3. Output

The coefficient actually applied is written at every output time:
`μ_dsgs_H-it<n>.png` (the three slots carry the same $\nu$, so one panel;
`:plot_dsgs_vars => ["H"]`), or the field `mu_dsgs` in the `.pvtu` files
with `:outformat => "vtk"`.

## 4. Deck keys of the shallow-water kernel

```julia
:visc_model    => DSGS_SW(),
:μ             => [1.0, 1.0, 1.0],   # per-equation multipliers of ν
:dsgs_CR       => 1.0,
:dsgs_Cmax     => 0.5,
:dsgs_Cmin     => 0.05,
:dsgs_norms    => "rank",            # scales over this rank's elements (default); "domain": the whole basin
:dsgs_sensor   => "residual",        # element-wise residual (the legacy sensor is also available)
:ldsgs_nodal   => false,             # element form; true = nodal form, with :dsgs_Cl
:dsgs_swe_g    => 9.81,              # g of the wave speed  (= _G_SWE of user_flux.jl)
:dsgs_swe_hmin => 1.0e-3,            # wet/dry threshold    (= _H_WET_SWE of user_flux.jl)
```
