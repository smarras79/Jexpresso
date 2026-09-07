# The 2D ideal GLM-MHD equations and the flux-emergence problem of Son, Jang & Magara (2025)

This document describes the system of equations solved by this case and the
test problem it is configured for. The test follows

> D. Son, Y. Jang, T. Magara,
> *A Comparative Analysis of High-resolution Shock-capturing Schemes for
> Two-dimensional Magnetohydrodynamic Simulation of Flux Emergence in the
> Solar Atmosphere*, ApJS **277**:46 (2025), Sections 2.1, 2.2 and 4.1.
> https://doi.org/10.3847/1538-4365/adb617 (open access; a copy is in
> `problems/MHD/Son_2025_ApJS_277_46.pdf`)

which is itself the classical two-dimensional emerging-flux model of Shibata
et al. (1989a, ApJ 345, 584): the nonlinear evolution of the Parker (undular
magnetic-buoyancy) instability of an isolated horizontal flux sheet in a
two-temperature (chromosphere + corona) stratified atmosphere, and its
self-similar expansion into the corona.

## 1. Nondimensional units

All quantities are normalized to their photospheric/chromospheric values
(paper Table 1): length by the pressure scale height $H_0 = k_BT_0/(mg_0)$,
velocity by the adiabatic sound speed $C_s = (\gamma k_BT_0/m)^{1/2}$, time
by $\tau_0 = H_0/C_s$, density by $\rho_0$, pressure by $p_0 = \rho_0C_s^2$,
temperature by $T_0 = mC_s^2/(\gamma k_B)$, magnetic field by
$B_0 = (\rho_0C_s^2)^{1/2}$ and gravity by $g_0 = C_s^2/(\gamma H_0)$. The
paper uses Heaviside–Lorentz units for $\mathbf{B}$ (magnetic pressure
$\tfrac12|\mathbf{B}|^2$, no $4\pi$), and so does this case. With
$H_0 = C_s = \rho_0 = 1$:

$$
p = \frac{\rho T}{\gamma},\qquad g_0 = \frac{1}{\gamma},\qquad \gamma = 1.05 .
$$

Jexpresso's $y$ is the paper's vertical coordinate $z$; the code's `v` is
the paper's $V_z$ and `By` its $B_z$.

## 2. Governing equations (paper Eqs. 6–8)

$$
\frac{\partial \mathbf{U}}{\partial t} + \nabla\cdot\mathbf{F} = \mathbf{S},
\qquad
\mathbf{U} = (\rho,\ \rho u,\ \rho v,\ E,\ \rho w,\ B_x,\ B_y,\ B_z,\ \psi)^T
$$

with the Dedner GLM fluxes (the same `user_flux.jl` as the Orszag–Tang and
Kelvin–Helmholtz cases — see their `EQUATIONS.md` for the full flux
vectors), the pressure

$$
p = (\gamma-1)\Big(E - \tfrac12\rho|\mathbf{v}|^2 - \tfrac12|\mathbf{B}|^2 - \tfrac12\psi^2\Big),
$$

and the source vector of the paper's Eq. (7),

$$
\mathbf{S} = \Big(0,\ 0,\ -\rho g_0,\ -\rho v g_0,\ 0,\ 0,\ 0,\ 0,\ -\frac{c_h^2}{c_p^2}\psi\Big)^T ,
$$

i.e. gravity $\mathbf{g} = (0, -g_0)$ in the vertical momentum and in the
energy ($\rho\,\mathbf{V}\cdot\mathbf{g}$), plus the parabolic damping of
the GLM field. Jexpresso's slot ordering puts $E$ in slot 4 and $\rho w$ in
slot 5 (the shared 2D kernels assume the energy lives in slot 4).

**γ = 1.05.** The paper deliberately uses a near-isothermal index: it raises
the growth rate of the undular mode and mimics the quasi-isothermal
chromosphere with limited heating of the corona. A consequence worth
knowing: with $\gamma - 1 = 0.05$ the gas pressure is a small residue of
the total energy wherever $\beta \ll 1$ (the emerged loop reaches
$\beta \sim 10^{-4}$, paper Fig. 6(e)), which makes the conservative
formulation prone to negative pressures. `user_flux.jl` floors $p$ inside
the flux at $10^{-9}$ (1% of the smallest pressure of the initial
atmosphere); the floor is inactive in a healthy run.

**GLM cleaning.** $c_h$ is the maximum $|\mathbf{v}| + c_f$ of the initial
condition (paper Eq. 11) — here the coronal sound speed
$\sqrt{T_{cor}/T_0} = 5\,C_s$ — kept constant. The damping is parametrized
as in the paper (Eq. 12, Mignone & Tzeferacos 2010) by
$\alpha_p = \Delta h\,c_h/c_p^2$, so $c_h^2/c_p^2 = \alpha_p c_h/\Delta h$;
the paper adopts $\alpha_p = 0.2$ for its WENO schemes, and $\Delta h$ is
the smallest LGL nodal spacing of the mesh (measured in `initialize.jl`).

## 3. Initial condition (paper Section 2.1)

Domain $[0, X_{max}]\times[0, Z_{max}] = [0, 80]\times[0, 35]$ (in $H_0$).

**Two-temperature atmosphere (Eq. 1)**

$$
T(z) = T_{ch} + \frac{T_{cor}-T_{ch}}{2}\left[\tanh\!\left(\frac{z-z_{cor}}{w_{tr}}\right)+1\right],
\qquad T_{ch} = T_0,\ T_{cor} = 25T_0,\ z_{cor} = 18H_0,\ w_{tr} = 0.6H_0 .
$$

**Magnetic flux sheet (Eqs. 2–3)**, parallel to $x$, inside the cold layer:

$$
B_x(z) = \left[\frac{2p(z)}{\beta(z)}\right]^{1/2},\qquad
\beta(z) = \frac{\beta_*}{f(z)},\qquad
f(z) = \frac14\left[1+\tanh\!\left(\frac{z-z_0}{w_0}\right)\right]\left[1-\tanh\!\left(\frac{z-z_1}{w_1}\right)\right],
$$

$z_0 = 4H_0$, $z_1 = z_0 + D = 8H_0$, $w_0 = w_1 = 0.5H_0$ (the paper writes
Eq. 2 in Gaussian units, $[8\pi p/\beta]^{1/2}$; the $2p/\beta$ above is
the same field in the Heaviside–Lorentz units of the code, i.e. the magnetic
pressure is $\tfrac12B_x^2 = p f/\beta_*$).

**Plasma beta $\beta_*$ at the sheet center: not printed by the paper.** It
was inferred from the paper's Fig. 1(b) by integrating the equilibrium below
for a range of values:

| $\beta_*$ | $\max B_x/B_0$ | $\log_{10}\rho$ at $z=8.5$ | at $z=20$ | at $z=35$ | $\log_{10}p$ at $z=35$ |
|---|---|---|---|---|---|
| 0.5 | 0.153 | −2.49 | −7.53 | −7.79 | −6.41 |
| **1.0** | **0.120** | **−2.81** | **−7.90** | **−8.16** | **−6.79** |
| 2.0 | 0.091 | −3.11 | −8.24 | −8.50 | −7.12 |
| paper Fig. 1(b) | ≈ 0.113 | ≈ −2.8 | ≈ −7.9 | ≈ −8.2 | ≈ −6.85 |

$\beta_* = 1$ reproduces the figure within reading accuracy and is the
standard case of Shibata et al. (1989a), whose expansion-law constants the
paper adopts (its Eqs. 40–43). It is a `Ref` (`fe_beta_star`) in
`initialize.jl`.

**Magnetostatic equilibrium (Eq. 4)**

$$
\frac{dp_{total}}{dz} = -\rho g_0,\qquad p_{total} = p + \tfrac12B_x^2 = p\left(1+\frac{f}{\beta_*}\right),\qquad \rho = \frac{\gamma p}{T}.
$$

With $P = p_{total}$ this is $dP/dz = -P/[T(z)(1+f/\beta_*)]$, integrated
numerically (trapezoids on a $10^{-4}H_0$ grid) from
$P(0) = p(0)(1+f(0)/\beta_*)$ with $p(0) = 1/\gamma$, i.e. $\rho(0) = 1$ as
in the paper's Fig. 1(b). Then $p = P/(1+f/\beta_*)$, $\rho = \gamma p/T$
and $B_x = (2pf/\beta_*)^{1/2}$. The resulting profiles: $\rho$ and $p$
fall by 8 and 6.8 decades from the photosphere to the top, the sheet's field
peaks at $B_x/B_0 = 0.12$ at $z = 4.3H_0$ (Alfvén speed $1.4\,C_s$), the
corona is unmagnetized with sound speed $5\,C_s$.

**Perturbation (Eq. 5)**, inside the sheet over the central region
$X_{max}/2 - \lambda/2 < x < X_{max}/2 + \lambda/2$:

$$
V_x = f(z)\,A\,C_s\sin\!\left[\frac{2\pi(x - X_{max}/2)}{\lambda}\right],\qquad A = 0.05,\ \lambda = 20H_0 .
$$

$V_z = w = B_z = \psi = 0$.

## 4. Boundary conditions

- **Horizontal**: periodic (mesh-level `periodicx` tags).
- **Bottom, $z = 0$**: symmetric. For the in-plane fields of a horizontal
  sheet this is the reflecting wall $V_z = B_z = 0$; implemented in
  `user_bc.jl` as the free-slip projection of the momentum and of
  $\mathbf{B}$.
- **Top, $z = 35H_0$**: the paper's boundary is "free with an absorbing
  layer (Machida & Matsumoto 2003)". The absorbing layer is a Rayleigh
  damping of the departure from the initial magnetostatic state,

  $$
  \mathbf{S} \mathrel{-}= \sigma(z)\,(\mathbf{U} - \mathbf{U}_e),\qquad
  \sigma(z) = \sigma_{max}\sin^2\!\left[\frac{\pi}{2}\frac{z - z_s}{Z_{max}-z_s}\right],\quad z > z_s ,
  $$

  with $z_s = 30H_0$ and $\sigma_{max} = 2/\tau_0$ (`user_source.jl`,
  both `Ref`s): a coronal wave ($5\,C_s$) crosses the $5H_0$ layer in one
  $\tau_0$ and is damped by $e^{-2}$ before it can reflect. Behind it the
  edge itself is a free-slip wall (normal momentum removed), the robust
  choice of Jexpresso's stratified-atmosphere cases; the strong-form CG
  discretization has no boundary flux that a "free" condition could
  prescribe.

## 5. Discretization

| | paper | this case |
|---|---|---|
| space | finite volume, 5th-order WENO-Z+M / WENO-NZ / IMWENO-P / TENO-LAD reconstructions, HLLD Riemann solver | continuous-Galerkin spectral elements, LGL, $N = 4$ |
| grid | $300^2$ to $2400^2$ uniform cells | $80\times35$ elements of $1H_0\times1H_0$, $320\times140$ unique points (0.25 $H_0$ mean, 0.17 $H_0$ min spacing) |
| time | SSP-RK3, Courant number 0.23 | Carpenter–Kennedy 2N54 (low-storage RK4), $\Delta t = 2.5\times10^{-3}\tau_0$ |
| shocks | the Riemann solver's dissipation | Marras–Nazarov DynSGS (`DSGS_MHD`, per-element norms, nodal density), residual-based eddy viscosity on momentum and $\mathbf{B}$, none on the energy (see README) |
| $\nabla\cdot\mathbf{B}$ | GLM, $c_h$ per step, $\alpha_p = 0.2$ | GLM, $c_h$ of the IC, $\alpha_p = 0.2$ |
| final time | $54\tau_0$ (Fig. 2 at $51\tau_0$) | $54\tau_0$, output every $\tau_0$ |

The resolution is the coarsest the transition widths admit (see
`FE_80x35.geo`): a $0.5H_0$ tanh gets three LGL nodes. It is also two times
coarser vertically than the paper's coarsest grid, so expect the contact
surface at the top of the loop and the current sheets at its feet to be
more diffuse than in the paper's $300^2$ panels; `ny = 70` in the `.geo`
matches the paper's vertical spacing at twice the cost.

**Time step.** The initial maximum wave speed is the coronal sound speed,
$5\,C_s$ ($= c_h$); at late times the loop's Alfvén speed reaches
$4$–$7\,C_s$ and the lateral downflows $4$–$5\,C_s$ (paper Sec. 4.1), so
$|\mathbf{v}| + c_f \approx 10\,C_s$ and $\Delta t = 2.5\times10^{-3}$ is a
Courant number of 0.15 on the smallest LGL spacing. 21,600 steps to
$t = 54\tau_0$.

## 6. What to expect (paper Section 4.1, Figs. 2, 5, 6)

1. **Linear phase, $t \lesssim 30\tau_0$**: the $\lambda = 20H_0$
   perturbation grows into an undulation of the sheet; gas slides down the
   rising field lines into the magnetic valleys at $x \approx 25$ and
   $55H_0$.
2. **Emergence, $t \approx 33$–$47\tau_0$**: the crest breaks through the
   transition region at $z = 18H_0$ and expands into the corona, driven by
   the magnetic pressure gradient. Along the centerline $x = 40H_0$
   (paper Fig. 5) the rise velocity follows Shibata's expansion law
   $V_z/C_s = a_1\Delta z$, $a_1 = 0.062$, the Alfvén speed
   $V_A/C_s = a_2\Delta z$, $a_2 = 0.3$, and $\rho \propto \Delta z^{-4}$,
   $B_x \propto \Delta z^{-1}$, with $\Delta z = (z - z_0)/H_0$.
3. **$t = 51\tau_0$** (paper Fig. 2): a loop spanning $x \approx 15$–$65H_0$
   and reaching $z \approx 23H_0$, with the coronal density pushed up to
   $z \approx 27H_0$ ($\log_{10}\rho \approx -5$ inside the loop against
   $-8$ in the ambient corona), dense pockets ($\log_{10}\rho \approx -3$)
   at the loop feet at $z \approx 6$–$8H_0$, lateral downflows of
   $4$–$5\,C_s$ along the loop sides where fast/intermediate shocks form,
   and field-aligned downflows of $2$–$3\,C_s$ near the footpoints.
4. **Late time, $t = 52$–$54\tau_0$**: the top of the loop approaches the
   absorbing layer; $V_z$ peaks at $\approx 1.25\,C_s$ at $z \approx 27H_0$,
   $V_A$ at $\approx 4\,C_s$, $\beta$ inside the loop $\sim 10^{-4}$.

Departures to expect from a CG-SEM + DynSGS run at this resolution: the
contact surface at the loop top is smeared over $\sim 1H_0$; the current
sheets at the feet and the Gibbs-free shocks of the paper's IMWENO-P are
regularized by the eddy viscosity (which also acts on $\mathbf{B}$ as a
turbulent resistivity, so the footpoint field will diffuse slightly); the
Alfvén-speed peak is set by the density evacuation and should be
reproduced; the $\psi$ field should stay $O(10^{-3})$ and concentrated at
the loop feet.

## 7. Output

`:outformat => "png"`: the solver writes, at every $\tau_0$,

- `ρ-it<n>.png` — $\log_{10}(\rho/\rho_0)$ on the `jet` scale $[-8.5, 0]$
  of the paper's Fig. 2, with magnetic field lines (isocontours of the
  vector potential $A_y$, black) and velocity vectors (white, reference
  arrow $= 5\,C_s$);
- `v`, `vA`, `Bx`, `p`, `T`, `β` panels (`p` and `β` in $\log_{10}$);
- `profile-it<n>.png` — $V_z/C_s$, $V_A/C_s$, $\log_{10}(B_x/B_0)$,
  $\log_{10}(\rho/\rho_0)$ along $x = X_{max}/2$ on the axes of the paper's
  Fig. 5, with $z_{cor} = 18H_0$ marked.

`:outformat => "vtk"` writes the same output variables plus the
`mu_dsgs_<var>` DynSGS fields for ParaView.
