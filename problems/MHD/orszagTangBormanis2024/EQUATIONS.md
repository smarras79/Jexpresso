# The 2D ideal GLM-MHD equations and the Orszag–Tang vortex

This document describes the system of equations solved by this case and the
test problem it is configured for. The test follows

> A. Bormanis, C. A. Leon, A. Scheinker,
> *Solving the Orszag–Tang vortex magnetohydrodynamics problem with
> physics-constrained convolutional neural networks*,
> Phys. Plasmas **31**, 012101 (2024), Section I A (governing equations and
> initial condition) and Section III (simulated data).
> https://doi.org/10.1063/5.0172075

Implementation notes (file layout, variable ordering, how to run) are in
[README.md](README.md).

---

## 1. Governing equations

### 1.1 What the paper solves

Bormanis et al. solve **ideal MHD** — perfect conductivity, no magnetic
diffusivity, no viscosity — in the conservative form of their Appendix B:

$$
\frac{\partial \rho}{\partial t} = -\nabla\cdot(\rho\mathbf{v}),
$$

$$
\frac{\partial}{\partial t}(\rho\mathbf{v})
 = -\nabla\cdot\left(\rho\mathbf{v}\mathbf{v}^T + p\,\mathbf{I} - \mathbf{B}\mathbf{B}^T\right),
$$

$$
\frac{\partial}{\partial t}(\rho e)
 = -\nabla\cdot\left(\rho e\,\mathbf{v} + p\,\mathbf{v} - \mathbf{B}(\mathbf{v}\cdot\mathbf{B})\right),
$$

$$
\frac{\partial \mathbf{B}}{\partial t}
 = -\nabla\cdot\left(\mathbf{B}\mathbf{v}^T - \mathbf{v}\mathbf{B}^T\right),
\qquad \nabla\cdot\mathbf{B} = 0,
$$

with total specific energy
$e = u + \tfrac{1}{2}\lVert\mathbf{v}\rVert^2 + \tfrac{1}{2}\lVert\mathbf{B}\rVert^2$,
total pressure $p = p_{\text{fluid}} + \tfrac{1}{2}\lVert\mathbf{B}\rVert^2$
and the monatomic perfect-gas closure $p_{\text{fluid}} = (\gamma-1)\rho u$
with $\gamma = 5/3$ (their Section I A). Units are Gaussian-with-$\sqrt{4\pi}$
absorbed, i.e. the magnetic pressure is $\lVert\mathbf{B}\rVert^2/2$ and the
Alfvén speed is $\lVert\mathbf{B}\rVert/\sqrt{\rho}$.

The paper maintains $\nabla\cdot\mathbf{B}=0$ with **constrained transport**
on a staggered finite-volume mesh (their Section III A and Appendix B): the
face-averaged field is updated from corner electric fields so that a specific
discrete divergence is conserved to machine precision.

### 1.2 What this case solves

Constrained transport has no analogue in a collocated **continuous-Galerkin
spectral element** discretization: there is no staggered grid and no face
average to update. Jexpresso therefore enforces the divergence constraint the
way the DG/SEM MHD literature does — with Dedner's **generalized Lagrange
multiplier (GLM)** reformulation, which adds one scalar field $\psi$ that
carries divergence error out of the domain at a constant speed $c_h$. This is
the *same equation set already implemented for*
[`problems/MHD/kelvinHelmholtzChan2022`](../kelvinHelmholtzChan2022/EQUATIONS.md),
reused here unchanged; only the initial condition, the domain and the
run parameters are new.

**Notation warning.** From here on $p$ denotes the **gas** pressure (the
paper's $p_{\text{fluid}}$) and the total pressure is written $p_{tot}$. In
§1.1 above, following the paper, $p$ was the *total* pressure. The two
formulations describe the same fluxes.

The system is

$$
\frac{\partial \mathbf{u}}{\partial t} + \nabla\cdot\mathbf{f}(\mathbf{u}) = \mathbf{s},
\qquad
\mathbf{u} = \left(\rho,\; \rho\mathbf{v},\; E,\; \mathbf{B},\; \psi\right)^T,
$$

with $\mathbf{v} = (u,v,w)$ and $\mathbf{B} = (B_x,B_y,B_z)$ — nine unknowns.
Following the standard 2D reduction, the flux in $z$ is dropped but the third
components $w$ and $B_z$ are **kept**, because plasmas admit three-dimensional
electromagnetic interactions even in two-dimensional problems. (For
Orszag–Tang specifically $w \equiv B_z \equiv 0$ for all time; the slots are
carried for generality.)

With the total (gas + magnetic) pressure
$p_{tot} = p + \tfrac{1}{2}\lVert\mathbf{B}\rVert^2$, the Cartesian fluxes are

$$
\mathbf{F}_1 =
\begin{bmatrix}
\rho u \\
\rho u^2 + p_{tot} - B_x^2 \\
\rho u v - B_x B_y \\
u\left(\tfrac{1}{2}\rho\lVert\mathbf{v}\rVert^2 + \tfrac{\gamma p}{\gamma-1} + \lVert\mathbf{B}\rVert^2\right) + B_x\left(c_h\psi - \mathbf{v}\cdot\mathbf{B}\right) \\
\rho u w - B_x B_z \\
c_h \psi \\
u B_y - v B_x \\
u B_z - w B_x \\
c_h B_x
\end{bmatrix},
\qquad
\mathbf{F}_2 =
\begin{bmatrix}
\rho v \\
\rho v u - B_y B_x \\
\rho v^2 + p_{tot} - B_y^2 \\
v\left(\tfrac{1}{2}\rho\lVert\mathbf{v}\rVert^2 + \tfrac{\gamma p}{\gamma-1} + \lVert\mathbf{B}\rVert^2\right) + B_y\left(c_h\psi - \mathbf{v}\cdot\mathbf{B}\right) \\
\rho v w - B_y B_z \\
v B_x - u B_y \\
c_h \psi \\
v B_z - w B_y \\
c_h B_y
\end{bmatrix},
$$

written in the row order $(\rho, \rho\mathbf{v}, E, \mathbf{B}, \psi)$ used by
the MHD literature. (In the code the total energy is stored in slot 4 and
$\rho w$ in slot 5 — see [README.md](README.md) for why.)

Setting $\psi \equiv 0$ and $c_h = 0$ recovers exactly the paper's
Eqs. (B3)–(B6). The two GLM couplings — $c_h\psi\,\mathbf{I}$ in the induction
rows and $c_h\mathbf{B}$ in the $\psi$ row — make $(\nabla\cdot\mathbf{B},\psi)$
satisfy a linear wave system with speed $c_h$, so divergence errors are
radiated away instead of accumulating.

### 1.3 Equation of state

$$
p = (\gamma - 1)\left(E - \tfrac{1}{2}\rho\lVert\mathbf{v}\rVert^2
      - \tfrac{1}{2}\lVert\mathbf{B}\rVert^2 - \tfrac{1}{2}\psi^2\right),
\qquad \gamma = \tfrac{5}{3}.
$$

$\gamma$ is set in `user_flux.jl` (`γ_mhd`). Note that the initial condition
below is written **in terms of $\gamma$**, so changing it changes $\rho$ and
$p$ as well.

### 1.4 Divergence-cleaning speed $c_h$

$c_h$ is set to the maximum wave speed of the initial condition and kept
constant for the whole run (standard GLM practice):

$$
c_h = \max_{\Omega}\left(\lVert\mathbf{v}\rVert + c_f^{\max}\right),
\qquad
c_f^{\max} = \sqrt{\frac{\gamma p}{\rho} + \frac{\lVert\mathbf{B}\rVert^2}{\rho}},
$$

where $c_f^{\max}$ bounds the fast magnetosonic speed over all propagation
directions. It is computed over the mesh in `initialize.jl` (MPI-reduced) and
stored in the `c_h_mhd` reference read by `user_flux!`. For the Orszag–Tang
initial condition it evaluates to $c_h \approx 2.6031$.

**ψ damping (mixed cleaning).** Purely hyperbolic cleaning only *transports*
divergence errors; nothing destroys $\psi$. In a DG setting the interface
Riemann fluxes upwind-damp the $\psi$ waves, but a continuous Galerkin
discretization has no interface dissipation, so on a periodic domain the
$\psi$ waves keep bouncing and alias into grid-scale (LGL odd–even) noise. We
therefore use Dedner's *mixed* (hyperbolic–parabolic) GLM and damp $\psi$
with the source

$$
S_\psi = -\frac{c_h^2}{c_p^2}\,\psi = -\frac{c_h}{c_r}\,\psi,
\qquad c_p^2 = c_h c_r,\quad c_r = 0.18 \ \text{(Dedner's value)},
$$

implemented in `user_source.jl` (`glm_cr_mhd` sets $c_r$). Only $\psi$ is
damped: since $E$ carries $\tfrac{1}{2}\psi^2$, removing $\psi$ while leaving
$E$ untouched converts the cleaned $\psi$-energy into heat, which is the
entropy-consistent behavior (Derigs et al. 2018).

This is the **only source term** in the case: the ideal MHD equations of the
paper have none.

### 1.5 The non-conservative term (NOT included)

The full non-conservative GLM-MHD system adds

$$
\Upsilon = (\nabla\cdot\mathbf{B})
\left(0,\; \mathbf{B},\; \mathbf{v}\cdot\mathbf{B},\; \mathbf{v},\; 0\right)
+ \left(0,\; \mathbf{0},\; \psi\,(\mathbf{v}\cdot\nabla\psi),\; \mathbf{0},\; \mathbf{v}\cdot\nabla\psi\right),
$$

i.e. the Powell term plus a Galilean GLM correction. It vanishes for exactly
divergence-free fields and is required for **entropy stability** and Galilean
invariance of split-form / two-point-flux discretizations. Since this case
uses neither the entropy-stable nor the kinetic-energy-preserving machinery
(see below), $\Upsilon$ is **omitted**; the conservative GLM system above is
the standard Dedner-type formulation and is well-posed on its own.

## 2. Discretization and stabilization

- Continuous spectral elements, LGL nodes, polynomial degree `:nop => 4` on a
  $32\times32$ element grid — i.e. $32 \times 4 = 128$ unique points per
  direction, **exactly the $128\times128$ resolution of the reference data**
  of the paper. Inexact (collocated) quadrature.
- Explicit low-storage RK: `CarpenterKennedy2N54()` with
  $\Delta t = 5\times10^{-4}$. With $c_h \approx 2.603$ and the smallest LGL
  spacing $\Delta x_{\min} \approx 5.40\times10^{-3}$ this is
  $\mathrm{CFL} \approx 0.24$. The margin is deliberate: the vortex steepens
  into shocks and local wave speeds grow. (The paper's finite-volume reference
  used $\Delta t = 8\times10^{-4}$ on the same $128^2$ grid.)
- **No entropy-stable / KEP flux differencing** (`:lkep => false`,
  `:entropy_variables => false`): the inviscid terms are the plain weak-form
  divergence of the pointwise fluxes above.
- **Stabilization: Smagorinsky LES viscosity** (`:visc_model => SMAG()`),
  with turbulent viscosity

$$
\mu_t = \rho\, C_s^2\, \Delta^2\, |S|,
\qquad
|S| = \sqrt{2\, S_{ij} S_{ij}},
\quad
S_{ij} = \tfrac{1}{2}\left(\partial_j v_i + \partial_i v_j\right),
$$

  applied per equation through the per-equation multipliers
  `:μ = [0, 8, 8, 8, 8, 8, 8, 8, 8]` (no mass diffusion):

  | equation | SGS term |
  |---|---|
  | $\rho u,\ \rho v$ | full deviatoric stress $\nabla\cdot\boldsymbol{\tau}$, $\ \boldsymbol{\tau} = 2\mu_e \mathbf{S} - \tfrac{2}{3}\mu_e (\nabla\cdot\mathbf{v})\mathbf{I}$ |
  | $E$ | heat flux $\nabla\cdot(\kappa_e \nabla T)$ + viscous work $\nabla\cdot(\boldsymbol{\tau}\,\mathbf{v})$, $\ T = p/\rho$ |
  | $\rho w$ | scalar diffusion $\nabla\cdot(\kappa_e \nabla w)$ |
  | $B_x, B_y, B_z$ | scalar diffusion of each component — an effective **turbulent resistivity** |
  | $\psi$ | scalar diffusion |

  with $\mu_e = \mu_{mol} + \mu_t$ and the eddy diffusivities
  $\kappa_e \sim \mu_t/(\rho\,\mathrm{Pr}_t)$ (energy) and
  $\mu_t/(\rho\,\mathrm{Sc}_t)$ (scalars). No Richardson correction
  (`:lrichardson => false`; there is no gravity in this problem) and no
  spectral filter (`:lfilter => false`).

  **This is the one substantive numerical departure from the paper.** The
  reference is an ideal (inviscid, non-resistive) finite-volume computation
  whose dissipation comes from the Riemann solver at cell faces; here the
  scheme is dissipation-free by construction and Smagorinsky supplies the
  regularization that the Orszag–Tang shocks require. Smagorinsky is an
  eddy-viscosity closure, not a shock-capturing scheme, so shocks are
  smeared over a few points rather than resolved sharply. See §4.

  The 8× multiplier is **not** a physically calibrated LES constant. On this
  grid the bare $\mu_t = \rho C_s^2\Delta^2|S|$ is $O(10^{-6}$–$10^{-4})$,
  far too small to keep the gas pressure positive through the shocks that
  form at $t\approx 0.5$: with `:μ = 1` (or 2) the run aborts at
  $t \approx 0.55$ when $p = (\gamma-1)(E - \frac12\rho|v|^2 -
  \frac12|B|^2)$ goes negative. `:μ = 8` corresponds to an effective
  $C_s \approx \sqrt{8}\,(0.23) \approx 0.65$. The measured stability
  boundary and the alternatives (including why the Boyd–Vandeven filter is
  *not* used) are tabulated in
  [README.md](README.md#stabilization-what-was-actually-tested).

## 3. Test: 2D Orszag–Tang vortex

The Orszag–Tang vortex (Orszag & Tang, JFM 90:129, 1979) is the standard 2D
MHD benchmark: a smooth, doubly periodic initial state whose nonlinear
evolution generates interacting shocks and a magnetic-reconnection-driven
transition to MHD turbulence.

**Domain.** $\Omega = [0,1]^2$, doubly periodic. Final time $T_{final} = 1$
(the interval over which the paper generated its data set).

**Initial condition** (paper Eqs. 6–9), with $\gamma = 5/3$:

$$
\rho(\mathbf{r},0) = \frac{\gamma^2}{4\pi},
\qquad
p(\mathbf{r},0) = \frac{\gamma}{4\pi},
$$

$$
\mathbf{v}(\mathbf{r},0) = \bigl(-\sin(2\pi y),\; \sin(2\pi x)\bigr),
$$

$$
\mathbf{B}(\mathbf{r},0) = \frac{1}{2\pi\sqrt{4\pi}}\,\nabla_{\!\perp}
\left(\frac{\cos(4\pi x)}{2} + \cos(2\pi y)\right),
\qquad
\nabla_{\!\perp} = \left(\frac{\partial}{\partial y},\; -\frac{\partial}{\partial x}\right).
$$

Carrying out the curl gives the familiar closed form used in `initialize.jl`:

$$
B_x = -B_0\sin(2\pi y), \qquad B_y = B_0\sin(4\pi x),
\qquad B_0 = \frac{1}{\sqrt{4\pi}},
$$

which is $\mathbf{B} = \nabla\times(A_z\hat{\mathbf{z}})$ with
$A_z = B_0\left(\frac{\cos(4\pi x)}{4\pi} + \frac{\cos(2\pi y)}{2\pi}\right)$.
Because $B_x$ depends only on $y$ and $B_y$ only on $x$, the initial field is
**exactly** divergence free pointwise — not merely to truncation — so
$\psi$ starts at, and stays near, zero.

Numerically:

| quantity | value |
|---|---|
| $\rho = 25/(36\pi)$ | $0.2210485$ |
| $p = 5/(12\pi)$ | $0.1326291$ |
| $B_0 = 1/\sqrt{4\pi}$ | $0.2820948$ |
| sound speed $a=\sqrt{\gamma p/\rho}$ | exactly $1$ |
| $\max\lVert\mathbf{v}\rVert$ | $\sqrt{2}$ |
| $\max\lVert\mathbf{B}\rVert$ | $\sqrt{2}\,B_0 = 0.3989$ |
| minimum plasma $\beta = p/(\tfrac{1}{2}\lVert\mathbf{B}\rVert^2)$ | $5/3$ |
| $c_h$ | $2.6031$ |

The remaining fields are $w = B_z = \psi = 0$, and

$$
E = \frac{p}{\gamma-1} + \frac{1}{2}\rho\lVert\mathbf{v}\rVert^2
  + \frac{1}{2}\lVert\mathbf{B}\rVert^2 + \frac{1}{2}\psi^2 .
$$

Note that although the paper calls $p$ the *total* pressure in its Section
I A, Eq. (7) is the classical Orszag–Tang value $5/(12\pi)$ of the **gas**
pressure (Picone & Dahlburg, its Ref. 16; Mocz et al., its Ref. 17), and that
is how it is used here — consistently with $a = 1$ and $\beta_{\min} = 5/3$,
the standard values quoted for this benchmark.

**What to expect.** The velocity field is a single large-scale cellular
vortex; the magnetic field has twice the periodicity in $x$, so the flow and
the field are misaligned from the start. The sequence over $t \in [0,1]$:

- $t \lesssim 0.2$: smooth compression. Density and pressure develop a
  four-lobe pattern; the field lines are wound up by the vortex.
- $t \approx 0.3$–$0.5$: the compressions steepen into shocks that cross the
  domain diagonally, producing the characteristic X-shaped density pattern.
  Current sheets form along the field reversals at $x = 1/4, 3/4$.
- $t \approx 0.5$–$1.0$: shock–shock and shock–current-sheet interactions;
  the central current sheet tears and reconnects, forming the magnetic island
  at the domain center that is the classic signature of this test. Small-scale
  structure fills the domain — the paper describes $t \gtrsim 0.75$ as "a
  turbulent part of the Orszag–Tang vortex". Compare the density snapshots in
  Figure 3 and the $B$ snapshots in Figures 5–8 of the paper.

A quantitative check: the "True" row of the paper's Figure 3 shows $\rho$ at
$t = 0.6992,\ 0.7992,\ 0.8992,\ 0.9992$ on a color scale spanning
$\approx 0.1$ to $0.4$ — i.e. the density stays within roughly a factor of two
either side of its uniform initial value $0.2210$ throughout the run. A run
that drifts far outside that band, or that produces negative $\rho$ or $p$, is
under-resolved or under-damped.

Nothing in the 2D system couples $w$ and $B_z$ once they vanish, so those two
output fields should stay at machine zero for the whole run — a free
consistency check.

**Diagnostics written** (VTK, every $0.05$ time units by default):
$\rho, u, v, w, p, B_x, B_y, B_z, \psi, T$.

## 4. Expected agreement with the paper — and its limits

This case reproduces the paper's **physical setup exactly**: same equations
(up to the GLM reformulation of the divergence constraint), same $\gamma$,
same initial condition, same domain, same periodicity, same resolution, same
time interval. It does **not** reproduce the paper's *numerics*: the
reference is a staggered-mesh finite-volume scheme with constrained transport,
this is a collocated CG spectral element method with GLM cleaning and
Smagorinsky regularization.

Consequences to expect, in decreasing order of importance:

1. **Shocks are smeared, not captured.** Smagorinsky viscosity is an
   eddy-viscosity closure. It keeps the run stable but spreads the
   Orszag–Tang shocks over several points and leaves some Gibbs ringing at
   the steepest fronts. Large-scale structure (the X pattern, the central
   island, the overall $\rho$ and $B$ morphology) should match the paper
   well; peak values across shocks will not match to the digit.
2. **Small-scale structure differs in character, not obviously in amount.**
   Comparing the computed density at $t = 0.70,\ 0.80,\ 0.90,\ 1.00$ against
   the "True" row of the paper's Figure 3 (same $[0.1, 0.4]$ color scale —
   see `assets/MHD_OT_rho.png`): the large-scale features agree well, with
   the same diagonal banding, the same dark high-density blobs near the
   center at $t \approx 0.8$, and the same overall amplitude range. But this
   solution is *more* filamentary than the reference, not less — the
   4th-order SEM operator resolves thin density filaments that the reference
   finite-volume scheme smooths into rounder blobs, despite the heavy SGS
   viscosity. Do not expect this case to look uniformly smoother than the
   paper's figures; expect it to look sharper in the filaments and smeared
   across the shocks.
3. **$\nabla\cdot\mathbf{B}$ is controlled, not zero.** Constrained transport
   holds a discrete divergence at machine precision; GLM cleaning instead
   holds it small and bounded. Expect $\psi$ to stay near zero away from the
   current sheets and to show bounded, saturating speckle at them, where the
   CG discretization genuinely produces grid-scale $\nabla\cdot\mathbf{B}$
   error. What must **not** happen is domain-wide $\psi$ striping that grows
   in time — that is the undamped (`:lsource => false`) failure mode.

If the run needs more (or less) regularization, the knob is `:μ` — the
per-equation SGS multipliers. `:μ = 4` is the least dissipation measured to
be stable here and retains more magnetic-field structure; below that the run
aborts at the shocks. Prefer `:μ` over `:lfilter`: the Boyd–Vandeven filter
acts on the conservative variables independently and can itself destroy
pressure positivity (see
[README.md](README.md#stabilization-what-was-actually-tested)). `:nop` (lower
order = more robust, less accurate) and `:Δt` are the remaining knobs, but
the failure mode here is Gibbs overshoot at shocks, not a CFL violation — the
acoustic CFL is ≈ 0.24 and steady right up to the abort.

### 4.1 Verified end state

The shipped configuration was run to $t = 1$ and produces, at $t = 1$:

| field | range |
|---|---|
| $\rho$ | $[0.084,\ 0.364]$ (paper's Fig. 3 colorbar: $\approx[0.1, 0.4]$) |
| $p$ | $[0.047,\ 0.362]$ — positive everywhere |
| $\psi$ | $[-3.2\times10^{-3},\ 2.5\times10^{-3}]$ |
| $w,\ B_z$ | exactly $0$ |

No non-finite values in any of the 21 snapshots. The $w = B_z = 0$ result is
the free consistency check noted above, and it holds to the bit.

The computed density at the paper's comparison times is reproduced in
`assets/MHD_OT_rho.png` (repository root), on the same $[0.1, 0.4]$ color
scale as the "True" row of its Figure 3:

![Orszag-Tang density](../../../assets/MHD_OT_rho.png)

This is a **qualitative** morphology check — same large-scale features, same
amplitude band — not a quantitative validation. No pointwise or norm-based
comparison against the reference data has been performed, and the paper's
data set is not redistributed with Jexpresso.

## 5. Run

```bash
julia --project=. src/Jexpresso.jl MHD orszagTangBormanis2024
```

The mesh (doubly periodic $[0,1]^2$, $32\times32$ quads) ships with the case.
It can be regenerated with

```bash
gmsh -2 problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.geo \
     -o problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.msh
```

## References

1. A. Bormanis, C. A. Leon, A. Scheinker, *Solving the Orszag–Tang vortex
   magnetohydrodynamics problem with physics-constrained convolutional neural
   networks*, Phys. Plasmas 31, 012101 (2024). doi:10.1063/5.0172075
2. S. A. Orszag, C.-M. Tang, *Small-scale structure of two-dimensional
   magnetohydrodynamic turbulence*, J. Fluid Mech. 90 (1979) 129–143.
3. J. M. Picone, R. B. Dahlburg, *Evolution of the Orszag–Tang vortex system
   in a compressible medium. II. Supersonic flow*, Phys. Fluids B 3 (1991)
   29–44. (The paper's Ref. 16.)
4. P. Mocz, M. Vogelsberger, L. Hernquist, *A constrained transport scheme for
   MHD on unstructured static and moving meshes*, Mon. Not. R. Astron. Soc.
   442 (2014) 43–55. (The paper's Ref. 17 — the source of its constrained
   transport scheme and of the γ = 5/3 Orszag–Tang parameters.)
5. A. Dedner, F. Kemm, D. Kröner, C.-D. Munz, T. Schnitzer, M. Wesenberg,
   *Hyperbolic divergence cleaning for the MHD equations*,
   J. Comput. Phys. 175 (2002) 645–673.
6. D. Derigs, A. R. Winters, G. J. Gassner, S. Walch, M. Bohm,
   *Ideal GLM-MHD: About the entropy consistent nine-wave magnetic field
   divergence diminishing ideal magnetohydrodynamics equations*,
   J. Comput. Phys. 364 (2018) 420–467.
7. J. Chan, H. Ranocha, A. M. Rueda-Ramírez, G. J. Gassner, T. Warburton,
   Front. Phys. 10:898028 (2022). doi:10.3389/fphy.2022.898028 — the source of
   the GLM-MHD equation set implemented in Jexpresso.
8. J. Smagorinsky, *General circulation experiments with the primitive
   equations: I. The basic experiment*, Mon. Weather Rev. 91 (1963) 99–164.
