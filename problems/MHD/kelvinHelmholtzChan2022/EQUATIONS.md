# The 2D ideal GLM-MHD equations and the magnetized Kelvin-Helmholtz test

This document describes the system of equations solved by this case and the
test problem it is configured for. Both follow

> J. Chan, H. Ranocha, A. M. Rueda-Ramírez, G. J. Gassner, T. Warburton,
> *On the entropy projection and the robustness of high order entropy stable
> discontinuous Galerkin schemes for under-resolved flows*,
> Frontiers in Physics 10:898028 (2022), Section 3.2 (equations) and
> Section 3.2.1 (test). https://doi.org/10.3389/fphy.2022.898028

Implementation notes (file layout, variable ordering, how to run) are in
[README.md](README.md).

---

## 1. Governing equations

The ideal GLM-MHD (generalized Lagrange multiplier magnetohydrodynamics)
equations augment ideal MHD with a scalar field $\psi$ that propagates
divergence errors of the magnetic field out of the domain at a constant
hyperbolic cleaning speed $c_h$, so that the solution evolves toward
$\nabla\cdot\mathbf{B}=0$.

The system solved here is the conservative part of the paper's Eq. (6):

$$
\frac{\partial \mathbf{u}}{\partial t} + \nabla\cdot\mathbf{f}(\mathbf{u}) = \mathbf{0},
$$

with state vector

$$
\mathbf{u} = \left(\rho,\; \rho\mathbf{v},\; E,\; \mathbf{B},\; \psi\right)^T,
\qquad
\mathbf{v} = (u, v, w), \quad \mathbf{B} = (B_x, B_y, B_z),
$$

where $\rho$ is the density, $\rho\mathbf{v}$ the momentum, $E$ the total
energy density, $\mathbf{B}$ the magnetic field, and $\psi$ the
divergence-correcting field — nine unknowns in total.

**Two-dimensionality.** Following Section 3.2 of the paper, the 2D system is
obtained by dropping the flux in $z$ while **keeping the third components of
velocity and magnetic field** ($w$, $B_z$), because plasmas admit
three-dimensional electromagnetic interactions even in two-dimensional
problems.

### 1.1 Fluxes

With the total (gas + magnetic) pressure

$$
p_{tot} = p + \tfrac{1}{2}\lVert\mathbf{B}\rVert^2,
$$

the Cartesian fluxes are

$$
\mathbf{F}_1 =
\begin{pmatrix}
\rho u \\
\rho u^2 + p_{tot} - B_x^2 \\
\rho u v - B_x B_y \\
u\left(\tfrac{1}{2}\rho\lVert\mathbf{v}\rVert^2 + \tfrac{\gamma p}{\gamma-1} + \lVert\mathbf{B}\rVert^2\right) + B_x\left(c_h\psi - \mathbf{v}\cdot\mathbf{B}\right) \\
\rho u w - B_x B_z \\
c_h \psi \\
u B_y - v B_x \\
u B_z - w B_x \\
c_h B_x
\end{pmatrix},
\qquad
\mathbf{F}_2 =
\begin{pmatrix}
\rho v \\
\rho v u - B_y B_x \\
\rho v^2 + p_{tot} - B_y^2 \\
v\left(\tfrac{1}{2}\rho\lVert\mathbf{v}\rVert^2 + \tfrac{\gamma p}{\gamma-1} + \lVert\mathbf{B}\rVert^2\right) + B_y\left(c_h\psi - \mathbf{v}\cdot\mathbf{B}\right) \\
\rho v w - B_y B_z \\
v B_x - u B_y \\
c_h \psi \\
v B_z - w B_y \\
c_h B_y
\end{pmatrix},
$$

written here in the row order used by the paper's state
$(\rho, \rho\mathbf{v}, E, \mathbf{B}, \psi)$ with
$\rho\mathbf{v}=(\rho u,\rho v,\rho w)$. (In the code the total energy is
stored in slot 4 and $\rho w$ in slot 5 — see README.md for why.)

The blocks are the familiar ideal-MHD ones — mass advection; momentum with
total pressure and Maxwell stress $-\mathbf{B}\mathbf{B}^T$; total energy
with pressure work, Poynting-like transport $-\mathbf{B}(\mathbf{v}\cdot\mathbf{B})$
and the GLM contribution $c_h\psi\,\mathbf{B}$; the induction equation
$\mathbf{v}\mathbf{B}^T - \mathbf{B}\mathbf{v}^T$ — plus the two GLM cleaning
couplings $c_h\psi\,\mathbf{I}$ in the induction rows and $c_h\mathbf{B}$ in
the $\psi$ row. Together the last two make $(\nabla\cdot\mathbf{B}, \psi)$
satisfy a linear wave system with speed $c_h$: divergence errors are radiated
away instead of accumulating.

### 1.2 Equation of state

The gas pressure closes the system:

$$
p = (\gamma - 1)\left(E - \tfrac{1}{2}\rho\lVert\mathbf{v}\rVert^2
      - \tfrac{1}{2}\lVert\mathbf{B}\rVert^2 - \tfrac{1}{2}\psi^2\right).
$$

The heat-capacity ratio is set in `user_flux.jl` (`γ_mhd`); the paper does
not state a numerical value, and $\gamma = 5/3$ (ideal monatomic plasma, the
standard choice of the GLM-MHD literature the paper builds on) is used here.
Set it to $1.4$ to match the companion Euler case
`problems/CompEuler/kelvinHelmholtzChan2022`.

### 1.3 Divergence-cleaning speed $c_h$

Following Section 3.2.1 of the paper, $c_h$ is set to the **maximum wave
speed of the initial condition and kept constant** for the whole run:

$$
c_h = \max_{\Omega}\left(\lVert\mathbf{v}\rVert + c_f^{\max}\right),
\qquad
c_f^{\max} = \sqrt{\frac{\gamma p}{\rho} + \frac{\lVert\mathbf{B}\rVert^2}{\rho}},
$$

where $c_f^{\max}$ bounds the fast magnetosonic speed over all propagation
directions. This is computed over the mesh in `initialize.jl` (MPI-reduced)
and stored in the `c_h_mhd` reference read by `user_flux!`. For the initial
condition below and $\gamma=5/3$ it evaluates to $c_h \approx 2.344$. The
paper notes that smaller $c_h$ degrades robustness while larger $c_h$
stiffens the system.

**ψ damping (mixed cleaning).** Purely hyperbolic cleaning only *transports*
divergence errors; nothing destroys $\psi$. In the paper's DG setting the
interface Riemann fluxes upwind-damp the $\psi$ waves, but a continuous
Galerkin discretization has no interface dissipation, so on a periodic
domain the $\psi$ waves keep bouncing and alias into grid-scale (LGL
odd–even) noise. We therefore use Dedner's *mixed*
(hyperbolic–parabolic) GLM and damp $\psi$ with the source

$$
S_\psi = -\frac{c_h^2}{c_p^2}\,\psi = -\frac{c_h}{c_r}\,\psi,
\qquad c_p^2 = c_h c_r,\quad c_r = 0.18 \ \text{(Dedner's value)},
$$

implemented in `user_source.jl` (`glm_cr_mhd` sets $c_r$). Only $\psi$ is
damped: since $E$ carries $\tfrac{1}{2}\psi^2$, removing $\psi$ while
leaving $E$ untouched converts the cleaned $\psi$-energy into heat, which is
the entropy-consistent behavior (Derigs et al. 2018).

### 1.4 The non-conservative term (NOT included)

The full non-conservative GLM-MHD system of the paper (Eq. 6–7) adds a term

$$
\Upsilon = (\nabla\cdot\mathbf{B})
\left(0,\; \mathbf{B},\; \mathbf{v}\cdot\mathbf{B},\; \mathbf{v},\; 0\right)
+ \left(0,\; \mathbf{0},\; \psi\,(\mathbf{v}\cdot\nabla\psi),\; \mathbf{0},\; \mathbf{v}\cdot\nabla\psi\right),
$$

i.e. the Powell term plus a Galilean GLM correction. It vanishes for exactly
divergence-free fields and is required to achieve **entropy stability** and
Galilean invariance of the split-form/two-point-flux discretizations studied
in the paper. Since this case deliberately uses neither the entropy-stable
nor the kinetic-energy-preserving machinery (see below), $\Upsilon$ is
**omitted for now**; the conservative GLM system above is the standard
Dedner-type formulation and is well-posed on its own.

## 2. Discretization and stabilization

- Continuous spectral elements, LGL nodes, polynomial degree
  `:nop => 7` on a $32\times 32$ element grid (the paper's $N=7$
  configuration of Figure 4), inexact (collocated) quadrature.
- Explicit low-storage RK: `CarpenterKennedy2N54()` with
  $\Delta t = 2\times 10^{-3}$.
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
  `:μ = [0, 1, 1, 1, 1, 1, 1, 1, 1]` (no mass diffusion):

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

## 3. Test: 2D magnetized Kelvin-Helmholtz instability (Sec. 3.2.1)

A magnetized modification of the Euler KHI of the paper's Section 3.1.1:
a double shear layer with a sinusoidal transverse perturbation, threaded by
a uniform **vertical** magnetic field.

**Domain.** $\Omega = [-1, 1]^2$, doubly periodic. Final time
$T_{final} = 15$.

**Initial condition** (paper Eq. 10, with the smoothed step of Eq. 5):

$$
\tilde{B}(x,y) = \tanh(15y + 7.5) - \tanh(15y - 7.5),
$$

$$
\begin{aligned}
\rho &= \tfrac{1}{2} + \tfrac{3}{4}\tilde{B}, &\quad p &= 1, &\quad \psi &= 0,\\
u &= \tfrac{1}{2}\left(\tilde{B} - 1\right), &\quad v &= \tfrac{1}{10}\sin(2\pi x), &\quad w &= 0,\\
B_x &= 0, &\quad B_y &= 0.125, &\quad B_z &= 0,
\end{aligned}
$$

and $E = \frac{p}{\gamma-1} + \frac{1}{2}\rho\lVert\mathbf{v}\rVert^2 +
\frac{1}{2}\lVert\mathbf{B}\rVert^2 + \frac{1}{2}\psi^2$.

$\tilde{B}\in(0,2)$ switches between the outer stream
($\rho = 1/2$, $u = -1/2$) and the inner stream ($\rho = 2$, $u = +1/2$)
across thin $\tanh$ layers at $y = \pm 1/2$; the $v$-perturbation seeds the
roll-up.

**What to expect.** The shear rolls up into KH vortices; the vertical field
is stretched by the shear into a strong horizontal component (in the runs
used to validate this case, $B_x$ grows from $0$ to $\approx \pm 0.4$ by
$t = 0.5$), and magnetic tension progressively suppresses the small-scale
vortical structures while extending the flow features in $y$ — compare the
paper's Figure 4 (MHD) against its Figure 1 (Euler). MHD turbulence develops
after $t = 10$; the paper uses this regime as a robustness stress test
(Table 3). Since this setup replaces the paper's entropy-stable Gauss DG
with Smagorinsky-regularized collocation, results are expected to be more
dissipative than the reference and are not meant to reproduce it to the
digit.

**Diagnostics written** (VTK, every $0.5$ time units by default):
$\rho, u, v, w, p, B_x, B_y, B_z, \psi, T$.

## 4. Run

```bash
julia --project=. src/Jexpresso.jl MHD kelvinHelmholtzChan2022
```

The mesh must be a doubly periodic $[-1,1]^2$ quad grid; it can be
regenerated with

```bash
gmsh -2 problems/MHD/kelvinHelmholtzChan2022/KHI_32x32_periodic.geo \
     -o meshes/gmsh_grids/hexa_TFI_32x32_unitsquare.msh
```

## References

1. J. Chan, H. Ranocha, A. M. Rueda-Ramírez, G. J. Gassner, T. Warburton,
   Front. Phys. 10:898028 (2022). doi:10.3389/fphy.2022.898028
2. A. Dedner, F. Kemm, D. Kröner, C.-D. Munz, T. Schnitzer, M. Wesenberg,
   *Hyperbolic divergence cleaning for the MHD equations*,
   J. Comput. Phys. 175 (2002) 645–673.
3. D. Derigs, A. R. Winters, G. J. Gassner, S. Walch, M. Bohm,
   *Ideal GLM-MHD: About the entropy consistent nine-wave magnetic field
   divergence diminishing ideal magnetohydrodynamics equations*,
   J. Comput. Phys. 364 (2018) 420–467.
4. A. M. Rueda-Ramírez, S. Hennemann, F. J. Hindenlang, A. R. Winters,
   G. J. Gassner, *An entropy stable nodal discontinuous Galerkin method for
   the resistive MHD equations. Part II: Subcell finite volume shock
   capturing*, J. Comput. Phys. 444 (2021) 110580.
5. J. Smagorinsky, *General circulation experiments with the primitive
   equations: I. The basic experiment*, Mon. Weather Rev. 91 (1963) 99–164.
