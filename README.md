# <img src="./assets/logo-ext2.png" width="500" title="JEXPRESSO logo">

| **Documentation** |
|:------------ |
 [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://smarras79.github.io/Jexpresso/dev/) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://smarras79.github.io/Jexpresso/dev/) |
|**Build Status** |
| [![CI](https://github.com/smarras79/Jexpresso/actions/workflows/CI.yml/badge.svg?branch=master&event=push)](https://github.com/smarras79/Jexpresso/actions/workflows/CI.yml?query=branch%3Amaster) [![Documentation](https://github.com/smarras79/Jexpresso/actions/workflows/Documentation.yml/badge.svg?branch=master&event=push)](https://github.com/smarras79/Jexpresso/actions/workflows/Documentation.yml?query=branch%3Amaster)
| **Contacts**  |
| [![Simone Marras](https://img.shields.io/badge/Simone%20Marras-smarras%40njit.edu-8e7cc3)](mailto:smarras@njit.edu) |
| [![Yassine Tissaoui](https://img.shields.io/badge/Yassine%20Tissaoui-tissaoui%40wisc.edu-8e7cc3)](mailto:tissaoui@wisc.edu) |
| [![Hang Wang](https://img.shields.io/badge/Hang%20Wang-hang.wang%40njit.edu-8e7cc3)](mailto:hang.wang@njit.edu) |
| **Citation** |
| [![DOI](https://img.shields.io/badge/article-arXiv:2401.05624-green)](https://doi.org/10.48550/arXiv.2401.05624) |

# JEXPRESSO
A CPU and GPU research software for the numerical solution of a system of arbitrary conservation laws using **continuous spectral elements** and finite differences in **1D, 2D, 3D**. DISCLAIMER: this will always be WIP! Contact us to join the team of developers!

Suggested Julia version: 1.11.9

# A note about the use of AI
Jexpresso has been developed by humans since 2021 and continues to be so. Since Spring 2026, AI has been assisting the developers for new problems additions, debugging, and code's documentation. As AI becomes more reliable, we foresee an increased use of it for code development under the direct supervision of a human expert. 
The Jexpresso core team uses Claude whereas some external developers have been successfully using OpenAI's Codex for their own implementations.

# Table of Contents

- [Installation](#installation)
- [Equations](#equations)
  1. [1D wave equation](#1-1d-wave-equation)
  2. [1D shallow water](#2-1d-shallow-water)
  3. [2D Helmholtz](#3-2d-helmholtz)
  4. [2D scalar advection-diffusion](#4-2d-scalar-advection-diffusion)
  5. [2D Euler equations of compressible flows with gravity and passive chemicals](#5-2d-euler-equations-of-compressible-flows-with-gravity-and-passive-chemicals)
  6. [3D Euler equations of compressible flows with gravity](#6-3d-euler-equations-of-compressible-flows-with-gravity)
- [Showcase](#turbulent-abl)
  - [Planet Saturn](#planet-saturn)
  - [Turbulent ABL](#turbulent-abl)
  - [Shallow cumuli](#shallow-cumuli)
- [Examples available in this branch](#examples-available-in-this-branch)
  - [1D shock tube with dynamic SGS (DynSGS) for shock capturing](#1d-shock-tube-with-dynamic-sgs-dynsgs-for-shock-capturing)
  - [1D acoustic wave](#1d-acoustic-wave)
  - [Flow at Mach 3 with forward-facing step](#flow-at-mach-3-with-forward-facing-step)
  - [Kelvin-Helmholtz instability](#kelvin-helmholtz-instability)
  - [Solid elasticicy](#Solid-elasticity)
  - [MHD: magnetized Kelvin-Helmholtz instability](#magneto-hydrodynamics-mhd-magnetized-kelvin-helmholtz-instability)
  - [MHD: Orszag-Tang vortex](#magneto-hydrodynamics-mhd-orszag-tang-vortex)
  - [Cloud simulation: shallow cumuli with BOMEX conditions](#cloud-simulation-shallow-cumuli-with-bomex-conditions)
  - [Shallow water on a spherical shell](#shallow-water-on-a-spherical-shell)
  - [2D Euler equations with buoyancy and two passive tracers](#2d-euler-equations-with-buoyancy-and-two-passive-tracers)
  - [3D Euler equations with buoyancy](#3d-euler-equations-with-buoyancy)
  - [Spectral convergence of the SEM: doubly periodic Poisson problem](#spectral-convergence-of-the-sem-doubly-periodic-poisson-problem) (SEM direct/AMG, static condensation, pseudo-spectral, FFT)
  - [Laguerre semi-infinite element test suite](#laguerre-semi-infinite-element-test-suite)
    - [Test 1: 1D wave equation with Laguerre absorbing layers](#test-1-1d-wave-equation-with-laguerre-semi-infinite-element-absorbing-layers)
    - [Test 2: 1D wave train for linearized shallow water equations](#test-2-1d-wave-train-for-linearized-shallow-water-equations)
    - [Test 3: 2D advection-diffusion equation](#test-3-2d-advection-diffusion-equation)
    - [Test 4: 2D Helmholtz equation](#test-4-2d-helmholtz-equation)
  - [Rising thermal bubble with semi-infinite Laguerre elements for outflows](#rising-thermal-bubble-with-semi-infinite-laguerre-elements-for-outflows)
  - [Hydrostatic linear mountain waves with semi-infinite Laguerre elements for outflows](#hydrostatic-linear-mountain-waves-with-semi-infinite-laguerre-elements-for-outflows)
  - [Non-hydrostatic mountain waves: comparison against WRF](#non-hydrostatic-mountain-waves-comparison-against-wrf)

# Installation:
Follow the instructins in [INSTALL.md](INSTALL.md)

Run into trouble? Check the [FAQ.md](FAQ.md) for common installation and run errors.

If you use Jexpresso please drop us a line to let us know. We'd like to add a link to your paper or work on this page.

Please cite Jexpresso using:

```
@article{tissaoui2024,
  author = {Y. Tissaoui and J. F. Kelly and S. Marras}
  title = {Efficient Spectral Element Method for the Euler Equations on Unbounded Domains},
  volume ={487},
  pages={129080},
  year = {2024},
  journal = {App. Math. Comput.},
}

@inproceedings{marrasJexpresso,
  author    = {S. Marras and Y. Tissaoui and H. Wang and S. Stechmann}
  title     = {Jexpresso V0.1.0: a Julia-language, user-friendly, multi-physics parallel solver for the solution of conservation laws on CPUs and GPUs},
  booktitle = {Proceedings of the 36th Parallel CFD international conference 2025},
  year      = {2025},
  address   = {Merida, Yucatan, Mexico},
  month     = {November},
  organization = {UNAM},
}
```

# Equations:
Jexpresso uses arbitrarily high-order (3rd and above) **continuous spectral elements** to solve

$$\delta\frac{\partial \bf q}{\partial t} + \sum_{i=1}^{nd}\nabla\cdot{{\bf F}_i({\bf q})} = \mu\nabla^2{\bf q} + {\bf S}({\bf q}) + ~{\rm b.c.}$$

where $\delta = 0,1$ simply indicates time-independent equations the vectors ${\bf q}$, ${\bf F}$, and ${\bf S}$ are problem-dependent as shown below,
and are taken to be zero vectors of the appropriate size when not explicitly stated otherwise.

The Julia package [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) is used for time discretization and stepping.

In order, we provide tests and results for the following equations:

### 1. 1D wave equation

$${\bf q}=\begin{bmatrix}
u \\
v
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
v\\
u
\end{bmatrix}$$

### 2. 1D shallow water

$${\bf q}=\begin{bmatrix}
h \\
u
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
Uh + Hu\\
gh + Uu
\end{bmatrix},$$

where $H$ and $U$ are a reference height and velocity, respectively.

### 3. 2D Helmholtz

$${\bf S}=\begin{bmatrix}
\alpha^2 u + f(x,z)
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
u_{xx} + u_{zz}
\end{bmatrix},$$

for a constant value of $\alpha$ and $\mu$, which are case-dependent.

### 4. 2D scalar advection-diffusion

$${\bf q}=\begin{bmatrix}
q\\
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
qu\\
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
qv\\
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
q_{xx} + q_{zz}
\end{bmatrix},$$

### 5. 2D Euler equations of compressible flows with gravity and passive chemicals

With $N$ passive chemicals $c_i, \forall i=1,...,N$:

$${\bf q}=\begin{bmatrix}
\rho \\
\rho u\\
\rho v\\
\rho \theta\\
\rho c1\\
...\\
\rho cN
\end{bmatrix}\quad {\bf F1}=\begin{bmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
\rho u \theta\\
\rho u c1\\
...\\
\rho u cN
\end{bmatrix}\quad {\bf F2}=\begin{bmatrix}
\rho v\\
\rho v u\\
\rho v^2 + p\\
\rho v \theta\\
\rho v c1\\
...\\
\rho v cN
\end{bmatrix}\quad {\bf S}=\begin{bmatrix}
0\\
0\\
-\rho g\\
0\\
0\\
...\\
0
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
0\\
u_{xx} + u_{zz}\\
v_{xx} + v_{zz}\\
\theta_{xx} + \theta_{zz}\\
c1_{xx} + c1_{zz}\\
...\\
cN_{xx} + cN_{zz}
\end{bmatrix}.$$

### 6. 3D Euler equations of compressible flows with gravity

$${\bf q}=\begin{bmatrix}
\rho \\
\rho u\\
\rho v\\
\rho w\\
\rho \theta\\
\end{bmatrix}\quad {\bf F}1=\begin{bmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
\rho u w\\
\rho u \theta\\
\end{bmatrix}\quad {\bf F}2=\begin{bmatrix}
\rho v\\
\rho v u\\
\rho v^2 + p\\
\rho v w\\
\rho v \theta\\
\end{bmatrix}\quad {\bf F3}=\begin{bmatrix}
\rho w\\
\rho w u\\
\rho w v\\
\rho w^2 + p\\
\rho w \theta\\
\end{bmatrix}\quad {\bf S}=\begin{bmatrix}
0\\
0\\
0\\
-\rho g\\
0\\
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
0\\
u_{xx} + u_{yy} + u_{zz}\\
v_{xx} + v_{yy} + v_{zz}\\
w_{xx} + w_{yy} + w_{zz}\\
\theta_{xx} + \theta_{yy} + \theta_{zz}\\
\end{bmatrix}.$$


If you are interested in contributing, please get in touch:
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:tissaoui@wisc.edu), [Hang Wang](mailto:hang.wang@njit.edu)

## Planet Saturn:
Example of a relatively coarse simulation of the atmosphere of planet Saturn during 500 days.

Original test described by Scott and Polvani "Forced-Dissipative Shallow-Water Turbulence on the Sphere and the Atmospheric Circulation of the Giant Planets" J. Atmos. Sci. Vol 64, 2007

<img src="assets/vorticity-grey-scale.gif"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />

## Turbulent ABL:
Example of coarse simulation of the turbulent atmospheric boundary layer. Domain size: 10240m X 10240m X 3000m using 64x64x24 spectral elements of order 4.
Surface and SGS: Monin-Obukhov Similarity Theory model with Richardson-corrected Smagorinsky.
<img src="assets/ABLfullDomain.gif"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />


## Shallow cumuli:
Example of shallow cumuli simulations (right) for the type of Barbados clouds shown on the left: (picture taken from [P. Blossey webpage](https://www.atmos.washington.edu/~bloss/) from U. Washington)

<img src="assets/barbados.jpg"
     alt="Markdown icon"
     style="float: left; margin-right: 3.5px;" />

# Examples available in this branch:
Below are just a few pre-packaged examples available in Jexpresso.
To add your own new problem, see [ADD_A_NEW_TEST.md](ADD_A_NEW_TEST.md).


## 1D shock tube with dynamic SGS (DynSGS) for shock capturing:
Classical Sod's tube with shock and expansion.
The DynSGS SGS model by Marras et al. 2015 and later is used to capture the shock.

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "sod1d")
```

<img src="assets/sod1d.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

## 1D acoustic wave:
```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "case1")
```

<img src="assets/1dacoustic.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

## Flow at Mach 3 with forward-facing step
Classical flow at Mach 3 with DynSGS shock capturing
```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "ffs_step")
```

<img src="assets/shock-MrhoSchielern.jpg"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

## Flow at Mach 3 with airfoil
Mach 3 flow over a (non-supersonic) airfoil with exact geometry (i.e. curved elements). DynSGS shock capturing.
```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "naca64A210")
```

<img src="assets/NACA64A210-TWOFIGS.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />


## Kelvin-Helmholtz instability
Classical shear-triggered instability test.

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "kelvinHelmholtzChan2022")
```

<img src="assets/KH-EC-SGSsmag.jpg"
     alt="Markdown icon"
     style="float: left; margin-right: 3.5px;" />


## Solid elasticity
Timoshenko's model of elasticity

```julia
using Jexpresso
Jexpresso.run_case("Elasticity", "beam2d")
```

https://github.com/user-attachments/assets/78872c85-e8b5-494f-b95a-4a94aeb5f07f

<img src="assets/beam2d.jpeg"
     alt="Markdown icon"
     style="float: left; margin-right: 3.5px;" />



## Magneto-Hydrodynamics (MHD), magnetized Kelvin-Helmholtz instability:

The problem is defined in [`problems/equations/MHD/kelvinHelmholtzChan2022`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/MHD/kelvinHelmholtzChan2022).

```julia
using Jexpresso
Jexpresso.run_case("MHD", "kelvinHelmholtzChan2022")
```

<img src="assets/MHD_By.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

## Magneto-Hydrodynamics (MHD), Orszag-Tang vortex:

The classical 2D MHD benchmark, with the setup of Bormanis, Leon & Scheinker,
*Phys. Plasmas* **31**, 012101 (2024): doubly periodic unit square, γ = 5/3,
128×128 points, t ∈ [0, 1]. The problem is defined in
[`problems/equations/MHD/orszagTangBormanis2024`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/MHD/orszagTangBormanis2024)
(see its `README.md` and `EQUATIONS.md`).

```julia
using Jexpresso
Jexpresso.run_case("MHD", "orszagTangBormanis2024")
```

Stabilized with **DynSGS** — the residual-based, parameter-free dynamic SGS
model of Marras, Nazarov & Giraldo (see [DSGS.md](DSGS.md)).

<img src="assets/MHD-OT-4plots.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Top row: density at t = 0.5 s (left) and t = 1.0 s (right).
Bottom row: residual viscosity. Simulation using 120 × 120 4th-order spectral elements in
a unit square.

The solver writes VTK, not PNG. These three figures are rendered from a
finished run by `julia --project=. tools/plot_orszag_tang.jl`; see the
[case README](problems/MHD/orszagTangBormanis2024/README.md#figures).

## Magneto-Hydrodynamics (MHD), flux emergence in the solar atmosphere:

The two-dimensional emergence of a horizontal magnetic flux sheet through a
two-temperature (chromosphere + corona) stratified atmosphere — the nonlinear
Parker instability of Shibata et al. (1989) — with the setup of Son, Jang &
Magara, *ApJS* **277**:46 (2025): γ = 1.05, [0, 80 H₀] × [0, 35 H₀],
t ∈ [0, 54 τ₀]. The problem is defined in
[`problems/MHD/fluxEmergenceSon2025`](problems/MHD/fluxEmergenceSon2025)
(see its `README.md` and `EQUATIONS.md`).

```bash
mpiexec -n 10 julia --project=. src/Jexpresso.jl MHD fluxEmergenceSon2025
```

Stabilized with **DynSGS**, integrated with Carpenter–Kennedy 2N54 on the
coarsest N = 4 grid the problem admits (80×35 elements). The solver writes PNGs
styled after the paper's figures (log₁₀ density on the paper's `jet` scale
with magnetic field lines and velocity vectors; centerline profiles of the
rise velocity, Alfvén speed, field and density on the axes of its Fig. 5)
directly, gathered on one rank under MPI.

[`problems/MHD/brioWu1d`](problems/MHD/brioWu1d) is the 1D Brio–Wu MHD
shock tube (Dao & Nazarov 2022, §5.2), the MHD counterpart of `CompEuler/sod1d`:
DynSGS in its conserved form with the 1D MHD kernel, 600 LGL points, the
reference solution overlaid at t = 0.2.

```bash
julia --project=. src/Jexpresso.jl MHD brioWu1d
```

[`problems/ShallowWater/SoliWaveIslandDSGS`](problems/ShallowWater/SoliWaveIslandDSGS)
is the solitary wave on a conical island of Marras et al. (2018, §5.5)
stabilized by DynSGS for the shallow-water system (`DSGS_SW()`) instead of the
constant viscosity of `ShallowWater/SoliWaveIsland`.

```bash
julia --project=. src/Jexpresso.jl ShallowWater SoliWaveIslandDSGS
```

[`problems/MHD/fluxEmergenceSon2025DSGS`](problems/MHD/fluxEmergenceSon2025DSGS)
is the same problem with **DynSGS alone** keeping the solution admissible:
no positivity limiter, the dissipation acting on the relative departure from
the magnetostatic reference state (`:dsgs_ref_weight`, DSGS.md §4.5).

```bash
mpiexec -n 10 julia --project=. src/Jexpresso.jl MHD fluxEmergenceSon2025DSGS
```

## Cloud simulation: shallow cumuli with BOMEX conditions:

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "3d_bomex")
```
<img src="assets/bomex.png"
     alt="Markdown icon"
     style="float: left; margin-right: 3.5px;" />

## Shallow water on a spherical shell
Benchmark: classical Galewki and Polvani's barotropic jet
```julia
using Jexpresso
Jexpresso.run_case("ShallowWater", "SWsphere")
```
<img src="assets/SWsphere-Galewki-visc1e5-36x36.jpg"
     alt="Markdown icon"
     style="float: left; margin-right: 3.5px;" />


This case also ships with **Proper Orthogonal Decomposition** switched on: at
the end of the run the code extracts the energy-ranked modes of the flow, draws
them on an equirectangular map together with the energy spectrum and the
temporal coefficients, and writes the basis out for a reduced-order model.

POD is a property of the framework rather than of this case: **any** problem
turns it on with `:lpod => true` in its deck and supplies nothing else, in 1D,
2D, 3D or on a manifold. `problems/AdvDiff/PODbenchmark` is the reference
benchmark, a problem whose POD is known in closed form. See
[`docs/POD.md`](docs/POD.md).



## 2D Euler equations with buoyancy and two passive tracers
The problem is defined in `problems/equations/CompEuler/thetaTracers`. To run it you would do the following:
```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "thetaTracers")
```

<img src="assets/thetaTracersMeshUnstr.png"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />


## 3D Euler equations with buoyancy
The problem is defined in `problems/equations/CompEuler/3d`. To run it you would do the following:
```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "3d")
```

<img src="assets/rtb3d.png"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />

## Spectral convergence of the SEM: doubly periodic Poisson problem
The problem is defined in `problems/Elliptic/poisson_periodic_sem`: $-\nabla^2 u = f$ on $[0,2\pi]^2$, periodic in $x$ and $y$, with the exact solution $u = A\,(p(x)\,p(y) - 1/(c^2-1))$, $p(s) = 1/(c-\cos s)$, the product of two periodic Poisson kernels ($c = (r+1/r)/2$ with $r = 0.8$, and $A$ scales the peak to $u(0,0)=1$; see `user_source.jl`). Its Fourier coefficients decay geometrically, like $r^{|k_x|+|k_y|}$, so $u$ is not band-limited: the Fourier solvers have a grid-dependent error too, not just round-off. The same deck solves it six ways:

| solver | deck flags | discretisation | solve |
|---|---|---|---|
| SEM direct | (default) | 16×16 spectral elements of order N (`:nop`) | sparse direct (Cholesky) on the full periodic system |
| SEM AMG | `:linsolve_amg => true` | same | AMG-preconditioned conjugate gradients on the full system |
| SC direct | `:lstatic_condensation => true` | same, statically condensed (below) | sparse Cholesky on the (symmetrised) skeleton system |
| SC AMG | `:lstatic_condensation => true, :EL_skeleton_solver => "amg"` | same, statically condensed | AMG-preconditioned conjugate gradients on the skeleton system |
| pseudo-spectral | `:lpseudospectral => true` | Fourier collocation on a uniform `:fft_N`² grid, Kopriva's derivative matrix (`FourierDerivativeMatrix`) | dense matrix diagonalisation, O(N³) |
| FFT | `:lfft => true` | Fourier spectral on a uniform `:fft_N`² grid | FFTW, O(N² log N) |

**Static condensation (SC)** is the algorithm of element learning (`elementLearning_Axb!`), used here as a solver: the interior unknowns of every element are eliminated with the local operators T^ie = (A_{vo,vo})⁻¹ A_{vo,vb}, computed from the element blocks of the SEM matrix (element learning replaces exactly these with a trained network), which leaves a Schur-complement system on the element skeleton only — 3 840 of the 16 384 unknowns at N = 8. The skeleton system is solved, and the interiors are recovered element by element. Nothing is approximated: SC reproduces the full SEM solution to round-off. On this periodic problem there is no Dirichlet boundary (Γ = ∅), so the skeleton system is singular like the full one; one unknown is pinned and the result is shifted to zero mean, as in every periodic solve. **AMG** is smoothed aggregation (AlgebraicMultigrid.jl) as the preconditioner of conjugate gradients (Krylov.jl), to a relative residual of 10⁻¹²; `:amg_method => "rs"` selects Ruge–Stüben. The same `:EL_skeleton_solver` option applies to the element-learning inference and to the Dirichlet decks (see `problems/Elliptic/poisson_dirichlet_sc`).

The skeleton matrix B is symmetric in exact arithmetic, but the subtraction that forms it leaves a round-off asymmetry (≈1e-16 relative). That was enough for Julia's `factorize`, which tests for exact symmetry, to fall back to UMFPACK LU. `el_skeleton_solve` now symmetrises B as (B + Bᵀ)/2 and factorises it with CHOLMOD Cholesky, like the full system. The table and figures below use Cholesky.

To run it you would do the following:
```julia
using Jexpresso
Jexpresso.run_case("Elliptic", "poisson_periodic_sem")
```

### Benchmark: error and time-to-solution of the six solvers
`tools/periodic_poisson_benchmark/pipeline.jl` runs all six solvers at the same number of unknowns, the four SEM solves at orders N = 2…8 and the two Fourier solvers on 16N × 16N grids, and writes the table below and the figures. Run it from the REPL:
```julia
julia --project=.
julia> using Jexpresso
julia> include("tools/periodic_poisson_benchmark/pipeline.jl")
julia> rows = run_periodic_poisson_benchmark()
```
For higher resolutions, `levels = 0:3` refines the 16×16 mesh uniformly through `:linitial_refine`, up to (16·2^L)² elements. The pipeline can cap the size per solver, resume an interrupted sweep, and draw h-refinement figures of time against unknowns. See [tools/periodic_poisson_benchmark/README.md](tools/periodic_poisson_benchmark/README.md) for running it on a laptop.

**Timing protocol.** Every configuration runs twice in the same Julia session and only the **second run** is recorded, so compilation never enters a number. Each time is a single wall-clock measurement of that run (no repetition, no minimum over samples). The mesh and SEM preprocess caches are switched off, so the SEM infrastructure is built, not loaded from disk, and output files are switched off.

**What the columns mean.**
- **solve**: the solve step alone (triangular solves for SEM direct; the CG iterations for SEM AMG; the skeleton solve plus the interior recovery for SC; four dense N×N products for pseudo-spectral; rfft/scale/brfft for the FFT).
- **setup**: the solver's own infrastructure (periodic reduction, then the Cholesky factorisation or the AMG hierarchy; for SC also the element blocks and the Schur complement; 1-D eigen-decompositions; FFTW plan).
- **solved for / CG its**: the unknowns of the system actually solved (the skeleton for SC) and the conjugate-gradient iterations of the AMG solves.
- **SEM infrastructure**: the mesh read and the SEM setup (basis, metrics, mass and Laplacian assembly). Only the SEM needs it.
- **time-to-solution**: everything the method needs: SEM infrastructure + RHS + setup + solve for the four SEM solves; RHS + setup + solve for the Fourier solvers.
- **run_case wall-clock**: the whole second `run_case` call. The Fourier solvers read no mesh: the driver dispatches to them before any SEM infrastructure is built. The FFT is planned with `FFTW.ESTIMATE` (`:fft_plan => "estimate"`), so its plan is built, and charged to its setup, in every run.

**Reading the results.** The four SEM solves compute the same discrete solution, so their error curves coincide (only the last one drawn is visible). What separates them is cost. AMG on the condensed skeleton system needs far fewer iterations than on the full system, and the gap grows with the order (41 against 138 CG iterations at N = 8): the condensation removes the element-interior modes that make the high-order SEM system hard for AMG, and the skeleton system is also 4× smaller. At N = 8 the SC AMG solve step is 29× faster than the full-system AMG one. The sparse direct solves stay the fastest solve steps at these sizes. For SEM direct, SC direct and SC AMG the time-to-solution is dominated by the SEM infrastructure (1.3–2.4 s, mostly reading and building the mesh); full-system AMG is the exception at high order (3.5 s of solver cost at N = 8).

The pseudo-spectral and FFT errors fall geometrically with the grid, from 5.6e-1 on the 32² grid to 5.0e-7 on the 128² grid, as the r^(N_g/2) decay of the Fourier coefficients predicts. The two compute the same discrete solution, and their errors agree to 6 digits. The SEM error falls exponentially with the order, from 7.3e-1 at N = 2 to 5.0e-5 at N = 8. At the same number of unknowns the Fourier error is about 100× smaller at N = 8: the solution is analytic and periodic, which suits a global Fourier basis, while the SEM gains its accuracy element by element. Likewise the four SEM solvers (direct, AMG, SC direct, SC AMG) compute the same discrete SEM solution, so their errors agree to 3–4 digits (see the table below); the two error plots draw them as a single SEM curve, and the time plots separate them.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_error_vs_order-dark.svg">
  <img src="assets/ppb_error_vs_order.svg" width="680" alt="L-infinity error versus SEM order N: the SEM (one curve for its four solvers) falls from 7.3e-1 at N = 2 to 5.0e-5 at N = 8; the pseudo-spectral and FFT solvers at the same number of unknowns fall from 5.6e-1 to 5.0e-7.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_error_vs_dofs-dark.svg">
  <img src="assets/ppb_error_vs_dofs.svg" width="680" alt="L-infinity error versus number of unknowns, 1024 to 16384, for the SEM (one curve for its four solvers), pseudo-spectral and FFT solvers.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_error_vs_solve_time-dark.svg">
  <img src="assets/ppb_error_vs_solve_time.svg" width="680" alt="L-infinity error versus the wall-clock of the solve step alone, second run, for the SEM, pseudo-spectral and FFT solvers.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_error_vs_total_time-dark.svg">
  <img src="assets/ppb_error_vs_total_time.svg" width="680" alt="L-infinity error versus time-to-solution including all the infrastructure each method needs, second run: the SEM solvers between 1.3 and 6 s, dominated by the mesh read and SEM setup except full-system AMG at high order; the pseudo-spectral solver from 2 to 31 ms; the FFT from 0.2 to 2 ms.">
</picture>

| method | elements | SEM order N | unknowns (grid) | solved for | CG its | ‖e‖∞ | relative ‖e‖₂ | solve | setup | RHS | SEM infrastructure | time-to-solution | run_case wall-clock |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SEM direct | 16² | 2 | 1 024 (32²) | 1 024 | — | 7.3e-01 | 3.2e+00 | 152 µs | 3.59 ms | 149 µs | 1.32 s | **1.32 s** | 1.76 s |
| SEM direct | 16² | 3 | 2 304 (48²) | 2 304 | — | 2.2e-01 | 5.1e-01 | 559 µs | 12.9 ms | 468 µs | 1.51 s | **1.53 s** | 2.01 s |
| SEM direct | 16² | 4 | 4 096 (64²) | 4 096 | — | 4.7e-02 | 1.3e-01 | 625 µs | 17 ms | 502 µs | 1.49 s | **1.51 s** | 1.97 s |
| SEM direct | 16² | 5 | 6 400 (80²) | 6 400 | — | 2.5e-03 | 2.0e-03 | 1.11 ms | 34.3 ms | 646 µs | 1.62 s | **1.65 s** | 2.13 s |
| SEM direct | 16² | 6 | 9 216 (96²) | 9 216 | — | 1.2e-03 | 2.8e-03 | 1.81 ms | 72.5 ms | 886 µs | 1.65 s | **1.73 s** | 2.23 s |
| SEM direct | 16² | 7 | 12 544 (112²) | 12 544 | — | 1.7e-04 | 2.5e-04 | 4.69 ms | 155 ms | 1.28 ms | 1.95 s | **2.11 s** | 2.57 s |
| SEM direct | 16² | 8 | 16 384 (128²) | 16 384 | — | 5.0e-05 | 4.0e-05 | 5.84 ms | 256 ms | 1.48 ms | 2.35 s | **2.61 s** | 3.12 s |
| SEM AMG | 16² | 2 | 1 024 (32²) | 1 024 | 31 | 7.3e-01 | 3.2e+00 | 5.66 ms | 2.47 ms | 142 µs | 1.37 s | **1.38 s** | 1.86 s |
| SEM AMG | 16² | 3 | 2 304 (48²) | 2 304 | 46 | 2.2e-01 | 5.1e-01 | 28.9 ms | 7.01 ms | 281 µs | 1.47 s | **1.51 s** | 2 s |
| SEM AMG | 16² | 4 | 4 096 (64²) | 4 096 | 64 | 4.7e-02 | 1.3e-01 | 87.4 ms | 13.5 ms | 450 µs | 1.45 s | **1.56 s** | 2.06 s |
| SEM AMG | 16² | 5 | 6 400 (80²) | 6 400 | 82 | 2.5e-03 | 2.0e-03 | 227 ms | 27 ms | 656 µs | 1.52 s | **1.78 s** | 2.27 s |
| SEM AMG | 16² | 6 | 9 216 (96²) | 9 216 | 102 | 1.2e-03 | 2.8e-03 | 513 ms | 60.9 ms | 1.06 ms | 1.73 s | **2.31 s** | 2.84 s |
| SEM AMG | 16² | 7 | 12 544 (112²) | 12 544 | 117 | 1.7e-04 | 2.5e-04 | 1.43 s | 139 ms | 1.2 ms | 1.77 s | **3.34 s** | 3.81 s |
| SEM AMG | 16² | 8 | 16 384 (128²) | 16 384 | 138 | 5.0e-05 | 4.0e-05 | 3.28 s | 227 ms | 1.49 ms | 2.35 s | **5.86 s** | 6.34 s |
| SC direct | 16² | 2 | 1 024 (32²) | 768 | — | 7.3e-01 | 3.2e+00 | 667 µs | 6.03 ms | 144 µs | 1.49 s | **1.5 s** | 2.02 s |
| SC direct | 16² | 3 | 2 304 (48²) | 1 280 | — | 2.2e-01 | 5.1e-01 | 1.4 ms | 13 ms | 263 µs | 1.52 s | **1.54 s** | 2.03 s |
| SC direct | 16² | 4 | 4 096 (64²) | 1 792 | — | 4.7e-02 | 1.3e-01 | 2.2 ms | 23.9 ms | 443 µs | 1.47 s | **1.5 s** | 2 s |
| SC direct | 16² | 5 | 6 400 (80²) | 2 304 | — | 2.5e-03 | 2.0e-03 | 6.38 ms | 64.2 ms | 874 µs | 1.64 s | **1.71 s** | 2.24 s |
| SC direct | 16² | 6 | 9 216 (96²) | 2 816 | — | 1.2e-03 | 2.8e-03 | 6.83 ms | 71 ms | 900 µs | 1.57 s | **1.65 s** | 2.15 s |
| SC direct | 16² | 7 | 12 544 (112²) | 3 328 | — | 1.7e-04 | 2.5e-04 | 11.1 ms | 120 ms | 1.15 ms | 1.84 s | **1.97 s** | 2.46 s |
| SC direct | 16² | 8 | 16 384 (128²) | 3 840 | — | 5.0e-05 | 4.0e-05 | 37.2 ms | 190 ms | 1.53 ms | 2.09 s | **2.32 s** | 2.79 s |
| SC AMG | 16² | 2 | 1 024 (32²) | 768 | 23 | 7.3e-01 | 3.2e+00 | 3.38 ms | 4.97 ms | 144 µs | 1.38 s | **1.39 s** | 1.97 s |
| SC AMG | 16² | 3 | 2 304 (48²) | 1 280 | 26 | 2.2e-01 | 5.1e-01 | 9 ms | 12.3 ms | 276 µs | 1.57 s | **1.59 s** | 2.06 s |
| SC AMG | 16² | 4 | 4 096 (64²) | 1 792 | 30 | 4.7e-02 | 1.3e-01 | 17.6 ms | 21.1 ms | 413 µs | 1.54 s | **1.58 s** | 2.11 s |
| SC AMG | 16² | 5 | 6 400 (80²) | 2 304 | 33 | 2.5e-03 | 2.0e-03 | 29.1 ms | 37.8 ms | 661 µs | 1.48 s | **1.55 s** | 2.03 s |
| SC AMG | 16² | 6 | 9 216 (96²) | 2 816 | 36 | 1.2e-03 | 2.8e-03 | 44.4 ms | 66.7 ms | 881 µs | 1.69 s | **1.81 s** | 2.3 s |
| SC AMG | 16² | 7 | 12 544 (112²) | 3 328 | 39 | 1.7e-04 | 2.5e-04 | 68.7 ms | 113 ms | 1.17 ms | 1.81 s | **1.99 s** | 2.44 s |
| SC AMG | 16² | 8 | 16 384 (128²) | 3 840 | 41 | 5.0e-05 | 4.0e-05 | 112 ms | 176 ms | 1.48 ms | 2.31 s | **2.6 s** | 3.08 s |
| pseudo-spectral | — | — | 1 024 (32²) | 1 024 | — | 5.6e-01 | 1.5e+00 | 14.3 µs | 1.71 ms | 129 µs | — | **1.85 ms** | 480 ms |
| pseudo-spectral | — | — | 2 304 (48²) | 2 304 | — | 4.6e-02 | 9.9e-02 | 36.2 µs | 3.56 ms | 233 µs | — | **3.83 ms** | 443 ms |
| pseudo-spectral | — | — | 4 096 (64²) | 4 096 | — | 3.5e-03 | 5.0e-03 | 127 µs | 6.25 ms | 384 µs | — | **6.76 ms** | 499 ms |
| pseudo-spectral | — | — | 6 400 (80²) | 6 400 | — | 3.0e-04 | 2.4e-04 | 167 µs | 9.55 ms | 664 µs | — | **10.4 ms** | 496 ms |
| pseudo-spectral | — | — | 9 216 (96²) | 9 216 | — | 3.1e-05 | 1.5e-05 | 422 µs | 16.2 ms | 1.59 ms | — | **18.2 ms** | 551 ms |
| pseudo-spectral | — | — | 12 544 (112²) | 12 544 | — | 3.8e-06 | 1.7e-06 | 503 µs | 19.1 ms | 1.2 ms | — | **20.8 ms** | 476 ms |
| pseudo-spectral | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 631 µs | 28.7 ms | 1.54 ms | — | **30.9 ms** | 482 ms |
| FFT | — | — | 1 024 (32²) | 1 024 | — | 5.6e-01 | 1.5e+00 | 12.1 µs | 97.4 µs | 110 µs | — | **220 µs** | 440 ms |
| FFT | — | — | 2 304 (48²) | 2 304 | — | 4.6e-02 | 9.9e-02 | 39.4 µs | 211 µs | 297 µs | — | **547 µs** | 461 ms |
| FFT | — | — | 4 096 (64²) | 4 096 | — | 3.5e-03 | 5.0e-03 | 36 µs | 125 µs | 370 µs | — | **531 µs** | 473 ms |
| FFT | — | — | 6 400 (80²) | 6 400 | — | 3.0e-04 | 2.4e-04 | 93 µs | 269 µs | 626 µs | — | **988 µs** | 453 ms |
| FFT | — | — | 9 216 (96²) | 9 216 | — | 3.1e-05 | 1.5e-05 | 196 µs | 482 µs | 1.56 ms | — | **2.24 ms** | 509 ms |
| FFT | — | — | 12 544 (112²) | 12 544 | — | 3.8e-06 | 1.7e-06 | 160 µs | 299 µs | 1.19 ms | — | **1.65 ms** | 466 ms |
| FFT | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 214 µs | 172 µs | 1.51 ms | — | **1.9 ms** | 438 ms |

To regenerate the figures from an existing `results.csv`: `python3 tools/periodic_poisson_benchmark/plot.py`.

<!-- ppb_highres_v2:begin — generated by tools/periodic_poisson_benchmark/readme_figures.py; rerun it instead of editing -->

### Benchmark with mesh refinement, 4 levels (16×16 to 128×128 elements)

The same six solvers, with the 16×16 mesh refined uniformly through `:linitial_refine`: level 0 = 16×16 elements, level 1 = 32×32 elements, level 2 = 64×64 elements, level 3 = 128×128 elements, SEM orders N = 2…8. The Fourier grids have the same number of unknowns, (16·2^L·N)². Timings follow the protocol above (second run in one Julia session). Produced by `pipeline.jl` with `outdir = "ppb_highres_v2"`; see [tools/periodic_poisson_benchmark/README.md](tools/periodic_poisson_benchmark/README.md).

<details>
<summary>Results table (all configurations)</summary>

| method | elements | SEM order N | unknowns (grid) | solved for | CG its | ‖e‖∞ | relative ‖e‖₂ | solve | setup | RHS | SEM infrastructure | time-to-solution | run_case wall-clock |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| SEM direct | 16² | 2 | 1 024 (32²) | 1 024 | — | 7.3e-01 | 3.2e+00 | 82.3 µs | 1.58 ms | 71.9 µs | 328 ms | **330 ms** | 491 ms |
| SEM direct | 16² | 3 | 2 304 (48²) | 2 304 | — | 2.2e-01 | 5.1e-01 | 121 µs | 3.41 ms | 113 µs | 333 ms | **337 ms** | 440 ms |
| SEM direct | 16² | 4 | 4 096 (64²) | 4 096 | — | 4.7e-02 | 1.3e-01 | 215 µs | 11.7 ms | 204 µs | 357 ms | **369 ms** | 484 ms |
| SEM direct | 16² | 5 | 6 400 (80²) | 6 400 | — | 2.5e-03 | 2.0e-03 | 484 µs | 22.2 ms | 294 µs | 405 ms | **428 ms** | 564 ms |
| SEM direct | 16² | 6 | 9 216 (96²) | 9 216 | — | 1.2e-03 | 2.8e-03 | 686 µs | 50.7 ms | 454 µs | 708 ms | **760 ms** | 928 ms |
| SEM direct | 16² | 7 | 12 544 (112²) | 12 544 | — | 1.7e-04 | 2.5e-04 | 778 µs | 44.8 ms | 528 µs | 513 ms | **559 ms** | 668 ms |
| SEM direct | 16² | 8 | 16 384 (128²) | 16 384 | — | 5.0e-05 | 4.0e-05 | 4.78 ms | 117 ms | 722 µs | 695 ms | **817 ms** | 987 ms |
| SEM direct | 32² | 2 | 4 096 (64²) | 4 096 | — | 1.9e-01 | 5.1e-01 | 274 µs | 6.99 ms | 193 µs | 376 ms | **384 ms** | 511 ms |
| SEM direct | 32² | 3 | 9 216 (96²) | 9 216 | — | 1.6e-02 | 4.2e-02 | 504 µs | 19.9 ms | 428 µs | 437 ms | **458 ms** | 592 ms |
| SEM direct | 32² | 4 | 16 384 (128²) | 16 384 | — | 8.0e-04 | 1.3e-03 | 823 µs | 35.2 ms | 678 µs | 498 ms | **534 ms** | 673 ms |
| SEM direct | 32² | 5 | 25 600 (160²) | 25 600 | — | 6.4e-05 | 4.3e-05 | 1.46 ms | 60.7 ms | 1.07 ms | 597 ms | **660 ms** | 825 ms |
| SEM direct | 32² | 6 | 36 864 (192²) | 36 864 | — | 1.3e-05 | 8.7e-06 | 2.79 ms | 168 ms | 1.5 ms | 885 ms | **1.06 s** | 1.22 s |
| SEM direct | 32² | 7 | 50 176 (224²) | 50 176 | — | 1.1e-06 | 7.7e-07 | 4.17 ms | 313 ms | 2.03 ms | 1.47 s | **1.79 s** | 1.98 s |
| SEM direct | 32² | 8 | 65 536 (256²) | 65 536 | — | 2.0e-07 | 9.6e-08 | 6.95 ms | 635 ms | 2.71 ms | 2.04 s | **2.68 s** | 2.89 s |
| SEM direct | 64² | 2 | 16 384 (128²) | 16 384 | — | 3.3e-03 | 1.9e-03 | 1.44 ms | 27.4 ms | 706 µs | 1.41 s | **1.44 s** | 1.99 s |
| SEM direct | 64² | 3 | 36 864 (192²) | 36 864 | — | 2.2e-04 | 1.7e-04 | 3.38 ms | 350 ms | 1.55 ms | 1.64 s | **2 s** | 2.76 s |
| SEM direct | 64² | 4 | 65 536 (256²) | 65 536 | — | 8.5e-06 | 7.0e-06 | 5.03 ms | 206 ms | 2.61 ms | 2.1 s | **2.32 s** | 3.04 s |
| SEM direct | 64² | 5 | 102 400 (320²) | 102 400 | — | 1.1e-06 | 4.7e-07 | 9.25 ms | 512 ms | 4.19 ms | 3.12 s | **3.64 s** | 4.3 s |
| SEM direct | 64² | 6 | 147 456 (384²) | 147 456 | — | 3.5e-08 | 1.5e-08 | 18.3 ms | 1.78 s | 5.89 ms | 4.12 s | **5.93 s** | 6.73 s |
| SEM direct | 64² | 7 | 200 704 (448²) | 200 704 | — | 6.5e-09 | 2.6e-09 | 24.5 ms | 2.35 s | 8.05 ms | 5.1 s | **7.48 s** | 8.32 s |
| SEM direct | 64² | 8 | 262 144 (512²) | 262 144 | — | 4.5e-10 | 2.0e-10 | 65.9 ms | 2.95 s | 10.9 ms | 8.07 s | **11.1 s** | 12.6 s |
| SEM direct | 128² | 2 | 65 536 (256²) | 65 536 | — | 1.5e-04 | 6.2e-05 | 4.72 ms | 139 ms | 2.61 ms | 2.06 s | **2.21 s** | 2.82 s |
| SEM direct | 128² | 3 | 147 456 (384²) | 147 456 | — | 4.5e-06 | 2.1e-06 | 11.1 ms | 498 ms | 5.82 ms | 2.74 s | **3.26 s** | 3.87 s |
| SEM direct | 128² | 4 | 262 144 (512²) | 262 144 | — | 1.8e-07 | 6.6e-08 | 23.2 ms | 1.4 s | 10.6 ms | 4.19 s | **5.63 s** | 6.31 s |
| SEM direct | 128² | 5 | 409 600 (640²) | 409 600 | — | 7.3e-09 | 2.5e-09 | 58.8 ms | 4.15 s | 16.9 ms | 6.96 s | **11.2 s** | 12.2 s |
| SEM direct | 128² | 6 | 589 824 (768²) | 589 824 | — | 2.6e-10 | 1.1e-10 | 136 ms | 8.45 s | 25.4 ms | 11.2 s | **19.8 s** | 21.6 s |
| SEM direct | 128² | 7 | 802 816 (896²) | 802 816 | — | 2.4e-11 | 1.6e-11 | 388 ms | 12.7 s | 37 ms | 21 s | **34.2 s** | 37.9 s |
| SEM AMG | 16² | 2 | 1 024 (32²) | 1 024 | 31 | 7.3e-01 | 3.2e+00 | 2.67 ms | 1.48 ms | 73.2 µs | 347 ms | **351 ms** | 508 ms |
| SEM AMG | 16² | 3 | 2 304 (48²) | 2 304 | 46 | 2.2e-01 | 5.1e-01 | 11.4 ms | 2.84 ms | 117 µs | 342 ms | **357 ms** | 460 ms |
| SEM AMG | 16² | 4 | 4 096 (64²) | 4 096 | 64 | 4.7e-02 | 1.3e-01 | 38.7 ms | 5.24 ms | 190 µs | 326 ms | **370 ms** | 471 ms |
| SEM AMG | 16² | 5 | 6 400 (80²) | 6 400 | 82 | 2.5e-03 | 2.0e-03 | 104 ms | 14.6 ms | 289 µs | 406 ms | **525 ms** | 677 ms |
| SEM AMG | 16² | 6 | 9 216 (96²) | 9 216 | 102 | 1.2e-03 | 2.8e-03 | 238 ms | 19.8 ms | 425 µs | 451 ms | **709 ms** | 858 ms |
| SEM AMG | 16² | 7 | 12 544 (112²) | 12 544 | 117 | 1.7e-04 | 2.5e-04 | 467 ms | 32 ms | 522 µs | 482 ms | **982 ms** | 1.09 s |
| SEM AMG | 16² | 8 | 16 384 (128²) | 16 384 | 138 | 5.0e-05 | 4.0e-05 | 887 ms | 89.3 ms | 731 µs | 679 ms | **1.66 s** | 1.78 s |
| SEM AMG | 32² | 2 | 4 096 (64²) | 4 096 | 34 | 1.9e-01 | 5.1e-01 | 12 ms | 4.43 ms | 203 µs | 489 ms | **505 ms** | 647 ms |
| SEM AMG | 32² | 3 | 9 216 (96²) | 9 216 | 50 | 1.6e-02 | 4.2e-02 | 52.5 ms | 16 ms | 481 µs | 519 ms | **588 ms** | 758 ms |
| SEM AMG | 32² | 4 | 16 384 (128²) | 16 384 | 68 | 8.0e-04 | 1.3e-03 | 167 ms | 34.1 ms | 741 µs | 630 ms | **833 ms** | 989 ms |
| SEM AMG | 32² | 5 | 25 600 (160²) | 25 600 | 83 | 6.4e-05 | 4.3e-05 | 421 ms | 109 ms | 1.09 ms | 706 ms | **1.24 s** | 1.41 s |
| SEM AMG | 32² | 6 | 36 864 (192²) | 36 864 | 102 | 1.3e-05 | 8.7e-06 | 952 ms | 126 ms | 1.57 ms | 1.04 s | **2.12 s** | 2.3 s |
| SEM AMG | 32² | 7 | 50 176 (224²) | 50 176 | 121 | 1.1e-06 | 7.7e-07 | 1.93 s | 242 ms | 2.06 ms | 1.34 s | **3.51 s** | 3.69 s |
| SEM AMG | 32² | 8 | 65 536 (256²) | 65 536 | 145 | 2.0e-07 | 9.6e-08 | 3.73 s | 564 ms | 2.68 ms | 2.01 s | **6.31 s** | 6.51 s |
| SEM AMG | 64² | 2 | 16 384 (128²) | 16 384 | 40 | 3.3e-03 | 1.9e-03 | 52.9 ms | 16.6 ms | 695 µs | 500 ms | **571 ms** | 675 ms |
| SEM AMG | 64² | 3 | 36 864 (192²) | 36 864 | 57 | 2.2e-04 | 1.7e-04 | 237 ms | 86.3 ms | 1.55 ms | 1.01 s | **1.34 s** | 1.53 s |
| SEM AMG | 64² | 4 | 65 536 (256²) | 65 536 | 74 | 8.5e-06 | 7.0e-06 | 732 ms | 153 ms | 2.75 ms | 1.16 s | **2.05 s** | 2.25 s |
| SEM AMG | 64² | 5 | 102 400 (320²) | 102 400 | 96 | 1.1e-06 | 4.7e-07 | 1.94 s | 280 ms | 4.28 ms | 1.89 s | **4.11 s** | 4.32 s |
| SEM AMG | 64² | 6 | 147 456 (384²) | 147 456 | 115 | 3.5e-08 | 1.5e-08 | 4.59 s | 2.37 s | 6.1 ms | 3.27 s | **10.2 s** | 10.6 s |
| SEM AMG | 64² | 7 | 200 704 (448²) | 200 704 | 129 | 6.5e-09 | 2.6e-09 | 8.25 s | 2.04 s | 8.22 ms | 4.23 s | **14.5 s** | 15 s |
| SEM AMG | 64² | 8 | 262 144 (512²) | 262 144 | 147 | 4.5e-10 | 1.9e-10 | 15.3 s | 3.31 s | 10.4 ms | 7.48 s | **26.1 s** | 26.5 s |
| SEM AMG | 128² | 2 | 65 536 (256²) | 65 536 | 46 | 1.5e-04 | 6.2e-05 | 246 ms | 61.2 ms | 2.59 ms | 1.15 s | **1.46 s** | 1.61 s |
| SEM AMG | 128² | 3 | 147 456 (384²) | 147 456 | 64 | 4.5e-06 | 2.1e-06 | 1.05 s | 455 ms | 5.83 ms | 1.86 s | **3.37 s** | 3.54 s |
| SEM AMG | 128² | 4 | 262 144 (512²) | 262 144 | 85 | 1.8e-07 | 6.6e-08 | 3.37 s | 515 ms | 10.4 ms | 3.1 s | **6.99 s** | 7.51 s |
| SC direct | 16² | 2 | 1 024 (32²) | 768 | — | 7.3e-01 | 3.2e+00 | 209 µs | 3.19 ms | 90.8 µs | 531 ms | **534 ms** | 683 ms |
| SC direct | 16² | 3 | 2 304 (48²) | 1 280 | — | 2.2e-01 | 5.1e-01 | 339 µs | 4.34 ms | 112 µs | 331 ms | **336 ms** | 434 ms |
| SC direct | 16² | 4 | 4 096 (64²) | 1 792 | — | 4.7e-02 | 1.3e-01 | 615 µs | 9.39 ms | 201 µs | 359 ms | **369 ms** | 481 ms |
| SC direct | 16² | 5 | 6 400 (80²) | 2 304 | — | 2.5e-03 | 2.0e-03 | 1.53 ms | 25.1 ms | 364 µs | 484 ms | **511 ms** | 650 ms |
| SC direct | 16² | 6 | 9 216 (96²) | 2 816 | — | 1.2e-03 | 2.8e-03 | 2.59 ms | 37.3 ms | 400 µs | 411 ms | **451 ms** | 553 ms |
| SC direct | 16² | 7 | 12 544 (112²) | 3 328 | — | 1.7e-04 | 2.5e-04 | 4.64 ms | 56.4 ms | 574 µs | 495 ms | **557 ms** | 674 ms |
| SC direct | 16² | 8 | 16 384 (128²) | 3 840 | — | 5.0e-05 | 4.0e-05 | 8.41 ms | 158 ms | 718 µs | 871 ms | **1.04 s** | 1.21 s |
| SC direct | 32² | 2 | 4 096 (64²) | 3 072 | — | 1.9e-01 | 5.1e-01 | 859 µs | 17.7 ms | 278 µs | 755 ms | **774 ms** | 887 ms |
| SC direct | 32² | 3 | 9 216 (96²) | 5 120 | — | 1.6e-02 | 4.2e-02 | 1.42 ms | 32.3 ms | 418 µs | 439 ms | **473 ms** | 596 ms |
| SC direct | 32² | 4 | 16 384 (128²) | 7 168 | — | 8.0e-04 | 1.3e-03 | 2.79 ms | 64.7 ms | 789 µs | 557 ms | **625 ms** | 768 ms |
| SC direct | 32² | 5 | 25 600 (160²) | 9 216 | — | 6.4e-05 | 4.3e-05 | 5.57 ms | 213 ms | 1.25 ms | 758 ms | **978 ms** | 1.14 s |
| SC direct | 32² | 6 | 36 864 (192²) | 11 264 | — | 1.3e-05 | 8.7e-06 | 18 ms | 284 ms | 2 ms | 2.28 s | **2.59 s** | 2.86 s |
| SC direct | 32² | 7 | 50 176 (224²) | 13 312 | — | 1.1e-06 | 7.7e-07 | 18.8 ms | 306 ms | 2.05 ms | 1.45 s | **1.78 s** | 1.96 s |
| SC direct | 32² | 8 | 65 536 (256²) | 15 360 | — | 2.0e-07 | 9.6e-08 | 34.8 ms | 529 ms | 2.73 ms | 2.03 s | **2.6 s** | 2.82 s |
| SC direct | 64² | 2 | 16 384 (128²) | 12 288 | — | 3.3e-03 | 1.9e-03 | 2.45 ms | 55.5 ms | 677 µs | 500 ms | **559 ms** | 990 ms |
| SC direct | 64² | 3 | 36 864 (192²) | 20 480 | — | 2.2e-04 | 1.7e-04 | 6.99 ms | 320 ms | 1.62 ms | 756 ms | **1.08 s** | 1.64 s |
| SC direct | 64² | 4 | 65 536 (256²) | 28 672 | — | 8.5e-06 | 7.0e-06 | 17 ms | 561 ms | 2.74 ms | 1.19 s | **1.77 s** | 2.48 s |
| SC direct | 64² | 5 | 102 400 (320²) | 36 864 | — | 1.1e-06 | 4.7e-07 | 55.9 ms | 1.5 s | 4.33 ms | 1.61 s | **3.17 s** | 3.77 s |
| SC direct | 64² | 6 | 147 456 (384²) | 45 056 | — | 3.5e-08 | 1.5e-08 | 72.9 ms | 1.99 s | 6.15 ms | 2.92 s | **4.98 s** | 5.63 s |
| SC direct | 64² | 7 | 200 704 (448²) | 53 248 | — | 6.5e-09 | 2.6e-09 | 134 ms | 3.03 s | 7.97 ms | 3.82 s | **6.99 s** | 7.67 s |
| SC direct | 64² | 8 | 262 144 (512²) | 61 440 | — | 4.5e-10 | 1.9e-10 | 169 ms | 3.65 s | 10.8 ms | 6.94 s | **10.8 s** | 12.2 s |
| SC direct | 128² | 2 | 65 536 (256²) | 49 152 | — | 1.5e-04 | 6.2e-05 | 10.6 ms | 187 ms | 2.63 ms | 1.14 s | **1.34 s** | 1.82 s |
| SC direct | 128² | 3 | 147 456 (384²) | 81 920 | — | 4.5e-06 | 2.1e-06 | 26.3 ms | 688 ms | 5.78 ms | 1.85 s | **2.57 s** | 3.15 s |
| SC direct | 128² | 4 | 262 144 (512²) | 114 688 | — | 1.8e-07 | 6.6e-08 | 74.3 ms | 1.69 s | 270 ms | 3.12 s | **5.16 s** | 5.89 s |
| SC direct | 128² | 5 | 409 600 (640²) | 147 456 | — | 7.3e-09 | 2.5e-09 | 178 ms | 3.74 s | 16.4 ms | 5.68 s | **9.61 s** | 11 s |
| SC direct | 128² | 6 | 589 824 (768²) | 180 224 | — | 2.5e-10 | 1.1e-10 | 254 ms | 6.37 s | 27.1 ms | 10.5 s | **17.2 s** | 18.9 s |
| SC direct | 128² | 7 | 802 816 (896²) | 212 992 | — | 3.0e-11 | 2.7e-11 | 10.2 s | 12.5 s | 45 ms | 20.6 s | **43.3 s** | 46.2 s |
| SC AMG | 16² | 2 | 1 024 (32²) | 768 | 23 | 7.3e-01 | 3.2e+00 | 1.74 ms | 7.7 ms | 383 µs | 565 ms | **574 ms** | 738 ms |
| SC AMG | 16² | 3 | 2 304 (48²) | 1 280 | 26 | 2.2e-01 | 5.1e-01 | 3.8 ms | 3.88 ms | 126 µs | 325 ms | **333 ms** | 437 ms |
| SC AMG | 16² | 4 | 4 096 (64²) | 1 792 | 30 | 4.7e-02 | 1.3e-01 | 7.55 ms | 10.8 ms | 186 µs | 337 ms | **356 ms** | 460 ms |
| SC AMG | 16² | 5 | 6 400 (80²) | 2 304 | 33 | 2.5e-03 | 2.0e-03 | 16 ms | 37.6 ms | 338 µs | 468 ms | **522 ms** | 644 ms |
| SC AMG | 16² | 6 | 9 216 (96²) | 2 816 | 36 | 1.2e-03 | 2.8e-03 | 19.8 ms | 28.2 ms | 394 µs | 414 ms | **462 ms** | 573 ms |
| SC AMG | 16² | 7 | 12 544 (112²) | 3 328 | 39 | 1.7e-04 | 2.5e-04 | 31.4 ms | 49 ms | 531 µs | 509 ms | **590 ms** | 703 ms |
| SC AMG | 16² | 8 | 16 384 (128²) | 3 840 | 41 | 5.0e-05 | 4.0e-05 | 43.6 ms | 115 ms | 784 µs | 762 ms | **922 ms** | 1.08 s |
| SC AMG | 32² | 2 | 4 096 (64²) | 3 072 | 25 | 1.9e-01 | 5.1e-01 | 7.07 ms | 6.71 ms | 194 µs | 476 ms | **490 ms** | 639 ms |
| SC AMG | 32² | 3 | 9 216 (96²) | 5 120 | 32 | 1.6e-02 | 4.2e-02 | 18.6 ms | 14.7 ms | 405 µs | 395 ms | **429 ms** | 534 ms |
| SC AMG | 32² | 4 | 16 384 (128²) | 7 168 | 36 | 8.0e-04 | 1.3e-03 | 35.4 ms | 29.6 ms | 677 µs | 485 ms | **550 ms** | 667 ms |
| SC AMG | 32² | 5 | 25 600 (160²) | 9 216 | 39 | 6.4e-05 | 4.3e-05 | 59.4 ms | 141 ms | 1.12 ms | 762 ms | **964 ms** | 1.14 s |
| SC AMG | 32² | 6 | 36 864 (192²) | 11 264 | 43 | 1.3e-05 | 8.7e-06 | 96.4 ms | 178 ms | 1.59 ms | 1.35 s | **1.62 s** | 1.89 s |
| SC AMG | 32² | 7 | 50 176 (224²) | 13 312 | 45 | 1.1e-06 | 7.7e-07 | 301 ms | 277 ms | 2.24 ms | 1.6 s | **2.18 s** | 2.52 s |
| SC AMG | 32² | 8 | 65 536 (256²) | 15 360 | 48 | 2.0e-07 | 9.6e-08 | 205 ms | 581 ms | 2.77 ms | 1.74 s | **2.53 s** | 2.71 s |
| SC AMG | 64² | 2 | 16 384 (128²) | 12 288 | 29 | 3.3e-03 | 1.9e-03 | 30.9 ms | 26.2 ms | 673 µs | 501 ms | **559 ms** | 664 ms |
| SC AMG | 64² | 3 | 36 864 (192²) | 20 480 | 36 | 2.2e-04 | 1.7e-04 | 89.5 ms | 238 ms | 1.61 ms | 1.28 s | **1.6 s** | 1.77 s |
| SC AMG | 64² | 4 | 65 536 (256²) | 28 672 | 41 | 8.5e-06 | 7.0e-06 | 161 ms | 215 ms | 2.61 ms | 1.31 s | **1.68 s** | 1.91 s |
| SC AMG | 64² | 5 | 102 400 (320²) | 36 864 | 45 | 1.1e-06 | 4.7e-07 | 787 ms | 752 ms | 4.19 ms | 1.65 s | **3.19 s** | 3.6 s |
| SC AMG | 64² | 6 | 147 456 (384²) | 45 056 | 49 | 3.5e-08 | 1.5e-08 | 443 ms | 1.48 s | 5.85 ms | 2.59 s | **4.52 s** | 4.74 s |
| SC AMG | 64² | 7 | 200 704 (448²) | 53 248 | 52 | 6.5e-09 | 2.6e-09 | 662 ms | 2.78 s | 8.23 ms | 3.92 s | **7.37 s** | 7.64 s |
| SC AMG | 64² | 8 | 262 144 (512²) | 61 440 | 55 | 4.4e-10 | 1.9e-10 | 1.14 s | 2.95 s | 10.7 ms | 6.78 s | **10.9 s** | 11.3 s |
| SC AMG | 128² | 2 | 65 536 (256²) | 49 152 | 36 | 1.5e-04 | 6.2e-05 | 153 ms | 106 ms | 2.61 ms | 1.14 s | **1.4 s** | 1.55 s |
| SC AMG | 128² | 3 | 147 456 (384²) | 81 920 | 42 | 4.5e-06 | 2.1e-06 | 383 ms | 563 ms | 5.78 ms | 1.84 s | **2.79 s** | 2.97 s |
| SC AMG | 128² | 4 | 262 144 (512²) | 114 688 | 48 | 1.8e-07 | 6.6e-08 | 823 ms | 2.64 s | 10.9 ms | 3.53 s | **7 s** | 7.33 s |
| SC AMG | 128² | 5 | 409 600 (640²) | 147 456 | 52 | 7.3e-09 | 2.5e-09 | 1.36 s | 4.04 s | 17 ms | 5.81 s | **11.2 s** | 11.7 s |
| SC AMG | 128² | 6 | 589 824 (768²) | 180 224 | 56 | 2.5e-10 | 1.1e-10 | 2.14 s | 5.75 s | 24.1 ms | 11.1 s | **19 s** | 19.9 s |
| SC AMG | 128² | 7 | 802 816 (896²) | 212 992 | 60 | 2.8e-11 | 2.4e-11 | 9.87 s | 11.2 s | 37.8 ms | 25.3 s | **46.4 s** | 48.6 s |
| pseudo-spectral | — | — | 1 024 (32²) | 1 024 | — | 5.6e-01 | 1.5e+00 | 10.3 µs | 863 µs | 67.2 µs | — | **941 µs** | 143 ms |
| pseudo-spectral | — | — | 2 304 (48²) | 2 304 | — | 4.6e-02 | 9.9e-02 | 23.3 µs | 1.55 ms | 104 µs | — | **1.67 ms** | 101 ms |
| pseudo-spectral | — | — | 4 096 (64²) | 4 096 | — | 3.5e-03 | 5.0e-03 | 96.8 µs | 2.94 ms | 172 µs | — | **3.21 ms** | 108 ms |
| pseudo-spectral | — | — | 6 400 (80²) | 6 400 | — | 3.0e-04 | 2.4e-04 | 111 µs | 5.61 ms | 300 µs | — | **6.02 ms** | 129 ms |
| pseudo-spectral | — | — | 9 216 (96²) | 9 216 | — | 3.1e-05 | 1.5e-05 | 173 µs | 6.92 ms | 390 µs | — | **7.48 ms** | 106 ms |
| pseudo-spectral | — | — | 12 544 (112²) | 12 544 | — | 3.8e-06 | 1.7e-06 | 251 µs | 8.9 ms | 505 µs | — | **9.66 ms** | 107 ms |
| pseudo-spectral | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 507 µs | 15.3 ms | 650 µs | — | **16.5 ms** | 169 ms |
| pseudo-spectral | — | — | 4 096 (64²) | 4 096 | — | 3.5e-03 | 5.0e-03 | 96.5 µs | 2.93 ms | 198 µs | — | **3.23 ms** | 112 ms |
| pseudo-spectral | — | — | 9 216 (96²) | 9 216 | — | 3.1e-05 | 1.5e-05 | 168 µs | 6.45 ms | 405 µs | — | **7.03 ms** | 104 ms |
| pseudo-spectral | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 358 µs | 12.2 ms | 647 µs | — | **13.2 ms** | 119 ms |
| pseudo-spectral | — | — | 25 600 (160²) | 25 600 | — | 1.0e-08 | 4.9e-09 | 684 µs | 19 ms | 1.16 ms | — | **20.8 ms** | 131 ms |
| pseudo-spectral | — | — | 36 864 (192²) | 36 864 | — | 2.3e-10 | 1.1e-10 | 1.31 ms | 27.9 ms | 1.55 ms | — | **30.7 ms** | 139 ms |
| pseudo-spectral | — | — | 50 176 (224²) | 50 176 | — | 5.5e-12 | 3.2e-12 | 1.88 ms | 38.6 ms | 2.14 ms | — | **42.6 ms** | 156 ms |
| pseudo-spectral | — | — | 65 536 (256²) | 65 536 | — | 1.3e-12 | 3.9e-12 | 3.9 ms | 55.6 ms | 2.69 ms | — | **62.2 ms** | 184 ms |
| pseudo-spectral | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 360 µs | 11.9 ms | 648 µs | — | **12.9 ms** | 107 ms |
| pseudo-spectral | — | — | 36 864 (192²) | 36 864 | — | 2.3e-10 | 1.1e-10 | 1.41 ms | 29.3 ms | 1.57 ms | — | **32.3 ms** | 190 ms |
| pseudo-spectral | — | — | 65 536 (256²) | 65 536 | — | 1.3e-12 | 3.9e-12 | 3.7 ms | 60.2 ms | 3.78 ms | — | **67.7 ms** | 232 ms |
| pseudo-spectral | — | — | 102 400 (320²) | 102 400 | — | 2.1e-12 | 4.2e-12 | 5.71 ms | 90.3 ms | 4.04 ms | — | **100 ms** | 263 ms |
| pseudo-spectral | — | — | 147 456 (384²) | 147 456 | — | 7.2e-13 | 1.2e-12 | 9.14 ms | 122 ms | 6.05 ms | — | **137 ms** | 249 ms |
| pseudo-spectral | — | — | 200 704 (448²) | 200 704 | — | 4.4e-13 | 1.2e-12 | 14.8 ms | 172 ms | 7.85 ms | — | **195 ms** | 298 ms |
| pseudo-spectral | — | — | 262 144 (512²) | 262 144 | — | 9.1e-13 | 2.1e-12 | 21.2 ms | 231 ms | 10.8 ms | — | **263 ms** | 376 ms |
| pseudo-spectral | — | — | 65 536 (256²) | 65 536 | — | 1.3e-12 | 3.9e-12 | 2.8 ms | 50.1 ms | 2.55 ms | — | **55.4 ms** | 154 ms |
| pseudo-spectral | — | — | 147 456 (384²) | 147 456 | — | 7.2e-13 | 1.2e-12 | 9.06 ms | 118 ms | 5.69 ms | — | **133 ms** | 232 ms |
| pseudo-spectral | — | — | 262 144 (512²) | 262 144 | — | 9.1e-13 | 2.1e-12 | 22.1 ms | 243 ms | 11 ms | — | **276 ms** | 385 ms |
| pseudo-spectral | — | — | 409 600 (640²) | 409 600 | — | 2.0e-12 | 4.7e-12 | 42.2 ms | 405 ms | 16.9 ms | — | **464 ms** | 603 ms |
| pseudo-spectral | — | — | 589 824 (768²) | 589 824 | — | 4.4e-12 | 1.6e-11 | 70.9 ms | 668 ms | 24.2 ms | — | **763 ms** | 965 ms |
| pseudo-spectral | — | — | 802 816 (896²) | 802 816 | — | 1.3e-11 | 4.2e-11 | 115 ms | 1.22 s | 39.7 ms | — | **1.37 s** | 1.62 s |
| FFT | — | — | 1 024 (32²) | 1 024 | — | 5.6e-01 | 1.5e+00 | 5.58 µs | 33.8 µs | 53.7 µs | — | **93 µs** | 97.5 ms |
| FFT | — | — | 2 304 (48²) | 2 304 | — | 4.6e-02 | 9.9e-02 | 11.5 µs | 72.7 µs | 103 µs | — | **187 µs** | 95 ms |
| FFT | — | — | 4 096 (64²) | 4 096 | — | 3.5e-03 | 5.0e-03 | 17 µs | 34.2 µs | 173 µs | — | **225 µs** | 102 ms |
| FFT | — | — | 6 400 (80²) | 6 400 | — | 3.0e-04 | 2.4e-04 | 52.8 µs | 119 µs | 275 µs | — | **446 µs** | 108 ms |
| FFT | — | — | 9 216 (96²) | 9 216 | — | 3.1e-05 | 1.5e-05 | 44.9 µs | 84.6 µs | 412 µs | — | **541 µs** | 107 ms |
| FFT | — | — | 12 544 (112²) | 12 544 | — | 3.8e-06 | 1.7e-06 | 60.2 µs | 78.8 µs | 504 µs | — | **644 µs** | 95.4 ms |
| FFT | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 69.8 µs | 61.5 µs | 648 µs | — | **779 µs** | 97.1 ms |
| FFT | — | — | 4 096 (64²) | 4 096 | — | 3.5e-03 | 5.0e-03 | 16.2 µs | 34.6 µs | 181 µs | — | **232 µs** | 98.2 ms |
| FFT | — | — | 9 216 (96²) | 9 216 | — | 3.1e-05 | 1.5e-05 | 42.5 µs | 75.2 µs | 378 µs | — | **496 µs** | 95.4 ms |
| FFT | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 69 µs | 57.4 µs | 659 µs | — | **786 µs** | 96.5 ms |
| FFT | — | — | 25 600 (160²) | 25 600 | — | 1.0e-08 | 4.9e-09 | 176 µs | 197 µs | 1.08 ms | — | **1.45 ms** | 98.3 ms |
| FFT | — | — | 36 864 (192²) | 36 864 | — | 2.3e-10 | 1.1e-10 | 167 µs | 136 µs | 1.44 ms | — | **1.75 ms** | 102 ms |
| FFT | — | — | 50 176 (224²) | 50 176 | — | 5.4e-12 | 2.6e-12 | 287 µs | 120 µs | 1.97 ms | — | **2.37 ms** | 99.3 ms |
| FFT | — | — | 65 536 (256²) | 65 536 | — | 1.3e-13 | 6.3e-14 | 474 µs | 697 µs | 2.94 ms | — | **4.11 ms** | 133 ms |
| FFT | — | — | 16 384 (128²) | 16 384 | — | 5.0e-07 | 2.3e-07 | 685 µs | 46.2 µs | 670 µs | — | **1.4 ms** | 95.2 ms |
| FFT | — | — | 36 864 (192²) | 36 864 | — | 2.3e-10 | 1.1e-10 | 227 µs | 167 µs | 1.5 ms | — | **1.9 ms** | 113 ms |
| FFT | — | — | 65 536 (256²) | 65 536 | — | 1.3e-13 | 6.3e-14 | 1.12 ms | 3.34 ms | 2.62 ms | — | **7.09 ms** | 126 ms |
| FFT | — | — | 102 400 (320²) | 102 400 | — | 6.6e-15 | 9.3e-15 | 1.17 ms | 2.65 ms | 4.12 ms | — | **7.94 ms** | 152 ms |
| FFT | — | — | 147 456 (384²) | 147 456 | — | 6.1e-15 | 6.1e-15 | 1.01 ms | 2.19 ms | 5.71 ms | — | **8.92 ms** | 111 ms |
| FFT | — | — | 200 704 (448²) | 200 704 | — | 6.3e-15 | 9.0e-15 | 1.37 ms | 1.36 ms | 7.75 ms | — | **10.5 ms** | 113 ms |
| FFT | — | — | 262 144 (512²) | 262 144 | — | 7.1e-15 | 1.2e-14 | 1.57 ms | 1.9 ms | 10.1 ms | — | **13.6 ms** | 119 ms |
| FFT | — | — | 65 536 (256²) | 65 536 | — | 1.3e-13 | 6.3e-14 | 304 µs | 117 µs | 2.55 ms | — | **2.97 ms** | 101 ms |
| FFT | — | — | 147 456 (384²) | 147 456 | — | 6.1e-15 | 6.1e-15 | 1.75 ms | 212 µs | 5.76 ms | — | **7.72 ms** | 105 ms |
| FFT | — | — | 262 144 (512²) | 262 144 | — | 7.1e-15 | 1.2e-14 | 1.72 ms | 911 µs | 10.4 ms | — | **13 ms** | 123 ms |
| FFT | — | — | 409 600 (640²) | 409 600 | — | 9.0e-15 | 1.7e-14 | 2.31 ms | 3.15 ms | 16.2 ms | — | **21.7 ms** | 132 ms |
| FFT | — | — | 589 824 (768²) | 589 824 | — | 6.7e-15 | 9.1e-15 | 3.85 ms | 2.8 ms | 23 ms | — | **29.6 ms** | 154 ms |
| FFT | — | — | 802 816 (896²) | 802 816 | — | 6.6e-15 | 9.3e-15 | 5.63 ms | 18.7 ms | 32.5 ms | — | **56.8 ms** | 192 ms |

</details>

#### Time versus unknowns under mesh refinement

At fixed SEM order N, one curve per solver as the mesh is refined: where the AMG curves cross the direct-solve curves, if they do. *Cost* is the solver's own setup plus solve, without the SEM infrastructure; *total* is the time-to-solution including it.

**SEM order N = 2**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_solve_N2-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_solve_N2.svg" width="680" alt="Wall-clock of the solve step versus number of unknowns under uniform mesh refinement at SEM order 2, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_cost_N2-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_cost_N2.svg" width="680" alt="Wall-clock of setup plus solve (solver cost) versus number of unknowns under uniform mesh refinement at SEM order 2, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_total_N2-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_total_N2.svg" width="680" alt="Wall-clock of the time-to-solution including all infrastructure versus number of unknowns under uniform mesh refinement at SEM order 2, one curve per solver.">
</picture>

**SEM order N = 3**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_solve_N3-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_solve_N3.svg" width="680" alt="Wall-clock of the solve step versus number of unknowns under uniform mesh refinement at SEM order 3, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_cost_N3-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_cost_N3.svg" width="680" alt="Wall-clock of setup plus solve (solver cost) versus number of unknowns under uniform mesh refinement at SEM order 3, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_total_N3-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_total_N3.svg" width="680" alt="Wall-clock of the time-to-solution including all infrastructure versus number of unknowns under uniform mesh refinement at SEM order 3, one curve per solver.">
</picture>

**SEM order N = 4**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_solve_N4-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_solve_N4.svg" width="680" alt="Wall-clock of the solve step versus number of unknowns under uniform mesh refinement at SEM order 4, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_cost_N4-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_cost_N4.svg" width="680" alt="Wall-clock of setup plus solve (solver cost) versus number of unknowns under uniform mesh refinement at SEM order 4, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_total_N4-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_total_N4.svg" width="680" alt="Wall-clock of the time-to-solution including all infrastructure versus number of unknowns under uniform mesh refinement at SEM order 4, one curve per solver.">
</picture>

**SEM order N = 5**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_solve_N5-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_solve_N5.svg" width="680" alt="Wall-clock of the solve step versus number of unknowns under uniform mesh refinement at SEM order 5, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_cost_N5-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_cost_N5.svg" width="680" alt="Wall-clock of setup plus solve (solver cost) versus number of unknowns under uniform mesh refinement at SEM order 5, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_total_N5-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_total_N5.svg" width="680" alt="Wall-clock of the time-to-solution including all infrastructure versus number of unknowns under uniform mesh refinement at SEM order 5, one curve per solver.">
</picture>

**SEM order N = 6**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_solve_N6-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_solve_N6.svg" width="680" alt="Wall-clock of the solve step versus number of unknowns under uniform mesh refinement at SEM order 6, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_cost_N6-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_cost_N6.svg" width="680" alt="Wall-clock of setup plus solve (solver cost) versus number of unknowns under uniform mesh refinement at SEM order 6, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_total_N6-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_total_N6.svg" width="680" alt="Wall-clock of the time-to-solution including all infrastructure versus number of unknowns under uniform mesh refinement at SEM order 6, one curve per solver.">
</picture>

**SEM order N = 7**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_solve_N7-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_solve_N7.svg" width="680" alt="Wall-clock of the solve step versus number of unknowns under uniform mesh refinement at SEM order 7, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_cost_N7-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_cost_N7.svg" width="680" alt="Wall-clock of setup plus solve (solver cost) versus number of unknowns under uniform mesh refinement at SEM order 7, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_total_N7-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_total_N7.svg" width="680" alt="Wall-clock of the time-to-solution including all infrastructure versus number of unknowns under uniform mesh refinement at SEM order 7, one curve per solver.">
</picture>

**SEM order N = 8**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_solve_N8-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_solve_N8.svg" width="680" alt="Wall-clock of the solve step versus number of unknowns under uniform mesh refinement at SEM order 8, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_cost_N8-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_cost_N8.svg" width="680" alt="Wall-clock of setup plus solve (solver cost) versus number of unknowns under uniform mesh refinement at SEM order 8, one curve per solver.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_hrefine_total_N8-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_hrefine_total_N8.svg" width="680" alt="Wall-clock of the time-to-solution including all infrastructure versus number of unknowns under uniform mesh refinement at SEM order 8, one curve per solver.">
</picture>

#### Error and time at each mesh level

**Level 0: 16×16 elements**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_order_L0-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_order_L0.svg" width="680" alt="L-infinity error versus SEM order N on the 16×16 mesh (one SEM curve for its four solvers), with the pseudo-spectral and FFT solvers at the same number of unknowns.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_dofs_L0-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_dofs_L0.svg" width="680" alt="L-infinity error versus number of unknowns on the 16×16 mesh.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_solve_time_L0-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_solve_time_L0.svg" width="680" alt="L-infinity error versus the wall-clock of the solve step on the 16×16 mesh, for the six solvers.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_total_time_L0-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_total_time_L0.svg" width="680" alt="L-infinity error versus time-to-solution including all infrastructure on the 16×16 mesh, for the six solvers.">
</picture>

**Level 1: 32×32 elements**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_order_L1-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_order_L1.svg" width="680" alt="L-infinity error versus SEM order N on the 32×32 mesh (one SEM curve for its four solvers), with the pseudo-spectral and FFT solvers at the same number of unknowns.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_dofs_L1-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_dofs_L1.svg" width="680" alt="L-infinity error versus number of unknowns on the 32×32 mesh.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_solve_time_L1-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_solve_time_L1.svg" width="680" alt="L-infinity error versus the wall-clock of the solve step on the 32×32 mesh, for the six solvers.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_total_time_L1-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_total_time_L1.svg" width="680" alt="L-infinity error versus time-to-solution including all infrastructure on the 32×32 mesh, for the six solvers.">
</picture>

**Level 2: 64×64 elements**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_order_L2-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_order_L2.svg" width="680" alt="L-infinity error versus SEM order N on the 64×64 mesh (one SEM curve for its four solvers), with the pseudo-spectral and FFT solvers at the same number of unknowns.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_dofs_L2-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_dofs_L2.svg" width="680" alt="L-infinity error versus number of unknowns on the 64×64 mesh.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_solve_time_L2-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_solve_time_L2.svg" width="680" alt="L-infinity error versus the wall-clock of the solve step on the 64×64 mesh, for the six solvers.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_total_time_L2-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_total_time_L2.svg" width="680" alt="L-infinity error versus time-to-solution including all infrastructure on the 64×64 mesh, for the six solvers.">
</picture>

**Level 3: 128×128 elements**

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_order_L3-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_order_L3.svg" width="680" alt="L-infinity error versus SEM order N on the 128×128 mesh (one SEM curve for its four solvers), with the pseudo-spectral and FFT solvers at the same number of unknowns.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_dofs_L3-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_dofs_L3.svg" width="680" alt="L-infinity error versus number of unknowns on the 128×128 mesh.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_solve_time_L3-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_solve_time_L3.svg" width="680" alt="L-infinity error versus the wall-clock of the solve step on the 128×128 mesh, for the six solvers.">
</picture>

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="assets/ppb_highres_v2/ppb_error_vs_total_time_L3-dark.svg">
  <img src="assets/ppb_highres_v2/ppb_error_vs_total_time_L3.svg" width="680" alt="L-infinity error versus time-to-solution including all infrastructure on the 128×128 mesh, for the six solvers.">
</picture>
<!-- ppb_highres_v2:end -->

## Laguerre semi-infinite element test suite
This section contains instructions to run all of the test cases presented in

```
@article{tissaoui2024,
  author = {Y. Tissaoui and J. F. Kelly and S. Marras}
  title = {Efficient Spectral Element Method for the Euler Equations on Unbounded Domains},
  volume ={487},
  pages={129080},
  year = {2024},
  journal = {App. Math. Comput.},
}
```

### Test 1: 1D wave equation with Laguerre semi-infinite element absorbing layers

The problem is defined in [`problems/CompEuler/wave1d_lag`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/CompEuler/wave1d_lag) and by default output will be written to `output/CompEuler/wave1d_lag`. To solve this problem run the following commands from the Julia command line:

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "wave1d_lag")
```

<img src="assets/wave_v_4.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

### Test 2: 1D wave train for linearized shallow water equations

The problem is defined in [`problems/equations/AdvDiff/Wave_Train`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/AdvDiff/Wave_Train) and by default output will be written to `output/AdvDiff/Wave_Train`. To solve this problem run the following commands from the Julia command line:

```julia
using Jexpresso
Jexpresso.run_case("AdvDiff", "Wave_Train")
```

<img src="assets/Wave_Train_final.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />


### Test 3: 2D advection-diffusion equation

The problem is defined in [`problems/equations/AdvDiff/2D_laguerre`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/AdvDiff/2d_Laguerre) and by default output will be written to `output/AdvDiff/2D_laguerre`. To solve this problem run the following commands from the Julia command line:

```julia
using Jexpresso
Jexpresso.run_case("AdvDiff", "2D_laguerre")
```

<img src="assets/ad2d-4s-line.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

### Test 4: 2D Helmholtz equation

The problem is defined in [`problems/equations/Helmholtz/case1_laguerre`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/Helmholtz/case1_laguerre) and by default output will be written to `output/Helmholtz/case1_laguerre`. To solve this problem run the following commands from the Julia command line:

```julia
using Jexpresso
Jexpresso.run_case("Helmholtz", "case1_laguerre")
```

<img src="assets/Helmholtz_from_jexpresso-line.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

## Rising thermal bubble with semi-infinite Laguerre elements for outflows

The problem is defined in [`problems/equations/CompEuler/theta_laguerre`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/CompEuler/theta_laguerre) and by default output will be written to `output/CompEuler/theta_laguerre`. To solve this problem run the following commands from the Julia command line:

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "theta_laguerre")
```

<img src="assets/48.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />


## Hydrostatic linear mountain waves with semi-infinite Laguerre elements for outflows

The problem is defined in [`problems/equations/CompEuler/HSmount_Lag`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/CompEuler/HSmount_Lag) and by default output will be written to `output/CompEuler/HSmount_Lag`. To solve this problem run the following commands from the Julia command line:

```bash      
using Jexpresso
Jexpresso.run_case("CompEuler", "HSmount_Lag")
```

<img src="assets/wvelo.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

## Non-hydrostatic mountain waves: comparison against WRF

<img src="assets/NHjexpVSwrf.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />
