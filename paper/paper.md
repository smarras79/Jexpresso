---
title: 'Jexpresso.jl: An Extensible Julia Software for the Solution of Partial Differential Equations on CPUs & GPUs with Minimal User Intervention'
tags:
  - Julia
  - spectral element method
  - conservation laws
  - computational fluid dynamics
  - atmospheric modeling
  - high-performance computing
  - GPU
authors:
  - name: Simone Marras
    orcid: 0000-0000-0000-0000
    corresponding: true
    affiliation: "1, 2"
  - name: Yassine Tissaoui
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Hang Wang
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
  - name: Department of Mechanical & Industrial Engineering, New Jersey Institute of Technology, Newark, NJ, USA
    index: 1
  - name: Department of Mathematical Sciences, New Jersey Institute of Technology, Newark, NJ, USA
    index: 2
  - name: Department of Mathematics, University of Wisconsin, Madison, WI, USA
    index: 3
date: 1 January 2026
bibliography: paper.bib
---

# Summary

Jexpresso is an open-source Julia package for solving systems of conservation
laws in one, two, and three spatial dimensions. Its guiding idea is that a user
should be able to add a new physical problem by writing down what they would
write on paper---a solution vector, its fluxes, its sources, and its boundary
conditions---and nothing else. Everything downstream, namely the spatial
discretization, stabilization, time integration, parallel decomposition, and the
choice of CPU or GPU execution, is supplied by the package and is deliberately
invisible to the problem definition.

The default spatial discretization is a continuous spectral element method (SEM)
on Legendre--Gauss--Lobatto nodes at arbitrary polynomial order, which is most
advantageous at fourth order and above, with a design that allows finite
difference, finite volume, or discontinuous Galerkin operators to be added
non-invasively. Time integration is delegated to the DifferentialEquations.jl
ecosystem [@rackauckas2017differentialequations], giving access to explicit and
implicit schemes without any solver code in Jexpresso itself. Distributed-memory
parallelism uses MPI.jl [@byrne2021mpi] through the Gridap.jl ecosystem
[@badia2020gridap], with meshes read via GridapGmsh.jl and adaptive refinement
and load balancing handled by p4est [@ranocha2023mpi]. Threading and GPU offload
are expressed once and dispatched to either target through KernelAbstractions.jl
and JACC.jl [@JACC].

# Statement of need

Systems of conservation laws underpin computational fluid dynamics, atmospheric
and ocean modeling, radiative transfer, and astrophysics. Spectral element
methods deliver high accuracy per degree of freedom, but their implementation
complexity has kept them out of reach of researchers who are domain experts
rather than numerical analysts. Mature frameworks such as Nek5000
[@fischer2008nek5000] and SPECFEM3D [@komatitsch2002spectral] are written in
Fortran or C++, which raises the cost of prototyping a new equation set; Julia
packages for related problems---ClimateMachine.jl [@clima2019] and TrixiAtmo.jl
[@ranocha2024trixiatmo]---target specific equation sets or use discontinuous
Galerkin rather than continuous spectral elements.

Jexpresso addresses the gap between "a framework that solves my equations" and
"a framework I can teach my equations". Adding a system requires no modification
of any solver: a self-contained problem directory supplies the initial
condition, fluxes, sources, boundary conditions, and the conserved-to-primitive
mapping, and that same directory then runs in 1D, 2D or 3D, in serial or under
MPI, on CPU or GPU, with any time integrator in the ecosystem. The repository
ships roughly one hundred such directories across nine equation families---from
advection--diffusion and Burgers to compressible Euler, elliptic and Helmholtz
problems, ideal MHD, radiative transfer, shallow water and acoustics---serving
both as regression tests and as templates.

Two capabilities are, to our knowledge, not available together elsewhere in the
Julia ecosystem. Laguerre-based semi-infinite spectral elements
[@tissaoui2024efficient] extend a bounded domain to infinity, removing spurious
reflections at artificial boundaries---essential for atmospheric and oceanic
problems. And a residual-based dynamic sub-grid-scale model [@marras2015dynsgs]
sets the artificial viscosity from the local residual of the governing equations
rather than from a tuned constant, so dissipation appears at shocks and
under-resolved features and vanishes where the solution is smooth; because it
reads only the residual it applies unchanged to whatever system the user
defined. Smagorinsky, Vreman and WALE closures, artificial viscosity and
Boyd--Vandeven spectral filtering are also available.

# Design: any system, written as it is on paper

Jexpresso is built around the general conservation law

$$
\left(\delta\,\partial_t \mathbf{q}\right)
  + \sum_{i=1}^{\mathrm{nsd}}\nabla\cdot\mathbf{F}_i(\mathbf{q})
  = \mu\nabla^2\hat{\mathbf{q}} + \mathbf{S}(\mathbf{q}) + \mathrm{b.c.},
$$

where $\nabla = [\partial_{x_i}]$, $\mathbf{q}$ is the solution vector,
$\mathbf{F}_i$ and $\mathbf{S}$ are the flux and source vectors, and
$\hat{\mathbf{q}}$ are the variables to be diffused (by default $\mathbf{q}$
itself). The time derivative is parenthesized because the system need not be
time dependent, which is how elliptic and Helmholtz problems are expressed.

**Any problem that can be written in this form is added to Jexpresso by
declaring $\mathbf{q}$, $\mathbf{F}_i$, $\mathbf{S}$, $\hat{\mathbf{q}}$ and the
boundary conditions, in the same order and with the same content a user would
write them on paper.** No solver, assembly routine, time integrator or
communication pattern is touched. Two-dimensional linear transport,
$\partial_t\rho + \partial_x(\rho u) + \partial_y(\rho v) = 0$, is the whole of
the listing below.

```julia
function user_flux!(F, G, q, args...)
    # q = [rho];  F = [u * rho];  G = [v * rho]
    u = 0.5   # x-velocity [m/s]
    v = 1.0   # y-velocity [m/s]
    F[1] = u * q[1]
    G[1] = v * q[1]
end
```

Scaling up changes only the length of $\mathbf{q}$ and the content of
$\mathbf{F}_i$ and $\mathbf{S}$. The compressible Euler equations with gravity
and two passive tracers,

$$
\partial_t\!
\begin{pmatrix}\rho\\ \rho u\\ \rho v\\ \rho\theta\\ \rho c_1\\ \rho c_2\end{pmatrix}
\!+\partial_x\!
\begin{pmatrix}\rho u\\ \rho u^2\!+\!p\\ \rho vu\\ \rho\theta u\\ \rho c_1u\\ \rho c_2u\end{pmatrix}
\!+\partial_y\!
\begin{pmatrix}\rho v\\ \rho uv\\ \rho v^2\!+\!p\\ \rho\theta v\\ \rho c_1v\\ \rho c_2v\end{pmatrix}
\!=\!
\begin{pmatrix}0\\ 0\\ -\rho g\\ 0\\ 0\\ 0\end{pmatrix},
$$

become a six-entry `user_flux!` and a `user_source!` whose only non-zero row is
$-\rho g$; adding a third tracer means adding one row. The same recipe covers
the quasi-1D Euler equations of a Laval nozzle, where the nozzle area appears in
both $\mathbf{F}$ and $\mathbf{S}$; the Helmholtz equation
$\mu(u_{xx}+u_{yy}) = \alpha^2u + f(x,y)$, which is the general form with
$\delta = 0$, a Laplacian on the left and the rest coded as a source; ideal
magnetohydrodynamics with hyperbolic divergence cleaning, a nine-field system;
and the five-dimensional radiative transfer equation, whose solution depends on
position $\mathbf{x}$ and direction $\mathbf{s}$ and whose scattering integral
enters as a source.

A problem lives entirely in its own directory: initial condition, fluxes,
sources, boundary conditions, the conserved-to-primitive mapping, and a
parameter file. Nothing there knows the polynomial order, the element type, the
rank count, or whether the run is on CPU or GPU, and nothing outside it is
edited to add new physics---the separation of concerns familiar from `OpenFOAM`,
except that Jexpresso asks for flux and source *components* rather than the
strong form.

# Key features

- **Generality:** a new system is added as fluxes, sources and boundary
  conditions only; the same definition runs in 1D, 2D and 3D.
- **High-order accuracy:** continuous spectral elements at arbitrary order, and
  Laguerre semi-infinite elements for unbounded domains.
- **Stabilization:** residual-based dynamic SGS, Smagorinsky, Vreman, WALE,
  artificial viscosity, and spectral filters.
- **HPC:** MPI domain decomposition, p4est adaptive refinement and load
  balancing, and one set of kernels for threaded CPU and GPU execution.
- **Interoperability:** Gridap.jl meshing, Gmsh input, VTK and PNG output, and a
  continuous-integration suite of analytical and benchmark tests.

# Problems shipped with the package

The repository ships roughly one hundred ready-to-run problem directories, which
double as regression tests and as templates. The most substantial are:

- **Atmospheric dynamics with clouds:** shallow cumulus and stratocumulus LES
  (BOMEX, DYCOMS, RICO-type), 2D/3D squall lines with bulk microphysics, moist
  and dry Rayleigh--Bénard convection, rising thermal bubbles.
- **Mountain and gravity waves:** hydrostatic, non-hydrostatic and Schär
  mountains, with Laguerre semi-infinite elements absorbing the vertically
  propagating waves in place of a sponge layer.
- **Turbulent boundary layers:** wall-modelled and resolved LES with
  Monin--Obukhov surface treatment, and urban-scale flow over buildings.
- **High-speed compressible flow:** Sod shock tube, quasi-1D Laval nozzle,
  isentropic vortex.
- **Magnetohydrodynamics:** ideal GLM-MHD with divergence cleaning, on the
  magnetized Kelvin--Helmholtz instability and the Orszag--Tang vortex.
- **Radiative transfer:** 2D/3D transport with scattering, coupled to a
  dynamical core for atmospheric radiation.
- **Others:** shallow water with wetting/drying, elliptic and Helmholtz
  problems, advection--diffusion, Burgers, acoustics.

These have supported published work on unbounded-domain spectral elements
[@tissaoui2024efficient], adaptive refinement for atmospheric radiation
[@tissaoui2025parcfd], and the solver itself [@marras2025jexpresso].

# Acknowledgements

We acknowledge the continuous support and input of Samuel Stechmann (University
of Wisconsin, Madison) on the extension of Jexpresso beyond fluid dynamics, and
discussions with Hendrik Ranocha (Johannes Gutenberg University Mainz) and the
other authors of Trixi.jl.

# References
