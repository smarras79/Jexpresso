# 2D ideal GLM-MHD: the Orszag–Tang vortex

This case runs the two-dimensional **Orszag–Tang vortex**, the standard MHD
benchmark, with the setup of

> A. Bormanis, C. A. Leon, A. Scheinker,
> *Solving the Orszag–Tang vortex magnetohydrodynamics problem with
> physics-constrained convolutional neural networks*,
> Phys. Plasmas **31**, 012101 (2024) — **Section I A, Eqs. (1)–(9)** and
> **Section III**. https://doi.org/10.1063/5.0172075

on the doubly periodic unit square $[0,1]^2$, $\gamma = 5/3$, over
$t \in [0, 1]$, at the paper's $128\times128$ resolution.

A full write-up of the governing equations, the GLM divergence cleaning, the
Smagorinsky stabilization, the initial condition and what to expect from the
run is in [EQUATIONS.md](EQUATIONS.md); this file covers the implementation
choices.

## Run

```bash
julia --project=. src/Jexpresso.jl MHD orszagTangBormanis2024
```

or, interactively,

```julia
julia> using Jexpresso
julia> Jexpresso.run_case("MHD", "orszagTangBormanis2024")
```

Output is VTK in `./output/`, written every 0.05 time units:
`ρ, u, v, w, p, Bx, By, Bz, ψ, T`.

## Initial condition (paper Eqs. 6–9)

With $\gamma = 5/3$ and $B_0 = 1/\sqrt{4\pi}$:

| field | value |
|---|---|
| $\rho$ | $\gamma^2/(4\pi) = 25/(36\pi) \approx 0.2210485$ |
| $p$ | $\gamma/(4\pi) = 5/(12\pi) \approx 0.1326291$ |
| $(u, v)$ | $(-\sin 2\pi y,\ \sin 2\pi x)$ |
| $(B_x, B_y)$ | $(-B_0 \sin 2\pi y,\ B_0 \sin 4\pi x)$ |
| $w,\ B_z,\ \psi$ | $0$ |

The magnetic field is Eq. (9),
$\mathbf{B} = \frac{1}{2\pi\sqrt{4\pi}}\nabla_\perp\!\left(\frac{\cos 4\pi x}{2} + \cos 2\pi y\right)$
with $\nabla_\perp = (\partial_y, -\partial_x)$, evaluated in closed form. Since
$B_x = B_x(y)$ and $B_y = B_y(x)$, the initial field is exactly divergence
free pointwise. Derived quantities: sound speed exactly $1$,
$\max\lVert\mathbf{v}\rVert = \sqrt{2}$, minimum plasma $\beta = 5/3$, and the
GLM cleaning speed $c_h \approx 2.6031$ (printed by `initialize.jl` at
startup).

## What is (and is not) implemented

- **Conservative GLM-MHD fluxes**: 9 unknowns
  `(ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)`. The 2D system keeps the third
  components of velocity and magnetic field and drops the z-flux. `ψ` is the
  GLM divergence-cleaning field; the cleaning speed `c_h` is set to the
  maximum wave speed of the initial condition and kept constant (computed in
  `initialize.jl`, stored in the `c_h_mhd` Ref). This is the **same equation
  set as `problems/MHD/kelvinHelmholtzChan2022`**, reused unchanged — only
  the initial condition, the domain and the run parameters are new.
- **Divergence control by GLM cleaning, not constrained transport.** The
  paper enforces $\nabla\cdot\mathbf{B}=0$ to machine precision with
  constrained transport on a staggered finite-volume mesh (its Sec. III A),
  which has no analogue in a collocated continuous-Galerkin SEM. GLM cleaning
  instead keeps the divergence small and bounded. `user_source.jl` adds
  Dedner's ψ-damping source `S_ψ = -(c_h/c_r) ψ` (`c_r = 0.18`, tunable via
  the `glm_cr_mhd` Ref); without it the undamped ψ waves have no dissipation
  mechanism in a CG discretization (no interface Riemann fluxes) and
  accumulate grid-scale noise on a periodic domain.
- **Stabilization: Smagorinsky SGS at 8× strength** (`:visc_model => SMAG()`,
  `:μ = [0, 8, 8, …]`), not the entropy-stable / KEP two-point fluxes
  (`:lkep => false`, `:entropy_variables => false`). The paper solves *ideal*
  (inviscid, non-resistive) MHD and gets its dissipation from a finite-volume
  Riemann solver; a collocated CG scheme has none, and the Orszag–Tang shocks
  need explicit regularization. Momentum gets the full deviatoric stress, the
  energy gets `κ∇T` plus the `τ·u` viscous work, and `ρw`, `B`, `ψ` are
  diffused as scalars (a turbulent resistivity for `B`). The 8× multiplier is
  not a tuned LES constant — it is the dissipation this discretization needs
  to get through the shocks; see the next section.
- **NOT implemented**: the non-conservative Powell / Galilean-GLM term Υ. It
  is only required for entropy stability of split-form discretizations, which
  are not used here.

Because Smagorinsky is an eddy-viscosity closure and not a shock-capturing
scheme, expect the large-scale morphology (the X-shaped density pattern, the
central magnetic island, overall `ρ` and `B` structure) to match the paper's
Figures 3 and 5–8 well, while peak values across shocks are smeared and the
finest late-time structures are damped. See §4 of
[EQUATIONS.md](EQUATIONS.md) for the full list of expected departures and the
knobs to tune.

## Stabilization: what was actually tested

The Orszag–Tang vortex steepens into shocks at $t \approx 0.5$, and that is
where an under-dissipated run dies — the gas pressure is squeezed toward zero
at a shock until `p = (γ-1)(ρE - ½ρ|v|² - ½|B|²)` goes negative and the
sound speed takes a square root of a negative number. Every configuration
below was run to $t = 1$ on the shipped mesh; "min p" is the minimum gas
pressure over the domain at $t = 0.55$, the snapshot where the failures
occur (initial pressure is $0.1326$, so it is also roughly the fraction of
the initial value that survives at the strongest shock).

| `:μ` | filter | result | min p @ t=0.55 | max\|Bx\| @ t=1 | max\|ψ\| |
|---|---|---|---|---|---|
| 1 | — | **aborts** at t ≈ 0.55 | — | — | — |
| 2 | — | **aborts** at t ≈ 0.55 | — | — | — |
| 4 | — | completes | 3.7e-2 | 0.460 | 1.2e-2 |
| **8** | **—** | **completes (shipped)** | **4.3e-2** | **0.448** | **3.2e-3** |
| 4 | erf 0.05 | completes | 3.8e-2 | 0.500 | 5.9e-3 |
| 1 | erf 0.10 | completes, but barely | 5.1e-4 | 0.618 | 1.3e-2 |
| 1 | erf 0.20 | **aborts** at t ≈ 0.55 | — | — | — |

Three things worth knowing before changing these:

1. **Plain Smagorinsky at `:μ = 1` does not survive this problem.** On this
   grid $\mu_t = \rho C_s^2 \Delta^2 |S|$ is $O(10^{-6}$–$10^{-4})$, which is
   essentially nothing against a shock. The 8× multiplier is an effective
   $C_s \approx \sqrt{8}\,(0.23) \approx 0.65$, far above the usual
   $0.1$–$0.23$. That is the honest cost of regularizing shocks with an
   eddy-viscosity closure in a scheme that has no Riemann solver.
2. **More filter is not monotonically safer.** The Boyd–Vandeven filter acts
   on the *conservative* variables independently, so filtering $\rho$,
   $\rho\mathbf{v}$ and $E$ separately can make $E - \frac12\rho|v|^2 -
   \frac12|B|^2$ negative where the unfiltered state was fine — which is
   exactly why `erf 0.2` fails while `erf 0.1` survives. Scaling `:μ` is
   better targeted: it adds dissipation in proportion to $|S|$, i.e. at the
   shocks and nowhere else.
3. **Raising `:μ` also improves divergence control**, because the same eddy
   viscosity acts on **B** as a turbulent resistivity: $\max|\psi|$ falls
   from $1.3\times10^{-2}$ (filter route) to $3.2\times10^{-3}$ at `:μ = 8`.

`:μ = 4` is the least dissipation measured to be stable and keeps noticeably
more magnetic-field structure; it is only 2× above the observed failure
point, so it is offered as the sharper alternative rather than the default.

## Notes

- **Variable ordering**: the MHD literature writes the state as
  `(ρ, ρv, E, B, ψ)`; here the total energy sits in **slot 4** and `ρw` in
  slot 5 because Jexpresso's shared 2D kernels (viscous τ·u augmentation,
  sound-speed/CFL diagnostic) assume slot 4 of a 2D system is the energy.
- **γ**: `5/3` (the paper's value, its Sec. I A) is set in `user_flux.jl`
  (`γ_mhd`). The initial condition is written *in terms of* γ, so changing it
  changes ρ and p too.
- **`w` and `Bz` must stay at machine zero** for this problem — nothing in
  the 2D system couples them once they vanish. They are a free consistency
  check on the output.
- **Resolution**: `:nop => 4` on the 32×32 element mesh gives 32·4 = 128
  unique points per direction, i.e. exactly the 128×128 grid of the paper's
  reference data.
- **Time step**: `:Δt => 5e-4` is CFL ≈ 0.24 against `c_h ≈ 2.603` and the
  smallest LGL spacing ≈ 5.4e-3. The margin is deliberate — the vortex
  steepens into shocks and local wave speeds grow. The paper's finite-volume
  reference used Δt = 8e-4 on the same grid.
- **Ignore the CFL printout.** The shared diagnostic in
  `src/kernel/physics/soundSpeed.jl` is written for CompEuler and is not
  meaningful for this case:
  - the *acoustic* CFL uses the Euler sound speed (γ = 1.4, no magnetic
    pressure) and the nominal spacing `Δs = h/(N+1)`, not the MHD fast speed
    and the smallest LGL spacing — so it reads ≈ 0.05 where the real MHD
    number is ≈ 0.24;
  - the *viscous* CFL uses `maximum(inputs[:μ])` — the per-equation
    **multiplier**, here 1.0 — in place of the actual eddy viscosity, so it
    prints ≈ 8.2. The true Smagorinsky $\mu_t \sim 10^{-6}$ here, i.e. a
    viscous CFL of order $10^{-5}$. Nothing is unstable.

  Both quirks are pre-existing and shared with the MHD KHI case.

## Mesh

The case ships its own mesh, `OT_32x32_periodic.msh` — a doubly periodic
`[0,1]²` grid of 32×32 quads — so that it runs out of the box. (Most
Jexpresso cases instead read from the separate `JexpressoMeshes` repository;
no periodic *unit-square* grid exists there, and the Orszag–Tang initial
condition has period 1 in both directions, so a `[-1,1]²` mesh would tile the
vortex four times.) `initialize.jl` warns if the mesh it is handed does not
span `[0,1]²`.

Regenerate it from the shipped `.geo` with

```bash
gmsh -2 problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.geo \
     -o problems/MHD/orszagTangBormanis2024/OT_32x32_periodic.msh
```

To move it into the standard mesh tree instead, write it to
`meshes/gmsh_grids/OT_32x32_periodic.msh` and switch to the commented-out
`:gmsh_filename` line in `user_inputs.jl`.

## Files

| file | contents |
|---|---|
| `user_inputs.jl` | solver parameters (Δt, tend, mesh, SGS, output) |
| `initialize.jl` | the Orszag–Tang initial condition and `c_h` |
| `user_flux.jl` | GLM-MHD fluxes `F`, `G`, `γ_mhd`, `pressure_mhd` |
| `user_source.jl` | Dedner ψ-damping source |
| `user_bc.jl` | periodic (mesh-level); free-slip fallback |
| `user_primitives.jl` | conserved → primitive/output mapping |
| `OT_32x32_periodic.geo` | gmsh geometry for the periodic `[0,1]²` grid |
| `OT_32x32_periodic.msh` | the generated mesh |
| `EQUATIONS.md` | equations, test description, expected results |
