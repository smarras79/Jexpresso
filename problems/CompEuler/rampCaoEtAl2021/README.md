# `CompEuler/rampCaoEtAl2021` — 2D Mach-7.7 compression ramp with laminar separation

```bash
julia --project=. src/Jexpresso.jl CompEuler rampCaoEtAl2021
```

The x–y flow of

> S. Cao, J. Hao, I. Klioutchnikov, H. Olivier, C.-Y. Wen, *Unsteady effects in
> a hypersonic compression ramp flow with laminar separation*, J. Fluid Mech.
> **912**, A3 (2021), [doi:10.1017/jfm.2020.1093](https://doi.org/10.1017/jfm.2020.1093),

at the paper's coarser grid, case **G1** of §2.2.

## What is and is not 2D about this case

The paper's DNS is three-dimensional, but the *problem* is not: the ramp is a
2D geometry, and the paper computes the 2D flow twice — as the field its 3D
runs are started from (§2.3, "initialised by duplicating the two-dimensional
converged solution in the spanwise direction") and as the base flow of the
global stability analysis of §3.2. That is what this deck solves. So:

**reproduced** — the laminar flat-plate boundary layer, δ = 1.38 mm at
separation; the separation bubble, 0.59 < x/L < 1.26; the leading-edge,
separation and reattachment shocks and their interaction; the surface pressure
and heat flux along the line of symmetry (figures 2a, 3c).

**not reproduced, by construction** — the streamwise heat-flux streaks, their
spanwise modulation, and the low-frequency unsteadiness of §§3–5. Those are a
*global instability of this base flow*: three-dimensional by definition
(spanwise wavenumber β ≠ 0; the most unstable mode is at λ_z/L = 0.066), so no
2D run can carry them. This case is the flow they grow on, not a cheap version
of them.

## Configuration

| | |
|---|---|
| geometry | sharp leading edge, flat plate L = 100 mm, ramp 15°, also 100 mm along the surface (§2.2). §2.3's 1 mm free-stream strip ahead of the leading edge is **not** included — see Mesh |
| free stream | M = 7.7, T = 125 K, p = 760 Pa, u = 1726 m/s, Re/m = 4.2 × 10⁶ (Table 1, shock tunnel TH2) |
| wall | no slip, isothermal, T_w = 293 K (T_w/T_∞ = 2.34) |
| gas | perfect, γ = 1.4, Pr = 0.71, Sutherland's law for μ |
| grid | `ramp15.msh`, 269 × 60 quads → 1077 × 241 LGL points at `:nop => 4` (paper G1: 1080 × 240), Δy_wall = 7.98 × 10⁻⁶ m (paper 8 × 10⁻⁶) |
| time step | Δt = 1 × 10⁻⁹ s, set by the **viscous** limit from the DynSGS cap, not by advection — derived in `user_inputs.jl` |

The Reynolds number is not an independent input: μ(125 K) = 8.656 × 10⁻⁶ Pa s
from the standard air Sutherland constants, with the paper's own p, T, u and L,
gives Re_L = 4.22 × 10⁵ against its 4.2 × 10⁵ — so those constants are the ones
the paper used, and changing `:sutherland_muref` changes the case.

## Numerics

Total energy (not ρθ — θ is not conserved across a shock), DynSGS residual-based
shock capturing, and the molecular viscosity on top of it
(`:lsutherland => true`, [DSGS.md §4.11](../../../DSGS.md)). The sensor is what
holds the shocks; Sutherland's law is what makes the boundary layer exist — a
DynSGS-only run of this geometry is inviscid flow over a ramp, with no bubble.
The 2D viscous assembly the two feed is the real Navier–Stokes operator:
deviatoric stress tensor, viscous work, Fourier conduction.

`:Δt => 5.0e-9` is CFL ≈ 0.3 on the wall-normal and leading-edge spacings; see
the derivation in `user_inputs.jl`. `:tend => 2.0e-3` is t·u_∞/L = 34.5, about
17 flow-throughs — 400,000 steps. This is a DNS grid; run it on several ranks.

## Mesh

`ramp15.geo` is the record of the grid; `generate_mesh.py` writes `ramp15.msh`
directly from the same numbers, so no gmsh installation is needed:

```bash
python3 generate_mesh.py            # -> ramp15.msh
gmsh -2 ramp15.geo -o ramp15.msh    # identical, if you have gmsh
```

Two conforming transfinite blocks (plate, ramp). The upper boundary is the wall contour shifted *vertically* by
H = 60 mm, so every streamwise grid line is vertical and the ramp block is a
uniform shear of a rectangle — no metric distortion. H puts the whole shock
system inside the domain: the separation shock crosses the outflow plane some
36 mm below the top, so the free-stream Dirichlet condition there is never
asked to swallow a discontinuity.

**No upstream strip.** §2.3 places 20 grid points in 1 mm ahead of the leading
edge; this grid does not, and spends all 269 streamwise elements on the body
(137 plate + 132 ramp). The strip was the single worst cell on the grid:
conforming blocks force it to carry the wall-clustered Δy = 8e-6 m, while its
lower boundary is a symmetry line, so it has no boundary layer and runs the full
1726 m/s at the finest wall-normal spacing in the domain. When DynSGS saturates
there, the viscous CFL it implies is ~3e-10 s — and measured, the case died at
step 304 at exactly that node. The cost is that the leading-edge singularity now
sits on the inflow plane, which is what the strip existed to avoid; `user_bc.jl`
gives that node to the wall.

### The unstretched diagnostic grid

`ramp15_uniform.msh` is the same geometry with **no wall-normal stretching** —
401 uniform elements over H, Δy_wall = 2.58e-5 m, **y⁺ = 1.00** (G1 is 0.3).
Point `:gmsh_filename` at it to run the case without the stretching.

```bash
python3 generate_mesh.py --py 1.0 --ny 401 -o ramp15_uniform.msh
```

It is a **diagnostic grid, not a substitute for G1**: at y⁺ = 1 the wall heat
flux is under-resolved relative to the paper, so Stanton-number comparisons
against figure 2(a) are not meaningful on it. What it is for is isolating the
stretching, which a five-rung ladder from `ffs_step` identified as the cause of
the blow-up at t = 2.53e-7 s:

| test | result |
|---|---|
| `ffs_step` as shipped | passes |
| ramp mesh + M 7.7, free-slip wall | fails |
| ramp mesh + M 7.7, no-slip adiabatic | fails |
| ramp geometry + M 7.7, **uniform** y-mesh | **passes** |
| `ffs_step` isotropic mesh, **M 7.7** conditions | **passes** |

Conditions, geometry and the wall treatment are all cleared; Δy_wall = 7.98e-6 m
is what is left. The failure is also independent of Δt (identical failure *time*
across a 5× range) and of MPI rank count (32 ranks and serial agree to 17
digits), so it is neither a CFL violation nor a decomposition artifact.

It costs 269 × 401 = 107,869 elements — 6.7× the production grid, because a
uniform mesh fine enough at the wall must carry that spacing to H = 60 mm.
`--ny 201` gives y⁺ ≈ 2 at half the cells if that is enough for a given test.

To refine, either raise `NX_PLATE` / `NX_RAMP` / `NY` in `generate_mesh.py`
(`NY = 80`, `NX_PLATE = 203`, `NX_RAMP = 197` is roughly case G2, 1600 × 320),
or set `:linitial_refine => true` with `:init_refine_lvl => 1`. Cut `:Δt` with
it — but **not** by the square of the element size: Δt here is set by the DynSGS
cap, whose ν_cap scales with Δelem while the stability limit scales with Δy², so
refinement costs one factor of element size in Δt, not two.

## Running it: cores and wall time

`ffs_step` is the same solver configuration at almost exactly the same size —
4032 elements at `:init_refine_lvl => 1` = 16,128 elements, `:nop => 4`, 4
equations, DynSGS with the legacy sensor — and its deck records 64,000 steps at
about 12 h on one core, i.e. **0.68 s/step**. This case is 16,140 elements with
the same everything, so the nominal run at Δt = 1e-9 (2,000,000 steps) is
**≈ 380 core-hours**, with maybe a factor 2 of slack in that "order 12 h".

| cores | elements/rank | wall time (~80% efficiency) |
|---|---|---|
| 16 | 1009 | ~30 h |
| **32** | **504** | **~15 h** |
| 64 | 252 | ~8 h |

**32 is the recommendation** — one node, and memory is trivial (a few hundred MB
in total). Don't go much past 64: ~250 elements/rank is where a 2D `nop = 4` SEM
starts being latency-bound, and this run makes 10 million RHS calls, each with a
halo exchange and — with `:dsgs_norms => "domain"`, now the default — two or
three `Allreduce`. Those collectives are only a minute or two of the run at
32–64 ranks, but they grow with rank count while the compute per rank shrinks.

Do **not** switch to `:dsgs_norms => "rank"` to dodge them. `ffs_step` records
what that costs: rank-local norms make the viscosity depend on the partition,
measured at 18× worse error on `smoothVortex` going from 2 to 8 ranks.

One refinement level (`:init_refine_lvl => 1`) is 4× the elements and 2× the
steps (see Mesh for why it is 2× and not 4×), so 8× the cost: ~3000 core-hours,
i.e. 128 cores (504 elements/rank again) for about a day.

## What to check the result against

| quantity | paper | where |
|---|---|---|
| separation | x/L = 0.59 | §3.1 |
| reattachment | x/L = 1.26 | §3.1 |
| bubble length | 0.67 L | §3.1 |
| C_p and St along the wall | figure 2(a), against the experiment of Roghelia et al. (2017b) | §2.5 |
| boundary layer profiles upstream of separation | figure 2(b), against the compressible similarity solution | §2.5 |

`:lschlieren => true` writes the numerical schlieren that makes all three shocks
and the shear layer visible at once — the direct counterpart of the
experimental schlieren the paper is validated against. Colour the `schlieren`
field with a **reversed** greyscale in ParaView.

## Not in CI

`test/generate_ci_ref.jl` has to *run* the case to store a reference, and at
400,000 steps this one is not a CI case as written. If you want it guarded,
copy the deck with `:tend` cut to a few hundred steps first — long enough for
the leading-edge shock to form, which is what would catch a regression in the
DynSGS or Sutherland paths.
