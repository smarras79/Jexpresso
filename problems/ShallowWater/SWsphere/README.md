# SWsphere — shallow water equations on a spherical shell

**Status: RUNS.** Grid → manifold metrics → initial condition → time
integration. The case solves the Galewsky et al. (2004) barotropically unstable
jet with the shallow water equations of Marras, Kopera & Giraldo (2015) on a
cubed-sphere shell, by continuous-Galerkin SEM in space and SSP-RK3 in time,
with the Lagrange multiplier keeping the flow on the shell.

Two early exits remain, for working on one stage at a time:
`:lgrid_only => true` stops after the grid, `:linit_only => true` after the
initial condition.

| flag | stops after |
|---|---|
| `:lgrid_only => true` | the grid. `initialize.jl` is never called. |
| `:linit_only => true` | the grid **and** the initial condition, without integrating. |

`:lgrid_only` wins if both are set. The shipped deck has both `false` — it runs.

## Run it

The cubed-sphere grid ships with the case (`cubed_sphere.msh`: 600 quads, 602
vertices, 10 elements per panel edge), so there is nothing to generate first.

```bash
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere")
```

To regenerate the grid at a different resolution, edit `n` in the `.geo` next
to this README and run gmsh:

```bash
gmsh -2 problems/ShallowWater/SWsphere/cubed_sphere.geo \
     -o problems/ShallowWater/SWsphere/cubed_sphere.msh
```

Any gmsh output format works — MSH 2.2 or 4.1, parametric nodes or not — because
the file is read by GridapGmsh, i.e. by gmsh itself, exactly like every other
Jexpresso case.

> **Run this in a fresh Julia session after pulling.** `run_case` re-includes
> `run.jl`, `problems/drivers.jl` and the case's `user_*.jl`, but *not* the rest
> of the module, which is only evaluated when `Jexpresso` itself is loaded. A
> session opened before you pulled will run the new `drivers.jl` against the old
> module and die with `UndefVarError`.

Output, in `problems/ShallowWater/SWsphere/output/`:

| file | what it is |
|---|---|
| `sphere_grid_ho.vtu` | the initial state |
| `sphere_0001.vtu` … | one per output time, `:ndiagnostics_outputs` of them |

Each file is the high-order grid — **(ngl-1)² sub-elements per spectral
element**, exactly as `write_vtk_grid_only` does for the flat cases, so every
LGL node is a corner of a sub-cell and the linear elements are never written on
their own — carrying the solution as point data:

* `phi, phiu, phiv, phiw` — the conservative state actually integrated
* `h, u, v, w` — the primitives, **with the velocity PROJECTED ONTO THE SHELL**:
  `u` is zonal (eastward), `v` meridional (northward), `w` radial. Plotting the
  raw Cartesian `uₓ` over a sphere shows a dipole straddling the prime meridian
  — an artefact of the frame, not the flow — where the zonal component shows the
  jet as the band it actually is. `w` is the velocity form of the constraint and
  must be ~0 everywhere.
* `vorticity` — relative vorticity `ζ = n̂·(∇ₛ×u)`. **This is the field the test
  is judged on**: `h` barely moves while the instability develops, so a height
  plot makes the roll-up nearly invisible.
* `velocity` — a true VTK vector for glyphs/streamlines. Built from the
  conservative state, so it is genuinely Cartesian; `u,v,w` above are tangent
  components and glyphing *those* would point every arrow nowhere.
* `momentum_normal` = `(φu)·x̂` — the constraint; **this is how you see whether
  the flow is leaving the shell**
* `ip`, `node_type`, `lon`, `lat`, `radius`, and the cell fields `iel`, `panel`

View with representation **"Surface With Edges"**: the edges drawn are the
sub-element boundaries, so the LGL point distribution is what you see.

## What was actually added

A spherical shell is a **2D manifold embedded in 3D**, and the existing 2D gmsh
path in `src/kernel/mesh/mesh.jl` cannot represent one: it stores only `(x,y)`
per node and interpolates only `(x,y)` when it adds the high-order points, which
would collapse the sphere onto its equatorial disc. The **grid** is therefore
read by the ordinary gmsh path, which recognises a 2D manifold embedded in 3D
from the model itself and keeps `z`; everything downstream of the grid is what
stays specific to the manifold:

* `src/kernel/mesh/mesh.jl` — `mod_mesh_read_gmsh!` sets `mesh.lmanifold`,
  interpolates the LGL nodes in (x,y,z) and snaps them onto the shell
  (`project_nodes_to_shell!`). No separate reader, no separate mesh type.
* `src/kernel/mesh/sphere_metrics.jl` — the **metric terms of the manifold**
  (contravariant basis, surface Jacobian) and the diagonal mass matrix.
* `src/kernel/operators/sphere_rhs.jl` — the **SEM right-hand side**: surface
  divergence + direct stiffness summation + `M⁻¹`, the artificial diffusion
  `δν∇ₛ²(φu)` in weak form, and the modal filter.
* `src/kernel/solvers/sphere_time_loop.jl` — **SSP-RK3**, with the Lagrange
  projection applied after every stage, CFL time step, diagnostics and output.
* `initialize.jl` — the Galewsky et al. (2004) jet (see below).
* `user_flux.jl` / `user_source.jl` — the equations (see below).
* `src/io/write_output.jl` — `write_vtk_sphere_grid`, next to the existing
  `write_vtk_grid_only` writers.
* `test/test_sphere_mesh.jl`, `test/test_sphere_visc.jl` — standalone tests
  (`julia --project=. test/test_sphere_mesh.jl`). The second one checks the
  surface Laplacian against the spherical-harmonic eigenvalues `−l(l+1)/R²`.

### The node numbering, which is the delicate part

The shell is **closed**: it has no boundary, so there is no "outer edge" to
anchor a numbering to, and every edge is an interior edge shared by exactly two
elements. Get this wrong and the panel seams are numbered twice: the grid still
*looks* watertight in ParaView, but CG/SEM assembly across every seam is
silently wrong.

None of it needs special handling. Gridap builds the element/edge topology from
the node ids in the `.msh`, so a **unique edge table is what the reader already
has**, and `add_high_order_nodes_edges!` creates the LGL points once per unique
edge and hands the same ids to both neighbouring elements. Continuity across a
seam is therefore *node identity*, exactly as it is in the interior of any flat
grid — there is nothing special about a seam.

Two consequences worth knowing:

1. **Your `.msh` must be watertight.** It is the grid file, not the reader, that
   decides whether the six panels share their seam nodes. The `.geo` shipped
   here builds all six surfaces on the *same* twelve arcs and eight corner
   points, so gmsh emits one node per seam location. A file built panel by panel
   as six independent surfaces is six disconnected patches, and Jexpresso will
   read it faithfully as six disconnected patches. Fix it in gmsh (share the
   curves) rather than stitching it back together at read time — which is what
   the old `:lmerge_coincident_nodes` did, and why it is gone.

2. **The node count is the check.** On a closed shell

   ```
   npoin = npoin_linear + nedges*(ngl-2) + nelem*(ngl-2)²
   ```

   holds *only* if every unique edge contributed its high-order nodes exactly
   once. For the shipped grid at `:nop => 5`: `602 + 1200·4 + 600·16 = 15002`,
   with `V - E + F = 602 - 1200 + 600 = 2`. `test/test_sphere_mesh.jl` asserts
   both, plus that no two nodes are coincident, at several polynomial orders and
   for both `MSH 2.2` and `MSH 4.1`-with-parametric-nodes files.

The metrics are verified separately at run time (`check_sphere_metrics`, switch
`:lcheck_grid`), and each test prints `PASS`/`FAIL`:

   | test | what it proves |
   |---|---|
   | M1 | metric identities `aⁱ·a_j = δⁱⱼ` |
   | M2 | the contravariant basis is tangential to the shell |
   | M3 | the outward unit normal satisfies `n̂·x̂ = 1` |
   | M4 | `Σ M[ip] = 4πR²` — the SEM quadrature integrates the sphere exactly |
   | M5 | curvature identity `∇ₛ·(J aⁱ)/J = -(2/R) n̂` |

   M4 is the sharpest of these: it is a single number that fails immediately if
   any high-order node is misplaced or any seam is numbered twice.
   `:lstop_on_bad_grid` (default `true`) aborts the run if any test fails.


### Geometry of the high-order points

Edge points are placed by **great-circle (slerp) interpolation** at the LGL
abscissae; element interiors by a **Coons transfinite patch** built from those
four great-circle edges, then projected radially onto the shell. Both are exactly
on the sphere, and the edge construction depends only on the edge's two
endpoints, so the two sides of a seam agree.

### Changing the cube-face → sphere map

`cubed_sphere.geo` builds the **classical gnomonic** cubed sphere: the six
panels are the central projection of the faces of an inscribed cube
(Sadourny 1972). That is the simplest map and the worst conditioned — all of
its distortion piles up at the eight cube corners.

One input slides the nodes onto a different map *after* the grid is read.
Connectivity, the panel decomposition and the panel boundaries are untouched;
only where the nodes sit **inside** each panel changes, so this is safe to do
once `mod_mesh_read_gmsh!` has already built the topology.

```julia
:cubed_sphere_map => :conformal,   # :none (default) | :gnomonic | :equiangular | :conformal
```

The map the nodes come **from** is always the gnomonic one, and there is
deliberately no input for it: reading a node's face coordinate back off the
sphere *is* the gnomonic inverse, and it is what `cubed_sphere.geo` emits. A
"source map" switch would be a claim about the `.msh` that nothing can check,
and getting it wrong would silently produce a grid that is neither map.

| value | map | what it buys |
|---|---|---|
| `:gnomonic` | equidistant central projection — Sadourny (1972) | nothing; this is what the file already is, so it is a no-op |
| `:equiangular` | face coordinate measured as an angle, `u = tan α`, `α ∈ [-π/4, π/4]` — Ronchi, Iacono & Paolucci (1996) | a more even node spacing along a panel EDGE than the gnomonic one |
| `:conformal` | Rančić, Purser & Mesinger (1996) | locally **orthogonal** everywhere except the eight cube corners, so grid lines meet at 90° right up to a corner instead of degenerating to 120° |

### Corner behaviour, which is what actually distinguishes them

Angle between the two families of grid lines, walking the face diagonal
`u = v = t` in to a cube corner. 90° is orthogonal:

| `t` | gnomonic | equiangular | conformal |
|---|---|---|---|
| 0.0 (face centre) | 90.0° | 90.0° | 90.0° |
| 0.9 | 116.6° | 114.9° | **90.000°** |
| 0.99 | 119.7° | 119.5° | **90.000°** |
| 0.9999 | 120.0° | 120.0° | **90.000°** |

Gnomonic and equiangular both collapse to a 120° corner — the three panels
meeting there each contribute 120°, and neither map fights it. Only the
conformal map holds 90° all the way in. Over the whole panel the worst
deviation from orthogonality is 29.3° (gnomonic), 29.0° (equiangular),
**0.000°** (conformal).

### What it costs, measured on the shipped grid

Element edge lengths over all 600 quads. The minimum is what sets the explicit
time step under `:lcfl_dt`:

| map | min edge | max edge | ratio | Δt vs gnomonic |
|---|---|---|---|---|
| gnomonic | 710.2 km | 999.7 km | 1.41 | 1.00× |
| equiangular | 561.7 km | 1365.7 km | 2.43 | 0.79× |
| conformal | **721.6 km** | 1322.3 km | 1.83 | **1.02×** |

So on THIS grid the conformal map is both the smoothest at the corners and
marginally the cheapest — `:equiangular` is the one that costs 21% of the time
step. Do not generalise the Δt column to finer grids: the conformal map's cell
scale falls off toward a corner (`|r_u|` is 0.55 of its face-centre value at
`t = 0.9`, 0.25 at `t = 0.99`), so the more elements per panel edge, the closer
the corner-adjacent nodes sit to the singular point and the more of that
shrinkage a refined grid will see. Re-measure if you refine.

The mechanics — the panel frames, the Rančić Table B1 coefficients and the
series reversion that gives the inverse — are in
`src/kernel/mesh/cubed_sphere_maps.jl`, and `test/test_cubed_sphere_maps.jl`
checks them (including that the conformal map really is conformal, and that
the two panels sharing a cube edge move a seam node to the same place, so the
shell stays watertight).

`:cubed_sphere_map` is part of the mesh-cache fingerprint, so changing it
re-reads and re-remaps the grid instead of reusing the previous run's node
positions out of `.jexpresso_cache/`.

### Inspecting the numbering in ParaView

Colour `sphere_grid_ho.vtu` by the point field **`ip`** (the global node index).
Across a panel seam the field must look *continuous* — a jump there means the two
sides carry different node numbers, i.e. the shell was torn open. `node_type`
(0 vertex / 1 edge / 2 interior), `lon`, `lat` and `radius` are written too, and
the cell field `panel` shows the gmsh physical tag of each panel.

## The initial condition

**Galewsky, Scott & Polvani (2004)**, *An initial-value problem for testing
numerical models of the global shallow-water equations*, Tellus **56A**, 429-440
— the barotropically unstable mid-latitude jet.

The jet is naturally written in the tangent basis,

```
u(φ) zonal (eastward, +λ),   v = 0
```

and then pushed to the conservative Cartesian state `q = [φ, φu, φv, φw]` that
the equations use, through `e_λ = (-sin λ, cos λ, 0)` and
`e_φ = (-sin φ cos λ, -sin φ sin λ, cos φ)`. Being a combination of `e_λ` and
`e_φ`, the result is tangential to the shell **by construction** —
`initialize.jl` asserts `max|(φu)·x̂|` is at round-off before returning, and
refuses to start otherwise.

The jet is confined to `φ ∈ [π/7, π/2 - π/7]` and is `C^∞`: every derivative
vanishes at both edges, so it joins the motionless fluid smoothly.

```
u(φ) = (u_max/e_n) exp[ 1/((φ-φ₀)(φ-φ₁)) ]   in the band, 0 outside
v    = 0
h(φ) = h₀ - (1/g) ∫ a u [ f + tan(φ) u/a ] dφ',   f = 2Ω sin φ
```

`h₀` is **not hard-coded**. It is fixed by requiring the global mean depth to be
10 km, and since `h` is zonally symmetric that mean is a 1-D integral whose order
of integration can be swapped, collapsing the double integral to a single
quadrature:

```
mean(h) = ½ ∫ h cos φ dφ = h₀ - (1/2g) ∫ a u [f + tan(φ)u/a] (1 - sin φ) dφ
```

With the published parameters this returns **10158.1861704546 m**, the constant
quoted in the literature — which is the sharpest available check that the
parameters and the quadrature agree with the paper, and `test/test_sphere_mesh.jl`
asserts it to 1e-6 m. The balance integral itself is verified there too, by
finite-differencing `h(φ)` and comparing against its own integrand: that is what
certifies `h` and `u` are a *balanced pair* rather than two independently
plausible profiles.

The instability is seeded by

```
h'(λ,φ) = ĥ cos φ exp[-(λ/α)²] exp[-((φ₂-φ)/β)²],  ĥ=120 m, α=1/3, β=1/15, φ₂=π/4
```

Set `:lgalewsky_perturbation => false` for the **balanced control run** — an
exact steady solution, and the natural first test of the scheme once it exists.
The unperturbed state is stored in `q.qe` either way, so the drift can be
measured directly.

Every parameter (`:galewsky_umax`, `:galewsky_phi0`, `:galewsky_hhat`, …) is
overridable from `user_inputs.jl`; `initialize.jl` falls back to the published
values. `:sphere_radius => 6.37122e6` in the deck pins the shell to the Galewsky
Earth radius whatever the `.msh` carries — the balance is computed on the radius
the grid actually has, and `initialize.jl` warns if that is not the Earth.

## The equations

**Marras, Kopera & Giraldo (2015)**, *Simulation of shallow-water jets with a
unified element-based continuous/discontinuous Galerkin model with grid
flexibility on the sphere*, QJRMS **141**: 1727-1739, Eq. (8):

```
∂φ/∂t    + ∇·(φu)   = 0
∂(φu)/∂t + ∇·(φu⊗u) = -φ∇φ - f(x × φu) + μx + δν∇²(φu)
```

`φ = gh` is the geopotential, `u = (u,v,w)` the **full Cartesian** velocity,
`x` the position vector, `f = 2ωz/r²`, and `μ` the Lagrange multiplier.

State vector: `q = [φ, φu, φv, φw]` — four equations. Not the tangent-plane pair
`(u,v)`: the multiplier exists precisely to remove the velocity component
**normal** to the shell, and in a two-component tangent basis that component
does not exist, so there would be nothing to constrain. The price is one
redundant degree of freedom, which is what `μ` and the projection deal with.

### Where each term lives

| term | where | why |
|---|---|---|
| `∇·(φu)`, `∇·(φu⊗u)` | `user_flux.jl` | under a divergence |
| `φ∇φ` | `user_flux.jl` | as `∇·((φ²/2)I)` — see below |
| `-f(x × φu)` | `user_source.jl` | pointwise |
| `μx` | `user_source.jl` | pointwise, closed form — see below |
| `δν∇²(φu)` | `:lvisc`, `:μ` | a second derivative, Jexpresso's viscous path |

The paper carries the pressure as the **source** `-φ∇φ`. Jexpresso's
`user_source!` is a pointwise callback — one node, no derivatives — so that form
cannot be evaluated there. The identity `∇·((φ²/2)I) = φ∇φ` moves it into the
flux instead, giving the standard conservative momentum flux tensor
`T_ij = φ u_i u_j + (φ²/2)δ_ij`. Algebraically the same equation, and it fits
Jexpresso's split exactly. On the manifold it stays correct because the RHS
differentiates the **Cartesian components with surface derivatives**:
`∂ⱼ(δᵢⱼφ²/2) = ∂ᵢ(φ²/2)` is then the surface gradient, tangential by
construction, so no spurious normal pressure force appears.

### The Lagrange multiplier — two mechanisms, both needed

The paper (section 3.2, after Coté 1988) presents `μ` as a **projection applied
to the state** after each step, Eqs (9)-(11):

```
(φu)ⁿ⁺¹_c = P (φu)ⁿ⁺¹_u ,   P = I - x xᵀ/r² ,   μ = -(φu)ⁿ⁺¹_u·x / r²
```

That is `project_momentum_to_sphere!` in `src/kernel/mesh/sphere_metrics.jl`.

But `μ` also has a **closed pointwise form**, which is what lets it live in a
pointwise source at all. Requiring `x·∂(φu)/∂t = 0` and taking `x·` of the
momentum equation kills the Coriolis term (`x·(x × φu) = 0`) and the pressure
term (a surface gradient is tangential), and leaves only the curvature part of
the flux divergence, `x·[∇ₛ·(φu⊗u)] = -φ|u|²`. Hence

```
μ = -φ|u|²/r²        μx = -(φ|u|²/r) n̂
```

— a **centripetal** force of magnitude `φ|u|²/r`, exactly what is needed to bend
a parcel of "mass" `φ` moving at speed `|u|` around a circle of radius `r`. That
physical reading is the quickest check on the sign: `μ` must be negative.

The two are complementary, not redundant. `μx` keeps the **continuous** equation
consistent, so the exact solution never leaves the shell and the discretization
is not fighting a wrong PDE. The projection removes the normal drift the
**discrete** operators accumulate anyway (`μ` above is exact only if the discrete
divergence is). `:llagrange_projection` controls the second; the first has no
switch because `user_source!` receives no `inputs` — comment the two `μ*` terms
out to run without it.

The VTK file carries `momentum_normal` = `(φu)·x̂`, the very quantity being held
at zero: plot it and you can see whether the flow is leaving the shell.

### How this was checked

The unperturbed Galewsky jet is an **exact steady solution**, so the whole
right-hand side must vanish on it — which only happens if the flux, the pressure
term, the Coriolis sign *and* the multiplier are all right at once.
`test/test_sphere_mesh.jl` evaluates `-∇ₛ·F(q) + S(q)` on the analytic balanced
state using the real `user_flux!` and `user_source!` with a finite-difference
surface divergence, and checks that the residual **converges to zero at second
order** under refinement (it cannot be exactly zero — it carries the stencil's
truncation error). Measured rate: 2.00, residual ~1e-4 of the size of the
individual terms that cancel. Each term is also checked separately against its
closed form, and the projection against its defining properties (kills the
normal component, leaves a tangential field untouched, idempotent).

## Space and time discretization

**Metrics** (`sphere_metrics.jl`). On a flat grid the map (ξ,η) → (x,y) is
square and `metric_terms.jl` stores its inverse. On a shell the map is
(ξ,η) → (x,y,z): a 3×2 Jacobian with no inverse. What replaces it is the
contravariant basis of the surface,

```
a_ξ = ∂x/∂ξ,  a_η = ∂x/∂η,  n = a_ξ × a_η,  J = |n|
a¹  = (a_η × n)/J²,          a² = (n × a_ξ)/J²
```

so `∇ₛf = a¹ ∂f/∂ξ + a² ∂f/∂η`, i.e. `a¹ = (dξdx, dξdy, dξdz)` — Jexpresso's
own naming with a third component. Both are tangential, which is why the
conservative pressure term produces no spurious normal force. The mass matrix is
diagonal (LGL/inexact integration) and assembled by direct stiffness summation.

**RHS** (`sphere_rhs.jl`). Continuous Galerkin, strong form, exactly the
convention of `_expansion_inviscid!` with the third flux component added:

```
rhs_el -= ω_i ω_j J [ (∂F/∂x + ∂G/∂y + ∂H/∂z) − S ]      then DSS, then M⁻¹
```

No boundary term: the shell is closed, so `∫_Γ` vanishes identically (the paper
notes this for the CG case). **All of it is SEM** — the differentiation is the
LGL differentiation matrix contracted with the metric terms; there are no finite
differences anywhere in the solver.

**Stabilization** (`sphere_rhs.jl`). Two independent mechanisms, and the run
needs **at least one** of them — the inviscid unfiltered high-order solution
grows grid-scale modes and blows up, which is exactly what the paper reports
for the cubed sphere. The time loop prints which are active at startup and
warns if both are off.

| deck | what it does |
|---|---|
| `:lfilter => true` | the exponential modal filter, applied once per step as the integrator's `step_limiter!`, followed by a mass-weighted DSS average, so it conserves `∫φ` exactly. The stand-in for the paper's Boyd-Vandeven filter; `:filter_alpha => 0.05` is its "reduce the highest modes by 5%". |
| `:lvisc => true`, `:μ => 1e5` | the artificial diffusion `δν∇²(φu)` of Eq. (8b), `ν = 1e5 m²/s` in the paper. `:ivisc_equations => [2,3,4]` puts it on the momentum only, where the paper has it. |

This deck keeps the **filter**, which is the paper's own choice for the
published test. The other branch is a case of its own —
[`../SWsphere_visc`](../SWsphere_visc/README.md), identical but for those two
switches — so both are covered by something runnable rather than by a comment.
Running with both on is legitimate too; it is simply more dissipative.

The diffusion is the **surface (Laplace-Beltrami) Laplacian**, built from the
same 3×2 manifold metrics as the divergence, and assembled in **weak** form —
integration by parts leaves no boundary term on a closed shell:

```
∫ ψ ν∇ₛ²q dΩ = −∫ ∇ₛψ·(ν ∇ₛq) dΩ ,     ∇ₛq = a¹ ∂q/∂ξ + a² ∂q/∂η
```

which needs only one derivative of the C⁰ solution and is symmetric negative
semi-definite by construction, so the term can only ever remove energy. This is
why it cannot go through the viscous path in `src/kernel/operators/rhs.jl`: that
one builds `∇q` from the flat 2×2 inverse Jacobian, which has no third metric
direction. Dropping the third component is not an approximation but a different
operator — it differentiates along the projection of the shell onto `z = 0`,
singular at the equator.

Diffusion is a second derivative, so it is explicit-stable only for
`Δt ~ Δ²/ν`; `sphere_cfl_dt` takes the min of that and the wave-speed limit. At
the shipped resolution the two are ~1e5 s and ~1e2 s, so `ν = 1e5` costs no time
step at all — the wave speed still sets `Δt`.

**Time** (`sphere_time_loop.jl`). SSP-RK3 (Shu-Osher), the paper's section 4.2.
The Lagrange projection runs **after every stage**, not only at the end of the
step: each stage is a full Euler update and can push momentum off the shell, and
projecting only at the end would feed an off-manifold state into the next
stage's flux evaluation. `:Δt` is optional — omit it and the step comes from
`Δt = :cfl · Δmin / max(|u| + √φ)`.

Diagnostics printed every `:ndiagnostics_prints` steps: mass, energy,
`max|(φu)·x̂|` (the constraint), `max|ζ|`, and `max|ζ-ζ₀|`.

The last one is the instability indicator. `max|ζ|` on its own is dominated by
the jet's own shear (~1e-4 across the band from t=0) and only rises ~30% over
six days, which understates what is happening; differencing against the initial
field isolates the growing perturbation.

The deck runs the full **144 h** of the published test. The perturbation is
linearly unstable but grows slowly: at one day the jet still looks laminar, and
the vortex train only appears around day 4–6. Shortening `:tend` hides the very
thing the test exists to show.

## How this was checked

Everything below was **run**, in Julia:

| what | result |
|---|---|
| metric identities `aⁱ·a_j = δⁱⱼ` | 6.7e-16 |
| contravariant basis tangential | 2.7e-15 |
| `Σ M[ip]` vs `4πR²` (SEM quadrature of the area) | 1.7e-12 relative |
| curvature identity `∇ₛ·(J aⁱ)/J = −(2/R)n̂` | 1.2e-5 at nop=4, 1.3e-8 at nop=6, 2.4e-11 at nop=8 |
| steady-state residual of the SEM operator | 6.2e-3 → 2.5e-3 → 2.7e-4 under refinement |
| mass conservation over a run | 6e-11 |
| `max\|(φu)·x̂\|` after integrating | 4e-10 against a momentum scale of 8e6 — **5e-17 relative** |
| error on the balanced (steady) jet | 60 m → 37 m → 3.3 m as resolution rises |
| filter mass conservation | 4e-13 |
| surface Laplacian vs `−l(l+1)/R²` for `l = 1…4` | 2.5e-7 → 2.9e-6 at nop=5; 1.5e-8 → 1.3e-7 at nop=6 |
| surface Laplacian annihilates constants | 1e-12 (times `R²`) |
| surface Laplacian symmetric / negative definite | 5e-15 asymmetry; `qᵀLq < 0` |

The **curvature identity** deserves a note. On a flat grid the free-stream
condition is `∇·(J aⁱ) = 0`. On a curved surface that is false: the surface
divergence of a constant vector leaves the mean-curvature term, `∇ₛ·c =
−H(c·n̂)` with `H = 2/R`. Asserting the flat identity flags correct metrics as
broken — it reads exactly `2/R`. This is not a formality: it is the *same*
curvature that produces `x·[∇ₛ·(φu⊗u)] = −φ|u|²`, which is what the multiplier
`μ = −φ|u|²/r²` cancels. If the discrete metrics got the curvature wrong, the
multiplier would not balance the discrete flux divergence and the flow would
leave the shell.

The unperturbed jet is an exact steady solution, so the residual of the SEM
operator on it must converge to zero — and it does, spectrally. That single test
covers the flux, the pressure term, the Coriolis sign, the multiplier *and* the
metrics at once.

## What comes next

1. Fold the shell into `sem_setup`, so it uses `params_setup` and Jexpresso's
   own SciML integrators (with the projection as a stage callback) instead of
   the self-contained SSP-RK3 here.
2. MPI: the shell path is serial today (the `:lspherical_shell` branch of `problems/drivers.jl` says so).
3. The reference diagnostic of the test is relative vorticity at day 6; that
   needs `∇ₛ×u` on the manifold, which the metrics now make straightforward.
