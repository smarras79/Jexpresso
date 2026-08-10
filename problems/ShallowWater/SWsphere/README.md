# SWsphere — shallow water equations on a spherical shell

**Status: GRID + INITIAL CONDITION.** This case builds the high-order
spectral-element grid on a closed spherical shell, verifies it, builds the
Galewsky et al. (2004) barotropically unstable jet on it, writes both to VTK, and
stops. There is no time integration yet — the equations and the metric terms on
the manifold are still missing.

Two stopping points, both in `user_inputs.jl`:

| flag | stops after |
|---|---|
| `:lgrid_only => true` | the grid. `initialize.jl` is never called. |
| `:linit_only => true` | the grid **and** the initial condition. This is what the shipped deck does. |

`:lgrid_only` wins if both are set. Flip both to `false` once the flux/source
kernels exist.

## Run it

```bash
# 1. get a closed, all-quadrilateral cubed-sphere grid (no gmsh needed)
julia --project=. tools/generate_cubed_sphere.jl 10 6.371e6 ./meshes/gmsh_grids/cubed_sphere.msh

#    ...or build it with gmsh from the .geo next to this README
# gmsh -2 problems/ShallowWater/SWsphere/cubed_sphere.geo \
#      -format msh2 -o ./meshes/gmsh_grids/cubed_sphere.msh

# 2. build the high-order grid
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphere")
```

> **Run this in a fresh Julia session after pulling.** `run_case` re-includes
> `run.jl`, `problems/drivers.jl` and the case's `user_*.jl`, but *not* the rest
> of the module: `src/kernel/mesh/sphere_mesh.jl` and `src/io/write_output.jl`
> are only evaluated when `Jexpresso` itself is loaded. A session opened before
> you pulled will run the new `drivers.jl` against the old module and die with
> `UndefVarError: mod_mesh_sphere_driver not defined in Jexpresso`.

Output, in `problems/ShallowWater/SWsphere/output/`:

| file | what it is |
|---|---|
| `sphere_grid_ho.vtu` | the high-order grid: **(ngl-1)² sub-elements per spectral element**, exactly as `write_vtk_grid_only` does for the flat cases. Every LGL node is a corner of a sub-cell — the linear elements are never written on their own. Carries the initial condition as point data when one has been built: `h`, `u`, `v`, plus `velocity`, the Cartesian form of `(u,v)` for glyphs and streamlines. |

One file per run, like the other cases. View it with representation **"Surface
With Edges"**: the edges drawn are the sub-element boundaries, so the LGL point
distribution is what you see.

## What was actually added

A spherical shell is a **2D manifold embedded in 3D**, and the existing 2D gmsh
path in `src/kernel/mesh/mesh.jl` cannot represent one: it stores only `(x,y)`
per node and interpolates only `(x,y)` when it adds the high-order points, which
would collapse the sphere onto its equatorial disc. So the shell gets its own
builder:

* `src/kernel/mesh/sphere_mesh.jl` — reads the linear quad shell (`MSH 2.2` and
  `MSH 4.1` ASCII), populates it with LGL points, verifies it.
* `initialize.jl` — the Galewsky et al. (2004) jet (see below).
* `src/io/write_output.jl` — `write_vtk_sphere_grid`, next to the existing
  `write_vtk_grid_only` writers.
* `tools/generate_cubed_sphere.jl` — equiangular gnomonic cubed-sphere generator.
* `test/test_sphere_mesh.jl` — standalone tests (`julia --project=. test/test_sphere_mesh.jl`).

### The node numbering, which is the delicate part

The shell is **closed**: it has no boundary, so there is no "outer edge" to
anchor a numbering to, and every edge is an interior edge shared by exactly two
elements. Four things are done about that.

1. **Coincident vertices are merged first.** A cubed-sphere `.msh` produced
   panel by panel very often carries *duplicate* nodes along the panel seams.
   Read naively, such a grid is six disconnected patches that merely look
   watertight, and every CG/SEM assembly across a seam is silently wrong. The
   merge (spatial hash, tolerance `:node_merge_tol` relative to the radius) is
   what turns it into one shell. If your file is already watertight the merge is
   a no-op and says so.

2. **High-order edge nodes are created once per _unique_ edge**, in a canonical
   direction, and the two neighbouring elements both point at those same nodes —
   reversing their traversal when they walk the edge the other way. Continuity
   across a seam is therefore *node identity*, not a floating-point coincidence.
   Creating the points per element instead would duplicate every seam node,
   which is the classic failure mode on a closed surface.

3. **Every quad is re-oriented outward** (normal away from the centre), because
   gmsh gives no global orientation guarantee on a multi-panel surface, and the
   metric terms and the tangent-plane velocity basis to come will need a
   right-handed `(i,j) → outward normal` frame in every element.

4. **Everything is then verified** (`check_sphere_mesh`, switch `:lcheck_grid`).
   Each test prints `PASS`/`FAIL`:

   | test | what it proves |
   |---|---|
   | T1 | every unique edge is used by exactly 2 elements → closed, manifold |
   | T2 | `V - E + F = 2` → a genus-0 closed surface, not a torn/duplicated one |
   | T3 | every edge traversed once in each direction → globally consistent orientation |
   | T4 | the element graph is one connected component → the panels really are joined |
   | T5 | every index in `1:npoin` is used, and `npoin` is exactly `npoin_linear + nedges*(ngl-2) + nelem*(ngl-2)²` |
   | T6 | no two distinct node indices sit at the same point → no seam numbered twice |
   | T7 | every node lies on the shell to round-off |
   | T8 | no folded/inverted sub-cell |
   | T9 | two elements sharing an edge list the **same** node ids along it |

   T5, T6 and T9 are the ones that catch a torn seam. `:lstop_on_bad_grid`
   (default `true`) aborts the run if any of them fails.

   The shell area is also reported against `4πR²`; it converges with `:nop`.

### Geometry of the high-order points

Edge points are placed by **great-circle (slerp) interpolation** at the LGL
abscissae; element interiors by a **Coons transfinite patch** built from those
four great-circle edges, then projected radially onto the shell. Both are exactly
on the sphere, and the edge construction depends only on the edge's two
endpoints, so the two sides of a seam agree.

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

State vector, in the **local tangent basis** of the shell:

```
q = [h, u, v]      h  depth [m]
                   u  zonal      (eastward, +λ) [m/s]
                   v  meridional (northward, +φ) [m/s]
```

with `e_λ = (-sin λ, cos λ, 0)` and `e_φ = (-sin φ cos λ, -sin φ sin λ, cos φ)`,
so the Cartesian velocity is `u·e_λ + v·e_φ` — written to the VTK file as
`velocity` because ParaView cannot interpret tangent components on its own.

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

## What comes next

The pieces still missing before this becomes a solver, in order:

1. metric terms of the manifold (covariant/contravariant basis, surface
   Jacobian, curvature) — the equivalent of `src/kernel/mesh/metric_terms.jl`
   for a 2D surface in 3D;
2. the tangent-plane velocity basis and the mass/differentiation matrices on it;
3. `user_flux.jl` / `user_source.jl` — the shallow water fluxes with the surface
   metric, the full Coriolis term `f = 2Ω sin φ`, and the curvature term that
   keeps the momentum tangent to the shell.

Then flip `:linit_only` to `false`.
