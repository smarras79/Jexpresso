# SWsphere — shallow water equations on a spherical shell

**Status: GRID ONLY.** This case builds the high-order spectral-element grid on a
closed spherical shell, verifies it, writes it to VTK, and stops. There is no
initial condition and no time integration yet — that is what `:lgrid_only => true`
in `user_inputs.jl` is for.

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
> `UndefVarError: write_vtk_sphere_… not defined in Jexpresso`.

Output, in `problems/ShallowWater/SWsphere/output/`:

| file | what it is |
|---|---|
| `sphere_grid_ho.vtu`        | the high-order grid: **(ngl-1)² sub-elements per spectral element**, exactly as `write_vtk_grid_only` does for the flat cases. Every LGL node is a corner of a sub-cell — the linear elements are never written on their own. |
| `sphere_grid_wireframe.vtu` | the same grid as explicit line segments: every LGL grid line of every element. `line_type == 1` are the spectral-element boundaries, `line_type == 0` the interior LGL lines. |
| `sphere_grid_points.vtu`    | the LGL nodes themselves, one vertex each. Representation "Points", colour by `node_type` (0 vertex / 1 edge / 2 interior). |

## What was actually added

A spherical shell is a **2D manifold embedded in 3D**, and the existing 2D gmsh
path in `src/kernel/mesh/mesh.jl` cannot represent one: it stores only `(x,y)`
per node and interpolates only `(x,y)` when it adds the high-order points, which
would collapse the sphere onto its equatorial disc. So the shell gets its own
builder:

* `src/kernel/mesh/sphere_mesh.jl` — reads the linear quad shell (`MSH 2.2` and
  `MSH 4.1` ASCII), populates it with LGL points, verifies it.
* `src/io/write_output.jl` — `write_vtk_sphere_grid` / `write_vtk_sphere_wireframe` /
  `write_vtk_sphere_points`, next to the existing `write_vtk_grid_only` writers.
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

## What comes next

The pieces still missing before this becomes a solver, in order:

1. metric terms of the manifold (covariant/contravariant basis, surface
   Jacobian, curvature) — the equivalent of `src/kernel/mesh/metric_terms.jl`
   for a 2D surface in 3D;
2. the tangent-plane velocity basis and the mass/differentiation matrices on it;
3. `initialize.jl` — the Williamson et al. (1992) test cases. The grid already
   carries `mesh.lon`, `mesh.lat` and `mesh.radius` per node, which is what they
   need;
4. `user_flux.jl` / `user_source.jl` — the shallow water fluxes with the surface
   metric, the full Coriolis term `f = 2Ω sin φ`, and the curvature term.

Then flip `:lgrid_only` to `false`.
