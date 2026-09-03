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

### In parallel

The case runs on any number of MPI ranks:

```bash
julia --project=. -e 'using MPI; run(`$(mpiexec()) -n 4 $(Base.julia_cmd()[1]) --project=. \
      -e "using Jexpresso; Jexpresso.run_case(\"ShallowWater\", \"SWsphere\")"`)'
```

or with whichever `mpiexec` your MPI.jl is configured against (see `INSTALL.md`
for the system-binary setup). Nothing in the deck changes.

WHAT IS AND IS NOT EXCHANGED. The grid is partitioned by the ordinary Jexpresso
reader — each rank owns a subset of the elements and holds, in addition, a
mirrored copy of every node those elements share with a neighbour. The metric
terms are element-local and need no communication at all. What does need it is
every quantity that is a SUM OVER THE ELEMENTS TOUCHING A NODE, because on a
partition seam each rank holds only its own share of that sum:

| assembled quantity | where |
|---|---|
| the mass matrix `M` | `build_sphere_metrics` |
| the RHS `∂q/∂t` | `sphere_rhs!`, every RK stage |
| the filter's mass-weighted average | `sphere_filter!`, every step |
| the relative vorticity | `sphere_relative_vorticity!`, every diagnostic |

Each of those completes its sum with `assemble_mpi!` — the cross-rank half of
the direct stiffness summation — BEFORE dividing by `M`. After it, a mirrored
node carries exactly the value its owner does, so the state stays identical on
every rank that holds it and no separate halo exchange is needed. The Lagrange
projection and the RK update are node-local and therefore consistent for free.

Everything global is reduced: `Δt` (a MAX over the wave speed, so all ranks take
the same step), the conserved integrals (a SUM over OWNED nodes only — a
mirrored node would otherwise be counted once per rank), the drift, `max|ζ|`,
and every residual in `check_sphere_metrics`, whose verdict is therefore the
same everywhere.

MEASURED, by running THIS deck — the full 10-day Galewsky integration, 11 521
steps at `nop = 5` — on 1, 2 and 4 ranks. Every diagnostic it prints agrees on
all three, to every digit printed:

| after 10 days | 1 rank | 2 ranks | 4 ranks |
|---|---|---|---|
| `δE/E` | -1.008e-03 | -1.008e-03 | -1.008e-03 |
| `max\|ζ\|` | 8.737e-05 | 8.737e-05 | 8.737e-05 |
| `max\|ζ-ζ₀\|` | 1.623e-04 | 1.623e-04 | 1.623e-04 |
| max drift removed | 1.361e+03 | 1.361e+03 | 1.361e+03 |
| `δmass/mass` | 9.555e-12 | 9.559e-12 | 9.561e-12 |

The mass residual is the only column that moves, in its last digits, because the
order of summation changes with the partition — that is round-off, not a
different answer. All seven metric checks pass at every rank count (M4, the
area, at 4.9e-14 / 5.7e-14 / 5.1e-14).

WHAT IT BUYS. Same runs, wall time of the time loop, on a 4-core sandbox — 600,
300 and 150 elements per rank:

| ranks | time loop | of which compiling | speed-up |
|---|---|---|---|
| 1 | 152.4 s | 5.0 % | — |
| 2 | 81.3 s | 9.6 % | **1.87×** |
| 4 | 46.9 s | 16.6 % | **3.25×** |

— 1.97× and 3.71× on the compute alone, i.e. close to ideal for a grid this
small. Do measure on a short run and you will see none of it: at one simulated
day the per-process JIT (~35 s a rank, most of it the first `with_mpi` call)
swamps the integration and four ranks come out no faster than one. The
integration is what scales; the startup does not.

THE PARTITION is the generic one (`_compute_xy_partition`, a box decomposition
of the projected x-y plane), so a part of a SPHERE wraps around the far side
rather than being a compact patch. It is balanced — 0.97 to 1.00 load ratio at 4
ranks — and the mirrored-node overhead is small (15 404 local nodes against
15 002 global, i.e. 2.7 %), so it has not been worth replacing with a
sphere-aware partitioner. Re-measure before assuming that still holds at high
rank counts.

Output, in `problems/ShallowWater/SWsphere/output/`:

| file | what it is |
|---|---|
| `sphere_grid_ho.vtu` | the initial state |
| `sphere_0001.vtu` … | one per output time, `:ndiagnostics_outputs` of them |

On more than one rank each of those becomes a `.pvtu` plus one `.vtu` piece per
rank in a directory of the same name — the convention `write_vtk_grid_only` uses
for the flat cases. Open the `.pvtu`. The `ip` point field is the GLOBAL node
number, so it stays continuous across a partition seam as well as across a panel
seam; the `part` cell field is the owning rank.

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

   **Winding is not your problem, though.** gmsh orients each panel's surface
   on its own, from the direction of the curve loop that bounds it, and the six
   panels do not always come out the same way: a 32×32 grid written by gmsh
   4.1 from this same `.geo` had its whole −z panel wound inward, which is a
   negative surface Jacobian and used to stop `build_sphere_metrics` on
   "element 1". `orient_shell_elements_outward!` (in `mesh.jl`, run right
   after the nodes are placed on the shell) now transposes the (ξ,η)
   numbering of every element whose `a_ξ × a_η` points inward, reports how
   many it flipped, and the metrics see a consistently outward grid whatever
   the file's winding. `test/test_sphere_orientation.jl` reads a copy of this
   grid with one panel reversed and checks the metrics come out identical.

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

`cubed_sphere.geo` does **not** build the gnomonic cubed sphere, despite what
this section used to say: gmsh spaces `Transfinite Line` points at equal *angle*
along a `Circle` arc, so the shipped `.msh` is already the **equiangular** map of
Ronchi, Iacono & Paolucci (1996). It reproduces
`tools/generate_cubed_sphere.jl`'s `:equiangular` output to 2.2e-16 of `R`, and
differs from the gnomonic grid by 535 km.

One input slides the nodes onto a different map *after* the grid is read.
Connectivity, the panel decomposition and the panel boundaries are untouched;
only where the nodes sit **inside** each panel changes, so this is safe to do
once `mod_mesh_read_gmsh!` has already built the topology.

```julia
:cubed_sphere_map => :none,   # :none (default) | :gnomonic | :equiangular | :conformal
```

The map the nodes come **from** is *measured*, by `detect_cubed_sphere_map`,
which pushes every linear vertex back through each candidate's inverse and
reports which one lands on the uniform panel lattice (residual 2e-16 for the
right map against 1e-1 for the wrong ones). There is deliberately no input for
it: a "source map" switch is a claim about the `.msh` that nothing can check,
and getting it wrong silently produces a grid that is neither map — which is
exactly what happened while the source was *assumed* gnomonic, since asking for
`:equiangular` then warped an already-warped grid and cut the minimum element
edge from 710 km to 562 km, 21% of the time step, in silence.

| value | map | what it buys |
|---|---|---|
| `:none` | leave the grid as read | what this deck ships with — the grid is already equiangular |
| `:gnomonic` | equidistant central projection — Sadourny (1972) | nothing here; it *un-warps* the grid to the coarser-cornered projection |
| `:equiangular` | face coordinate as an angle, `u = tan α`, `α ∈ [-π/4, π/4]` — Ronchi, Iacono & Paolucci (1996) | a no-op on this grid; the remap detects that and says so |
| `:conformal` | Rančić, Purser & Mesinger (1996) | grid lines at 90° into a cube corner — and **refused on this grid**, see below |

### Corner behaviour, which is what actually distinguishes them

Angle between the two families of grid lines, walking the face diagonal
`u = v = t` in to a cube corner. 90° is orthogonal:

| `t` | gnomonic | equiangular | conformal |
|---|---|---|---|
| 0.0 (face centre) | 90.000° | 90.000° | 90.000° |
| 0.9 | 116.584° | 114.947° | **90.000°** |
| 0.99 | 119.668° | 119.482° | **90.000°** |
| 0.9999 | 119.997° | 119.995° | **90.000°** |

Gnomonic and equiangular both open to a 120° corner; only the conformal map
holds 90° all the way in. Over the whole face the worst deviation from
orthogonality is 30.000° for both of the first two — attained *at* the corner —
and below the finite-difference floor (1e-6) for the conformal map.

### Why `:conformal` is refused here, and why no regularisation of it works

**The 120° is topological.** Three panels meet at a cube corner and 360° is all
there is to share, so each panel opens 120° there. Now let `r` be a face map that
is *differentiable* at the corner with a *non-singular* differential `A`. The two
grid lines through the corner are `u ↦ r(u,1)` and `v ↦ r(1,v)`, their tangents
are `A e₁` and `A e₂`, and those two curves are the two cube-edge arcs. Hence

> `∠(A e₁, A e₂) = 120°` for **every** differentiable, non-degenerate map.

So **30° is the smallest deviation from orthogonality any usable map can have at
a cube corner** — `:gnomonic` and `:equiangular` both attain it — and a map that
keeps the square's 90° there must have `det A = 0`. That is what the conformal
map does: near a corner it behaves like `ζ^(4/3)`, so its local scale `|∂r/∂u|`
falls off as `d^(1/3)`. Measured, as a fraction of the face-centre value:

| `d = 1-u` | 1e-1 | 1e-2 | 1e-3 | 1e-4 |
|---|---|---|---|---|
| conformal | 0.444 | 0.206 | 0.0957 | 0.0444 |
| equiangular | 0.925 | 0.940 | 0.943 | 0.943 |

Fitted exponent 0.333 over four decades — it goes to **zero** at the corner, and
the grid puts a mesh vertex on each of the eight cube corners (`Point(2)`…
`Point(9)`, shared by the panels). Measured on this grid remapped to the pure
map, the surface Jacobian at a corner node comes out *positive but collapsing* —
1/27, 1/51, 1/93 of the grid median at `nop = 3, 5, 8`, smaller the finer the
grid, because the nearest LGL node keeps closing on the singular point — so
`build_sphere_metrics` does not necessarily reject the element as degenerate. It
builds wrong metrics instead, and `check_sphere_metrics` says so: M6 = 0.10,
0.14, 0.17 (`:radial`) and 1.5, 1.8, 1.7 (`:cross_product`) at those same orders,
against a 5e-2 tolerance. A node landing exactly on the zero would give the
blunter failure,

```
ERROR sphere_metrics.jl: non-positive surface Jacobian in element 1
at node (1,1). The element is degenerate or wound inward.
```

Either way the run is worthless, so `remap_cubed_sphere_nodes!` refuses up front,
where the cause can be named. The map itself is correct —
`test/test_cubed_sphere_maps.jl` shows it is conformal to 4e-9 — it is the
combination of that map with a node ON the singular point that does not work.
This is why conformal cubed spheres are used by cell-centred finite-volume
codes (Rančić's own model, CCAM), where the singular point is a cell corner
nobody evaluates at, and not by nodal spectral elements.

**A corner-stretched `:conformal` was shipped briefly and has been removed.** It
reparametrised the face coordinate diagonally, `r(u,v) = C(φ(u), φ(v))` with
`φ(t) = sign(t)(1-(1-|t|)^(3/4))`, the exponent picked to cancel the `d^(1/3)`
collapse. It does make the Jacobian at the corner node finite, and because `φ`
acts on each coordinate separately the grid lines stay orthogonal wherever the
derivatives exist. Both of those are true and neither is enough:

* orthogonality *plus* a non-zero Jacobian at the corner is precisely what the
  argument above forbids, so the 30° of shear gets spread over the angular
  sector and the corner becomes a **cone point** — the map turns positively
  homogeneous of degree 1 instead of differentiable. Measured `|dr/ds|` along
  the ray `(1-u, 1-v) = t(cos θ, sin θ)`: 0.639, 0.673, 0.705, 0.717 at
  `θ = 0°, 15°, 30°, 45°`, flat in `t` over `1e-2 … 1e-6`. No linear map fits
  that, so no polynomial element can represent it;
* a **separable** stretch cannot be confined to the corner. `φ` acts on `u`
  alone, so `φ'(1) = ∞` hits the whole cube edge: at `v = 0.5`,
  `|r_u| = 0.523·(1-u)^(-1/4)` — 0.94, 1.65, 2.94, 5.23, 9.30 at
  `1-u = 1e-1 … 1e-5`, against a flat 0.783 for `:equiangular`. The map is not
  `C¹` on any of the twelve cube edges.

`check_sphere_metrics` measures both, on this grid remapped
equiangular → stretched-conformal (`:radial`, tolerances in `sphere_metrics.jl`):

| nop | 3 | 4 | 5 | 7 | equiangular, nop = 4 |
|---|---|---|---|---|---|
| M3 `(J n̂)·x̂` | 3.2e-04 | 3.0e-04 | 8.3e-05 | 3.3e-05 | 1.5e-10 |
| M4 area | 2.8e-05 | 8.5e-06 | 3.3e-06 | 8.0e-07 | 1.3e-11 |
| M6 curl normal | 0.414 | 0.408 | 0.403 | 0.399 | 2.7e-05 |
| M7 rigid rotation | 1.3e-07 | 1.1e-07 | 7.0e-08 | 4.4e-08 | 2.3e-12 |

M6 does not move with `nop`, and it does not move with `h` either — 0.403,
0.408, 0.410, 0.412 at `n = 5, 10, 20, 40` elements per panel edge — and all of
it sits in the 24 corner elements (edge elements 2.6e-03, interior 4.3e-05 at
`nop = 4`). An O(1) error that neither `p`- nor `h`-refinement touches is not a
tolerance that needs loosening: the same tolerances leave the equiangular grid
four to nine orders of margin at every `nop` from 3 to 8.

Blending the conformal face coordinate into the equiangular one near the corners
— which the argument above says is the only direction left — does pass all seven
checks (best case, weight `w = u²v²`: M3 1.6e-08, M4 1.0e-12, M6 2.0e-03, M7
7.4e-10 at `nop = 5`), and is still not worth shipping: M6 stays ~100× the
equiangular grid's and converges only first-order in `nop`, and the blend has to
be wide enough that the map is no longer conformal over most of the panel. What
is left at that point is a worse `:equiangular`. **Use `:equiangular`** — i.e.
this grid, as read.

### What it costs, measured on the shipped grid

Element edge lengths over a panel (all six are congruent), `n = 10`,
`R = 6.37122e6` m. The minimum is what sets the explicit time step under
`:lcfl_dt`:

| map | min edge | max edge | ratio | Δt vs the grid as shipped |
|---|---|---|---|---|
| equiangular (as shipped) | **710.2 km** | 999.8 km | **1.41** | 1.00× |
| gnomonic | 641.1 km | 1255.6 km | 1.96 | 0.90× |
| conformal | 476.0 km | 1033.0 km | 2.17 | 0.67× (and refused) |

So the grid as shipped is both the most homogeneous and the cheapest of the
three, and `:cubed_sphere_map` has nothing to add to it. Do not generalise the
Δt column to finer grids: the conformal map's cell scale falls off toward a
corner, so the more elements per panel edge, the closer the corner-adjacent
nodes sit to the singular point and the more shrinkage a refined grid sees.
Re-measure if you refine.

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
