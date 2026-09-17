# `CompEuler/ffs_step_M7_round` — Mach-7 forward-facing step, **filleted corner**

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "ffs_step_M7_round")
```

`CompEuler/ffs_step_M7` with the step's convex corner rounded, and **nothing
else changed**. The two case directories differ in exactly two things:

| | `ffs_step_M7` | here |
|---|---|---|
| mesh | `ffs_step_transfinite.msh` (structured) | `ffs_step_round.msh` (**unstructured quads**) |
| deck | — | `:exact_geometry => Dict("fillet" => (:circle, 0.65, 0.15, 0.05))` |
| `user_bc.jl` | skips the vertical-face projection at `(0.6, 0.2)` | no corner special case at all |

Free stream, fluxes, primitives, source, `:Δt`, `:tend`, `:μ`, the DynSGS
settings and every `JEXPRESSO_M7_*` sweep switch are identical, so a
side-by-side run isolates one variable: **is the corner singularity what kills
the Mach-7 step?**

## Why

At a sharp `(0.6, 0.2)` the fluid's interior angle is **270°** — a point
singularity of the Euler solution. The expansion there is a Prandtl–Meyer fan
centred on that single point, and what it sheds downstream is an *entropy
layer*: a contact-type feature that convects and never self-heals.

A residual-based viscosity cannot fix that, and not because it is tuned wrong —
**an expansion fan is a smooth solution, so the sensor correctly returns almost
nothing in it.** That is why every knob on `ffs_step_M7` bought a factor and
none of them fixed it, and why `JEXPRESSO_M7_FILTER=0.005` — a *complete*
removal of the top mode at P4 — moved the death only from 7.5e-4 to 1.1e-3
while the structures downstream of the corner survived it.

Woodward & Colella (1984) §IV say this outright for this exact problem: the
step corner is *"a singular point of the flow"*, and they reset the state in
the cells beside it — same entropy and stagnation enthalpy as upstream —
because the spurious entropy layer otherwise streams downstream and corrupts
the Mach stem. **At Mach 3.** At Mach 7 the fan is far stronger and expands
into far lower density.

This deck removes the singularity from the **geometry** instead of patching the
solution. The ramp's sharp leading edge is the same kind of point, so whatever
this shows transfers to `rampCaoEtAl2021`.

## The mesh

An arc of radius `r = 0.05 m` centred at `(0.65, 0.15)`, tangent to the step
face at `(0.6, 0.15)` and to the step top at `(0.65, 0.2)`.

```
   ...........                        ...........
             |                                  \
             |  <- 270° point           r        )  <- turn spread over an arc
   __________|     singularity          _________/
```

Two things have to be right, and the first version of this mesh got both wrong.

### 1. Unstructured, not a transfinite block

Forcing a structured block around the fillet — the step face *and* the arc on
one side, a straight line opposite — makes the transfinite map shear the cells
at exactly the place the case is failing. Measured:

| mesh | minSICN (1 = perfect) |
|---|---|
| sharp, rectangular | 1.000 |
| **filleted, transfinite block (first attempt)** | **0.235** |
| filleted, unstructured quads (this mesh) | **0.676** |

A cell that distorted, right where the flow is failing, is worse than the sharp
corner it was meant to cure. gmsh's quasi-structured quad algorithm
(`Mesh.Algorithm = 11`) stays near-Cartesian away from the fillet and absorbs
the geometry locally. The algorithm/recombination/subdivision combination was
picked by measurement, not by taste — the `.geo` records why.

### 2. The arc must be *curved*, not a polygon

gmsh writes a **linear** grid: the `Circle` comes back as four straight
segments whose endpoints happen to sit on the circle. Filling those elements
with LGL nodes puts every high-order node on the **chord**, so however large
`:nop` is, the wall the solver sees is a polygon — and a free-slip wall then
generates spurious vorticity at every polygon corner. Rounding the corner and
discretizing it as a polygon just trades one corner for several.

That is what the separate `"fillet"` physical group is for:

```julia
:exact_geometry => Dict("fillet" => (:circle, 0.65, 0.15, 0.05)),
```

`src/kernel/mesh/exact_geometry.jl` snaps the high-order nodes of those edges
onto the true circle and blends the element interiors (Kopriva, *J. Sci.
Comput.* **26**(3):301–327, 2006, §3 — the linear-blending transfinite map,
which stays in `P^N` and therefore preserves the discrete metric identities and
the free stream exactly). Same mechanism `shock_circle` uses for its cylinder.

The **explicit** `(:circle, xc, yc, r)` form, not the `:circle` shorthand: the
shorthand fits centre and radius from the linear vertices, and a refined grid
puts new vertices at chord midpoints, which makes the fit ambiguous and it is
refused. Stating it keeps `JEXPRESSO_M7_REF` working.

The arc carries only **four** linear segments, deliberately: after the snap
those give 17 boundary nodes exactly on the circle at `:nop => 4`.
Over-refining the arc to chase the geometry would only cut `Δt` for nothing.
Fold margin (element thickness against the segment sagitta) is **17.5×**, so
`_check_curved_elements` has plenty of room.

### Numbers

| | sharp | filleted |
|---|---|---|
| elements | 4032 | **4150** |
| minSICN | 1.000 | 0.676 |
| inverted cells | 0 | **0** |
| edge length | uniform 0.025 | 0.0168 – 0.0380 |

`Δx_min` is **1.49× smaller**, so at the same `:Δt` the printed advective CFL
runs ~1.5× and the viscous one up to ~2.2× the sharp case's — which printed
0.036 and 0.12, so both still have room at `Δt = 5.0e-8`. For the comparison at
matched *CFL* rather than matched `Δt`, use `JEXPRESSO_M7_DT=3.4e-8`.

The `Mesh.*` options live inside the `.geo`, so this reproduces the committed
mesh exactly, and `rfac` changes the fillet radius:

```bash
gmsh -2 ffs_step_round.geo -o ffs_step_round.msh
```

## Why the corner BC special case is gone

On the sharp mesh, `(0.6, 0.2)` belongs to both the vertical step face
(`n = ±(1,0)`) and the horizontal step top (`n = ±(0,1)`). `build_custom_bcs_dirichlet!`
walks the boundary edge by edge, so that node is projected **twice** — zeroing
`u` *and* `v`, planting a no-slip stagnation point in the middle of an
expansion fan. `ffs_step_M7/user_bc.jl` skips the vertical-face projection
there to work around it.

The fillet makes the workaround unnecessary. There is no `(0.6, 0.2)` node any
more, and every node on the arc lies on a single wall segment with one
well-defined normal. The two *tangent* points are each shared by two segments,
but the arc is tangent to the face at one and to the step top at the other, so
the two normals differ only by the arc's turn over one segment (30° across the
three-segment arc). Two nearly-parallel projections are nearly idempotent —
they do not destroy the velocity the way two orthogonal ones do. **That is the
whole point of the fillet**, and it applies to the discrete boundary condition
just as much as to the continuous flow.

## The comparison to run

Same command, two cases:

```bash
J="julia --project=. src/Jexpresso.jl CompEuler"
JEXPRESSO_M7_TEND=1.5e-3 $J ffs_step_M7          # sharp:  dies at 7.5e-4
JEXPRESSO_M7_TEND=1.5e-3 $J ffs_step_M7_round    # filleted: ?
```

If the filleted run goes substantially further **and** the density field is
free of the bead chain streaming downstream of the corner, the singularity was
the cause and the next step is the ramp's leading edge. If it beads anyway, the
corner is exonerated and the problem is the Mach-7 expansion itself — at which
point `JEXPRESSO_M7_CMIN` (a background viscosity floor, which acts in an
expansion where the residual sensor by construction does not) is the remaining
lever, and a positivity floor is the one after that.

## Status

**Not yet run.** The mesh is generated and verified — 4150 quads, no inverted
cells, minSICN 0.676, physical groups `inflow`/`outflow`/`wall`/`fillet`
intact, four linear segments on the arc with a 17.5× fold margin. No Julia was
available in the environment this was written in, so the deck itself has not
been executed or parsed, and in particular the `:exact_geometry` snap has not
been seen to run on this mesh. Watch for the
`# SNAP HIGH-ORDER NODES ONTO EXACT GEOMETRY` banner and the reported wall
distance on the first run.
