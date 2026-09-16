# `CompEuler/ffs_step_M7_round` — Mach-7 forward-facing step, **filleted corner**

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "ffs_step_M7_round")
```

`CompEuler/ffs_step_M7` with the step's convex corner rounded, and **nothing
else changed**. The two case directories differ in exactly two things:

| | `ffs_step_M7` | here |
|---|---|---|
| mesh | `ffs_step_transfinite.msh` | `ffs_step_round.msh` |
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

An arc of radius `r = 0.05 m` (two elements) centred at `(0.65, 0.15)`, tangent
to the step face at `(0.6, 0.15)` and to the step top at `(0.65, 0.2)`.

```
   ...........                        ...........
             |                                  \
             |  <- 270° point           r        )  <- turn spread over an arc
   __________|     singularity          _________/
```

The three transfinite blocks survive: the block corner moves from `(0.6, 0.2)`
to `(0.6+r, 0.2)`, and block A's right-hand side becomes the step face *plus*
the arc (five curves on four sides, so `Transfinite Surface` names its corners
explicitly).

| | sharp | filleted |
|---|---|---|
| elements | 4032 | **4033** |
| edge length | uniform 0.025 | 0.0189 – 0.0260 (max/min **1.38**) |
| all quads, no inverted cells | yes | yes |

**This is not a stretched grid** — but `Δx_min` *is* 1.32× smaller, so at the
same `:Δt` the printed advective CFL is ~1.3× and the viscous one up to ~1.75×
the sharp case's. Both still have room at the default `Δt = 5.0e-8`. For a
comparison at matched *CFL* rather than matched `Δt`, use
`JEXPRESSO_M7_DT=3.8e-8`.

`n_arc = 3` was chosen over 4 deliberately: 4 gives max/min 1.70 instead of
1.38. Regenerate, or change the radius through `rfac`, with

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

**Not yet run.** The mesh is generated and verified (4033 quads, no inverted
cells, physical groups `inflow`/`outflow`/`wall` intact); no Julia was
available in the environment this was written in, so the deck itself has not
been executed or parsed.
