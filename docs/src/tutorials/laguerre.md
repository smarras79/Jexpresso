# Semi-infinite (Laguerre) elements

A boundary of the finite domain can be extended to infinity by a layer of
semi-infinite elements: Lagrange polynomials (`:nop`) along the boundary,
scaled Laguerre functions (`:nop_laguerre`) away from it. Waves leave the
finite domain through the layer instead of reflecting off an artificial
boundary.

## Turning it on

**1D.** `:llaguerre_1d_right => true` and/or `:llaguerre_1d_left => true`
attach one semi-infinite element to that end of the interval.

**2D.** Tag the boundary edges to extend with the physical name `"Laguerre"`
in the gmsh file. One semi-infinite element is built on each of them.

| Key | Default | Meaning |
|---|---|---|
| `:nop_laguerre` | 18 | polynomial order in the semi-infinite direction (`ngr = nop_laguerre + 1` nodes) |
| `:laguerre_beta` | 1.0 | scaling β of the Laguerre functions |
| `:xfac_laguerre`, `:yfac_laguerre` | 1.0 | stretching of the layer in x and y (1D uses `:yfac_laguerre`) |
| `:lperiodic_laguerre` | false | the layer is periodic where the finite domain is |
| `:llaguerre_bc` | false | impose `:laguerre_tag` on the far end of the layer |
| `:laguerre_tag` | `"none"` | the boundary tag `user_bc_dirichlet!` receives there |

## Boundary conditions on the layer

The semi-infinite elements use the same Dirichlet machinery as the rest of the
mesh: the values `user_bc_dirichlet!` sets are inferred once at setup and then
imposed on the state and on its tendency (no-penetration is recognised and
projected), see `src/kernel/boundaryconditions/bc_constraints.jl`. Which tag a
layer node gets:

* **Sides of the layer** (2D). Nothing is assumed about the direction of the
  layer (it may extend in x, in y, or along any straight face) or about the
  order of its elements: the sides are found from the mesh topology. A column
  side whose base node lies on a boundary of the finite domain — for example
  the lateral walls under a top layer, or the top and bottom boundaries beside
  a layer extending in x — continues that boundary and **inherits its tag and
  its condition**: a Dirichlet boundary stays Dirichlet, a wall stays a wall, a
  periodic boundary stays periodic (no constraint). Normal: perpendicular to
  the layer, outward. Sides shared by two semi-infinite elements are interior.
* **Far end** (the last Laguerre node of every column). Unconstrained unless
  `:llaguerre_bc => true`, in which case it gets `:laguerre_tag`, with the
  direction of the column as its normal. (The GPU kernel and the linear-solve
  path, e.g. Helmholtz, always impose it; the latter with a homogeneous value.)
* **1D.** An end with a semi-infinite element attached is not a boundary (no
  `"left"`/`"right"` condition there); its far node gets `:laguerre_tag` when
  `:llaguerre_bc => true`.

The same rule is used by every path that imposes boundary values: the inferred
constraints above, the legacy Dirichlet path, the GPU kernel (which receives the
inherited sides as per-element normals, since it cannot use the string tags) and
the linear-solve path (`laguerre_side_columns_2D` in `bc_constraints.jl`).

The semi-infinite elements are assembled together with the finite ones, before
the boundary constraints are applied, so a node shared by both is constrained
exactly like any other boundary node.
