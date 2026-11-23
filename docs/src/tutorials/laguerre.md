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

## TOTAL and PERT variables

Every scaled Laguerre function decays, ``e^{-\beta\xi/2}p(\xi)``, so the layer
can only represent fields that vanish far from the domain. A perturbation does;
a total field (the hydrostatic density and pressure, θ, the fluxes built from
them) does not, and its derivative in the layer would be wrong.

* **PERT**: nothing to do, the prognostic variables are the perturbation.
* **TOTAL**: inside the semi-infinite elements the right-hand side is computed
  on the departure from the reference state `qe`:
  ``\nabla\cdot(F(q)-F(q_e)) - (S(q)-S(q_e))`` for the inviscid part and the
  diffusion of `prim(q) - prim(qe)` for the viscous part. This is exactly the
  TOTAL system wherever `qe` is a steady state (``\nabla\cdot F(q_e) = S(q_e)``,
  e.g. hydrostatic balance) -- the same assumption PERT makes. The finite
  elements are unchanged.

## Filter

`:lfilter => true` filters the semi-infinite elements too, with `:mu_x` along
the base edge and `:mu_y` in the Laguerre direction (the reference directions
ξ and η of the element), for TOTAL and PERT. The Laguerre modes of the filter
use the same `:laguerre_beta` as the basis. With TOTAL the filter acts on
`q - qe`.

## Viscosity

With `:lvisc => true` the semi-infinite elements diffuse with the constant
coefficients `:μ` (one per equation), in 1D and 2D. Turbulence models
(`:visc_model`) are not used in the layer yet. The layer nodes are at
`x = xmax + yfac_laguerre*ξ` (1D, right end), so a smaller `:yfac_laguerre`
packs the nodes near the interface and, with explicit time stepping, needs a
smaller `:Δt` for the diffusion to stay stable.
`problems/AdvDiff/laguerre1d_visc` advects and diffuses a Gaussian into the
layer and stores the exact solution in `qe` for comparison.

## Boundary fluxes

`:bdy_fluxes => true` works on meshes with a Laguerre layer: fluxes are applied
on the boundary edges of the finite domain. They are not applied on the sides
of the layer.

## Linear solves (Helmholtz, Elliptic)

Both storage modes (`:lsparse`) impose Dirichlet rows on the boundary edges of
the finite domain (periodic and `"Laguerre"` interface edges excluded) and on
the layer: inherited sides as above, homogeneous far end.

## Accuracy

The layer converges spectrally in both `:nop` and `:nop_laguerre`. A
manufactured Helmholtz solution decaying into the layer
(``\sin(x/2)e^{-x/2}\cos y``) is recovered to ``10^{-13}``, and the 1D
advection-diffusion test above to ``10^{-13}`` (`:nop` 8), with
`:nop_laguerre` up to at least 63; nothing measurable is reflected at the
interface. Three practical points:

* Check the initial state against the inflow boundary values: a Gaussian
  whose tail is ``10^{-10}`` at a Dirichlet boundary held at 0 leaves an
  error of that size behind it, whatever the resolution.

* The layer must reach far enough for the solution to have decayed at its
  last node, which is homogeneous: the nodes are at ξ times the stretching
  factor (`:xfac_laguerre`/`:yfac_laguerre`) from the base, with
  ``\xi_{max} \approx 3.5\,`` `nop_laguerre`. A short layer (small factor,
  few nodes) truncates the solution and converges slowly in `:nop_laguerre`.
* The basis decays like ``e^{-\beta\xi/2}``; a solution decaying much more
  slowly than that in ``\xi`` needs more Laguerre nodes.
