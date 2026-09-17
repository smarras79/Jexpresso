# `CompEuler/shock_circle_M7` — Mach-7 laminar cylinder

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "shock_circle_M7")
```

`CompEuler/shock_circle` raised from Mach 3 to Mach 7 and made **viscous**, with
the grid clustered at the wall so the laminar boundary layer and the surface
heat flux are resolved. The curved-wall rung on the way to `rampCaoEtAl2021`.

## Free stream — deliberately not sea level

A Mach-7 cylinder at `p = 101325 Pa` would have `Re_D = 6e7`: the boundary layer
would be micrometres thick, turbulent, and unresolvable, and asking for a heat
flux would be meaningless. These are shock-tunnel conditions of the kind Cao
et al. (2021) use, with `p∞` chosen to put `Re_D` at 10⁴:

| | |
|---|---|
| `M∞`, `T∞`, `p∞` | 7, 125 K, 5 Pa |
| `ρ∞`, `c∞`, `u∞` | 1.3937e-4 kg/m³, 223.98 m/s, **1567.8 m/s** |
| `μ(T∞)` (Sutherland) | 8.656e-6 Pa·s |
| `Re_D` (D = 0.4 m) | **1.01e4** |
| `T₀` | 1345 K |
| `T_wall` | **300 K**, isothermal (`T_w/T₀ = 0.22`, a cold wall) |
| Knudsen | 1.0e-3 — comfortably continuum |

`T₀ = 1345 K` keeps calorically-perfect air defensible; at sea-level `T∞` the
same Mach number would demand 3150 K and real-gas chemistry. The gas is still
ideal: this is a **numerical** test at Mach 7.

## Grid — a modest cluster, not a y⁺ = 1 mesh

`δ ~ R/√(Re_R) = 2.8 mm` on a 200 mm radius. The size field puts **2.2 mm**
cells on the cylinder — about six LGL nodes across `δ` at `:nop => 4` — relaxing
to 60–77 mm in the far field.

| | |
|---|---|
| elements | 14296, **all quadrilateral** |
| minSICN | 0.609, **no inverted cells** |
| edge lengths | 0.00223 – 0.0770 m |
| `Δx_min` (LGL, nop 4) | 3.84e-4 m |

Near-square quads throughout — there is no anisotropic stretching anywhere.

### The circle is curved, not polygonal

gmsh writes a **linear** grid, so the circle arrives as 64 straight segments
whose endpoints merely happen to lie on it. Filling those elements with LGL
nodes would put every high-order node on a **chord**, and the wall the solver
sees would stay a polygon however large `:nop` is. On a **no-slip** wall that is
fatal twice over: the polygon corners shed spurious vorticity into the boundary
layer, and the wall-normal direction the temperature gradient — the heat flux —
is taken along is wrong by `O(h)` at every node.

So the deck carries

```julia
:exact_geometry => Dict("cylinder" => (:circle, 1.0, 0.0, 0.2)),
```

and `exact_geometry.jl` snaps the high-order nodes onto the true circle and
blends the element interiors (Kopriva, *J. Sci. Comput.* **26**(3):301–327,
2006, §3). On a curved wall this is the default expectation, not an
optimisation. Explicit centre and radius rather than the `:circle` shorthand,
which is refused on a refined grid.

**On the first run, watch for `# SNAP HIGH-ORDER NODES ONTO EXACT GEOMETRY` and
the wall distance it reports** — `shock_circle` gets 3e-16.

## The stabilising choice: `:lkep`, not more viscosity

Set in the deck, not left to a switch, because it is this deck's answer rather
than a sweep. It follows from what the `ffs_step_M7` experiments actually showed:

- filleting the step corner changed **nothing** (same death time to the digit) →
  not a geometric singularity;
- a *complete* removal of the top mode bought 1.47× → not the top mode;
- the finer grid died **earlier**, and the element-scale checkerboard was present
  in the **undisturbed free stream** upstream of the bow shock → not a physical
  feature, and not something a shock sensor can legitimately see.

What is left is the discretization. A collocation CG integrates the nonlinear
flux with the same LGL rule it interpolates on, so the flux is **aliased**;
energy goes into the grid-scale modes and CG has nothing to take it back out.
What Mach 7 changes is the *price*: `p = (γ-1)(ρE − ½ρ|u|²)` is a difference of
nearly equal numbers, so a relative error in `ρE` or `ρu` is amplified by
`1 + γ(γ-1)M²/2` — **3.5 at Mach 3, 14.65 at Mach 7**.

`:lkep => true` with `:volume_flux => ranocha()` assembles the inviscid RHS by
flux differencing from symmetric entropy-conservative two-point volume fluxes,
which removes that aliasing-driven transfer *by construction* — and, unlike a
filter or a viscosity floor, costs nothing inside the boundary layer, which is
what this case exists to resolve.

`central_euler()` is the honest control: if the case behaves identically under
it, the two-point machinery is not what is helping.

## Two settings that differ from the other decks, on purpose

**`:dsgs_sensor => "residual"`** (the default; `ffs_step`/`ramp` use `"legacy"`
for historical reasons). Two independent reasons here:

1. A cylinder's bow shock becomes **stationary**. `"legacy"` is `R ≈ |∂ₜq|`
   (`rhs.jl:205`) — a rate sensor that fires on a moving front and goes quiet on
   a standing one, which is exactly the steady state this case integrates
   towards. The element residual is `O(1/h)` at a discontinuity whether it moves
   or not.
2. A no-slip isothermal wall is a strong Dirichlet constraint on three of four
   slots. The residual path explicitly zeroes the residual at Dirichlet nodes
   (`_dsgs_boundary_pairs!`), because otherwise the constraint force reads as
   under-resolution and pins `ν` at its cap along the whole wall — measured on
   the rising bubble, where it blew the run up. The legacy branch returns before
   that correction.

**`:dsgs_Cmax => 0.1`** and **`:μ = [1,1,1,1]`**, following `rampCaoEtAl2021`:
the cap must be low enough that residual viscosity cannot smear the laminar
boundary layer beside it. `ffs_step`'s ×4 was measured on an inviscid run with
no boundary layer to protect.

**No `:dsgs_Cmin`.** A background floor is the obvious lever for a free-stream
checkerboard, but it applies viscosity everywhere *including inside the boundary
layer* and would corrupt the wall heat flux. `:lkep` is meant to make it
unnecessary. If the free stream still quilts, `:dsgs_Cmin => 0.01` is the first
thing to try — and the heat flux must then be re-checked.

## Boundary conditions

| tag | condition |
|---|---|
| `inflow`, `top`, `bottom` | free stream (all characteristics enter) |
| `outflow` | supersonic — nothing imposed |
| `cylinder` | **no slip, isothermal**: `u = v = 0`, `ρE = ρ cv T_w` |

`top`/`bottom` carry the free stream rather than a slip wall because the bow
shock asymptotes to the Mach angle `asin(1/7) = 8.2°` and reaches only
`|y| ≈ 0.49` at the outflow — it never touches `y = ±1`. Prescribing the
undisturbed state there is exact and cannot reflect.

Density is left unconstrained at the wall — it is the one variable a wall does
not fix, and continuity supplies it. `cv` is taken as `Rair/(γ-1)` from
`PhysConst` so the wall energy and the pressure use the same gas.

## Cost and what to look at

`1.0e-3 / 1.5e-8` = **66 700 steps**. `tend` is ~8 body times `R/u∞ = 1.28e-4 s`,
enough for a settled shock layer and a settled wall heat flux, and about half a
box flow-through. Advective CFL ≈ 0.08; the viscous number is ≈ 0.006 from the
molecular viscosity plus about as much again from DynSGS at `C_max = 0.1`.

Expect the bow shock to stand off `0.212 R = 42 mm`. `schlieren` shows it best.
For the heating, `T` is in `qoutvars`: the wall heat flux is
`q_w = k ∂T/∂n` along the now-exact wall normal.

## Startup

A uniform free stream is **not** a legal starting field for a no-slip isothermal
wall. At the wall the BC sets `u = v = 0` and `ρE = ρ cv T_w`, so `ρE` drops
183.85 → 30.02 while the node one 2.2 mm cell away is still at 183.85 —
dominated by the kinetic energy going (−171.3) rather than the temperature
(+17.5). That transient drove `p` to **−0.296 Pa** against `p∞ = 5 Pa` in the
first steps.

`rampCaoEtAl2021` already learned this and starts from a compressible laminar
boundary-layer profile instead of a uniform stream. Same idea here, adapted: a
cylinder has no similarity profile to lay down before the bow shock even
exists, so this is not a boundary layer — it is a **boundary-condition-consistent
field**. Velocity → 0 and `T` → `T_w` over `δ₀ = 5 mm` (≈2 wall cells, ≈2δ)
with the ramp's Pohlhausen blend `su = 2ζ − 2ζ³ + ζ⁴`, pressure held at `p∞`
across the layer, density following from `p∞` and the blended `T`. The wall node
then starts at `ρE = 12.55`, which is what the BC would set it to anyway — no
jump. The real boundary layer grows out of it.

Separately, `user_flux.jl` now floors `ρ` and `p` at `FLUXAUX_FLOOR = 1e-14`
before `ranocha()` takes their logarithms. A single node at `p ≤ 0` was not
degrading the answer, it was throwing `DomainError` out of `log` and killing the
run from inside the RHS. **If that floor engages anywhere but the first few
steps the run is wrong** — it is a guard, not a fix.

## Status

**Not yet run.** The mesh is generated and verified; no Julia was available in
the environment this was written in, so the deck has not been executed or
parsed. Two things to check on the first run: the `SNAP HIGH-ORDER NODES` banner
and its reported wall distance, and that the first CFL line's viscous number
matches the estimate above.
