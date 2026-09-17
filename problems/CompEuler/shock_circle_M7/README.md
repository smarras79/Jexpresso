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

## Run log

| run | config | died at | steps | flow travel |
|---|---|---|---|---|
| 1 | hold 0, Cmax 0.1, Δt 1.5e-8 | 1.33e-5 | 887 | 20.9 mm |
| 2 | hold 2, Cmax 0.03, Δt 1.5e-8 | 5.97e-6 | 398 | 9.4 mm |
| 3 | hold 0, Cmax 0.1, **Δt 7.5e-9** | 2.98e-5 | 3969 | 46.7 mm |
| 4 | run 3 + **`:ldsgs_nodal`** | — | — | — |

Run 3 confirmed `Δt` was the right lever: halving it took the failure from
1.33e-5 to 2.98e-5 s and from 887 to 3969 steps, and the shock layer finally
got past its 42 mm standoff. But the state at failure is the *same disease*
`ffs_step_M7` has:

| reported | ceiling from `h₀` | ratio |
|---|---|---|
| `max\|u\| = 2839 m/s` | `√(2h₀) = 1646` | **1.72×** |
| `max(\|u\|+c) = 2839` | 1803 | 1.57× |
| `max c = 989` → T = 2439 K | `T₀ = 1349 K` | **1.81×** |
| `p_min ≈ −150 Pa` | 0 | — |

The advective and acoustic CFLs are *identical*, both 2839, which means at the
worst node `c` was clamped to zero by `sqrt(max(γp/ρ, 0))` — i.e. `p ≤ 0`
there. (`ρ_max = 3.1e-3` is **not** an overshoot: `p_w/(R·T_w) = 460/(287·300)
= 5.3e-3`, so the cold-wall boundary-layer density is physical.)

**`:lkep => true` with `ranocha()` was active for all of this.** So
entropy-conservative flux differencing does not remove the element-scale
speckle — a clean negative for the aliasing hypothesis as that flux tests it.

### What the `mu_dsgs` field says instead

It is blocky and salt-and-pepper: adjacent elements at ~0.8 and ~0. That is
what the default coefficient *is* — element-wise constant, then spread to nodes
by `max` over the sharing elements (`SGS.jl:3017`), a dilation rather than a
smoothing. A CG method cannot absorb a μ that jumps across element interfaces:
`∇·(μ∇q)` then carries a spurious interface forcing proportional to the jump,
which makes element-scale oscillation, which the sensor reads as
under-resolution, which speckles μ further. That loop explains a speckled μ on
top of a speckled velocity far better than anything about the shock does.

Run 4 sets **`:ldsgs_nodal => true`** — the nodal Dao & Nazarov form
(`compute_dsgs_viscosity_nodal!`, reached for `DSGS()` 2D at `rhs.jl:1681`):
ν built at every node from the mass-weighted average of the residual over the
elements containing it, interpolated by the element loop. No element-wise
constant, no interface jump, no broadcast. Same model, in a form CG can carry.

The shock standoff is 42 mm, so **both runs died with the shock layer only
half formed.** The whole difficulty is the *formation* of the normal shock,
not any developed state.

## Where it fails: the stagnation streamline

Run 2 reported per-rank first-bad nodes that cluster tightly — every one within
**±13° of the stagnation streamline**, spanning the wall out through the shock
standoff:

| x | y | wall distance | angle from stagnation |
|---|---|---|---|
| 0.8000 | 0.0000 | **0.0 mm** (the stagnation point itself) | 0.0° |
| 0.8049 | 0.0438 | 0.0 mm | 12.7° |
| 0.7921 | −0.0109 | 8.2 mm | −3.0° |
| 0.7610 | −0.0189 | 39.8 mm (≈ the 42 mm standoff) | −4.5° |
| 0.7570 | 0.0534 | 48.8 mm | 12.4° |
| 0.7452 | 0.0252 | 56.0 mm | 5.6° |

That is a *localised, physical* failure, unlike the `ffs_step_M7` runs where the
whole field went at once. Four things coincide at that point and nowhere else:
the shock is **normal** (its strongest, ρ₂/ρ₁ = 5.46, p₂/p₁ = 57), the
temperature is highest (T₀ = 1345 K), the boundary layer is thinnest, and the
wall is coldest (300 K) — a 1045 K drop across ~2.8 mm.

## Run 1 diagnostics — two findings

**The printed `Viscous CFL` is a diagnostic artifact on this mesh.** `computeCFL`
forms it as `max(ν)` over the whole mesh × `Δt` / `min(Δx)²` over the whole
mesh, and on a graded grid those are at opposite ends:

| cell `h` | `ν_cap = Cmax·(h/5)·ρ∞·(|u|+c)` | local `ν Δt/Δx²` |
|---|---|---|
| 2.2 mm (wall) | 0.080 | **8.1e-3** |
| 20 mm | 0.72 | 9.0e-4 |
| 77 mm (far field) | **2.76** | 2.4e-4 |

`max ν = 3.07` reported by the run is the far-field cap; `min Δx = 3.85e-4` is
the wall cell. `3.0749 × 1.5e-8 / (3.85e-4)² = 0.312`, exactly what was printed
— against a true per-cell maximum of 8.1e-3. **A factor of 38 of artifact.**
The number is right on the near-uniform meshes it was written against
(`ffs_step`) and meaningless here.

**`max ν = 3.07` is the cap in the 77 mm cells**, i.e. the coefficient at its
ceiling in an undisturbed free stream — which looked like the documented
symptom of the missing startup hold, so run 2 turned the hold on *and* dropped
`:dsgs_Cmax` to 0.03. **Both were wrong, and changing two things at once made
the result harder to read.** Both cut dissipation; the run died sooner. Both
reverted.

The reasoning error is worth keeping: the hold protects against a sensor
misreading a smooth *field*, but what kills this case is a violent first few
*steps*, and those are violent however smooth the field is — a 1568 m/s stream
is standing on a no-slip wall at `t = 0` and a normal shock has to form in
front of it. That is `ffs_step`'s own argument for `hold => 0` ("holding ν at
zero there integrates the most violent steps of the run with no dissipation at
all"), and it applies here for the same reason. I was looking at the initial
condition instead of the first steps.

## Status

**Not yet run.** The mesh is generated and verified; no Julia was available in
the environment this was written in, so the deck has not been executed or
parsed. Two things to check on the first run: the `SNAP HIGH-ORDER NODES` banner
and its reported wall distance, and that the first CFL line's viscous number
matches the estimate above.
