# rampCaoEtAl2021_M7

Mach-7.7 flow over a 15-degree compression ramp with laminar separation —
Cao, Hao, Klioutchnikov, Olivier & Wen, *J. Fluid Mech.* **912**, A3 (2021).

This is a fork of `problems/CompEuler/rampCaoEtAl2021` carrying the three
things that came out of the Mach-7 debugging campaign (`ffs_step_M7`,
`ffs_step_M7_round`, `shock_circle_M7`). The original is untouched and is
still the paper-faithful reference deck; run both if you want the
comparison — that is the point of the fork.

## What is different, and nothing else

| | change | where |
|---|---|---|
| **(a)** | `:lkep => true`, `:volume_flux => ranocha()` — entropy-conservative two-point flux differencing replaces the aliased collocation volume term | `user_inputs.jl`, and `user_fluxaux!` / `flux_turbo` / `FLUXAUX_FLOOR` in `user_flux.jl` |
| **(b)** | `:lpositivity => true` with floors at 1e-6 of the free stream, plus its self-audit | `user_inputs.jl` |
| **(c)** | `const RAMP_MACH = 7.7`, read by both the free stream and the recovery temperature of the starting profile (it used to be written twice) | `initialize.jl` |

Plus one deprecation fix: `initialize.jl` reads `mesh.coords[dim,ip]`
instead of `mesh.x[ip]` / `mesh.y[ip]`.

**Every DynSGS setting is unchanged** — `:dsgs_sensor => "legacy"`,
`:dsgs_Cmax => 0.1`, `:μ => [1,1,1,1]`, `:ldsgs_nodal => false`. Changing
the sensor at the same time as adding (a) and (b) would make the outcome
unattributable. If this deck survives longer than the original, (a) is why.

## Why (a) is the one aimed at the actual mechanism

`p = (γ-1)(ρE - ρ|u|²/2)` is a difference of two nearly equal numbers, so a
relative error in `ρE` or `ρu` comes out of that subtraction amplified by
`(γ-1)ρE/p = 1 + γ(γ-1)M²/2` — **3.5 at M = 3, 17.5 at M = 7.7**. The same
aliasing error buys five times the pressure error here; the pressure drives
the momentum flux; the loop closes. Flux differencing removes the
aliasing-driven transfer by construction. Artificial viscosity only damps
the symptom, and enough of it to hold a Mach-7.7 shock also drowns the
boundary layer this case exists to resolve.

## Why (b) is a floor and not a scheme

`src/kernel/positivity/` repairs `ρ` and the internal energy in `uaux` once
per RHS call. It is **not** a positivity-preserving scheme — Zhang–Shu
cannot be built on a continuous Galerkin space; the CG-native answer is
invariant-domain preservation with convex limiting (Guermond–Popov–Tomas).
Read `src/kernel/positivity/README.md` before trusting a run it engaged on.

Its value here is as much diagnostic as numerical. The MPI-collective
report on stdout says how many node-visits were repaired out of how many
RHS calls, how much mass and energy that injected, the smallest `ρ` and `p`
seen anywhere, and **where and when the first repair happened**:

```
# POSITIVITY REPAIR ENGAGED — repaired N node-visits in M RHS calls
  [ρ-floor a, momentum-scaled b, energy-RAISED c]
  injected: mass ..., energy ...
  min ρ seen ..., min p seen ...
  first at (x, y, t) = (..., ..., ...)
```

If repairs appear away from the first few startup steps, **the run is
wrong**. Do not raise the floors to silence it — the message is the finding.

## Starting-field fix (measured on a 1-rank run)

A 1-rank run named the first negative pressure in the **whole domain** at
`(x, y) = (7.103784e-4, 1.496259e-4)` — the mid-LGL node of streamwise
element 2, at the top of the first wall element, `y/δ = 0.988`. Two defects
met in that one cell, both in `initialize.jl`, neither in the scheme:

1. **`δ = max(δfloor, δref√(s/sref))` is a kink.** The branches cross at
   `s = 6.0723e-4`, which is **28% along that same element**. A C⁰-but-not-C¹
   field inside a spectral element rings, and at M = 7.7 the ringing reaches
   `p` amplified 17.5×. Now `δ = √(δfloor² + δref²·s/sref)` — same asymptotes,
   no kink.
2. **The floor was one element thick.** `δ = 1.514e-4` against an element
   height of `0.06/401 = 1.496e-4`, so the whole boundary layer was 5 LGL
   nodes. At x = 100 mm it spans 12 elements — the leading edge was **12×
   less resolved** than the rest of the plate. `δfloor` is now `6.0e-4`,
   four element heights, dominant over the first 11 mm (separation is at
   59 mm, untouched).

Still not fixed: the leading edge sits *on* the inflow plane because this
deck removed the paper's 1 mm upstream strip. That jump is inherent to the
no-strip choice.

## Running

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "rampCaoEtAl2021_M7")
```

Meshes are **not** duplicated into this directory; `:gmsh_filename` points at
`../rampCaoEtAl2021/ramp15_uniform.msh` (401 uniform wall-normal elements —
the deck's own measurements say the wall stretching of `ramp15.msh` is what
breaks this case). Regenerate them in the original directory with
`python3 generate_mesh.py`.

`:Δt => 1.0e-9` and `:tend => 2.0e-3` are inherited unchanged, i.e. 2,000,000
steps. Score (a) on **how far it gets**, not on reaching `tend` — the
original dies early, and the first useful number is whether this one dies
later, and whether the positivity report fires before it does.

## Stepping stone at M = 7.0

Set `RAMP_MACH = 7.0` in `initialize.jl`; nothing else needs to change and
the free-stream banner reports back what you chose. Note that everything
else in the case (`δ_ref = 1.38 mm` at `s = 59 mm`, `T_w = 293 K`, the mesh)
is calibrated on 7.7, so a lower Mach is a numerics experiment, not the
paper's flow.
