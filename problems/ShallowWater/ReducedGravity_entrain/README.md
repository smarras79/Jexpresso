# ReducedGravity_entrain

`problems/ShallowWater/ReducedGravity` with entrainment from the deep layer
(McCreary & Kundu 1988) in `user_source.jl`:

```
S[1] = w_e = (H_e − H)² / (t_e H_e)   for H < H_e = 10 m,  t_e = 1 day;   0 otherwise
```

with the entrained volume removed uniformly over the basin
(`:swe_conserve_volume => true`, so the layer volume is constant), and
`:swe_umax => 2.0` kept as a safety net. Everything else is `include`d from
`ReducedGravity`.

In the standalone clone, from the deck's own day-1370 state (the deck as-is
dies at day 1434 there), this configuration runs to the 6-year end with no
clipping and no cap events, the layer bottoming out at 4–7 m and max |u| at the
physical 0.7–0.8 m/s; the volume is constant to 1e-10. Without the
compensation the layer gains the entrained water, +4.2 % by day 2100.

## Why

For H₀ = 150 m, g′ = 0.015 m/s², τ₀ = 0.03 N/m² the Sverdrup interior of the
subpolar gyre has to outcrop in the north-west (Parsons-type: it would need an
eastern thickness of ~360 m not to). The Ekman suction there (≈ 1×10⁻⁶ m/s,
≈ 9 cm/day) empties the layer, which then:

- goes down to the 1 mm film, where its momentum is set by its thick neighbours
  (in the clone's term-by-term budget the modal filter is the largest source,
  then the pressure flux) and u = Hu/H blows up — the deck dies;
- with `:swe_umax` (`ReducedGravity_umax`) survives, but the 1 mm clip then
  keeps H ≥ 1 mm by **adding** the water the suction removes, every stage: at
  day 1835 of that run 5.8×10⁸ m³ per step, ≈ 5 % of the basin volume per year
  and growing with the outcrop area.

Entrainment adds that water physically and smoothly instead — it is the
upwelled deep water — and holds the layer near H_e, so the thin-film machinery
(clip, cap, velocity desingularisation) should never engage. The entrained
volume is removed again as a uniform detrainment w̄ = ∫w_e dA / ∫dA, applied by
the kernel after the global assembly of the RHS
(`src/kernel/physics/shallow_water_volume.jl`): ∫H dA is then constant to
round-off. Set `:swe_conserve_volume => false` to let the layer keep the
entrained water instead (about +2 % of its volume per year while it outcrops).

## Run

```bash
mpiexec -n 8 julia --project=. src/Jexpresso.jl ShallowWater ReducedGravity_entrain
```

## What to look for

- No `# SWE limiter diagnostics` (clip) and no `# SWE velocity cap` lines at all.
  If they appear, H_e is too small or t_e too long for this forcing.
- The minimum H settling near ~9 m in the NW subpolar gyre instead of reaching
  the 1 mm film.
- The total volume Σᵢ Mᵢ Hᵢ constant to round-off (with the compensation on):
  the run writes it to `mass.csv` in its output directory every model hour
  (`:lmass_diagnostic`); follow it live with

  ```bash
  python3 tools/plot_mass_live.py output/ShallowWater/ReducedGravity_entrain
  ```

McCreary, J. P., and P. K. Kundu, 1988: A numerical investigation of the Somali
Current during the Southwest Monsoon. *J. Mar. Res.*, 46, 25–58.
