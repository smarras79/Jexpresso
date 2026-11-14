# ReducedGravity_entrain

`problems/ShallowWater/ReducedGravity` with entrainment from the deep layer
(McCreary & Kundu 1988) in `user_source.jl`:

```
S[1] = w_e = (H_e − H)² / (t_e H_e)   for H < H_e = 10 m,  t_e = 1 day;   0 otherwise
```

and `:swe_umax => 2.0` kept as a safety net. Everything else is `include`d
from `ReducedGravity`. **Status: under validation** (clone test from day 1370
to the 6-year end running; the deck as-is dies at day ~1380–1470).

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
volume is not compensated: the upper layer gains the upwelled water, at about
the Ekman-suction rate integrated over the outcrop region (expected ~1 % of the
basin volume per year).

## Run

```bash
mpiexec -n 8 julia --project=. src/Jexpresso.jl ShallowWater ReducedGravity_entrain
```

## What to look for

- No `# SWE limiter diagnostics` (clip) and no `# SWE velocity cap` lines at all.
  If they appear, H_e is too small or t_e too long for this forcing.
- The minimum H settling near ~9 m in the NW subpolar gyre instead of reaching
  the 1 mm film.
- The total volume growth (the entrained water), e.g. from the H field of the
  yearly VTK snapshots.

McCreary, J. P., and P. K. Kundu, 1988: A numerical investigation of the Somali
Current during the Southwest Monsoon. *J. Mar. Res.*, 46, 25–58.
