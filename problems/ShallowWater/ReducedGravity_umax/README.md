# ReducedGravity_umax

`problems/ShallowWater/ReducedGravity` with one added deck line:

```julia
:swe_umax => 2.0,   # m/s
```

which makes the existing `:swe_limiter => "clip"` stage limiter also rescale
the momentum wherever |(Hu, Hv)| > 2 m/s × H, after the 1 mm depth floor. H is
never changed by the cap, so it conserves mass (the clip itself does not). All
other files are `include`d from `ReducedGravity`.

## Why

In the standalone clone of this deck (`tools/reduced_gravity_clone`) the run
does not die from H < 0: the clip keeps H ≥ 1 mm. It dies from u = Hu/H on
layers just above the floor once the NW subpolar gyre outcrops (day ≈ 1377):
intermittent spikes of 15 m/s on a 7.5 mm layer (day 1405) and 129 m/s on a
2.5 mm layer (day 1450), then a non-finite state at day 1470.8. The CFL velocity
of the mesh is Δx_min/Δt = 13.8 km / 600 s = 23 m/s. The cap holds every node
below 2 m/s, about 2.5× the fastest physical flow in the run (≈ 0.8 m/s in the
western boundary current), so it should act only on the thin outcrop layers.

## Run

```bash
mpiexec -n 8 julia --project=. src/Jexpresso.jl ShallowWater ReducedGravity_umax
```

## What to look for in the log

- `# SWE velocity cap: step=… capped node-visits=… max|u| removed=…` — when it
  starts (expected near day 1375–1400, when the outcrop forms), how often it
  fires and how large the removed velocity is. If it fires well before the
  outcrop, or in the open gyre, 2 m/s is too low for this flow.
- `# SWE limiter diagnostics: …` — the 1 mm clip, with the mass it adds.
- Whether the run passes day ≈ 1470 (the clone deck's blow-up) and reaches the
  6-year end time.
