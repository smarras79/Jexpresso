# ReducedGravity_EkmanWind

> **Result: this does NOT fix the blow-up — kept as a negative control.**
> In the standalone clone (`tools/reduced_gravity_clone`) the deck dies at
> day 1470.8 and the deck with the wind ramped out below 20 m dies *earlier*,
> at day 1445.4, with the same thin-layer velocity spikes (17 m/s on a 4 mm
> layer at day 1440, where the local wind is 2×10⁻⁴ of full). The momentum on
> the thin layers does not come from the local wind; see the clone README for
> the diagnosis in progress. Do not use this case as a fix.

A one-variable experiment on the `ReducedGravity` double gyre: the deck, mesh,
time step, artificial viscosity, filter, clipping limiter and thin-film
relaxation are those of `problems/ShallowWater/ReducedGravity` (the files are
`include`d from there); only `user_source.jl` differs, in the thickness scale
over which the wind stress is handed to the active layer:

| | `ReducedGravity` | `ReducedGravity_EkmanWind` |
|---|---|---|
| wind stress fraction | smoothstep in H from 1 mm to **0.1 m** | smoothstep in H from 1 mm to **H_E = 20 m** |
| thin-film momentum relaxation | ramps in below 0.1 m | ramps in below 0.1 m (unchanged) |

## Why this is the hypothesis for the thin-layer blow-up

The wind stress is a force per unit area; applied to a layer of thickness H it
accelerates it at τ/(ρH), unbounded as the layer outcrops. In `ReducedGravity`
a 10 cm layer already carries the full Ekman transport τ₀/(ρf₀) = 0.36 m²/s,
i.e. 3.6 m/s (Fr ≈ 90), and because the same upwelling keeps thinning the layer
at roughly fixed transport, u = Hu/H exceeds Δx_min/Δt = 23 m/s (the CFL = 1
velocity of the 25×25, N = 4 mesh at Δt = 600 s) near H ≈ 1.5 cm. That is an
order of magnitude above the 1 mm clip of the `:swe_limiter`, which is why
clipping the depth does not prevent the blow-up. None of this happens in
`SoliWaveIsland(DSGS)`: there every momentum source vanishes with the layer
(pressure gH∇H, bathymetry g(H−He)∇Hb), so the velocity stays bounded by the
shallow-water characteristics however thin the run-up tongue gets.

With H_E = 20 m the acceleration is bounded by τ₀/(ρ H_E) and the velocity of a
thin layer by τ₀/(ρ f₀ H_E) = 1.8 cm/s (Fr ≈ 0.03): the outcrop then behaves
like an ordinary wet/dry front, which the 1 mm film, the desingularized
velocity and the clip already handle in the island case. Physically, H_E is
the depth of the Ekman layer that carries the stress; where the active layer is
thinner than that, the remainder of the stress goes into the deep, motionless
layer.

## Run

```bash
julia --project=. src/Jexpresso.jl ShallowWater ReducedGravity_EkmanWind
mpiexec -n 8 julia --project=. src/Jexpresso.jl ShallowWater ReducedGravity_EkmanWind
```

Run `ReducedGravity` alongside it with the same settings. Both print the
`:lmin_h_diagnostic` line when H goes negative and the `:swe_limiter_diagnostic`
summary whenever the clip fires; the first clipping events mark the onset of the
outcrop. The expectation is that `ReducedGravity` dies within a few steps of the
first clipping events, with |u| growing past Δx_min/Δt at nodes with H of a few
centimetres, while this case clips a thin film at rest and continues.

A standalone, serial clone of the discretisation that reproduces the deck
ingredient by ingredient (and lets each be switched off) is in
`tools/reduced_gravity_clone/`.
