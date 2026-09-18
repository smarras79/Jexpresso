# ReducedGravityDSGS

This case uses the wind-driven reduced-gravity DoubleGyre configuration from
`problems/ShallowWater/ReducedGravity`, replacing its constant artificial
viscosity with the `DSGS_SW()` stabilization configured in
`problems/ShallowWater/SoliWaveIslandDSGS`.

The physics, initial condition, boundary conditions, primitive variables, and
mesh are loaded directly from the ReducedGravity sibling. Only
`user_inputs.jl` is case-specific.

The DynSGS wave-speed cap uses the reduced gravity `g' = 0.015 m/s²`, and the
wet/dry regularization uses `Hmin = 10⁻³ m`, matching the ReducedGravity flux.

Run from the repository root with:

```bash
julia --project=. src/Jexpresso.jl ShallowWater ReducedGravityDSGS
```

The production deck integrates six 360-day model years with a 600-second
timestep, so use a shortened copy of the deck for a quick smoke test.
