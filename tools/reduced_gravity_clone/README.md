# Standalone clone of the ReducedGravity CG-SEM discretisation

`rg2d_clone.jl` is a plain-Julia (stdlib only) re-implementation of what the
Jexpresso CG path does for `problems/ShallowWater/ReducedGravity`, written to
dissect the thin-layer blow-up with each ingredient switchable:

* tensor-product LGL nodes (N = 4), 25 × 25 uniform quads on [0, 2000 km]²
* strong-form inviscid RHS with lumped (inexact) mass and DSS
  (`rhs_el -= ω_i ω_j J (∂F/∂x + ∂G/∂y − S)`, as `_expansion_inviscid!`)
* weak-form Laplacian AV on (H − He, Hu, Hv) with μ = [0, 300, 300] (as `_expansion_visc!`, `AV()`)
* Boyd–Vandeven "erf" modal filter (μ = 0.05) applied in place to q − qe at the
  start of every RHS evaluation (port of `filter.jl`, including `vandeven_modal`)
* free-slip walls as a projection of the state at the start of every RHS
  evaluation (as `build_custom_bcs_dirichlet!`)
* the case's `user_flux.jl` / `user_source.jl` verbatim (desingularised
  velocity, perturbation pressure, Coriolis, smoothstep wind, drag, (1 − wet) σ)
* SSPRK54 with the `"clip"` stage limiter (coefficients checked against
  OrdinaryDiffEq's `SSPRK54ConstantCache`)

It matches an independent Python implementation to all printed digits over the
first two model days. It is serial: ~10 ms per step (≈ 2 model days per
second), so a six-year deck run takes about an hour; run the variants
concurrently.

## Run matrix

```bash
cd tools/reduced_gravity_clone
julia -O3 rg2d_clone.jl tag=deck      dump=deck      > deck.log      &   # the deck as it is
julia -O3 rg2d_clone.jl tag=nofilter  dump=nofilter  filter=0  > nofilter.log &   # :lfilter => false
julia -O3 rg2d_clone.jl tag=noclip    dump=noclip    limiter=0 > noclip.log   &   # no :swe_limiter
julia -O3 rg2d_clone.jl tag=wind20    dump=wind20    hwind=20  > wind20.log   &   # = ReducedGravity_EkmanWind
julia -O3 rg2d_clone.jl tag=wind20lin dump=wind20lin wind=linear hwind_lin=20 > wind20lin.log &
wait
python3 plot_rg2d.py --logs deck.log nofilter.log noclip.log wind20.log wind20lin.log \
                     --dumps deck nofilter noclip wind20 wind20lin --out figs
```

Switches (`key=value`): `filter`, `mu_filter`, `limiter`, `hmin`, `mu`,
`visc=cons|vel` (Laplacian of Hu, or of u weighted by H), `wind=smoothstep|linear|off`,
`hwind` (upper end of the smoothstep wind ramp, 0.1 m in the deck), `hwind_lin`,
`hrelax` (upper end of the thin-film relaxation ramp, 0.1 m), `sigma`,
`sigma_mode=smooth|hard`, `dt`, `nel`, `N`, `tend_days`, `report_days`,
`umax_stop`, `dump`.

Every `report_days` the log prints min H and where it is, max |u| and the
thickness there, max |Hu|, the number of nodes with H < 1 m, the cumulative
number of clipped nodes and the relative mass change. The run stops when the
state is non-finite or max |u| exceeds `umax_stop`, and dumps the last finite
state (`<dump>_prev.dat`) and the stop state (`<dump>_last.dat`) as
`x y H Hu Hv u v` columns; `plot_rg2d.py` maps them and plots |u| against H for
every node against the Ekman line u = τ₀/(ρ f₀ H) and the CFL = 1 line Δx_min/Δt.
