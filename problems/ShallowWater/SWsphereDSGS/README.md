# SWsphereDSGS — shallow water on the sphere, stabilized by DynSGS

**Status: RUNS.** Same equations, same cubed-sphere grid and the same Galewsky
et al. (2004) initial condition as [`SWsphere`](../SWsphere/README.md). What is
different is the **stabilization**, and that is the entire point of the case:

| | `SWsphere` | `SWsphereDSGS` |
|---|---|---|
| modal filter | **on** (α = 0.05, every step) | **off** |
| viscosity | constant ν = 2·10⁵ m²/s, momentum only | **ν = ν(x, t)**, all four equations |
| tuned constants | ν, α, filter order, cutoff | none — `C1 = 1`, `C2 = 0.5` are the paper's |

DynSGS measures how badly the discrete solution fails to satisfy the PDE and
puts viscosity exactly there. Per element,

```
μ_res|e = C1 Δ²  max_i ‖R_i‖∞,e / ‖q_i − ⟨q_i⟩‖∞,Ω        residual viscosity
μ_max|e = C2 Δ   (‖u‖ + √φ)∞,e                            first-order-upwind cap
ν|e     = max(0, min(μ_max, μ_res))
R_i     = (3qⁿ_i − 4qⁿ⁻¹_i + qⁿ⁻²_i)/(2Δt) − M⁻¹·RHS_i     BDF2 residual
```

with `Δ = Δ_elem/(N+1)`. For an exact solution `R ≡ 0`; for a discrete one it is
small where the solution is resolved and spikes by orders of magnitude at
fronts. The model is *parameter-free* in the sense that matters: the sensor, not
a constant, decides where the dissipation goes.

The model itself, its formulation and the other three equation sets that use it
are documented in [`DSGS.md`](../../../DSGS.md) at the repository root. This
README is about what it is on the **shell**, and what it buys on this test.

## Run it

The cubed-sphere grid ships with the case (`cubed_sphere.msh`: 600 quads, 602
vertices, 10 elements per panel edge), so there is nothing to generate first.

```bash
julia --project=.
julia> using Jexpresso
julia> Jexpresso.run_case("ShallowWater", "SWsphereDSGS")
```

> **Run this in a fresh Julia session after pulling.** `run_case` re-includes
> `run.jl`, `problems/drivers.jl` and the case's `user_*.jl`, but *not* the rest
> of the module. A session opened before you pulled will run the new
> `drivers.jl` against the old module and die with `UndefVarError`.

Output goes to `problems/ShallowWater/SWsphereDSGS/output/`, one `.vtu` per
output time, carrying everything the `SWsphere` files carry —
`phi, phiu, phiv, phiw`, the shell-projected primitives `h, u, v, w`,
`vorticity`, `velocity`, `momentum_normal` — plus, new here:

| field | what it is |
|---|---|
| `mu_dsgs_phi` | ν the model chose for the continuity equation |
| `mu_dsgs_phiu`, `mu_dsgs_phiv`, `mu_dsgs_phiw` | ν for the three momentum components |

These are **piecewise constant per element** by construction (a shared node
takes whichever element writes it last), and they are the coefficients of the
last RK stage of the last completed step. Look at them to answer the question
the model exists to answer: *is the viscosity sitting on the fronts, or is it
smeared over the whole shell?* Do not differentiate them.

The run also prints, every `:ndiagnostics_prints` steps:

```
 #     DynSGS ν(momentum): mean = 7.888e+04  max = 4.145e+05 m²/s  (cap = 2.874e+07)
```

`ν` far below the cap is the healthy regime — it means `μ_res` is governing and
`min(μ_max, μ_res)` is not saturating. **`max ≈ cap` is the failure signal**:
the model has degraded to uniform first-order-upwind dissipation and is no
longer discriminating, and on this grid that also puts ν within a few percent of
the explicit diffusive stability limit. See *The frame trap* below for the one
way that is easy to cause by accident.

## Measured, 6 days, everything else equal

Both decks, same grid (nop = 5, 10 elements per panel edge), same `:cfl => 0.35`,
same Δt = 74.2 s. `SWsphere` numbers are read off its own (10-day) run. The
table is at **t = 5.90 d**, the last diagnostics step both decks print at the
same time; the DynSGS run's own final line at 6.000 d is given underneath.

| at t = 5.90 days | `SWsphere` (filter + ν = 2·10⁵) | **`SWsphereDSGS`** |
|---|---|---|
| max\|ζ\| | 6.79·10⁻⁵ | **9.25·10⁻⁵** |
| max\|ζ − ζ₀\| | 8.88·10⁻⁵ | **1.50·10⁻⁴** |
| δE/E | −8.84·10⁻⁴ | −8.14·10⁻⁴ |
| δmass/mass | −5.6·10⁻¹³ | −1.2·10⁻¹¹ |
| max\|(φu)·x̂\| | 3.5·10⁻¹⁰ | 3.5·10⁻¹⁰ |
| ν (momentum) | 2·10⁵ everywhere, always | mean 8.6·10⁴, max 4.5·10⁵ (1.6% of cap) |

`SWsphereDSGS` at its final step, t = 6.000 d: max\|ζ\| = 9.70·10⁻⁵,
max\|ζ − ζ₀\| = 1.47·10⁻⁴, δE/E = −8.25·10⁻⁴, ν mean 7.9·10⁴ / max 4.1·10⁵.

**Relative vorticity is the field this test is judged on** — the height field
barely moves while the barotropic instability develops, so `h` makes the roll-up
nearly invisible. The initial jet carries max\|ζ₀\| = 1.12·10⁻⁴.

Read the first two rows together. The filtered constant-ν deck has destroyed
39% of the jet's vorticity by day 6 (1.12 → 0.68·10⁻⁴); DynSGS has kept it at
0.93·10⁻⁴, a 17% loss. And the *perturbation* vorticity max\|ζ − ζ₀\|, which is
what the instability actually is, is **68% larger** — the roll-up develops,
rather than being damped while it tries to.

That is the whole argument for a residual sensor. The two decks are not carrying
different *amounts* of dissipation — DynSGS's mean ν is 7.9·10⁴, well **below**
the constant 2·10⁵ it replaces, and its energy decay is slightly smaller too.
They are carrying it in different *places*: the constant-ν deck spends it
uniformly, including on the smooth 90% of the shell where nothing needs it,
where DynSGS concentrates it on the elements whose residual says they need it
and leaves the rest alone.

Mass is conserved an order of magnitude less tightly, and that is expected
rather than a defect: this deck diffuses the continuity equation
(`:ivisc_equations => [1,2,3,4]`, see below) where `SWsphere` does not, so φ
passes through one more weak-form operator every stage. 10⁻¹¹ over 6913 steps is
still accumulated round-off, not a leak — the surface Laplacian integrates to
zero exactly on a closed manifold.

### Cost

Wall clock, serial, on the same machine: **≈ 20 ms/step** against ≈ 15 ms/step
for filter + constant ν, i.e. DynSGS costs roughly a third more per step. The
extra work is one full pass over the domain to form the residual, two reduction
passes for ⟨q⟩ and ‖q − ⟨q⟩‖∞, and a second element loop for the viscous term
(see *Why two element loops* below). Δt is unchanged: the diffusive branch of
the CFL condition is bounded by the model's own `C2·Δ·(‖u‖+c)` cap, which on
this grid leaves it slacker than the gravity-wave limit that actually sets the
step. The banner prints both so it is checkable rather than assumed.

## What is specific to the shell

DynSGS already existed in Jexpresso for 1D CompEuler, 2D CompEuler θ-form and 2D
GLM-MHD (`DSGS.md`). `DSGS_SW()` is a fourth tag rather than a reuse of
`DSGS()`, for four reasons that are not stylistic:

**1. There is no equation of state.** The fast wave of the shallow water system
is the *gravity* wave, `c = √φ = √(gh)` — the same speed `sphere_cfl_dt` uses —
so the cap is `C2·Δ·(‖u‖ + √φ)`. The Euler paths' `c = √(γp/ρ)` has nothing to
evaluate here.

**2. No ρ factor.** ν from the model is **kinematic**, and on the shell it is
applied to the **conservative** variables: `ν∇ₛ²(φu)`, which is literally the
`δν∇²(φu)` of Marras/Kopera/Giraldo Eq. (8b). The Euler-θ and MHD paths need an
extra ρ̄ because their `user_primitives!` hands the viscous operator `(u, v, θ)`
instead. Carrying that factor here would be a units error, and φ is not a
density anyway.

**3. The operator is Laplace-Beltrami, not the flat Laplacian.** The viscous
machinery in `rhs.jl` builds ∇q from the flat 2×2 inverse Jacobian, which has no
slot for the shell's third metric component; dropping it is not an approximation
but a different operator, singular at the equator. So the model's coefficient
feeds `_sphere_visc_el!` in `sphere_rhs.jl` — the surface operator, already
tested against its spherical-harmonic eigenvalues in `test/test_sphere_visc.jl`.

**4. The frame trap**, which is the one that bites.

### The frame trap

The three Cartesian momentum components on a sphere are a property of the
**frame**, not of the flow. A purely **zonal** jet — the Galewsky test, and most
classic shallow-water test cases with it — has `u = u·e_λ` with
`e_λ = (−sinλ, cosλ, 0)`, so its third component **φw is identically zero**,
everywhere and for all time.

Its residual is not. The φw equation carries a large Coriolis source,
`−f(x × φu)|_z ≈ −2Ω sinφ · r cosφ · φu`, which the flux divergence cancels only
to discretization error. Normalize *that* by φw's own spread — exactly zero, so
in practice by whatever floor the code puts under the denominator — and the φw
ratio is inflated by the ratio of the floor to the real momentum scale.
Measured here: **~260×**. `max_i` is then won by φw on every element, μ_res
saturates the cap everywhere, and the model becomes uniform first-order-upwind
dissipation — the exact opposite of what it is for. On this grid that also sits
right at the explicit diffusive stability limit, so the run does not merely
degrade, it **dies at t ≈ 1100 s**.

The fix is to normalize the momentum equation by the norm of the **vector** it
is a component of. The three momentum slots therefore share one denominator,

```
‖φu − ⟨φu⟩‖∞,Ω = max_x |(φu − ⟨φu⟩)(x)|      Euclidean, over the 3-vector
```

which is frame-independent: a zonal jet normalizes φw by the jet's own momentum
scale, the φw ratio drops to the same order as the other three, and `max_i` goes
back to being a statement about the flow. With it, ν sits at 1.4% of the cap and
the run completes. `test/test_sphere_dsgs.jl` keeps this as a regression test.

## The deck, and the two switches that mean something different here

```julia
:lvisc                => true,
:visc_model           => DSGS_SW(),
:ivisc_equations      => [1, 2, 3, 4],
:μ                    => 1.0,
:dsgs_C1              => 1.0,
:dsgs_C2              => 0.5,
:lfilter              => false,
```

**`:μ` is not a viscosity here.** With `DSGS_SW()` the model supplies ν itself,
per element and per RK stage; `:μ` is the **dimensionless per-equation
multiplier** of Marras eq. (10) — how much of the one element coefficient each
equation gets. `1.0` means "as the model says", and it *defaults* to 1.0 rather
than to 0.0. A constant-ν value left in `:μ` by accident would multiply the
model's output by 10⁵, so `build_sphere_viscosity` rejects anything above 100
with an error that says what to use instead.

**`:ivisc_equations => [1,2,3,4]` includes continuity**, which is *not* where the
paper puts its constant ν (it uses `[2,3,4]`, momentum only, so the scheme stays
strictly conservative). The reason is the filter. `SWsphere` runs with the modal
filter **on**, and the filter damps every field including φ. This deck has no
filter, and continuous Galerkin has no upwinding, so with `[2,3,4]` the
geopotential would carry no dissipation of any kind — the failure the error
message at the end of `sphere_time_loop.jl` warns about. To reproduce the
paper's placement exactly, set `[2,3,4]` **and** switch `:lfilter` back on.

`C1` and `C2` are the paper's own values and are deliberately not tuned per
case. Raising `C2` only lifts the ceiling and does nothing while μ_res governs,
which — per the diagnostics line — it does here.

`:dsgs_gamma` and `:dsgs_Prt` exist for the MHD path and are ignored: this
system has neither an equation of state nor an energy equation for them to act
on.

## Why two element loops

`sphere_rhs!` assembles the inviscid RHS over **all** elements, then sizes ν,
then runs a **second** element loop for the viscous term. The order is forced,
not a structuring choice: the residual is `∂q/∂t − M⁻¹·RHS_inviscid`, and RHS is
only the complete inviscid right-hand side once every element has contributed to
it through direct stiffness summation. Folding the viscous term into the first
loop — which is what the constant-ν path does, correctly, because its ν does not
depend on RHS — would size ν from a partially assembled RHS whose value at a
shared node depends on the element numbering.

RHS is still in **weak** form at that point; the `M⁻¹` is applied once at the
end. That is what `compute_dsgs_viscosity!` expects — it multiplies by `Minv`
itself when forming the residual — and the viscous term it then adds goes in in
weak form like everything else.

## The BDF2 history

`∂q/∂t` is not available inside an RK stage, so it comes from the second-order
backward difference on `qⁿ, qⁿ⁻¹, qⁿ⁻²`. Those must be **step** states: the
buffers are rolled once per time step by a gate in `_sphere_ode_rhs!` that fires
on the first stage of each step (`t` sweeps `t + cᵢΔt` within a step). Rolling
per *stage* instead would difference consecutive intermediate states over a full
Δt, which is not an approximation of `∂q/∂t` at all — the defect the flat paths
carried until they got their own `dsgs_qnm1/qnm2` pair (`DSGS.md` §4.4, §7.4).

Both buffers start at the initial state, so the first residual is `−M⁻¹·RHS(q⁰)`
— the actual `∂q/∂t` of the initial condition, which for the balanced Galewsky
jet is nearly zero — rather than the `3q⁰/(2Δt)` an empty history would give.

Δt is published to `sp.Δt[]` *after* the "land exactly on `tend`" adjustment: a
residual built on the requested step rather than the taken one is wrong by that
rounding on every step of the run. `sphere_rhs!` errors rather than proceeding
if `sp.Δt[]` is unset, because an unset Δt does not degrade the model, it
inverts it — R → ∞ everywhere and μ pins at the cap on every element.

## Tests

```bash
julia --project=. test/test_sphere_dsgs.jl
```

Checks, in order: the element length scale; that a discretely exact state gets
**no** viscosity; that the cap is never exceeded and is reached exactly under an
absurd residual; the formula `C1·Δ²·‖R‖/denom` term by term against a
hand-placed residual, including linearity in C1 and **locality** (an element
sharing no node with the residual sees nothing); the per-equation multipliers;
**selectivity** — a tanh front one-seventh of an element wide draws ~34× the
viscosity of a smooth l = 1 harmonic and keeps it *concentrated* (mean/max ≈
0.12, so it is not simply turning itself up everywhere); the zonal
non-saturation regression above; and the deck switches, including the guard on a
constant-ν `:μ`.

`test/test_sphere_visc.jl` covers the surface Laplacian this model feeds.

## References

- S. Marras, M. Nazarov, F. X. Giraldo, *Stabilized high-order Galerkin methods
  based on a parameter-free dynamic SGS model for LES*, J. Comput. Phys. **301**
  (2015) 77–101.
- M. Nazarov, J. Hoffman, *Residual-based artificial viscosity for simulation of
  turbulent compressible flow using adaptive finite element methods*, Int. J.
  Numer. Meth. Fluids **71** (2013) 339–357.
- S. Marras, M. A. Kopera, F. X. Giraldo, *Simulation of shallow-water jets with
  a unified element-based continuous/discontinuous Galerkin model with grid
  flexibility on the sphere*, Q. J. R. Meteorol. Soc. **141** (2015) 1727–1739.
- J. Galewsky, R. K. Scott, L. M. Polvani, *An initial-value problem for testing
  numerical models of the global shallow-water equations*, Tellus **56A** (2004)
  429–440.
