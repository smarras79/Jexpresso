# rampFreeStreamTest

**This is not a flow. It is a rung ladder of one-property unit tests of the
discretisation.** One property, one number, one flip at a time — each rung is
~500 steps, seconds of wall clock, and the result is a single number: the
**range of `dp_rel`** in ParaView's Information tab at the last frame.

## Rungs 0 and 1 — DONE, both PASSED

| rung | `dp_rel` range | per RHS evaluation | verdict |
|---|---|---|---|
| 0 — pointwise flux | −3.6e-13 … **+5.8e-13** | 2.3e-16 = **1.0 eps** | pure round-off |
| 1 — `:lkep` + `ranocha()` | −1.9e-11 … **+2.5e-11** | 1.0e-14 = **45 eps** | pure round-off |

500 steps × 5 stages = 2500 RHS evaluations. Rung 1's 43× is what flux
differencing *costs*: instead of one pointwise flux per node it evaluates a
two-point flux against every other node in the line, each carrying a
logarithmic mean with a series expansion. A longer arithmetic chain
accumulates more round-off; it does not manufacture a source. A broken FSP
would show as **structure that grows** — a pattern locked to the mesh, orders
of magnitude larger. Both pictures are uniform, with nothing at the outflow,
the kink or the block junction.

The gap still to be explained is **2.7e9×**: the production first repair was
at p = −52.1 Pa, i.e. `dp_rel` = −6.9e-2.

**The whole inviscid path is now cleared** — volume operator, metrics,
outflow, and flux differencing. Rung 1 lives in the base Dict from here on,
so every later rung carries it.

## Rung 0 in detail

| | |
|---|---|
| `dp_rel` range | **−3.6e-13 … +5.8e-13** |
| verdict | **machine zero**, uniform over the whole domain including the outflow plane |

A uniform free stream is an exact Euler solution — every flux is constant, so
`∂F/∂x + ∂G/∂y = 0` *pointwise* and a correct scheme must return machine zero
at every node for ever. It does. So, established:

- The **metric terms are exact** across the 15° kink and the two-block junction.
- The **"impose nothing" outflow manufactures nothing.** The outflow plane is
  as clean as the interior — the missing-boundary-term hypothesis is dead.
- The inviscid volume operator holds a constant state at nop = 4 for 500 steps.

**What rung 0 cannot see, and it matters:** a uniform field is invariant under
*any permutation of the nodes*, so it cannot detect a halo-exchange or
node-indexing bug. Rung 0 clears the operator and the metrics. It does **not**
clear the parallel plumbing — see *The MPI check* below.

> Why there is no shock in the rung-0 picture even though the ramp is there:
> `FSP_WALL = :freestream` **removes the wall**. The plate and the ramp are
> prescribed free stream like every other boundary, so the ramp is only a
> *shape of the domain* and its surface injects uniform flow rather than
> deflecting it. Nothing turns the flow, so nothing compresses. That is the
> point — every piece of physics is removed so that only the operator is
> being measured.

## What is being hunted

In `rampCaoEtAl2021_M7` the positivity repair named the first negative
pressure in the run:

```
first at (x, y) = (0.19659258262890678, 0.08585606812651451)
```

`x = 0.1 + 0.1·cos15°` is the outflow plane to sixteen digits, and that node
sits **60 mm above the wall**. It was corrupted by **step 200**.

| | |
|---|---|
| fastest signal speed | `u + c` = 1725 + 224 = **1949 m/s** |
| time to cross 60 mm | 3.1e-5 s = step **31,000** |
| when it actually happened | step **≤ 200** |

**150× too early for anything physical to have reached it.** So from rung 2
on — where a real wall makes `dp_rel` legitimately large — the question is
never *how big* it is. It is **where it is**. In 500 steps the fastest signal
travels 0.97 mm; anything lighting up at the top of the domain did not travel
there, and whatever put it there is the bug.

**Read the location, not the magnitude.**

## The ladder

Rungs are **cumulative**: rung N means every block up to and including N is
uncommented. Two live in `user_bc.jl` (they cannot be deck keys —
`user_bc_dirichlet!` is not handed `inputs`); the rest are `inputs[:key] =
value` blocks at the bottom of `user_inputs.jl`, so uncommenting can never
collide with the base Dict.

| rung | what it adds | where | result |
|---|---|---|---|
| **0** | bare inviscid operator, no wall | base deck | **PASS 5.8e-13** |
| **1** | `:lkep` + ranocha flux differencing | base deck | **PASS 2.5e-11** |
| **2** | the no-slip isothermal wall | `user_bc.jl` | ? |
| **3** | Sutherland molecular viscosity | `user_inputs.jl` | ? |
| **4** | DynSGS, `:dsgs_norms => "domain"` | `user_inputs.jl` | ? |
| **5** | DynSGS, `:dsgs_norms => "rank"` | `user_inputs.jl` | ? |
| (+) | positivity, for the coordinate | `user_inputs.jl` | optional |

**Rung 1** is the largest untested difference between this case and
production: `rampCaoEtAl2021_M7` runs with `:lkep => true` and rung 0 ran with
it off. Flux differencing does **not** inherit free-stream preservation from
the pointwise scheme — the two-point form preserves a constant state only if
the discrete metric identities hold in the same form it uses them. Rung 0
proved the metrics are exact for the standard operator and proves nothing
about this one.

**Rung 4** is where the causality argument points. `:dsgs_norms => "domain"`
is an `MPI.Allreduce` over the whole domain — a *literal* non-local channel,
the one mechanism in the deck that can carry a number from the boundary layer
to a node 60 mm away in zero time. Run it only with rung 2 active: with the
wall removed the field stays identically equal to `qe`, the departure norm in
the denominator is exactly 0, and you would be measuring a 0/0.

**Rung 5** is the discriminator. If rung 4 corrupts the far field and rung 5
corrupts it *differently* — different place, different time, or not at all —
the normalisation is carrying the corruption.

## The MPI check — no code change, and not optional

A uniform field cannot detect a node-indexing bug, so this has to be done on
a non-uniform field, i.e. on the production case. Run `rampCaoEtAl2021_M7` on
**1 (or 2) ranks and on 32**, and compare:

- the printed CFL numbers at the same `t`
- the positivity counts and its first `(x, y)`
- the step at which it dies

`:dsgs_norms => "domain"` exists precisely so the answer does *not* depend on
rank count. If those numbers differ, the parallel path is implicated and every
Mach-7 result in this campaign is suspect.

## What stays off, and why each is load-bearing

- `:lvisc => false` (rungs 0–2) — viscosity would *damp the error being
  measured*.
- `:lpositivity => false` through rung 5 — the repair exists to hide negative
  pressure; here negative pressure is the signal. Turn it on only as a
  **second** run of a rung that already failed, to harvest the first-repair
  coordinate.
- `:lfilter => false`, `:lsource => false` — same argument.

Mesh, `:nop`, state and `:Δt` are the production ones unchanged, so a result
here transfers directly to `rampCaoEtAl2021_M7`.

## Running

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "rampFreeStreamTest")
```

Colour by **`dp_rel`**; its range on the Information tab is the measurement.
