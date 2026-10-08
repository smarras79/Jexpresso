# `src/kernel/positivity/` — realizability repair for the compressible Euler state

Two files, deliberately split:

| file | knows about |
|---|---|
| `Positivity.jl` | nothing but arrays and floats — a self-contained module, no Jexpresso types, no MPI, no mesh |
| `positivity_driver.jl` | `params`, `inputs`, MPI, the solution layout — the adapter, and the only Jexpresso-aware part |

**Off unless a deck sets `:lpositivity => true`.** No existing case changes.

## What this is not

It is **not** a positivity-*preserving* scheme. A preserving scheme carries a
proof that the update cannot leave the realizable set under a stated CFL
condition. Two things are worth being explicit about:

- **Zhang–Shu does not apply to CG.** It rescales each element's polynomial
  about its own cell average, which needs elements that can be modified
  independently and a cell average kept positive by a positivity-preserving
  Riemann flux. CG-SEM has neither: nodes are shared, so scaling one element
  breaks C⁰ continuity, and there is no interface flux to make positive.
- **The CG-native answer is invariant-domain preservation with convex
  limiting** — Guermond, Popov & Tomas; the `ryujin` code of Maier &
  Kronbichler. A low-order update with graph (Rusanov-type) viscosity along
  every node–neighbour edge, which is monotone, provably positive and
  carbuncle-resistant without a Riemann solver; then a node-by-node convex
  combination with the high-order update, limited to enforce ρ > 0, internal
  energy > 0 and a local minimum-entropy principle. One precondition is already
  met here: **Jexpresso uses a lumped LGL mass matrix**, which is what that
  limiting is formulated on.

That is the real answer and it is a project, not a file.

## What this is

The bound enforcement without the invariant-domain guarantee — a **local,
bounded, audited repair** that sits under the eventual limiter. It promises
three things:

1. the RHS is never evaluated on a non-realizable state, so one bad node stops
   producing NaN fluxes, NaN sound speeds and `log(p < 0)`;
2. the repair is local and bounded, so a defect cannot spread by arithmetic;
3. **every intervention is counted and reported**, including where the first one
   happened, so it can never quietly rescue a run whose answer is meaningless.

(3) is the point. The open question on the Mach-7 cases is whether the failure
is a handful of nodes at the bow shock that then poison the field, or a field
that is globally garbage. If this engages at ten node-visits and the run
continues sensibly, the first is true and convex limiting is worth building. If
it engages at ten thousand, the second is true and no limiter would have helped.
Either answer is cheaper to buy here than after writing the limiter.

## The repair

Per node, in order:

1. `ρ < ρ_min` → `ρ = ρ_min`. Injects mass; recorded.
2. `p = (γ−1)(ρE − ke) < p_min`, with `ke = |ρu|²/(2ρ)`:
   - **(a)** if `ρE > e_min = p_min/(γ−1)` and `ke > 0`, scale momentum by
     `θ = √((ρE − e_min)/ke) ∈ [0,1)`, which makes `p = p_min` exactly.
     **`ρE` is untouched, so total energy is conserved exactly** — the repair
     converts kinetic energy into internal energy and nothing else. That is
     what a viscous term does: dissipative, entropy-increasing, the right sign.
     A wrong answer here is locally *over-damped*, never locally energised.
   - **(b)** only if `ρE ≤ e_min`: zero the momentum and raise `ρE`. This one
     **does** inject energy, and is counted separately for that reason. If this
     branch fires, the state is badly broken, not marginally so.

**NaN is left alone, deliberately.** `NaN < ρ_min` is false, so a NaN node
passes through and the solver's own non-finite check still fires. Repairing NaN
would turn a detectable failure into a silent one.

## Where it runs

Inside `rhs!`, immediately after the stage state is loaded into `params.uaux`
and before any flux is evaluated — so `user_flux!`, `user_fluxaux!` (which takes
`log(p)` for `ranocha`), the DSGS sensor and the sound speed cannot be handed a
non-realizable state. It makes the `FLUXAUX_FLOOR` guard in the Mach-7 decks
redundant; that guard was this idea, done badly and in the wrong place.

It edits a low-storage RK's own registers. `:lfilter` already does exactly that,
so the pattern is the house one — but it is worth knowing.

Coordinates come from **`mesh.coords[dim, ip]`** (the 3 × npoin array), not the
deprecated per-axis fields, matching the rest of the kernel (`rhs.jl:1131`).
They are used only to report *where* the first repair happened, and the lookup
is size-guarded so a missing coordinate can never cost a repair.

## Scope

Two state layouts, each with its own function in `Positivity.jl`, both requiring
`TOTAL()`, `:energy_equation => "energy"` and the CPU backend:

| layout | | |
|---|---|---|
| 2D/3D CompEuler, `neqs == nsd + 2` exactly | `positivity_limit!` | |
| 2D ideal GLM-MHD, `(ρ, ρu, ρv, ρE, ρw, Bx, By, Bz, ψ)` | `positivity_limit_mhd!` | γ from `:dsgs_gamma` |

The θ-form is still a clear error on the first call rather than a silent wrong
repair: it carries ρθ, which is positive for a different reason.

### The GLM-MHD branch

Recognised from the case's **own `qvars`**, not guessed from `neqs == 9`: a
nine-equation system that is not this one must not be repaired as if it were.
Two things differ structurally, which is why it is a separate function and not a
flag:

* **The momentum slot map is not contiguous.** These cases put ρE in slot 4 and
  the out-of-plane momentum ρw in slot 5, so momentum is `(2, 3, 5)` — the Euler
  loop's `for k = 2:(ien-1)` would scale `Bx` as a momentum component.
* **The internal energy owes the field**: `e = ρE − ke − ½|B|² − ½ψ²`. And
  `½|B|²` is **not reducible** — rescaling `B` would break the discrete
  `∇·B = 0` that the GLM cleaning and the initial condition maintain, which is a
  worse defect than the one being repaired. So the magnetic energy is a *fixed
  charge* against ρE, and that changes which branch is reachable: on a low-β
  problem `½|B|²` can exceed `ρE − e_min` on its own, and then no momentum
  scaling can restore `p` and branch 2b is the only option. On the magnetized jet
  (`β_a = 1e-2`, `½|B|² = 100` against an ambient `ρE` of 102.5) that is the
  normal regime, so **2b firing there is a statement about the field, not
  necessarily about a broken momentum.**

Branch 2a still conserves total energy *exactly* (`ρE`, `B` and `ψ` are all
untouched; only the momentum is scaled), so `∇·B` is unchanged by the repair.

`γ` comes from `:dsgs_gamma`, **not** from `PhysicalConst`: the MHD cases in this
tree run γ = 1.4, 5/3 and 1.05, and air's 1.4 would be silently wrong for two of
the three.

### How small `p_min` can usefully be

`p` is recovered by cancellation against `ρE`, so no repair can place it more
accurately than the spacing of `ρE` itself. The achievable *absolute* accuracy on
`p` is `(γ−1)·eps(ρE)`; the *relative* accuracy on `p_min` is `eps(ρE)/e_min`.
Measured on the magnetized jet, where `ρE = 4.48e5` in the beam:

```
(γ−1)·eps(ρE) = 2.3e-11     <- p cannot be resolved below this at all
p_min         = 1.0e-6      <- 4.3e4 above it: safe
p lands on p_min to 2.3e-5 relative, not to machine precision
```

So keep `p_min` several orders above `(γ−1)·eps(ρE_max)`, and do not expect
`p == p_min` afterwards to better than `eps(ρE)/e_min`. This is a property of the
state, not of the repair: the same limit binds any scheme that carries `ρE` and
recovers `p` from it.

## Settings

| input | default | |
|---|---|---|
| `:lpositivity` | `false` | |
| `:positivity_rho_min` | 0.0 | absolute; **must** be set > 0 when enabled |
| `:positivity_p_min` | 0.0 | absolute; **must** be set > 0 when enabled |
| `:positivity_report` | `true` | first engagement, then once per decade |

The floors are absolute and have **no safe default** — a deck must state them
from its own scales (e.g. 1e-6 of the free-stream ρ and p), so the repair
engages only outside the realizable set and never inside the solution.
`positivity_validate` errors rather than guessing.

## Reading the report

A few node-visits near a shock is the repair doing its job. Engagement growing
without bound, or **anywhere inside the boundary layer**, means the answer is
wrong and the repair is only hiding it — check the reported first-engagement
coordinates against the wall before trusting any wall heat flux from that run.
