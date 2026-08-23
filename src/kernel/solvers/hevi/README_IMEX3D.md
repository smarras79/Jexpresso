# IMEX3D — fully three-dimensional implicit acoustics

A time integrator in which the **entire** linear acoustic-gravity subsystem is
implicit — in all three directions — and advection, diffusion, the sources and
the closures stay explicit. Δt stops being limited by sound and becomes limited
by advection, which on an atmospheric mesh is slower by a factor of `1/Mach`.

It is the third of the three schemes in this directory and shares almost all of
its machinery with the first:

| scheme | what is implicit | Δt bounded by | stage solve |
|---|---|---|---|
| `HEVI_ARK` | vertical acoustics | horizontal sound | banded LU per column, direct |
| `SPLIT_EXPLICIT` | nothing (substeps) | advection (outer) | none; inner acoustic steps |
| **`IMEX_ARK`** | **all acoustics, 3D** | **advection** | **preconditioned GMRES** |

The third row is the only one whose stage solve is iterative, and that is the
whole of its cost and the whole of its risk.

The HEVI README names this as one of the two ways past HEVI's ~2× ceiling:

> **a 3D implicit acoustic solve.** Removes the constraint entirely, at the cost
> of a global elliptic solve per stage. A different project.

This is that project. It turned out not to be a *different* one, because the
preconditioner the global solve needs was already sitting in this directory —
HEVI's own column solve, see *The preconditioner* below.

Read *What it costs* before adopting it. The scheme does what it says: on the
`rtb` benchmark it runs stably at 20× the explicit step and 8× HEVI's, where
both of those diverge. Whether that translates into wall-clock depends on how
expensive your `rhs!` is relative to a linear acoustic operator, and on the
`rtb` case — where `rhs!` is a bare inviscid RHS — it does not.

---

## Turning it on

```julia
:ode_solver => IMEX_ARK(:ARS343),   # was CarpenterKennedy2N54() or HEVI_ARK(...)
:Δt         => 1.0,                 # limited by ADVECTION now
:imex_umax  => 20.0,                # the flow speed you expect, in m/s
```

Nothing else in the deck changes. The RHS, the SGS closure, the wall model, the
sponge, the boundary conditions, the diagnostics, restart and the LES statistics
all behave exactly as under an explicit integrator, because the explicit half is
computed as `rhs!(u) − f_imp(u)` rather than re-derived. No source term or
closure can go missing in the split, whatever the deck turns on.

`:imex_umax` is the one setting worth thinking about before the first run — see
*The one number you have to supply* below.

| key | default | meaning |
|---|---|---|
| `:imex_verify` | `true` | setup self-check, including one full Krylov solve |
| `:imex_rtol` | `1e-8` | relative tolerance of the stage solve |
| `:imex_atol` | `1e-30` | absolute floor, for a zero right-hand side |
| `:imex_restart` | `20` | GMRES(m) restart length |
| `:imex_maxiter` | `200` | total iterations before giving up on a stage |
| `:imex_precond` | `:column` | `:column` or `:none` (`:none` is for measuring) |
| `:imex_umax` | measured | largest expected `|v|`, m/s |
| `:imex_mach` | derived | wedge opening; override only to be deliberately conservative |
| `:imex_lateral_walls` | `:auto` | `:auto` `:bc` `:box` `:none` |
| `:imex_wall_flux` | `true` | zero the implicit vertical mass flux at floor and lid |
| `:imex_linearization` | `:RS` | `:RS` (freeze at `qe`) or `:PS` (refresh from the solution) |
| `:imex_update_freq` | `5` | steps between refreshes under `:PS` |
| `:imex_stability_guard` | `:warn` | `:warn` `:error` `:off` |
| `:imex_monitor` | `false` | periodic line reporting the Krylov iteration count |
| `:lcfl_report` | `false` | print the direction-wise stability table at startup |
| `:lcfl_report_every` | `0` | re-print it every N steps (0 = off). Collective; `JEXPRESSO_CFL_REPORT_EVERY` overrides. |
| `:implicit_vdiff` | `false` | make the **vertical** SGS diffusion implicit too — see [below](#implicit-vertical-diffusion). Needs `:imex_linearization => :PS` under a dynamic closure. |

`:ode_adaptive_solver => true` is refused, for the same reason `HEVI_ARK`
refuses it: `rhs!` contains MPI collectives, and an adaptive controller runs
independently on each rank and will disagree about the step size — so the ranks
reach those collectives a different number of times and the run **deadlocks**
rather than failing.

---

## Why ARS343 here, when HEVI refuses it

`ark.jl` says ARS343 is unusable. `imex3d.jl` makes it the default. Both are
right, and the difference is not a matter of taste — it is which set of
`(zE, zI)` pairs the mesh can actually produce.

An ARK pair's amplification factor `R(zE, zI)` takes two arguments: what the
explicit half sees and what the implicit half sees, at the same time, for the
same mode.

**Under HEVI** the acoustics are *split between the halves*: `zE = Δt·c·kx` is
horizontal sound, `zI = Δt·c·kz` is vertical sound, and `kx` and `kz` are
independent. Every corner of the rectangle `[0,zEmax] × [0,zImax]` is a real
mode. ARS343's unstable region sits at `zI ≈ 1` with `zE ≳ 0.1`, squarely
inside it, and it amplifies by ~1e-4 per step — slow enough to present as a
blow-up at a fixed model time that barely moves with Δt and shifts with the rank
count, i.e. as anything but a tableau problem. Hence `ark_joint_amplification`,
and hence ARS232.

**Under IMEX3D** the acoustics are entirely implicit. The explicit half sees
advection, `zE = Δt·|v|·k`, and *the same wavenumber* `k` puts
`zI = Δt·c·k` on the implicit half. So

```
zE / zI = |v| / c = Mach
```

and the reachable set is a **wedge** of opening `Mach`, not a rectangle. Its one
detached piece is the `zI = 0` edge — the entropy and vorticity modes, advected
with no acoustic partner. Reaching `zE = 1` costs `zI = 1/Mach ≥ 10`, and by
there the L-stable implicit half has damped the mode by an order of magnitude.

`ark_wedge_amplification` measures this, and it reverses the ranking:

```
tableau   RHS/step   rectangle              wedge, Mach ≤ 0.3   per RHS
ARS232        3      neutral                zE ≤ 1.15            0.38
ARS443        4      neutral                zE ≤ 1.57            0.39
ARS343        4      amplifies for zE>0.07  zE ≤ 2.83            0.71
```

ARS343 is not marginally better here — it is nearly twice the step per unit
work, and its wedge limit is its *full explicit imaginary radius*: at
atmospheric Mach numbers the acoustic coupling costs it nothing at all.

`test/imex3d/test_wedge_stability.jl` recomputes that table from the tableaux
and asserts the ordering, so changing either default without redoing the
analysis fails a test rather than passing review.

---

## The one number you have to supply

Every other input to the stability estimate is measured at setup. This one
cannot be.

HEVI's explicit half is dominated by horizontal **sound**, which is a property
of the mesh and the base state and is therefore fully known at `t = 0`.
IMEX3D's explicit half is dominated by **advection** — and on a rising bubble,
a spin-up, or any case that starts from rest, `|v| = 0` at `t = 0`. An estimate
made from the initial state would report that any Δt whatsoever is safe, and it
would be wrong by however much the flow subsequently accelerates.

So `:imex_umax` is the deck's statement of the largest flow speed it expects.
The setup report says which of the two it used:

```
│      |v| used: 20.0 m/s (:imex_umax; the state at t=0 shows only 0.0);  Mach 0.0576
```

Set it from the physics, not from the initial condition. Over-estimating costs
step size linearly; under-estimating costs the run.

---

## How it works

### The split

Implicit — the complete linear acoustic-gravity operator, linearised about the
reference state `qe`:

```
∂(δρ)/∂t   = −∇·(δρu, δρv, δρw)
∂(δρu)/∂t  = −∂/∂x [ β · δρθ ]
∂(δρv)/∂t  = −∂/∂y [ β · δρθ ]
∂(δρw)/∂t  = −∂/∂z [ β · δρθ ] − g · δρ            β = ∂p/∂(ρθ) = c²/θ
∂(δρθ)/∂t  = −∇·(θ̄ δρu, θ̄ δρv, θ̄ δρw)
```

Explicit — everything else, as `rhs!(u) − f_imp(u)`.

All five equations are implicit, because the horizontal pressure gradient and
the horizontal mass flux live in the two HEVI leaves out. Equations beyond the
fifth (moisture, tracers) are untouched by the implicit half and pass through
the stage solve unchanged.

The operator is **not written twice**: it is `build_hevi_operator(...; full =
true)`, the same object `acoustic.jl` built for the substepping scheme, with the
ξ and η derivative sweeps put back alongside ζ. Same linearised fluxes, same
reference state, same floor/lid treatment, same DSS and inverse mass matrix.

It acts on `u − qe`, so `f_imp(qe) = 0` **exactly**. The reference state is an
exact steady state of the implicit part, and whatever discrete hydrostatic
imbalance the code already has stays entirely in the explicit part. The split
cannot create or destroy balance.

### Lateral walls are part of the operator

At a free-slip boundary the normal mass flux vanishes. HEVI imposed that at the
floor and the lid (it carries no horizontal flux to zero); this operator carries
all three, so it needs the lateral version too.

Leaving it out is not a rounding error on a walled domain. The implicit half
would let sound run out through the side walls, the difference would land in
`f_exp = rhs! − f_imp` as a stiff residual at exactly the nodes where `rhs!`
projects the normal momentum out, and the horizontal acoustic step-size limit
this whole scheme exists to remove would come back along the boundary with
nothing in the output to say so.

`:imex_lateral_walls => :bc` (what `:auto` picks when the mesh carries boundary
faces) walks exactly the faces `apply_boundary_conditions_dirichlet!` walks,
skipping the same `periodic*` tags, and reads the same face normals — so the
implicit operator's idea of where a wall is cannot drift from the one `rhs!`
enforces. A laterally periodic case gets no lateral walls automatically, because
a periodic face is not a boundary face.

Note that this is equivalent to applying the free-slip projection to the
operator's *input*: `δρu` appears in the continuity and θ fluxes and nowhere
else, so zeroing the normal component there is exactly `A ∘ P`, which is the
correct linearisation of `u ↦ rhs!(P u)`.

### The stage solve

`I − γΔt A` is no longer block-banded per column, so there is no direct
factorisation and no "factorise once". Every stage is a **right-preconditioned
restarted GMRES** over the whole distributed domain (`krylov.jl`).

Right-preconditioned, not left: the iteration runs on `A M⁻¹ y = b` and the
residual it minimises is the residual of the *original* system, so `:imex_rtol`
means what it says. With a preconditioner as strong as a column solve,
`‖M⁻¹(b − Ax)‖` is a poor proxy for the norm the time integrator cares about.

Two things about it are not the textbook choice, both deliberate:

**The inner product is multiplicity-weighted.** A continuous-Galerkin field
keeps a private copy of every node on a rank interface, so a rank-local `dot`
summed and reduced counts those nodes once per holder. That is a bilinear form,
but not the inner product of the space the operator is self-consistent in, and
Arnoldi built on it produces a basis that is not orthogonal in any norm the
residual is measured in. Each node is weighted by `1/multiplicity`, and the
multiplicity is *measured* — `assemble_mpi!` of a field of ones — rather than
derived from `gip2owner`, which would get local periodic duplicates wrong.

**Orthogonalisation is classical Gram-Schmidt, run twice.** Modified
Gram-Schmidt is the textbook answer and costs `j` separate `Allreduce`s at
Arnoldi step `j` — each a full latency round trip, and the reason a Krylov
method that looks fine on one rank falls over at scale. Classical
Gram-Schmidt batches all `j` projections into one reduction, and run twice
("CGS2") it is as stable as MGS to working precision. Three reductions per step,
independent of `j`.

**Every rank runs the same number of iterations.** The loop is driven entirely
by quantities that come out of `MPI.Allreduce`, which MPI guarantees identical
on every process, so the loop conditions are bit-identical across ranks. This is
not a nicety: the operator contains `assemble_mpi!`, so a rank that stopped one
iteration early would leave the others blocked inside a halo exchange forever.
Anything added to that loop that branches on a rank-local quantity reintroduces
the deadlock.

### The preconditioner

This is what makes the scheme affordable, and it was already written.

HEVI's `I − γΔt A_vertical` is exactly a **line-relaxation smoother in the
vertical** — the direction in which an atmospheric mesh is stiffest. It is
already factorised once per γΔt into one banded LU per column, and applying it
is two triangular solves. What is left for the Krylov iteration is the
horizontal acoustic coupling, whose CFL number is smaller by the mesh's
anisotropy ratio.

So the iteration count scales with the **horizontal** acoustic Courant number,
not the vertical one — and the mesh property that makes HEVI weak (little
anisotropy) is the one that makes this preconditioner strong. The two schemes
are complements, not alternatives.

It preconditions `(ρ, ρw, ρθ)` and leaves `ρu, ρv` alone. On a mesh whose ζ
lines are vertical that is not an approximation: `dζ/dx` and `dζ/dy` vanish, no
horizontal momentum enters the ζ flux divergence, and those two rows of the
vertical operator are exactly the identity — inverting them would be inverting
`I`. Dropping them is `(3/5)² = 0.36` of the banded arithmetic and of the factor
storage, and `3/5` of the bytes the gather/scatter moves for split columns. On a
terrain-following mesh they are not the identity, `hevi_choose_vars` says so
from the metrics, and the set widens to all five with no second code path.

Measured (`test/imex3d/runtests_standalone.jl`, identical on 1–5 ranks):

```
isotropic mesh (Δx = Δy = Δz)     7 iterations preconditioned, 8 without
8:1 anisotropic mesh              3 iterations preconditioned, 14 without
```

and on the two real meshes, the iteration count tracks `c/h_x` and is
independent of `h_z` — which is the statement "the vertical is gone" in the
only form that can be measured. See *What it costs*.

The isotropic figure is the honest one and is *why the test uses the
anisotropic mesh for the assertion*: on an isotropic mesh the vertical acoustic
term the preconditioner inverts is no stiffer than the two horizontal ones it
leaves behind, so it saves one iteration in eight — a true number and a useless
test, one that would stay green if the preconditioner were doing nothing at all.

Run with `:imex_precond => :none` to measure what it buys on your own case; the
setup report prints the iteration count either way.

---

## What it costs — read this before choosing Δt

The iteration count is the entire price of the scheme, and it is **not a
constant**. Measured on two different meshes, over an 8× range of Δt:

```
rtb test mesh (h_x = 21.6 m)          rtb_imex production mesh (h_x = 172.7 m)
 Δt     γΔt     CFL_h  iters           Δt     γΔt     CFL_h  iters
 0.05   0.0218   0.35    12            0.35   0.153    0.31    11
 0.10   0.0436   0.70    19            0.50   0.218    0.44    13
 0.20   0.0872   1.40    35            1.00   0.436    0.88    22
 0.40   0.1743   2.80    71            2.00   0.872    1.75    42
                                       4.00   1.744    3.50    86
```

One law fits both, to within a couple of iterations:

```
iterations  ≈  25 · CFL_h ,      CFL_h = γΔt · c / h_x
```

`CFL_h` is the **horizontal** acoustic Courant number — the part the vertical
preconditioner does not remove. That is the design working as intended, and it
is also the bill.

**Why it is linear, and why the restart length does not help.** The operator is
skew, so the eigenvalues of `I − γΔt A` lie on a segment of the vertical line
`Re = 1`, of half-length `CFL_h`. GMRES on a spectrum spread along a line needs
iterations proportional to its length, and there is no clustering for a
restarted Krylov polynomial to exploit — so restarting costs almost nothing.
Measured: 71 iterations at `m = 20` against 68 at `m = 240`. This is the
opposite of the familiar elliptic case, where GMRES(m) stagnates on restart and
a bigger `m` is the fix. **Do not raise `:imex_restart` to chase iterations
here** — it buys ~4% and costs memory linearly. The default of 20 is right.

**The consequence for Δt.** Cost per step is linear in Δt, so cost per unit
*simulated* time approaches a floor instead of falling:

```
cost to reach T  =  T · [ (4·T_rhs + 4·T_op)/Δt  +  3·25·γ·(c/h_x)·(T_op + T_col) ]
                          ^-- falls with Δt          ^-- independent of Δt
```

Raising Δt past the point where the first term stops dominating buys nothing.
On the production mesh that crossover is around Δt = 1 s, which is why the
`rtb_imex` deck sits there rather than at the 1.89 s the stability limit would
allow.

**What one iteration actually costs.** Measured on the `rtb_imex` production
mesh (41 205 points), in units of one HEVI vertical-operator application:

```
A_full   (3D acoustic operator, 5 equations)   2.70 × A_vert
column solve (preconditioner, 3 equations)     1.04 × A_vert
------------------------------------------------------------
one Krylov iteration                           3.74 × A_vert
```

The HEVI profiler measured `A_vert = 0.47 × rhs!` on the `rtb_hevi` case, so on
a case with *that* RHS one Krylov iteration costs **1.8 full `rhs!`
evaluations**. At Δt = 1.0 s the step is then

```
4 rhs! + 4 A_full + 3 × 22 iterations  ≈  125 rhs!-equivalents per step
                                       ≈  125 per simulated second
```

against `CarpenterKennedy2N54`'s 25 per simulated second at Δt = 0.2, and
HEVI's ≈ 17 at Δt = 0.35.

**So on `rtb_imex` this scheme buys step size, not wall-clock, and that is not
a defect to be tuned away — it is what the arithmetic says.** It removes the
acoustic restriction completely: it runs stably at 20× the explicit step and 8×
HEVI's, where both of those diverge, which is what it was built to demonstrate.
It is not the fastest way to integrate a 0.5 K bubble.

**The break-even is a single inequality**, and it is worth checking against your
own case rather than against this one. The cost per simulated second is

```
IMEX3D    ≈  4·T_rhs/Δt  +  75·γ·(c/h_x)·T_iter        T_iter = 3.74 · A_vert
explicit  ≈  5·T_rhs/Δt_explicit
```

`A_vert` depends only on the mesh; `T_rhs` depends on the physics. So IMEX3D
wins exactly when **`rhs!` is expensive relative to the linear acoustic
operator**, and on `rtb_imex` it is not — `rhs!` there is a bare inviscid RHS
with constant-coefficient viscosity, only 2.1 × `A_vert`. Break-even against
explicit on this mesh needs `T_rhs ≳ 16 · A_vert`, i.e. an `rhs!` about seven
times heavier than this case's. A full LES `rhs!` — SGS closure, wall model,
sponge, moisture — is plausibly in that range, but this repository has not
measured it, so treat it as a condition to check and not as a claim.

`JEXPRESSO_HEVI_PROFILE=1` prints the per-step breakdown for any case, which is
the measurement that settles it.

Two other things move the balance:

* **A coarse horizontal mesh.** The floor goes as `c/h_x`. The production
  mesh's `h_x` is 8× the `rtb` test mesh's, and its floor is 8× lower.
* **Sound genuinely binding Δt.** Run `:lcfl_report => true` first. If the
  report names SGS diffusion as the dominant term, none of this helps.

The route to removing the linear growth in iterations altogether is a reduced
(Helmholtz) stage solve — see *Not implemented*.

Per step, against `CarpenterKennedy2N54`'s 5 `rhs!` and no solve:

```
ARS343      4 rhs!
          + 4 applications of the 3D acoustic operator   (f_imp)
          + 3 stage solves, each ~N Krylov iterations, and each iteration is
              1 operator application + 1 column solve
```

`T_op` is roughly three times the vertical operator's (three derivative sweeps
instead of one). `T_col` is two triangular solves per column over `(ρ, ρw, ρθ)`
— the preconditioner deliberately drops `ρu` and `ρv`, whose rows in the
vertical operator are exactly the identity on a mesh with vertical ζ lines;
that is `(3/5)² = 0.36` of the banded arithmetic and of the factor storage,
for free.

Memory, per rank: the Krylov basis is `(m+1)` fields of `npoin × 5` doubles
plus three more for scratch — at `m = 20` and 100 000 local points, about
84 MB. The banded factors are HEVI's own, unchanged.

## The test case

`problems/CompEuler/rtb_imex/` is the rising thermal bubble, and it is the
`rtb_hevi` deck with one switch at the top:

```bash
julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler","rtb_imex")'

DBG_SCHEME=hevi     julia --project=. -e '...'    # the same case under HEVI
DBG_SCHEME=explicit julia --project=. -e '...'    # and fully explicit
```

Same domain, same bubble, same mesh (byte-identical to `rtb_hevi`'s, so the
comparison is a comparison), same viscosity, same output cadence, same `tend` —
each scheme at its own largest safe step, into its own output tree.

**Do not time them at a common Δt.** That measures cost per step, which is not
the question: every implicit scheme here costs more per step and pays for it by
taking a bigger one. Only wall-clock to a fixed physical time answers it, and
`:lstep_heartbeat => true` is already on, reporting `s/step` measured over the
interval since the previous heartbeat — a steady-state rate with the JIT
excluded, which total wall time is not.

```
wall to t = tend  =  (s/step) × tend / Δt
```

To put it in CI, use the repo's own tooling rather than editing the registry by
hand:

```bash
julia --project=. test/generate_ci_ref.jl CompEuler/rtb_imex
```

---

## Testing

Four suites, none of which needs the Jexpresso dependency stack — only MPI,
LinearAlgebra, Printf and OrdinaryDiffEq, and the last one not even those (it
is pure arithmetic on the Butcher tableaux).

```bash
julia --project=. test/imex3d/test_load.jl               # it loads at all
julia --project=. test/imex3d/runtests_standalone.jl     # the machinery
julia --project=. test/hevi/test_vdiffusion.jl           # implicit vertical diffusion
julia --project=. test/imex3d/test_rtb.jl                # the physics
julia --project=. test/imex3d/test_wedge_stability.jl    # the tableau choice
mpiexecjl -n 3 julia --project=. test/imex3d/runtests_standalone.jl
```

**Run the machinery suite on several rank counts.** Two of the things it tests
exist only in parallel and are no-ops on one rank: the multiplicity-weighted
inner product (on one rank every weight is 1) and the preconditioner's
gather/scatter for columns the partition cuts through.

`test_load.jl` earns its place for a second reason beyond the constructor trap
`test/hevi/test_load.jl` documents: IMEX3D *adds methods* to functions the HEVI
files declare — `ark_hooks`, `ark_relinearize!`. A new method is fine; one whose
signature happened to match HEVI's would silently replace it and every HEVI run
would relinearise through the 3D path, and nothing else in either suite would
notice.

`test_rtb.jl` runs the same benchmark, on the same mesh, as
`test/hevi/test_rtb.jl`, and builds **both** schemes so the step sizes come out
of one process on one grid. HEVI is in it deliberately and is not a straw man —
it is stable at Δt = 0.05 on that very mesh, and what kills it at 0.40 is the
horizontal acoustic CFL, which is precisely the term this scheme takes
implicitly and HEVI does not.

It also pins the iteration count. That assertion is not decoration: the
iteration count *is* the cost of the scheme, and a regression in the
preconditioner would leave every other test in this directory green while the
run quietly got several times slower.

Measured, on 1, 3 and 4 ranks (bit-identical on all three):

```
runtests_standalone.jl   25 assertions
test_rtb.jl              14 assertions
test_wedge_stability.jl  14 assertions
test_load.jl              3 assertions
```

and `test/hevi/runtests_standalone.jl` still passes unchanged — the ARK stepper
is shared, so a change here that broke HEVI would show up there.

---

---

## Implicit vertical diffusion

`:implicit_vdiff => true`. Shared with HEVI — same operator, same banded LU,
one implementation. The full account is in
[`README.md`](README.md#implicit-vertical-diffusion); what matters *here* is why
it is the natural companion to this scheme rather than an extra.

**IMEX3D removes the acoustic limit completely, which is exactly what makes the
viscous one visible.** With sound gone in all three directions the explicit
budget is advection plus SGS diffusion, and on a mesh refined in `z` the second
of those goes as `ν/Δz²` — the term the acoustic split cannot touch. Removing it
too is the only combination that leaves nothing vertical in the explicit budget.
`cfl_report` prints that row (`IMEX3D + implicit vertical diffusion`) alongside
the others so the two can be compared before either is switched on.

It is also the case where the extra cost lands hardest, and the report says so:
the `pvars` optimisation — preconditioning only `(ρ, ρw, ρθ)` because the
vertical operator's `ρu`, `ρv` rows are exactly the identity — is *precisely* what
implicit diffusion invalidates. Those rows stop being the identity, `pvars`
widens to all five, and the banded solve goes back to `(5/3)² ≈ 2.8×` the
arithmetic of the three-equation one. That is a real per-iteration cost, paid
against a `Δt` that is no longer capped by the viscous rate.

It requires `:imex_linearization => :PS` under a dynamic SGS closure, and is
refused at setup otherwise — the coefficient is read from the closure, which
returns its molecular floor at `t = 0`, so `:RS` would freeze it there and the
operator would carry no diffusion at all. In the decks `DBG_VDIFF=1` switches
the linearisation with it.

## Not implemented

* **GPU.** The column preconditioner is a LAPACK banded LU on the host and the
  gather/scatter moves host arrays. Refused at setup.
* **`:SOL_VARS_TYPE => PERT()`.** Refused at setup, same as HEVI: the operator
  acts on `u − qe`, which is the deviation only under `TOTAL()`.
* **AMR (`:ladapt`).** Refused at setup: adaptation invalidates the column
  topology the preconditioner is built on, its factorisations and its
  gather/scatter plan.
* **Implicit HORIZONTAL diffusion.** The vertical part is available
  (`:implicit_vdiff => true`, below); the horizontal and cross terms stay
  explicit, and on an LES mesh they are `Δz/Δx` of the vertical one, so they
  are not what binds.
* **An embedded error estimate**, hence no adaptive stepping. See above for why
  that would deadlock rather than merely be inaccurate.
* **A Helmholtz (reduced) stage solve.** This inverts the full
  five-equation acoustic system directly. Eliminating the momenta analytically
  leaves one scalar elliptic equation for `δρθ`, which is what operational NWP
  semi-implicit solvers do, and it is the route to an iteration count that does
  not grow linearly with Δt — an elliptic problem admits a genuinely scalable
  preconditioner (multigrid, or additive Schwarz with this same column solve as
  the vertical smoother, which is what `asm_preconditioner.jl` /
  `Axb_rad_mpi.jl` in this repo are for). That is the single change that would
  turn the step-size win measured here into a wall-clock win on cases where
  `rhs!` is cheap.
* **Nonlinear (Newton) implicit acoustics.** The implicit part here is linear by
  construction. Giraldo et al. (arXiv:2311.11425) time the nonlinear variants
  against the linear one in NUMA and find linear HEVI fastest by 5× over
  NHEVI-LU and 10× over NHEVI-GMRES; that result is an argument for the
  formulation used here, not against it.
