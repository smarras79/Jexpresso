# HEVI — horizontally-explicit, vertically-implicit time integration

A time integrator for compressible LES on grids whose vertical spacing is
finer than their horizontal spacing. It takes the vertical acoustic terms
implicitly — a direct banded solve inside each vertical column — and leaves
everything else explicit.

---

## Before you switch anything on: find out what actually limits Δt

Add to the deck:

```julia
:lcfl_report => true,
```

and run one step. You get a table of the per-direction stability limits,
measured from the real LGL node spacing and the current state — not from the
`Δ/N²` rule of thumb, which for `nop = 4` is pessimistic by a factor of 2.8
(the smallest LGL gap is `0.173·h`, not `h/16`).

The last block of the report names the dominant term and states what each
scheme would buy. Four outcomes, three of which are *not* "use HEVI":

| dominant term | what to do |
|---|---|
| vertical acoustics | HEVI. This is what it is for. |
| SGS diffusion | Set **`:implicit_vdiff => true`** — the acoustics are not what binds. The parabolic limit goes as `Δz²`, so on a mesh refined to resolve a surface layer it easily beats the acoustic one, and HEVI's acoustic split alone buys almost nothing. See [Implicit vertical diffusion](#implicit-vertical-diffusion) below. |
| horizontal acoustics | HEVI cannot help. Acoustic substepping (`substep.jl`) or the fully 3D implicit solve (`README_IMEX3D.md`) removes it. |
| advection | Already at the floor. No time-integration trick helps; the mesh is the only lever. |

The report costs one RHS evaluation, and that evaluation is also what fills
the SGS viscosity cache — without it `ν_t` is zero everywhere and the viscous
limits read as infinite, which is the one answer guaranteed to be wrong.

---

## What it is worth on the LESICP2 target mesh

Be clear-eyed about this. HEVI removes exactly one term from the explicit
eigenvalue budget — the vertical acoustic one — so the gain is bounded by the
grid's acoustic anisotropy, and nothing else.

For 128 × 128 × 120 elements over 10240 × 10240 × 5000 m at `nop = 4`:

```
element size      Δx = Δy = 80 m        Δz = 41.7 m
LGL min gap       Δx_eff = 13.8 m       Δz_eff = 7.2 m
anisotropy                              1.92
```

Same polynomial order in all three directions, so the LGL clustering factor
cancels and the ratio is just `Δx/Δz`. Removing the vertical acoustic term
halves the acoustic rate — and it is not the 10–30× the same split buys on a
mesh with `Δx ≫ Δz`; that ratio does not exist here.

**But halving the rate does not double Δt, because the tableau gives most of it
back.** This section previously claimed "about two on Δt" and "2.0–2.4× end to
end". Those figures were computed for `ARS343`, whose explicit half has an
imaginary radius of 2.83. The joint-stability analysis later ruled ARS343 out
for HEVI (see *Choosing a tableau* in `ark.jl`) and the default became
`ARS232`, whose radius is 1.732 — and the step size went with it. Measured on
this mesh with `ark_joint_dt_max`:

```
                      rate seen    imag. radius   RHS/step   Δt
CarpenterKennedy2N54   145.6 /s        3.34          5        0.0229 s
HEVI, ARS232            76.3 /s        1.732         3        0.0202 s
HEVI, ARS443            76.3 /s        1.570         4        0.0206 s
```

**HEVI's step on this mesh is about 0.9× the explicit one, not 2×.** Halving
the rate and roughly halving the stability radius very nearly cancel.

What is left is the cheaper step — 3 RHS evaluations plus 3 operator
applications and 2 column solves, against `CarpenterKennedy2N54`'s 5 RHS — and
that is where the whole gain now lives:

```
cost per simulated second, in rhs!-equivalents, a = A_vert/T_rhs
  a       explicit    HEVI     ratio
  0.47      223        268     0.83x   (rtb-like: a bare inviscid rhs!)
  0.25      223        212     1.05x
  0.12      223        179     1.25x   (plausible for SMAG + wall model)
  0.06      223        164     1.36x
```

So **roughly 1.0–1.35× end to end**, not 2.0–2.4×, and it depends on how
expensive your `rhs!` is relative to the vertical acoustic operator. HEVI is
*slower than explicit* on a case with a cheap RHS. `JEXPRESSO_HEVI_PROFILE=1`
prints `t_rhs` and `t_fimp` per call, so `a` is one run away.

> **This number is a property of the MESH, not of HEVI.** The contrast is
> `problems/CompEuler/LESICP2-coarse-hevi`, whose setup report measures
>
> ```
> node spacing    xi 138.1 m   eta 138.1 m   zeta 6.907 m     (20:1)
> joint IMEX      explicit half 5.1 1/s, implicit half 34.6 1/s
>                 Δt_max = 0.2410 s, recommended 0.1687 (70%)
> ```
>
> — HEVI at Δt = 0.16 against the explicit deck's 0.01, i.e. **16×**. Same
> split, same tableau, same code; the difference is that its elements are
> 800 m wide against the target mesh's 80 m, so its anisotropy is 20:1 rather
> than 1.92:1. **HEVI's gain is the anisotropy and nothing else**, so refining
> horizontally destroys it while refining vertically creates it. Read the
> `dominant term` line of your own CFL report before believing any figure in
> this section.
>
> Note also that `κ` — the ratio of the measured operator spectrum to
> `c/Δz_min`, which converts node spacings into rates — is itself
> mesh-dependent: 1.44 on `rtb_hevi` (uniform), 0.69 on that stretched coarse
> mesh, where only a few nodes carry the smallest gap. The table above uses
> the uniform-mesh value, which is the right one for the unstretched
> 128×128×120 target but not for a stretched grid. The setup report measures
> it for whatever mesh you hand it.

> **A caveat on the explicit baseline.** `CarpenterKennedy2N54`'s imaginary
> radius is quoted as 3.34 in `ark.jl` and 3.95 in `docs/hevi/hevi.md`, and the
> ratio above moves from 1.25× to 1.03× between them. The LESICP2 deck runs
> Δt = 0.02 against a computed rate of 145.6 /s, i.e. an effective 2.9 — so it
> has margin and the true baseline is somewhere in that band. Treat the HEVI
> gain on an isotropic-ish mesh as "order 1, measure it", not as a number.

### Getting past 2× on this grid

The remaining fast term is horizontal sound, and only two things remove it:

* **acoustic substepping on top of HEVI** (Klemp–Wilhelmson / MIS). Outer RK3
  step at the advective limit, inner substeps at the *horizontal* acoustic
  limit with the vertical still implicit. Estimated 8–10× on this mesh. The
  column solver in this directory is exactly the piece that makes the
  vertically-implicit half of such a substep cheap; the substepping loop
  itself is not written yet.

  Note that substepping *without* HEVI buys much less: the substep would then
  be limited by the vertical acoustic CFL, i.e. by the same anisotropy HEVI
  removes. The two compose; neither replaces the other.

* **a 3D implicit acoustic solve.** Removes the constraint entirely, at the
  cost of a global elliptic solve per stage. **This now exists** — see
  `README_IMEX3D.md` and `imex3d.jl`. It turned out not to be a different
  project, because the column solver in this directory is exactly the
  preconditioner that global solve needs: the Krylov iteration count then
  scales with the *horizontal* acoustic Courant number, and the two schemes
  compose rather than compete.

Cheaper things worth knowing about, not implemented here: reduced sound speed
(RSST) scales the continuity equation's time derivative by `1/ξ²` and buys `ξ`
directly for a few lines of code — a good stopgap and a good way to confirm
acoustics really are the problem, at the cost of distorting the acoustic
response; and vertical grid stretching or a lower polynomial order in `z`
attacks the same ratio from the mesh side and compounds with HEVI.

---

## Turning it on

```julia
:ode_solver           => HEVI_ARK(:ARS232),   # was CarpenterKennedy2N54()
:Δt                   => 0.02,                # measure it -- see above
:ode_adaptive_solver  => false,               # already the default
```

Nothing else in the deck changes. The RHS, the SGS closure, the wall model,
the sponge, the boundary conditions, the diagnostics, restart and the LES
statistics all behave exactly as they do under an explicit integrator,
because the explicit half is computed as `rhs!(u) − f_imp(u)` rather than
re-derived. No source term or closure can go missing in the split, whatever
the deck turns on.

Optional knobs:

| key | default | meaning |
|---|---|---|
| `:hevi_verify` | `true` | setup self-check; cheap, leave it on |
| `:hevi_wall_flux` | `true` | zero the implicit vertical mass flux at floor and lid |
| `:hevi_vars` | auto | override the implicit variable set |
| `:lcfl_report` | `false` | print the stability table at startup |
| `:lcfl_report_every` | `0` | re-print it every N steps (0 = off). The startup table is taken on a laminar sounding, where `ν_t ≈ 0`; this is how the viscous row is watched as the boundary layer spins up. Collective. `JEXPRESSO_CFL_REPORT_EVERY` overrides. |
| `:implicit_vdiff` | `false` | make the **vertical** SGS diffusion implicit too — see [below](#implicit-vertical-diffusion). Needs `:hevi_linearization => :PS` under a dynamic closure. |
| `:lprecompile_warmup` | `true` | one throw-away step before the real solve |

`:ode_adaptive_solver => true` is refused, and not only because the tableaux
carry no embedded error estimate. `rhs!` contains MPI collectives, so an
adaptive controller — which runs independently on each rank and will disagree
about the step size — reaches those collectives a different number of times
per rank and the run deadlocks rather than failing. (This is a property of
the code as it stands, not of HEVI: the same is true of the explicit solvers,
which is why `:ode_adaptive_solver` already defaults to `false`.)

Tableaux: `:ARS111` `:ARS121` `:ARS122` `:ARS222` `:ARS232` (default)
`:ARS343` `:ARS443`.

**`:ARS343` is refused for HEVI** by the startup guard — see *Choosing a
tableau* in `ark.jl`. It is the default for `IMEX_ARK`, which is not a
contradiction: the two splits reach different sets of `(zE, zI)`. See
`README_IMEX3D.md`.

**`:ARS222` is a trap.** It is the cheapest per step — 2 RHS — and its
explicit half has the stability polynomial of Heun's method, which is
unstable everywhere on the imaginary axis. Fine for a diffusion-dominated
test, wrong for advection. `ark_imaginary_radius` measures this for any
tableau and the setup report prints it.

---

## "HEVI refuses to start: tableau ARS232 amplifies by ..."

The guard is doing its job, but the message names the tableau and the tableau
is usually not the cause. **Read the `explicit half = acoustic + advective +
viscous` line in the setup report before changing anything.** The three terms
add into one rate and want three different remedies:

| binding term | what fixes it | what does NOT |
|---|---|---|
| acoustic | smaller Δt, or the 3D IMEX (`:limex`) | the closure |
| advective | smaller Δt | anything else |
| **viscous** | the LES filter width, or `:implicit_vdiff` | **smaller Δt, another tableau** |

The viscous one is the trap, because `2ν/h_z²` does not scale with Δt at all —
halving Δt to buy stability costs half the throughput and leaves the term
exactly where it was.

### The case this was written from

`LESICP2-coarse-imex`, 10 × 1 × 60 elements over 6400 × 1600 × 5000 m, first z
element 40 m, ARS232 at Δt = 0.2. The guard reported `max|R| = 49.9941` and
`Δt_max = 0.047 s`. Inverting those two numbers for the rates gives

```
rate_imp = 35.0 1/s      against c/Δz_min = 340/6.91 = 49.2   ->  κ = 0.71
rate_exp = 32.3 1/s      of which  acoustic 3.1 + viscous 29.2
```

κ = 0.71 is an ordinary spectral correction, so the vertical grading was *not*
inflating the estimate — that was the first theory and it was wrong. The
horizontal acoustic rate is 3.1 1/s, i.e. a horizontal acoustic Courant number
of 0.6 at that Δt: the acoustics were never close to the limit. **90 % of the
refusal was SGS diffusion.**

And that diffusion was an artefact. `ν_t` goes as the filter width **squared**,
and `:les_filter_width => :max` (the default) takes the width from the
*coarsest* element edge. On a 10 × 1 × 60 mesh the elements are 640 × 1600 ×
40 m and the coarsest edge is the **one element across y** — a direction the
case never meant to resolve. 1600 m against a 40 m vertical gap is 40× on the
width and **1600× on ν_t**, which is what turned a ~0.02 1/s term into 29 1/s.

`compute_element_size_driver` warns about this at setup (`LES filter width :max
gives ... a factor N on the filter and N² on the eddy viscosity`). On a
quasi-2D mesh set `:les_filter_width => :geometric`. `:max` is for *genuine*
anisotropy — on the 64 × 64 × 60 production mesh it is 1.59× on the width and
2.5× on ν_t, which is a modelling choice someone can defend.

Two related notes:

* These rates are **not** taken on a zero-viscosity state. `hevi_verify_physics`
  calls `rhs!` a few lines before the guard runs, so `sgs.μ_turb` is populated
  by the time `cfl_limits` reads it. (The CFL *report* is a different call and
  does carry the "taken at t = 0" caveat.)
* With `:implicit_vdiff => true` the vertical viscous term is dropped from the
  explicit half — `d/dz(μ d/dz)` is in the column operator, so charging it
  would refuse a Δt the scheme can take. It is a real (decaying) eigenvalue, so
  it is not added to the implicit half either: that solve takes it
  unconditionally.

A misspelt key here is a silent no-op with a squared cost — `:les_filter_widthc
=> :max` sat in two decks in this repo reading exactly like a setting and doing
nothing. `compute_element_size_driver` now rejects near-misses of the name.

## When it looks like it hangs (usually it doesn't)

**Check this first.** The integrator warm-up prints

```
# Integrator warm-up with real callbacks (PATIENCE: ONLY DONE ON 1st RUN!) .........
```

and then **never prints anything when it finishes** — the `DONE` line is
commented out in `TimeIntegrators.jl`. So the run goes silent at that message
even though it has moved on to the production solve. If the case's first
diagnostic output is far away (at Δt = 0.05, `t = 100 s` is 2000 steps), that
silence can last many minutes and looks exactly like a stall.

`:lstep_heartbeat => true` prints per-step progress from step 1 and removes the
ambiguity entirely. The `rtb_hevi` deck sets it for this reason.

Interrupting with Ctrl-C and reading the top frame of the backtrace is the
fastest test. If it names something under `_build_rhs!` — `_expansion_visc!`,
`inviscid_rhs_el!` — the run is executing normally and none of it is HEVI; an
explicit integrator would be sitting in exactly the same place.



The first `solve` is also where Julia compiles `perform_step!`, the whole
implicit chain, and the integrator specialised on the real `CallbackSet` and
on Jexpresso's ~90-field `params`. On a full-size case that takes minutes and
is indistinguishable from a deadlock. Set

```bash
JEXPRESSO_HEVI_TRACE=1 julia --project=. ...
```

and every phase announces itself with a wall-clock stamp:

```
[hevi +  11.23 s] initialize!: first rhs! evaluation (this is where the JIT ...)
[hevi +  11.29 s] step 1 begin: t=0.0 dt=0.05
[hevi +  11.35 s]   stage 2: entering column solve, gdt=0.0217933
[hevi +  11.35 s]     gather / 165 banded solves / scatter
[hevi +  11.35 s]   stage 2: rhs! done
```

Read it like this:

* **stalls before `step 1 begin`** — compiling, or stuck in `rhs!` itself
  (which is not HEVI: an explicit run would stall in the same place);
* **stalls between `entering column solve` and `column solve done`** — the
  implicit solve, i.e. genuinely this code;
* **`refactorising` on every step** — γΔt is not matching what setup computed,
  so every step pays a full rebuild. That is a bug, not a tuning problem;
* **steps completing but slowly** — it works; it is just a big problem.

Two ways to bisect further: `:lprecompile_warmup => false` skips the
throw-away warm-up step and goes straight to the production solve, and
`:lstep_heartbeat => true` prints per-step progress once the loop starts.
Interrupting with Ctrl-C and reading the backtrace also settles it in one
shot — the top frame names the function.

## How it works

### The split

Implicit: the terms carrying the vertical sound wave, linearised about the
reference state `qe`.

```
∂(δρ)/∂t   = −∂/∂z [ δ(ρw) ]
∂(δρw)/∂t  = −∂/∂z [ β · δ(ρθ) ] − g · δρ         β = ∂p/∂(ρθ) = c²/θ
∂(δρθ)/∂t  = −∂/∂z [ θ̄ · δ(ρw) ]
```

Explicit: everything else, as `rhs!(u) − f_imp(u)`.

Vertical *advection* stays explicit — `|w|/Δz` is about 1/800 of `c/Δz` here,
so making it implicit would buy nothing and would cost a refactorisation
every step.

The operator acts on `u − qe`, so `f_imp(qe) = 0` exactly. The reference state
is an exact steady state of the implicit part, and whatever discrete
hydrostatic imbalance the code already has stays entirely in the explicit
part, untouched. The split cannot create or destroy balance.

Because the operator is linear, each stage is **one direct banded solve** — no
Newton iteration, no convergence tuning, no `W` assembly, and **no
preconditioner**: there is no Krylov iteration to precondition. That is the
point of restricting the implicit part to the column-local linear acoustic
operator rather than taking the full vertical flux divergence implicitly.
`I − γΔt A` is block-banded per column, it is factorised once, and every stage
of every step is two triangular solves.

The repo's `asm_preconditioner.jl` / `Axb_rad_mpi.jl` become relevant one step
further out: a **3D** implicit acoustic solve — the other route past the ~2×
ceiling — is a genuine elliptic problem, needs a Krylov method, and needs
preconditioning. The column solver built here is exactly the line-relaxation
smoother such a preconditioner wants in the vertical, with additive Schwarz
handling the horizontal. Nothing in this directory needs it yet. That is also why
`KenCarp4` and friends out of OrdinaryDiffEq are the wrong tool: they assume a
nonlinear implicit part and would have to be taught about a Jacobian that is
block-banded per column and distributed over MPI.

### Files

Only the first six are HEVI proper. `acoustic.jl` / `substep.jl` are the
acoustic-substepping scheme and `krylov.jl` / `imex3d.jl` are the fully 3D
implicit one (`README_IMEX3D.md`); both are built on the operator and column
machinery below, which is why they live here.

| file | what it does |
|---|---|
| `cfl_diagnostics.jl` | direction-wise stability limits and what each scheme would buy |
| `columns.jl` | global column identification, and the MPI gather/scatter for columns that straddle ranks |
| `operator.jl` | the linear vertical acoustic-gravity operator |
| `factorize.jl` | one banded LU per column; the matrix is *probed out of the operator*, not hand-derived |
| `ark.jl` | ARS IMEX tableaux and the ARK stepper, shared with IMEX3D |
| `hevi.jl` | setup, self-check, report |
| `acoustic.jl` | the full 3D acoustic operator, shared with IMEX3D |
| `substep.jl` | the acoustic-substepping integrator |
| `krylov.jl` | distributed GMRES (IMEX3D) |
| `imex3d.jl` | fully 3D implicit acoustics: setup, stage solve, report |

### The matrix is extracted, not derived

`A` couples a node only to nodes within `±(ngl−1)` levels of it in the same
column. Colouring the levels modulo `S = 2(ngl−1)+1` makes a probe vector hit
exactly one band entry per row, so `S · nimp` applications of the operator
recover the matrix exactly.

This is not a shortcut, it removes a failure mode. A hand-derived matrix that
disagreed with the operator — in the ζ collocation weights, a metric term, the
DSS sum, the inverse mass matrix, the wall flux — would **not** break any
convergence test, because the split stays consistent whatever the implicit
matrix is. It would quietly degrade the implicit solve into a mediocre
preconditioner, and the vertical acoustic CFL would come back with nothing in
the output to say why. Probing makes the matrix *be* the operator, by
construction, and the setup self-check reports the agreement as a number.

### One factorisation, mostly

`W = I − γΔt A` is built and factorised once at setup, and every stage of every
step is then two triangular solves. It is rebuilt when γΔt changes, and that
does happen in a normal run: OrdinaryDiffEq shortens the step that lands on a
`tstop`, and Jexpresso passes every diagnostic and LES-statistics time as a
tstop — on the LESICP2 deck that is about 180 of them, each costing two
changes of γΔt (onto the short step and back).

Each rebuild is `S · nimp` applications of the cheap operator plus one banded
LU per column, together roughly the cost of one time step, and it writes into
the storage already allocated. About 360 rebuilds in a 200,000-step run is
~0.2% overhead. Solving with a stale `W` instead would be cheaper and wrong,
and would present as an unexplained stability limit rather than as a mistake.

### Footprint at production scale

For 128 × 128 × 120 elements at `nop = 4`, with the implicit set `(ρ, ρw, ρθ)`:

```
levels per column      481
global node columns    512 × 512 = 262 144
band half-width        kl = ku = nimp·ngl − 1 = 14
per column             43 × 1443 doubles ≈ 496 KB
at 1728 ranks          152 columns/rank ≈ 77 MB/rank ≈ 4.9 GB/node at 64 ranks
```

Comfortable against the 7 GB/rank the runs already request. The five global
arrays indexed by column id (`zbot`, `ztop`, `score`, `owner`) are ~2 MB each
and are reduced once, at setup.

### Why HEVI still works under a p4est partition

The column solve is only cheap if a column lives on one rank. Under p4est it
often does not: leaves are ordered by (tree, Morton-within-tree) and the
partition cuts contiguous ranges of that order, so a rank's subdomain is a run
of whole trees plus a Morton sub-cube at each end — and a sub-cube spans only
part of a tree's vertical extent.

Switching to `:lxy_partition` is not an option for the large runs: that path
reads the full 1.97M-element mesh on rank 0, which is what ran the nodes out
of memory and is the reason these cases go through p4est at all.

So this implementation does not require column-local partitioning. Columns
that are rank-local are solved in place with no communication; columns that
straddle ranks are gathered onto an owner (the rank holding the most levels of
it), solved, and scattered back. The setup report states what fraction had to
be split and how many bytes move per stage.

If you *want* full locality, `tools/pick_nranks.jl --columns` works out which
rank counts give it. For the LESICP2 target the answer is a mesh change:

```
coarse mesh   refine   column-local rank counts   elements/rank
16×16×15         3     ≤ 256                      7680 at 256   (too heavy)
32×32×30         2     ≤ 1024                     1920 at 1024  ← usable
64×64×60         1     ≤ 4096                      960 at 2048
```

The current 16×16×15 coarse mesh caps column-local running at 256 ranks with
7680 elements each. A 32×32×30 coarse mesh refined twice reaches the same
128×128×120 grid, gives 1024 fully column-local ranks at 1920 elements each,
and is still only 30,720 elements to read serially. That is close to the
1728-rank / 1137-element point the production run was aimed at.

This is guidance, not a guarantee — the coarse mesh must also be ordered
column-major (`tools/reorder_msh_columns.jl`) for tree columns to be
contiguous. The ground truth is the split percentage the setup prints.

---

## Implicit vertical diffusion

`:implicit_vdiff => true` moves `∂/∂z(μ ∂/∂z)` on `u`, `v`, `w` and `θ` into the
same column operator that already carries the vertical acoustics. It serves HEVI
and the 3D IMEX from one implementation, because both build the same operator
and both invert it with the same banded LU.

### Why it exists

Removing the acoustic limit exposes whatever was second, and on a mesh refined
in `z` that is not advection — it is the SGS diffusion of the vertical
derivative:

    λ_visc ~ ν / Δz²

The trap is that this is **invisible at startup**. `:lcfl_report` is taken on the
initial sounding, where a laminar profile gives `ν_t ≈ 0`, so the viscous row
says diffusion never limits `Δt` — and says it truthfully, at `t = 0`. As a
convective boundary layer spins up, `ν_t` climbs into the tens of m²/s over a few
hundred seconds and that row crosses the line, with nothing printing it again.
A run that dies at a fixed *model* time, having been perfectly stable up to it,
looks like a bug in the scheme and is not one. `:lcfl_report_every => N` re-prints
the table during the run so the crossing is visible rather than inferred.

The signature to check first: **a blow-up whose time moves with `Δt` is a
stability limit, not a physical instability.** The same deck on HEVI (small `Δt`)
and on IMEX3D (large `Δt`) blowing up at different model times is that signature.

### What is implicit, and what is deliberately not

On a mesh whose ζ lines are vertical, the physical viscous flux in `z` is

    τ_xz = μ (∂u/∂z + ∂w/∂x)          ← second term horizontal
    τ_yz = μ (∂v/∂z + ∂w/∂y)          ← second term horizontal
    τ_zz = (4/3) μ ∂w/∂z − (2/3) μ (∂u/∂x + ∂v/∂y)
    q_z  = κ ∂θ/∂z                    ← entirely column-local

and the operator takes the first term of each line. The cross terms stay
explicit; they scale as `μ/(Δx Δz)`, which is `Δz/Δx` of what is removed — 1/8
on the LESICP2 mesh. They could not be made implicit here in any case: a
horizontal flux is not something one column can see.

The `4/3` on `w` is the `2 − 2/3` that survives when `∂w/∂z` is collected from
both terms of `τ_zz`.

Because the column-local stress couples each variable only to itself, the
operator is block-diagonal in the variable index: mass is untouched, and `u`,
`v`, `w`, `θ` each get their own scalar vertical diffusion. There is no
cross-variable coupling to get wrong.

### Why it is nearly free to add

Two properties of the existing machinery do all the work.

1. **The split is a subtraction.** `f_exp = rhs!(u) − f_imp(u)`, so the scheme is
   consistent for *any* implicit operator. An approximate `D` cannot lose or
   double-count physics — it can only fail to remove stiffness. That is what
   makes dropping the cross terms safe.
2. **The band is assembled by probing.** `assemble_column_band` extracts the
   matrix *of the operator that is actually applied*, so HEVI's column LU and the
   3D IMEX's preconditioner absorb the diffusion with no code of their own.

What it does cost: the implicit variable set widens from `(ρ, ρw, ρθ)` to all
five, because `u` and `v` diffuse at the same `μ/Δz²` as `θ` and leaving them out
would remove none of the limit. The banded LU is `O(n·b²)` with both `n` and `b`
proportional to the variable count, so the factorisation goes as `nimp³` and the
triangular solves as `nimp²`. For the 3D IMEX this also retires the `pvars`
optimisation, which existed precisely because the acoustic operator's `ρu`, `ρv`
rows were the identity — with diffusion they are not.

### It needs `:PS`

The coefficient is read from the SGS closure, which returns its molecular floor
on a laminar initial state. `:RS` would freeze it there for the whole run: the
operator would carry no diffusion at all while still paying for the wider band,
and the only symptom would be the blow-up arriving on schedule and looking like
a bug somewhere else. Setting `:implicit_vdiff => true` under a dynamic closure
with `:RS` is refused at setup rather than allowed to do nothing. A deck-constant
`AV()` viscosity does not move, so `:RS` is fine there.

In the decks, `DBG_VDIFF=1` switches the linearisation with it; `DBG_LIN` still
overrides.

### What was verified

`test/hevi/test_vdiffusion.jl` (27 assertions, 1 and 3 ranks). The reference
implementations in it are transcribed *verbatim* from `_expansion_visc!` in
`rhs.jl` — all three sweeps, the full metric contraction, the same `dψ` index
order — because a cancellation between a weak form and anything else is not a
cancellation.

| check | result |
|---|---|
| `D` on `ρ̄ sin(kz)` vs `−μk² sin(kz)`, p=3, nelz 6 → 12 → 24 | `1.6e-2 → 1.6e-3 → 2.1e-4` (O(hᵖ)); p=5 gives `4.6e-6` |
| self-adjointness `⟨Dx,y⟩` vs `⟨Dy,x⟩` | `2.2e-12` relative |
| negative semidefinite `⟨Dx,x⟩` | `−1.5e6` |
| AV path vs the explicit RHS, all equations, all nodes | `9.0e-13` |
| SGS path vs the explicit RHS, interior, eq 2/3/4/5 | `3.7e-16 / 4.9e-16 / 7.2e-14 / 5.6e-13` |
| the same with the `4/3` removed (eq 4) | `1.0e-1` — the constant is load-bearing |
| probed band vs operator, with diffusion | `3.9e-16`; banded solve residual `7.8e-16` |

The momentum equations agree **in the interior** and not on the lateral faces,
which is the design and not a defect: the stress tensor's horizontal flux
components (`−(2/3)μ ∂w/∂z` in the `u` equation, `μ ∂u/∂z` in the `w` equation)
do not vanish even on a horizontally uniform field. The explicit ξ sweep
contracts them to `Σ_k ψ'_i(ξ_k) w_k = ψ_i(1) − ψ_i(−1)`, a term carried by the
two ξ-endpoint test functions alone, which DSS cancels between neighbouring
elements and around a periodic seam and cannot cancel where the domain ends. On
a laterally periodic LES the cancellation in `f_exp` is therefore complete to
round-off.

---

## Testing

Two suites, neither of which needs the Jexpresso dependency stack — only MPI,
LinearAlgebra, Printf and OrdinaryDiffEq.

```bash
julia --project=. test/hevi/test_load.jl               # it loads at all
julia --project=. test/hevi/runtests_standalone.jl     # the machinery
julia --project=. test/hevi/test_vdiffusion.jl         # implicit vertical diffusion
julia --project=. test/hevi/test_rtb.jl                # the physics
mpiexecjl -n 3 julia --project=. test/hevi/runtests_standalone.jl
```

### `test_load.jl` — the one that catches what the others cannot

Julia does **not** warn about method overwriting by default. During module
precompilation it both warns and treats it as a hard error. So a clash of
signatures is invisible to any test that loads the sources with `include` —
which every other test here does — and fatal the moment the package is
precompiled.

That is not hypothetical: the first version of this code shipped with

```julia
HEVI_ARK(name::Symbol = :ARS343) = HEVI_ARK(ark_tableau(name))
HEVI_ARK(; tableau::Symbol = :ARS343) = HEVI_ARK(ark_tableau(tableau))
```

Keyword arguments do not participate in dispatch, so **both** of those define
the zero-argument positional method `HEVI_ARK()` and the second silently
replaced the first. Every test passed; `using Jexpresso` failed outright. (The
repo had already been bitten by exactly this once — see the comment above
`assemble_mpi!` in `src/kernel/mpi/mpi_communications.jl`.)

`test_load.jl` loads all six HEVI sources into a module in a subprocess with
`--warn-overwrite=yes` and fails on any overwrite warning, without needing the
Jexpresso dependency stack. It reproduces the failure above at the exact line
numbers Julia reports.

### The rising thermal bubble

`test/hevi/test_rtb.jl` runs the dry nonhydrostatic benchmark (Robert 1993;
Giraldo & Restelli, JCP 227, 2008) — a warm cosine bubble, θ' = 0.5 K over a
250 m radius, released into a neutrally stratified atmosphere at rest — under
the full nonlinear compressible Euler system with gravity.

It is the right case for this integrator because it exercises exactly what the
split touches and nothing it does not: the vertical acoustic terms that move
to the implicit side, the buoyancy coupling `-g·dρ` inside the implicit
operator, and a base state at rest, so an imbalance introduced by the split
shows up as the whole domain drifting rather than as a subtle error in the
bubble.

The mesh is deliberately anisotropic — 8 × 1 × 24 elements over a 1000 m box,
so Δz is three times finer than Δx. On an *isotropic* bubble mesh HEVI buys
nothing, because there is no vertical acoustic term to remove that the
horizontal one does not match. Choosing an anisotropic mesh is what makes the
step-size claim testable rather than asserted.

What it asserts, and what came out:

| check | measured |
|---|---|
| discrete hydrostatic residual is small | 1.2e-4 relative in `dp/dz` (6.9% of the buoyant acceleration) |
| `f_imp` on the bubble-free state is **exactly** 0 | `0.0` |
| HEVI and explicit agree on `w` at a shared Δt | 1.9e-4 |
| HEVI at the big step agrees with itself at the small one | 4.7e-4 |
| bubble rises, θ' peak preserved, `\|v\|` at round-off | +7.5 m in 60 s, θ' = 0.4997 K, `\|v\|` = 1e-12, `w` max 0.46 m/s |
| explicit at that Δt diverges | it does |

`w` is the comparison variable because it starts at exactly zero and is what
buoyancy drives, so it has no large background to hide a discrepancy in.

A note on writing tests like this, learned the hard way here: **never put an
MPI collective behind a short-circuiting operator.** The natural way to write
the divergence check is

```julia
blew = !all(isfinite, u) || MPI.Allreduce(maximum(abs, u .- u0), MPI.MAX, comm) > 1e3
```

and it deadlocks — `||` skips the `Allreduce` on the first rank whose local
half is already true, while every other rank sits inside it. It passes in
serial and hangs on 3 ranks *after all the assertions have gone green*, which
is close to the worst way for a test to be wrong. Reduce unconditionally, then
combine locally. For the same reason the fixture's equation of state clamps
`(ρθ)^γ` rather than letting it throw: an exception raised inside the RHS
unwinds one rank out of the halo exchange the others are waiting in.

The step sizes are **measured, not estimated**. A sweep of both schemes on this
mesh puts the explicit limit between 0.02 and 0.04 s and the HEVI limit between
0.06 and 0.08 s — a factor of **2.5×** on a mesh with 3:1 vertical anisotropy.
A first-principles estimate from the LGL node gaps put both limits ~1.8× higher
(the spectral radius of the SEM derivative grows like `N²/h`, not `1/Δ_min`);
the factor is common to both schemes and cancels out of the ratio.

Two of these are worth dwelling on. `f_imp = 0` exactly on the bubble-free
state is the well-balancedness of the split, and it is exact rather than small
because the operator acts on `u − qe` and that state *is* `qe`. And the
hydrostatic residual is a property of the *fixture's* strong-form
discretisation in TOTAL variables, not of HEVI: it lives in `f_full`, so both
schemes inherit exactly the same one and it cancels out of every comparison
between them. (Jexpresso's own answer to it is `:SOL_VARS_TYPE => PERT()`,
which subtracts the base state analytically; `CompEuler/theta` uses that.)

### Timing HEVI against the explicit scheme

One deck, one switch at the top of `user_inputs()`:

```julia
lexplicit = false          # HEVI_ARK(:ARS232)      -> ./output_hevi/
lexplicit = true           # CarpenterKennedy2N54() -> ./output_explicit/
```

```bash
julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler","rtb_hevi")'
```

**Do not time them at the same Δt.** That measures cost per step, which is not
the question. HEVI costs *more* per step — 4 `rhs!` evaluations, 4 applications
of the implicit operator and 3 column solves, against `CarpenterKennedy2N54`'s
5 `rhs!` and no solve — and pays for it by taking a bigger one. Only wall-clock
to a fixed physical time answers it. Each scheme therefore runs at its own
largest safe step (`Δt` is set from `lexplicit` a few lines below the switch).

**What to read.** `:lstep_heartbeat => true` is already on, and it prints
`s/step` measured over the interval since the *previous* heartbeat — so after
the first few steps it is a steady-state rate with the JIT excluded, which
total wall time is not. The figure of merit is

```
wall to t = tend   =   (s/step)  ×  tend / Δt
```

**Finding each scheme's actual limit.** The defaults are conservative and the
honest comparison uses the largest step each can really take. Sweep `Δt` with a
short `tend` (30 s is plenty) and watch for the pressure going negative —
divergence shows up as a `DomainError` out of `perfectGasLaw_ρθtoP`. The CFL
table printed at startup reports the CFL each run is actually sitting at, so a
run over its limit is visible before it diverges.

**If HEVI comes out slower, profile it before tuning anything.**

```bash
JEXPRESSO_HEVI_PROFILE=1 julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler","rtb_hevi")'
```

Every 50 steps it prints where the time went:

```
┌── HEVI cost breakdown over 100 steps
│ rhs!               4.0 /step    0.00328 s each    0.0131 s/step   48.2 %
│ f_imp              4.0 /step    0.00155 s each    0.0062 s/step   22.7 %
│ column solve       3.0 /step    0.00229 s each    0.0069 s/step   25.2 %
│ refactorise        0.0 /step        0.0 s each       0.0 s/step    0.0 %
│ accounted                                        0.0262 s/step
│ measured step                                    0.0273 s/step   gap 0.0011 s
└────────────────────────────────────────────────
```

That is the reference run — `test_rtb.jl`'s mesh (16 005 points, 165 columns,
97 levels), where HEVI costs **1.29× per step** and wins **1.9× per simulated
second** over `Tsit5` at its own limit. Note the `gap`: everything is
accounted for.

**The one that has already bitten:** `f_imp` dominating. On the `rtb_hevi` case
it measured 0.0929 s per application — 75% of the step and 3.3× the cost of a
full nonlinear `rhs!` — against 0.0016 s for the same code on the same mesh in
`test_rtb.jl`.

The cause was Jexpresso's own struct definitions. `St_metrics` and `St_mesh`
are `Base.@kwdef mutable struct`s with **no field type annotations**, so
`metrics.dζdz` and `mesh.connijk` are `::Any` and every index inside an element
loop is a boxed dynamic dispatch. Hoisting them into locals does not help — the
local is still `::Any`. Only a *typed function barrier* fixes it: pass the
arrays positionally so the callee specialises on their runtime types, which is
what `rhs.jl` does for `_inviscid_rhs_el_3d!` and what `hevi_apply_A!` now does
for `_hevi_A_elements!`.

Measured penalty on the same kernel, typed vs `::Any` structs: **14× without
the barrier, 1.0× with it.** The tests here could not see it, because the mock
structs are concretely typed — which is why `runtests_standalone.jl` now builds
a deliberately untyped variant and asserts the barrier holds.

The other two things to read:

* **`refactorise` non-zero per step.** γΔt is not matching what setup
  computed, so every step rebuilds every column matrix. A rebuild costs about
  one time step, so this roughly doubles the cost. A bug, not a tuning issue.
* **a large `gap`.** Then the cost is *not* in HEVI's own work — not `f_imp`,
  not the column solve — and the place to look is outside this directory:
  callbacks, saving, allocation, GC.

`f_imp` and the column solve depend only on the mesh, not on how expensive the
case's `rhs!` is. So on a case whose `rhs!` costs `T` per call, a HEVI step
should be about `4T` plus the small fixed cost the table above shows. If the
measured step is much more than that, the extra is not the split.

(Finding each scheme's actual limit: edit `Δt` and `tend` at the top of
`user_inputs()` — 30 s of simulated time is enough to see divergence — and
watch for a `DomainError` out of `perfectGasLaw_ρθtoP`, which is the pressure
going negative.)

### The Jexpresso case

`problems/CompEuler/rtb_hevi/` is the same benchmark as a real deck — same
domain, same bubble, same anisotropic mesh (committed, plus the `.geo` that
generates it). Switching between HEVI and the explicit reference is two
commented lines at the top of `user_inputs.jl`.

To put it in CI, use the repo's own tooling rather than editing the registry
by hand:

```bash
julia --project=. test/generate_ci_ref.jl CompEuler/rtb_hevi
```

Runs without the Jexpresso dependency stack — only MPI, LinearAlgebra, Printf
and OrdinaryDiffEq — against a synthetic Cartesian SEM mesh partitioned so
that columns are deliberately cut across ranks. Verified on 1, 2, 3, 4, 5 and
7 ranks; on 3 ranks, 28 of the local column-instances are split and the answer
is still bit-identical to serial.

Results, identical on every rank count tried:

```
column locality        leakage 0.0 outside the lit column
band == operator       1.9e-16 relative
column spectrum        max|Re|/max|Im| = 5.6e-17  (hyperbolic, as it should be)
factorised solve       4.5e-16 residual
stage solve residual   3.1e-11 relative
cross-rank consistency 0.0
ARS343 vs Tsit5        order 3.01, err 2.1e-8
at 8x the explicit acoustic limit   HEVI bounded at 2e-12, explicit NaN
```

The `ARS343 vs Tsit5` comparison is the one that matters most: it is the only
check that would catch `f_imp` and the stage solve disagreeing with *each
other*. Both self-checks validate one half against itself and would pass.

Note the tests are slow under MPI on a small machine — the mock assembler does
a global reduction per RHS call, which is fine for correctness and terrible
for throughput. That is the fixture, not the code under test.

---

## Not implemented

* **Acoustic substepping** — the ~8–10× on this grid. See above.
* **GPU** — the column solve is LAPACK banded LU on the host.
* **`:SOL_VARS_TYPE => PERT()`** — refused at setup. The operator acts on
  `u − qe`, which is the deviation only under `TOTAL()`; under `PERT()` the
  solution vector already *is* the deviation. Supporting it is a small,
  well-defined change (feed `u` straight in, and re-derive the buoyancy term
  against the case's perturbation-form `user_source!`) but it has to be
  verified per case. Worth doing: `PERT()` is the well-balanced formulation and
  would remove the hydrostatic residual the bubble test reports.
* **AMR** (`:ladapt`) — refused at setup: adaptation invalidates the column
  topology, the factorisations and the gather/scatter plan.
* **Terrain-following meshes** are handled (the implicit variable set widens to
  all five components when `dζ/dx`, `dζ/dy` do not vanish, and levels are
  identified by normalised height) but have not been exercised on a real
  warped case.
