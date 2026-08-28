# IMEX3D: where everything is

A location map for the fully-3D implicit-acoustics solver — what is built,
where, and what runs on every stage of every step. Line numbers are from the
commit that added this file; the function names are the stable handle.

For *why* the scheme is built this way, see `README_IMEX3D.md`. This file is
strictly *where*.

**New to this code? Read §1b first.** The single most common wrong assumption
about this solver is that there is an implicit matrix somewhere. There is not,
and §1b says what exists instead.

Line numbers drift; the function names do not. If a number here does not land
where it claims, `grep` the name.

---

## 1. The one-paragraph summary

An IMEX Runge–Kutta scheme (ARS343) splits the right-hand side into an
**implicit** half carrying all acoustic terms and an **explicit** half carrying
everything else. Each stage requires solving

```
(I − γΔt·A) U = B
```

where `A` is the linearised 3D acoustic operator over the five fields
`(ρ, ρu, ρv, ρw, ρθ)`. That solve is a **preconditioned GMRES**, and the
preconditioner is a **direct banded solve down each vertical column** — exact
for the vertical acoustics, silent about the horizontal, which is why the
iteration count is governed by the *horizontal* acoustic Courant number
`CFL_h`.

Optionally (`:imex_schur`) the five-field system is reduced algebraically to a
**single scalar Helmholtz equation** in `P = β·Θ` over `Np` unknowns instead of
`5·Np`, solved with its own scalar column preconditioner, and the five fields
are rebuilt pointwise afterwards.

---

## 1b. "Where is the implicit matrix assembled?" — it isn't

**The implicit operator is never assembled. There is no matrix to find.**

This is the first thing that surprises a reader coming from a code that builds
a sparse Jacobian, so it is stated here before the location tables: `(I − γΔt·A)`
exists only as a **function that multiplies a vector**. Nothing in the run ever
holds its entries. Searching for an assembly routine, a `SparseMatrixCSC`, or a
`spzeros` on the implicit path will find nothing, because there is nothing.

### What the Krylov solve actually sees

The whole "matrix" is this closure (`imex3d.jl:289`):

```julia
matvec! = let params = params, op = op, g = g
    (W, V) -> begin
        hevi_apply_A!(W, V, params, op)     # W = A·V, element by element
        @inbounds @. W = V - g * W          # W = (I − γΔt·A)·V
        return W
    end
end
```

`gmres_solve!` is handed that closure and never asks for anything else — no
entries, no sparsity pattern, no transpose. Everything below is where the
`A·V` on line 1 comes from.

| what | where |
|---|---|
| the matvec `(I − γΔt·A)·V`, five fields | `imex3d.jl:289` (closure) → `operator.jl:578` `hevi_apply_A!` |
| the matvec `H[P]`, scalar Schur system | `schur_stage.jl:194` (closure) → `schur.jl:369` `schur_H!` |
| **the element kernel that IS `A`** | `operator.jl:578` `hevi_apply_A!` |
| implicit vertical diffusion, added into `A` | `vdiffusion.jl:313` `_hevi_D_elements!` |
| the operator OBJECT — coefficients (β, θ̄), wall masks, scratch. Not entries. | `acoustic.jl:60` `build_hevi_fast_operator`, `operator.jl:189` `build_hevi_operator` |

### Where a real matrix does exist: the preconditioner, and only there

The one place entries are stored is the preconditioner, and even there they are
**recovered by probing the operator, never derived**: apply the vertical
operator to coloured unit vectors and read its columns off the result. The band
is therefore the operator *by construction*, and the two cannot drift apart —
which is also why there is no derivation to read anywhere.

| what | where |
|---|---|
| five-field band `W = I − γΔt·A` per column, by probing | `factorize.jl:80` `assemble_column_band` |
| → banded LU (build / refactorise) | `factorize.jl:171` and `:199` `LAPACK.gbtrf!` |
| → applied, one triangular solve pair per column | `factorize.jl:258` `LAPACK.gbtrs!`, via `factorize.jl:218` `hevi_column_solve!` |
| scalar Schur band, by probing `schur_H!` | `schur_precond.jl:108` `assemble_schur_column_band!` |
| → LU / solve | `schur_precond.jl:230` `gbtrf!`, `:286` `gbtrs!` |

These are **block-diagonal by vertical column**, in LAPACK general-band storage
— `nvar·nlev` square per column, bandwidth `nvar·ngl − 1`. They are not the 3D
operator; they are its vertical part, which is the whole point of the scheme
(see §7).

### If you want to look at the matrix anyway

Both builders take a hook that is handed the **unfactorised** band. That is the
only moment it can be inspected: `gbtrf!` overwrites `AB` with the LU factors
in place, so reading it afterwards diagonalises the factors rather than the
matrix.

```julia
build_column_factorization(params, op, cc, topo, gdt;
                           verify_hook = (AB, kl, ku, n) -> ...)   # factorize.jl:162
build_schur_column_precond(params, topo, comm, lam;
                           verify_hook = (AB, kl, ku, n) -> ...)   # schur_precond.jl:183
```

`imex3d.jl:1102` uses exactly that hook to take the column spectrum
(`hevi.jl:563` `hevi_column_spectrum`). For a **dense** `(I − γΔt·A)` or a dense
`H`, built by probing the full operator with unit vectors and compared against
the matrix-free form, see `test/hevi/test_schur_helmholtz.jl`,
`test_schur_solve.jl`, `test_schur_blocks.jl` and `test_schur_precond.jl` —
those are the reference implementations to copy if you need the entries.

### Why it is built this way

Assembling `(I − γΔt·A)` would mean storing a sparse matrix over `5·Np`
unknowns — 79.6 M degrees of freedom on the production mesh — and rebuilding it
whenever γΔt or the reference state moves. The matvec costs one element-loop
pass either way, so the assembly buys nothing the iteration needs, and it adds
a second statement of the operator that can disagree with the first. The
preconditioner needs entries and gets them by probing, which is the one route
that cannot disagree.

---

## 2. Build path — once, at setup

| step | location |
|---|---|
| gate: is IMEX3D on? | `imex3d.jl:331` `imex3d_enabled` — true if `:limex` or `:ode_solver isa IMEX_ARK` |
| dimension/backend guard | `params_setup.jl:706-714` — 3D and CPU only, and it errors rather than falls back |
| **entry point** | `params_setup.jl:715` → `params = merge(params, (imex = build_imex3d(params, inputs),))` |
| **the builder** | `imex3d.jl:720` `build_imex3d(params, inputs)` |

Inside `build_imex3d`, in order:

| what | location |
|---|---|
| vertical column topology (who owns which column, how many levels) | `imex3d.jl:820` → `columns.jl:618` `build_column_topology` |
| the **operator** `A` — full 3D, five fields | `imex3d.jl:836` → `acoustic.jl:60` `build_hevi_fast_operator` |
| which fields to precondition | `imex3d.jl:870` → `hevi_choose_vars` (reads the metrics, not a deck flag) |
| vertical-only operator for the preconditioner | `imex3d.jl:879` → `operator.jl:189` `build_hevi_operator(...; full=false)` |
| column ownership + gather/scatter plan | `imex3d.jl:882-883` → `columns.jl:793` `build_column_comm` |
| **assemble + LU-factorise the bands** | `imex3d.jl:914` → `factorize.jl:162` `build_column_factorization` |
| the preconditioner object | `imex3d.jl:915` `IMEX3DPrecond(opv, cc, fac, topo, pvars)` |
| **or the deck's own**, `:imex_precond => :custom` | `precond_api.jl` `build_custom_precond` |
| distributed inner product | `imex3d.jl:923` `build_distributed_inner` |
| GMRES workspace (`5·Np`) | `imex3d.jl:928` `GMRESWorkspace` |
| **Schur path**, if `:imex_schur` | `imex3d.jl:934` → `schur_stage.jl:81` `build_imex3d_schur` |

### How the band is built

`factorize.jl:80` `assemble_column_band` recovers the matrix **by probing the
operator**, not by deriving entries: it applies the vertical operator to
coloured unit vectors and reads the columns off the result. The matrix is the
operator by construction, so the two cannot drift apart.

`factorize.jl:171` `LAPACK.gbtrf!` then factorises each column's band in place.
`factorize.jl:156` notes the one consequence: `gbtrf!` **overwrites** `AB`, so
the verification hook is handed the band *before* factorisation — that is the
only moment the assembled matrix can be checked against the operator it claims
to be.

---

## 3. Apply path — every GMRES iteration of every stage

| step | location |
|---|---|
| ARK driver calls the stage solve | `ark.jl` → `params.imex.solve!` |
| **five-field stage solve** | `imex3d.jl:251` `imex3d_solve!` |
| — matvec closure `(I − γΔt·A)·V` | `imex3d.jl:274` → `operator.jl:578` `hevi_apply_A!` |
| — **preconditioner closure** | `imex3d.jl` → `imex_precond_apply!` → `imex3d.jl:156` `imex3d_precond!` |
| — the Krylov loop | `imex3d.jl:287` → `krylov.jl:272` `gmres_solve!` |
| — where GMRES actually calls it | `krylov.jl:326` (Arnoldi) and `krylov.jl:400` (right-preconditioned update) |

### What `imex3d_precond!` does

It delegates to `factorize.jl:218` `hevi_column_solve!`, which is three steps:

| step | location |
|---|---|
| gather each rank's column pieces onto the owner | `factorize.jl:254` `column_gather!` |
| one banded back-substitution per owned column | `factorize.jl:258` `LAPACK.gbtrs!` |
| scatter the solved columns back | `factorize.jl:263` `column_scatter!` |

A column may be **split across ranks** — z is never partitioned by
`:lxy_partition`, but a column's nodes can still live on more than one rank —
which is why the gather/scatter exists at all rather than a purely local solve.

---

## 4. The Schur path (`:imex_schur => true`)

Same shape, one field instead of five.

| step | location |
|---|---|
| build | `schur_stage.jl:81` `build_imex3d_schur` |
| — scalar column preconditioner | `schur_precond.jl:183` `build_schur_column_precond` |
| — its band, by probing `schur_H!` | `schur_precond.jl:108` `assemble_schur_column_band!` |
| **stage solve** | `schur_stage.jl:149` `imex3d_solve_schur!` |
| — five-field RHS → scalar RHS | `schur_stage.jl:` `schur_setup_rhs!` (per solve, instrumented separately) |
| — matvec `H[P]` | `schur_stage.jl:184` → `schur.jl:369` `schur_H!` |
| — **preconditioner** | `schur_stage.jl` → `imex_precond_apply!` → `schur_precond.jl:259` `schur_precond!` |
| — the same three steps | `schur_precond.jl:283` gather, `:286` `gbtrs!`, `:292` scatter |
| — rebuild the five fields | `schur_stage.jl:` `schur_recover!` |

`schur_H!` is assembled from two sweeps — a scalar gradient and a masked
divergence — in `schur_kernel.jl` (`:imex_schur_kernel`, default on). The
reference form, two full five-field `hevi_apply_A!` calls, is kept reachable in
`schur.jl` as the independent statement the kernel is checked against; they
agree to 1.9e-16.

**The Schur reduction requires the advective Θ row.** With the flux form the
elimination leaves a 2×2 system in `(ρ, P)` rather than one scalar, so
`:imex_schur` forces `theta_advective` on the operator rather than offering it
as a separate choice.

---

## 4b. Swapping the preconditioner — `:imex_precond => :custom`

The stage solve reaches every preconditioner — both built-ins and anything a
deck supplies — through three generic functions in
**`precond_api.jl`**. Nothing else in the scheme knows which one is in use:
`IMEX3DCache` is parametric in the preconditioner type, and GMRES asks only for
a callable `precon!(Z)` that overwrites `Z` with `M⁻¹Z`.

| what | location |
|---|---|
| the contract, in full | `precond_api.jl` (header) |
| **apply** `M⁻¹` in place — required | `precond_api.jl` `imex_precond_apply!` |
| refresh under `:PS` — optional, no-op by default | `precond_api.jl` `imex_precond_refresh!` |
| the setup report's one line — optional | `precond_api.jl` `imex_precond_describe` |
| what the builder is handed | `precond_api.jl` `imex_precond_context` |
| resolve `:imex_precond_build` and prove it applies | `precond_api.jl` `build_custom_precond` |
| the two failures that are otherwise silent | `precond_api.jl` `imex_precond_selfcheck` |
| built-in methods, five-field | `imex3d.jl`, just below `imex3d_precond!` |
| built-in methods, scalar | `schur_precond.jl`, at the foot of the file |
| build site, five-field | `imex3d.jl` `build_imex3d`, after the `:column` block |
| build site, scalar | `schur_stage.jl` `build_imex3d_schur` |
| refresh site | `imex3d.jl:~720` `ark_relinearize!` |

A user writes one builder and one `imex_precond_apply!` method, both in their
case's `user_inputs.jl` — that file is `include`d into the `Jexpresso` module,
so no `import` and no edit to `src/` is involved. The step-by-step version,
with a complete worked example, is **"Writing your own preconditioner"** in
`README_IMEX3D.md`.

`:imex_precond` steers this solve only. The radiative-transfer solver has an
unrelated preconditioner menu in `src/kernel/operators/asm_preconditioner.jl`.

---

## 5. Relinearisation — when the reference state moves

| step | location |
|---|---|
| driver hook | `ark.jl:600` `ark_relinearize!` (declaration) |
| implementation | `imex3d.jl:650` `ark_relinearize!(params, imex, u)` |
| refresh + re-probe the Schur band | `schur_precond.jl:243` `refactorize_schur!` |
| re-factorise the five-field bands | `factorize.jl:199` `gbtrf!` |

`:imex_linearization => :RS` freezes the coefficients at the reference state and
never calls this; `:PS` refreshes them every `:imex_update_freq` steps. Both
preconditioners must be refreshed together with the operator — a band assembled
from stale coefficients still *converges*, just to a worse iteration count, so
the failure is silent.

---

## 6. The deck keys that steer it

| key | effect | read at |
|---|---|---|
| `:limex` / `:ode_solver => IMEX_ARK` | turns IMEX3D on | `imex3d.jl:331` |
| `:imex_schur` | scalar Schur stage solve | `schur_stage.jl` |
| `:imex_schur_kernel` | fast matvec for `H` (default on) | `schur_stage.jl` |
| `:imex_precond` | `:column`, `:none`, or `:custom` | `imex3d.jl`, `schur_stage.jl` |
| `:imex_precond_build` | `f(ctx) -> pc`, required by `:custom` | `precond_api.jl` |
| `:implicit_vdiff` | implicit vertical SGS diffusion. **Refuses to coexist with `:imex_schur`** — the pair demotes `:imex_schur` and says so in orange | `imex3d.jl` `imex_resolve_schur_vdiff!` |
| `:imex_allow_schur_vdiff` | run that pair anyway, to measure it | `imex3d.jl` |
| `:imex_rtol`, `:imex_restart`, `:imex_maxiter` | Krylov tolerance, cycle length, cap | `imex3d.jl` |
| `:imex_linearization`, `:imex_update_freq` | `:RS` frozen / `:PS` refreshed | `imex3d.jl:650` |
| `:imex_verify` | setup self-check on the assembled operator and a full stage solve | `imex3d.jl` |
| `:imex_umax` | expected max flow speed, for the stability report | `imex3d.jl` |
| `:imex_lateral_walls`, `:imex_wall_flux` | free-slip treatment in the implicit operator | `imex3d.jl` |

---

## 7. What governs the cost

The column preconditioner removes the vertical acoustic coupling **exactly** and
does nothing for the horizontal, so the Krylov iteration count is set by

```
CFL_h = γΔt·c/h_x
```

and by nothing else — not `h_z`, not the total acoustic CFL, not the node count.
Computed at `imex3d.jl` `imex_horizontal_cfl`, printed in the setup report.

The corollary is that **this scheme's advantage is proportional to the grid's
acoustic anisotropy** `h_x/h_z`. At 1:1 there is none: neither this nor HEVI has
a cheap direction left to precondition.

Measured on `CompEuler/rtb3d_schur`, 25 ranks, 4:1 anisotropy, both arms back to
back on the same nodes:

| s/step | five-field | Schur |
|---|---|---|
| stage solve | 9.542 | 1.880 |
| — matvec | 4.955 | 1.092 |
| — preconditioner | 2.042 | 0.250 |
| — orthogonalise | 1.922 | 0.199 |
| iterations/step | 69.4 | 36.2 |
| **step** | **10.651** | **2.990** |

Run-to-run variation on that cluster is ~13%, so the per-term ratios are the
transferable numbers, not the third digit of a step time.

---

## 8. Tests

| what it pins | file |
|---|---|
| the operator's block structure | `test/hevi/test_schur_blocks.jl` |
| the scalar Helmholtz against a dense solve | `test/hevi/test_schur_helmholtz.jl`, `test_schur_solve.jl` |
| the column preconditioner's band vs a dense `H_v` | `test/hevi/test_schur_precond.jl` |
| the fast matvec vs the reference form | `test/hevi/test_schur_kernel.jl` |
| both production stage solves agree | `test/imex3d/test_schur_stage.jl` |
| `build_imex3d`, both arms | `test/imex3d/test_schur_build.jl` |
| the `:custom` preconditioner hook | `test/imex3d/test_custom_precond.jl` |
| the `:imex_schur` + `:implicit_vdiff` repair | `test/imex3d/test_schur_build.jl` |
| answers independent of rank count | `test/hevi/test_rank_independence.jl` |
