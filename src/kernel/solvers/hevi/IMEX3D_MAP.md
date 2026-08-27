# IMEX3D: where everything is

A location map for the fully-3D implicit-acoustics solver — what is built,
where, and what runs on every stage of every step. Line numbers are from the
commit that added this file; the function names are the stable handle.

For *why* the scheme is built this way, see `README_IMEX3D.md`. This file is
strictly *where*.

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
| — **preconditioner closure** | `imex3d.jl:281-283` → `imex3d.jl:156` `imex3d_precond!` |
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
| — **preconditioner** | `schur_stage.jl:190` → `schur_precond.jl:259` `schur_precond!` |
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
| `:imex_precond` | `:column` or `:none` | `imex3d.jl` |
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
| answers independent of rank count | `test/hevi/test_rank_independence.jl` |
