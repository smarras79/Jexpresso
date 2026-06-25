# IMEX / Burgers: from `sl/imex` to `sl/imex_newmaster`

This note is for Santolo. It explains **what changed** when the IMEX
(implicit–explicit) time integrator and the 2-D Burgers test case were moved
from the `sl/imex` branch onto a branch off `sm/newmaster`
(`sl/imex_newmaster`), and **why**.

The numerical method is the same. What changed is *how it is wired into
Jexpresso*: the standalone, copy-pasted entry chain is gone, and IMEX is now
selected exactly like every other time integrator — through
`inputs[:ode_solver]`.

---

## 1. TL;DR

| Topic | `sl/imex` (old) | `sl/imex_newmaster` (new) |
|---|---|---|
| How IMEX is launched | Separate entry chain: `IMEXJexpresso.jl` → `IMEXrun.jl` → `IMEXdrivers.jl` → `imex_time_loop!` | Normal entry: `Jexpresso.jl` → `run.jl` → `drivers.jl` → `time_loop!`, which dispatches to `imex_time_loop!` when `inputs[:ode_solver] isa IMEX` |
| Selecting IMEX in a case | comment `:ode_solver`, and *run a different `main`* (`IMEXJexpresso.jl`) | `:ode_solver => IMEX()` in `user_inputs.jl`, run the normal `Jexpresso` |
| Case directory | `problems/equations/Burgers/case2d_imex_sl/` | `problems/Burgers/case2d_imex_sl/` |
| Integrator file | `src/kernel/solvers/IMEXTimeIntegrators.jl` (+ `IMEXMatrixStorage.jl`, `MyPrec.jl`) | `src/kernel/solvers/IMEXTimeIntegrators.jl` only |
| Linear solve | `Krylov.fgmres` + `MyPrec`/AMG/ILU **or** direct `:lsolve` | user `:lsolve` (direct) + direct sparse fallback; Krylov/AMG stack dropped |
| Sparse SEM Laplacian | added by the branch to `element_matrices.jl` | **already present on `sm/newmaster`** — reused as-is |
| Extra dependencies | `DiffEqBase`, `DiffEqDevTools`, `AlgebraicMultigrid`, `SciMLOperators`, `IncompleteLU`, … | **none added** |

> **Note on the name `_sl`.** On `sl/imex`, the five non-`user_inputs.jl`
> files of `case2d_imex_sl` are byte-for-byte identical to `case2d_imex`.
> There is no semi-Lagrangian transport anywhere on the branch (no
> departure-point / trajectory / upstream-interpolation code). `_sl` only
> marks the case that routes through the *generic* `imex_time_loop!` instead
> of the self-contained `imex_ars232_time_loop!`. The new branch reproduces
> that behaviour; it does **not** add semi-Lagrangian advection. If genuine
> SL is wanted, it has to be written from scratch.

---

## 2. The core change: no more duplicated `main`

### Old (`sl/imex`)

To run IMEX you had to launch a **second copy of the program**:

```
src/IMEXJexpresso.jl     ≈  src/Jexpresso.jl      (different include list; ends with IMEXrun.jl)
src/IMEXrun.jl           ≈  src/run.jl            (one line changed: driver_file = IMEXdrivers.jl)
problems/equations/IMEXdrivers.jl
                         ≈  problems/equations/drivers.jl
                                                  (non-linsolve branch calls imex_time_loop!)
```

`IMEXJexpresso.jl` and `Jexpresso.jl` both define `module Jexpresso`, so the
two could not coexist in one session. Selecting IMEX meant *editing which
top-level file you include*, not editing the case.

### New (`sl/imex_newmaster`)

The three `IMEX*` entry files are **deleted**. IMEX is a dispatch tag:

```julia
# src/kernel/abstractTypes.jl
abstract type AbstractTimeIntegrator end
struct IMEX <: AbstractTimeIntegrator end
```

```julia
# src/kernel/solvers/TimeIntegrators.jl  (inside time_loop!)
if inputs[:ode_solver] isa IMEX
    return imex_time_loop!(inputs, params, u)
end
# ... otherwise the usual OrdinaryDiffEq solve(...) path
```

So a case opts in with one line — `:ode_solver => IMEX()` — and runs through
the *same* `Jexpresso.jl` / `run.jl` / `drivers.jl` as every explicit case.
This is the whole point of the migration: one driver, one entry point.

The driver-side SciML warm-up (`precompile_warmup_run!`) is also skipped for
IMEX, since `IMEX()` is not an OrdinaryDiffEq algorithm.

---

## 3. File-by-file mapping

| `sl/imex` | `sl/imex_newmaster` | Notes |
|---|---|---|
| `src/IMEXJexpresso.jl` | *(removed)* | replaced by dispatch in `time_loop!` |
| `src/IMEXrun.jl` | *(removed)* | normal `run.jl` is used |
| `problems/equations/IMEXdrivers.jl` | *(removed)* | normal `drivers.jl` is used |
| `src/kernel/solvers/IMEXTimeIntegrators.jl` | `src/kernel/solvers/IMEXTimeIntegrators.jl` | rewritten: `imex_time_loop!(inputs, params, u)` (3-arg, params-driven) |
| `src/kernel/solvers/IMEXMatrixStorage.jl` | *(not ported)* | mixed-precision `StoredIMEXMatrix` + cached preconditioner — unused by the direct-solve case |
| `src/kernel/solvers/MyPrec.jl` | *(not ported)* | AMG/ILU/Jacobi preconditioner wrappers + `ldiv!` overload — unused by the direct-solve case |
| `src/kernel/solvers/IMEX_ARS.jl` (`imex_ars232_time_loop!`) | *(not ported)* | the self-contained ARS path used by `case2d_imex` (no `_sl`); not what `case2d_imex_sl` uses |
| `element_matrices.jl` IMEX additions (`build_laplace_matrix`, `DSS_laplace_sparse`) | already on `sm/newmaster` | reused unchanged |
| `src/io/mod_inputs.jl` IMEX defaults | `src/io/mod_inputs.jl` IMEX defaults | re-added, trimmed to the lean integrator (and the `matrix_free` bareword bug is fixed → `:matrix_free`) |
| `problems/equations/Burgers/case2d_imex_sl/` | `problems/Burgers/case2d_imex_sl/` | physics files taken from newmaster `case2d`; `user_inputs.jl` carries the IMEX config |

---

## 4. The integrator (`imex_time_loop!`)

### What stayed the same (the method)

The ARS(2,3,2) additive Runge–Kutta logic is reproduced exactly:

* Per-stage RHS
  `rhs_i = uⁿ + Σ_{j<i} Δt[ A_ex[i,j] S(U_j) + (A_im[i,j] − A_ex[i,j]) L(U_j) ]`
  (`construct_rhs_rk!`).
* Per-stage implicit operator `L_curr = I − λ L`, `λ = Δt · A_im[i,i]`
  (`L_update`).
* A fixed-point nonlinear loop solving `L_curr U_i = rhs_i`.
* Final step update `uⁿ⁺¹ = uⁿ + Δt Σ_i b_ex[i] S(U_i)` — **explicit
  combination only**, exactly as in the original loop. (This is why
  `:lvisc => true` matters; see §6.)
* The `"multistep"` family is kept too, with the explicit-Euler warm-up.

`S(u)` is still the standard `rhs!`, and `L = −μ M⁻¹ K` is still assembled by
the user closure with `build_laplace_matrix` + `DSS_laplace_sparse`.

### What changed (the plumbing)

1. **Signature.** `imex_time_loop!(inputs, sem, qp, params, u)` →
   `imex_time_loop!(inputs, params, u)`. `mesh`, `qp`, `neqs`, `npoin`,
   `basis`, `ω`, `metrics`, `Minv` are all pulled from `params` (newmaster
   already carries them). The user `S_fun!`/`bcs_fun!` keep their legacy
   `(…, params, sem)` / `(…, sem, qp)` signatures — `sem` is passed as
   `nothing` because the Burgers closures don't use it.

2. **Linear solve.** Dropped the `Krylov.fgmres` + `MyPrec` (AMG/ILU) +
   `StoredIMEXMatrix` machinery. The integrator now uses `inputs[:lsolve]`
   when provided (the case supplies a direct `L_curr \ b`) and otherwise a
   direct sparse factorization fallback. Rationale: `case2d_imex_sl` already
   selected the direct solver (`:lsolve = (L,b) -> L \ b`), so the Krylov/AMG
   path was never exercised at solve time — keeping it would have pulled in
   five heavy dependencies for dead-at-runtime code. A custom `:lsolve`
   closure is the extension point if iterative/preconditioned solves are
   wanted later.

3. **Constant-operator factorization caching.** When the implicit operator is
   time-independent (`:upd_L => false`, the Burgers case), the per-stage system
   matrix `I − λᵢ L` never changes. It is therefore assembled **and
   LU-factorized once** and the factorization is reused for every step (the
   solve receives the cached `lu(...)` object, so `F \ b` does only
   forward/back substitution). Without this, a fresh sparse LU — with full 2-D
   fill-in — was built and discarded on every stage of every step
   (`k · nsteps` factorizations), which made this small case allocate tens of
   GiB. This mirrors how Jexpresso's linear-solve path assembles/factorizes
   once and solves many times.

   The constant-operator path also runs inside a **type-stable function
   barrier** (`_imex_rk_run_const!`). `imex_time_loop!` reads the time step,
   work buffers and user closures (`:S_fun`, `:L_fun`) out of the `inputs`
   Dict and runtime-typed allocations, so inside it they are inferred `Any`.
   Handing an `Any`-typed stage vector to `rhs!` makes the whole RHS
   evaluation type-unstable and box on every array op — that allocated several
   GiB even after the factorization fix. The barrier takes `u` and the
   closures as positional arguments (so Julia specializes on their concrete
   types) and allocates the stage vectors/buffers with `similar(u)`, keeping
   `rhs!` allocation-light and the per-step work in the low-MiB range. The
   solve is an in-place `ldiv!` with the cached LU.

4. **Output.** Writes through newmaster's `write_output(SD, u, uaux, t, iout,
   mesh, mp, connijk_original, …, qp.qvars, qp.qoutvars, outformat; nvar,
   qexact)` — the same call newmaster's own callbacks use. `uaux` is refreshed
   from `u` (`u2uaux!`) before each write. The diagnostic schedule
   (`:diagnostics_at_times`) is honoured with a round-off-robust "write at the
   first step that reaches the next requested time" test.

---

## 5. `user_inputs.jl`: old vs new

Same ARS(2,3,2) tableaux, same `S_fun!` / `L_fun!` / `build_L` / `_imex_L`
(`L = −μ M⁻¹ K`, cached in a closure). The differences:

| key | `sl/imex` case2d_imex_sl | new case2d_imex_sl | why |
|---|---|---|---|
| `:ode_solver` | commented out; selected by running `IMEXJexpresso.jl` | `IMEX()` | native dispatch |
| `:sp`, `:prec_sp` | AMG solver/preconditioner params | *removed* | Krylov/AMG path not used |
| `:matrix_storage` | `:assembled` | *removed* (defaulted) | lean integrator solves directly |
| `:linitial_refine`, `:init_refine_lvl` | present (refine 10×10 once) | *removed* | not a newmaster input; baseline `case2d` runs the 10×10 mesh directly |
| `:lvisc` | `true` | `true` | **kept** — see §6 |
| `:method`, `:delta`, `:k`, `:coeff`, `:S_fun`, `:L_fun`, `:build_L`, `:lsolve` | present | present | unchanged config |

`_imex_L` builds `-μ[1] * (Minv .* K)`, where `K` is the global Galerkin
stiffness. Crucially it does **not** re-assemble `K` with
`build_laplace_matrix` / `DSS_laplace_sparse` at run time — that assembly is
extremely allocation-heavy (tens of millions of node-pair allocations, ~8 GiB
for this small case). Instead it reuses `params.Lap_sparse`, the **same** `K`
that `matrix_wrapper` already assembled once during `sem_setup` (because
`:ldss_laplace => true`). `params_setup` exposes `sem.matrix.L` as
`params.Lap_sparse` for exactly this purpose. With this, the IMEX run's
allocation is dominated by `rhs!` and VTK output (tens of MiB), not operator
assembly.

---

## 6. Two subtleties worth remembering

1. **`:lvisc => true` is load-bearing.** The final step update adds only the
   *explicit* combination `Δt Σ b_ex[i] S(U_i)`. The implicit `L` enters the
   per-stage solves but not the final `u` update. So the diffusion that
   actually reaches `u` comes through the explicit `S` (which includes the
   viscous term when `:lvisc => true`); the implicit operator's role is to
   keep the stage solves stable. If you set `:lvisc => false`, `u` would get
   **no** diffusion and the Riemann fronts would ring. This matches the
   original `sl/imex` behaviour; it is intentionally preserved.

2. **`_sl` ≠ semi-Lagrangian** (repeated because it is easy to assume
   otherwise). It is the *generic Krylov-style* IMEX path, here run with a
   direct solve.

---

## 7. How to run

```bash
# from the repo root, on branch sl/imex_newmaster
julia --project=. src/Jexpresso.jl Burgers case2d_imex_sl
# or in the REPL:
#   using Jexpresso; Jexpresso.run_case("Burgers", "case2d_imex_sl")
```

Requirements unchanged from `case2d`: the doubly-periodic mesh
`./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh` (the `meshes/` tree is
git-ignored, so it is provided/generated locally, not committed).

---

## 8. Suggested next steps

* Run `case2d_imex_sl` and compare the output against the `sl/imex` run to
  confirm the migration is numerically faithful (the porting was verified by
  inspection against the newmaster APIs; it still needs a runtime check).
* If an iterative/preconditioned solver is desired, reintroduce it behind a
  custom `:lsolve` closure (and, if AMG is wanted, add `AlgebraicMultigrid` /
  `IncompleteLU` back to `Project.toml`) rather than as a second code path.
* Optionally add a `test/CI-runs/Burgers/case2d_imex_sl/` mirror so CI
  exercises the IMEX path.
