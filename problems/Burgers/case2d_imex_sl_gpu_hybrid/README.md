# Burgers `case2d_imex_sl_gpu_hybrid` — host pipeline + GPU implicit solve

Same 2-D viscous Burgers IMEX problem as `../case2d_imex_sl`, configured as a
**hybrid**: the whole Jexpresso pipeline (mesh, `sem_setup`, the explicit `rhs!`,
all vector algebra) runs on the **host**, and **only the per-stage implicit
solve** `(I - λL) x = b` is offloaded to the **GPU** via JACC.

Why hybrid: `:backend => CUDABackend()` turns every array into a `CuArray`, but
large parts of Jexpresso's setup/RHS are not yet GPU-ported (they do host-style
scalar indexing). The hybrid keeps `:backend => CPU()` so all of that runs on the
proven host path, while the one piece your request targeted — the IMEX implicit
solver — runs on the device.

| | mesh / `sem_setup` / `rhs!` | implicit solve `(I-λL)x=b` |
|---|---|---|
| Hybrid (this case) | **host (CPU)** | **GPU (JACC BiCGSTAB)** |
| `case2d_imex_sl` | host (CPU) | host (CPU JACC / LU) |
| `case2d_imex_sl_gpu` | GPU (WIP, not fully ported) | GPU |

How it works (`_imex_rk_run_const_jacc_offload!` in
`src/kernel/solvers/IMEXTimeIntegrators.jl`): the per-stage operators `I - λᵢL`
are assembled on the host and uploaded to the device as CSR **once**; each stage
copies its small `rhs` vector host→device, solves with `jacc_bicgstab!` on the
device, and copies the solution device→host. `S(U)` is the host `rhs!`; `L(U)` is
a host sparse mat-vec.

## Run

> **JACC's backend is a compile-time preference.** It is read when JACC loads, so
> setting it to CUDA requires a **Julia restart** — a runtime `set_backend` does
> not switch the current session. This is the #1 reason the solve "runs on the
> CPU" despite a GPU being present.

**One-time setup** (writes `LocalPreferences.toml`):

```julia
using CUDA, JACC
JACC.set_backend("cuda")
```
**Then exit and restart Julia.** In the fresh session:

```julia
using CUDA, JACC
using Jexpresso
Jexpresso.jacc_status()        # MUST show a CuArray type and "GPU? YES"
Jexpresso.run_case("Burgers", "case2d_imex_sl_gpu_hybrid")
```

No `backend = :cuda` kwarg is needed: this case keeps `:backend => CPU()`, so only
JACC needs to be on its CUDA backend. The run's banner prints the solve array type
(`CuArray` ⇒ GPU, `Vector` ⇒ CPU) and warns if it fell back to the CPU. For AMD
GPUs use `JACC.set_backend("amdgpu")` + AMDGPU.

Requirements: CUDA.jl in the project (`] add CUDA`), a functional NVIDIA GPU
(`using CUDA; CUDA.functional()`), and the mesh
`./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh`.

## Precision study (half / single / double)

Only the offloaded solve's precision is configurable (`:imex_jacc_solve_precision`,
default = host precision); the host pipeline always runs in `TFloat`. To sweep
half/single/double and collect diagnostics in one call:

```julia
using CUDA, JACC; using Jexpresso          # JACC already on its CUDA backend (restarted)
Jexpresso.run_imex_precision_study("Burgers", "case2d_imex_sl_gpu_hybrid")
```

It runs the case once per precision `(Float16, Float32, Float64)` and prints a
comparison table plus writes `imex_precision_study.csv`. Columns:

| column | meaning |
|---|---|
| `nsolve`, `avg_it` | # device solves and average BiCGSTAB iterations/solve |
| `mean_res`, `ncnv` | mean residual ‖b−Ax‖ and # solves that didn't reach tol |
| `solve_s`, `us/solve` | **device implicit-solve time** (host rhs!/assembly/output excluded) |
| `speedup` | (double solve time) / (this solve time) — the **gain** vs double |
| `relerr` | ‖u_prec − u_f64‖₂ / ‖u_f64‖₂ — the **accuracy cost** of reduced precision |

```
 # ===== IMEX/JACC implicit-solve precision study  (reference = Float64) =====
  prec      dev  nsolve  avg_it    mean_res  ncnv   solve_s  us/solve  speedup        relerr
  Float16   GPU    1500   40.00   2.000e-03  1500     2.100   1400.00    0.45x     3.529e-02
  Float32   GPU    1500    5.00   7.000e-07     0     0.450    300.00    2.11x     1.721e-06
  Float64   GPU    1500    5.00   4.000e-14     0     0.950    633.33    1.00x     0.000e+00
```
(numbers above are illustrative). Each run also prints its own diagnostics block.
Note: pure `Float16` BiCGSTAB usually will **not** reach a tight tolerance — the
large residual / `ncnv` count / `relerr` is itself the diagnostic, and it may be
*slower* than double because it iterates to `itmax`. `Float32` is typically the
useful operating point (large speedup, small `relerr`). To set the precision for
a single run instead of the sweep, put `:imex_jacc_solve_precision => Float32` in
`user_inputs.jl`.

## Performance note

At this problem size (~5k points) the GPU solve will **not** necessarily beat the
CPU's cached sparse LU — host↔device transfer per stage and a small SpMV dominate.
The offload pays off on large meshes / many DOF, where the linear solve is the
bottleneck. Correctness is identical to the host LU solve (verified to ~1e-9 in
`test/test_imex_jacc.jl`).
