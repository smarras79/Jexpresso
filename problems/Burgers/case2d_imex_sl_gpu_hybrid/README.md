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

It runs the case once per precision `(Float16, Float32, Float64)`, prints a
comparison table (BiCGSTAB iterations, residual norms, non-converged solves,
solve wall-time, final ‖u‖₂ and its drift from the double run) and writes
`imex_precision_study.csv`. Each run also prints its own diagnostics block at the
end. Note: pure `Float16` BiCGSTAB often will **not** reach a tight tolerance —
that (large residual / non-converged count) is itself the diagnostic.

## Performance note

At this problem size (~5k points) the GPU solve will **not** necessarily beat the
CPU's cached sparse LU — host↔device transfer per stage and a small SpMV dominate.
The offload pays off on large meshes / many DOF, where the linear solve is the
bottleneck. Correctness is identical to the host LU solve (verified to ~1e-9 in
`test/test_imex_jacc.jl`).
