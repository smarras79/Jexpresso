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

On a node with an NVIDIA GPU and CUDA.jl in the project:

```julia
using Jexpresso
Jexpresso.run_case("Burgers", "case2d_imex_sl_gpu_hybrid"; backend = :cuda)
```

`backend = :cuda` loads CUDA and points JACC at the GPU; the case itself keeps
`:backend => CPU()`, so the Jexpresso compute backend stays on the host. For AMD
GPUs use `backend = :amdgpu`. If CUDA is not loaded (JACC on its CPU backend) the
solve simply runs on the CPU — still correct.

Requirements: CUDA.jl in the project (`] add CUDA`), a functional NVIDIA GPU
(`using CUDA; CUDA.functional()`), and the mesh
`./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh`.

## Performance note

At this problem size (~5k points) the GPU solve will **not** necessarily beat the
CPU's cached sparse LU — host↔device transfer per stage and a small SpMV dominate.
The offload pays off on large meshes / many DOF, where the linear solve is the
bottleneck. Correctness is identical to the host LU solve (verified to ~1e-9 in
`test/test_imex_jacc.jl`).
