# Burgers `case2d_imex_sl_gpu` — GPU twin of `case2d_imex_sl`

Same 2-D viscous Burgers problem and the same native IMEX integrator as
`../case2d_imex_sl`, but configured to run on the GPU through JACC.jl.

The two directories are identical except for `user_inputs.jl`:

| | `case2d_imex_sl` (CPU) | `case2d_imex_sl_gpu` (GPU) |
|---|---|---|
| `:backend` | unset → `CPU()` | `CUDABackend()` |
| `:limex_jacc` | `true` | `true` |

With `:limex_jacc => true` the per-stage implicit solve `(I - λL) x = b` runs
on the device via the portable JACC BiCGSTAB + SpMV
(`src/kernel/solvers/imex_jacc.jl`) instead of a host sparse LU. Setting
`:backend => CUDABackend()` additionally puts `u`, the explicit `rhs!` and all
vector algebra on the GPU.

## Run

On a node with an NVIDIA GPU and CUDA.jl in the project, **load CUDA before the
run** with `Jexpresso.enable_cuda!()` (this is required — see the note below):

```julia
using Jexpresso
Jexpresso.enable_cuda!()                        # loads CUDABackend + JACC CUDA backend
Jexpresso.run_case("Burgers", "case2d_imex_sl_gpu")
```

### Why `enable_cuda!()` and not just `using CUDA`?

`CUDABackend` is defined by CUDA.jl, an optional dependency Jexpresso does not
load by default. A case file cannot load CUDA lazily and construct
`CUDABackend()` in the same call — Julia rejects calling a constructor that was
added to the method table *after* the running code started (a "world-age"
error). `enable_cuda!()` loads CUDA ahead of the run, into the `Jexpresso`
module, and selectively imports only `CUDABackend` (so it doesn't clash with
`Jexpresso.CG`). It also nudges JACC onto its CUDA backend.

Requirements:
* `CUDA.jl` in the active project — add it on the GPU machine with `] add CUDA`;
* a functional NVIDIA GPU (`using CUDA; CUDA.functional()`);
* the doubly-periodic mesh `./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh`
  (the `meshes/` tree is git-ignored, so it is provided/generated locally).

For AMD GPUs, replace `CUDABackend()` with `ROCBackend()` in `user_inputs.jl`,
add AMDGPU.jl, and call `Jexpresso.enable_amdgpu!()` instead.

## Kernel-level tests (no full run needed)

* CPU: `julia --project test/test_imex_jacc.jl` (also part of `test/runtests.jl`)
* GPU: `julia --project test/test_imex_jacc_gpu.jl`
