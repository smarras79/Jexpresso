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

```bash
# On a node with an NVIDIA GPU, CUDA.jl available, and JACC configured for CUDA:
julia --project=. -e 'using CUDA; using JACC; JACC.set_backend("cuda")' \
                  src/Jexpresso.jl Burgers case2d_imex_sl_gpu
# or in the REPL:
#   using CUDA, JACC; JACC.set_backend("cuda")
#   using Jexpresso; Jexpresso.run_case("Burgers", "case2d_imex_sl_gpu")
```

Requirements:
* `CUDA.jl` loaded in the session (so `CUDABackend` resolves) and a functional
  NVIDIA GPU;
* JACC configured for CUDA (a `LocalPreferences.toml` selecting the CUDA backend,
  or `JACC.set_backend("cuda")`);
* the doubly-periodic mesh `./meshes/gmsh_grids/hexa_TFI_10x10_burgers2d.msh`
  (the `meshes/` tree is git-ignored, so it is provided/generated locally).

For AMD GPUs, replace `CUDABackend()` with `ROCBackend()` in `user_inputs.jl`
and use `JACC.set_backend("amdgpu")` (AMDGPU.jl).

## Kernel-level tests (no full run needed)

* CPU: `julia --project test/test_imex_jacc.jl` (also part of `test/runtests.jl`)
* GPU: `julia --project test/test_imex_jacc_gpu.jl`
