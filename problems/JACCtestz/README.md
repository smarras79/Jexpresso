# JACC.jl in Jexpresso

One kernel body, five backends. `src/kernel/operators/rhs_jacc.jl` evaluates the
2D continuous-Galerkin spectral-element right-hand side through
[JACC.jl](https://github.com/JuliaORNL/JACC.jl), so the same source runs on

| backend   | hardware                    |
|-----------|-----------------------------|
| `threads` | multi-core CPU (the default)|
| `cuda`    | NVIDIA                      |
| `amdgpu`  | AMD                         |
| `oneapi`  | Intel                       |
| `metal`   | Apple                       |

The first case wired up is the 2D non-linear shallow water solitary wave running
up a conical island, `problems/ShallowWater/SoliWaveIsland` — the 2D warm-up for
the shallow water equations on the sphere.

---

## Running it

### 1. The standalone check (start here)

Needs JACC and nothing else — no MPI, no Gridap, no gmsh, no `using Jexpresso`.
Seconds to run, and the right first thing to try when bringing a backend up on a
new machine.

```
julia --project=. -t 8 problems/JACCtestz/JACC_2d_swe.jl
julia --project=. -t 8 problems/JACCtestz/JACC_2d_swe.jl 32 16 5   # nelx nely nop
```

It builds its own curvilinear SEM mesh, then checks that

1. the lake at rest is an exact steady state of the discrete operator
   (well-balancedness of the flux/source split — a statement about the
   discretisation, not just the PDE),
2. the JACC right-hand side equals an independently written serial evaluation of
   the same weak form, inviscid and viscous, **to the last bit**,
3. two evaluations are bitwise identical,

and reports throughput. The `@testset` version of the same checks, plus the
connectivity-map and free-slip-wall assertions, is

```
julia --project=. -t 4 test/test_jacc_rhs.jl
```

### 2. The real case

Flip one switch in `problems/ShallowWater/SoliWaveIsland/user_inputs.jl`:

```julia
:ljacc => true,
```

and run the case as usual, asking for threads on the command line:

```
julia --project=. -t 8 src/Jexpresso.jl ShallowWater SoliWaveIsland
```

Everything else is unchanged: same mesh, same `Float64` setup, same SSPRK54, same
output files.

Measured end to end on this deck — the case run twice in one session, serial CPU
RHS then JACC RHS, five steps of SSPRK54, HDF5 output diffed:

| field | max\|q\| | max\|Δ\| | relative |
|-------|---------|---------|----------|
| H     | 3.82e-01 | 1.67e-16 | 4.4e-16 |
| Hu    | 1.21e-01 | 9.37e-17 | 7.7e-16 |
| Hv    | 5.89e-10 (zero by symmetry) | 7.75e-17 | — |

i.e. agreement to floating-point round-off. Not bit-for-bit, and it cannot be:
the direct stiffness summation adds the same numbers in a different order (gather
per node instead of scatter per element), and floating-point addition is not
associative. What *is* bit-for-bit is any two JACC runs of the same deck — see
"No atomics" below.

### 3. On a GPU

The backend is **not** a deck setting. JACC reads it from a `Preferences` value,
so it is written once per project and needs a restart:

```
julia --project=. -e 'using JACC; JACC.set_backend("cuda")'   # amdgpu|oneapi|metal|threads
```

Then run exactly the same commands as above. `JACC.set_backend` adds the vendor
package (CUDA.jl, AMDGPU.jl, …) to `Project.toml` — think before committing that
change.

Leave `:backend` at its `CPU()` default. That key selects the *separate*
KernelAbstractions path in `rhs_gpu.jl`, which also flips `TFloat`/`TInt` to 32
bit and moves the entire setup onto the device; the two paths must not both be
on, and `rhs_jacc.jl` refuses to start if they are.

---

## How it is put together

```
u (host, from OrdinaryDiffEq)
  │  copyto!
  ▼
u ─► uaux                                    parallel_for(npoin)
  ├─► Dirichlet / free-slip boundary values  parallel_for(nbnode)
  ├─► F, G, S and the primitive variables    parallel_for(nelem, ngl, ngl)
  ├─► ν∇q in the contravariant basis         parallel_for(nelem, ngl, ngl)   [viscous only]
  ├─► element residual  -ωJ(∇·F - S) - ∫∇ψ·(ν∇q)
  │                                          parallel_for(nelem, ngl, ngl)
  └─► direct stiffness summation + M⁻¹       parallel_for(npoin)
  │  copyto!
  ▼
du (host)   ── cross-rank DSS on the host (a no-op on one rank)
```

Six launches, no host round trip in between.

### No atomics, on purpose

The KernelAbstractions kernels do the direct stiffness summation with
`@atomic RHS[ip,ieq] -= …` from every `(element, node)` thread at once. Correct,
but the floating-point summation order then depends on how threads happened to be
scheduled: two runs of the same deck do not agree bit for bit, and a GPU run
cannot be compared against the CPU reference except by eyeballing tolerances.

Here it is inverted. The element kernels write **only their own slot** of
`rhs_el[iel,i,j,ieq]` — no two threads share an output — and a second pass,
parallelised over global nodes, gathers through a point→(element,i,j) map built
once at setup:

```
RHS[ip,ieq] = Minv[ip] * Σ_{(iel,i,j) ∈ p2e[ip]} rhs_el[iel,i,j,ieq]
```

Race-free by construction, identical on every backend and at every thread count,
and it costs one extra `npoin × neqs` pass plus an integer CSR map. Two JACC runs
of the same deck therefore agree to the last bit, on CPU and GPU alike; against
the CPU path they agree to round-off, since the summation order differs.

The boundary condition is handled the same way: one thread per **distinct
boundary node**, applying each of that node's edges in a fixed order. That also
removes a race in the per-edge spelling, where a corner node belongs to two edges
and whichever normal landed last decided the corner's velocity.

### Matching the CPU boundary contract

Two details of `build_custom_bcs_dirichlet!` (BCs.jl) are not incidental, and
getting either wrong moves the answer by ~1e-3 within the first time step:

* the corrected field is pushed back into the integrator's state — the builder's
  last statement is `uaux2u!(u, uaux, neqs, npoin)`, which contradicts the
  `# WARNING: ... applied to uaux[:,:] and NOT u[:]!` comment at its top;
* whether to write a Dirichlet value at all is decided by `AlmostEqual`
  (Kopriva alg. 139) with **ε = 1e-6, absolute** — so any wall correction smaller
  than about a microsecond-scale momentum is declined. An exact `!=` here enforces
  a slightly different wall, and the two trajectories separate.

### What a case needs

The kernels call the *device* spelling of the user callbacks, the same ones the
KernelAbstractions path uses:

| callback                  | needed              |
|---------------------------|---------------------|
| `user_flux_gpu`           | always              |
| `user_source_gpu`         | when `:lsource`     |
| `user_bc_dirichlet_gpu`   | when the mesh has boundary edges |
| `user_primitives_gpu`     | when `:lvisc`       |

In 2D, `user_bc_dirichlet_gpu` takes `(q, qe, x, y, t, nx, ny, qbdy, lpert)` —
the two coordinates **separately**. `qbdy` arrives pre-filled with a sentinel;
returning it in a slot means "leave this component alone" (that is how the depth
`H` is left free at a free-slip wall).

---

## Coverage

Supported:

- 2D flat CG-SEM, any `neqs`, `TOTAL()` or `PERT()`
- volume flux divergence + user source, well-balanced
- constant artificial viscosity, `:visc_model => AV()`
- Dirichlet / free-slip boundaries
- MPI (the cross-rank DSS runs on the host, exactly as the KA path does)

Not supported — and refused **by name at setup**, loudly, rather than quietly
solving a different equation:

- 1D / 3D
- the spherical shell (`:lspherical_shell`)
- Laguerre semi-infinite elements
- the CG filter (`:lfilter`)
- microphysics (`:lmoist`)
- AMR (`:ladapt`, `:lamr`)
- viscosity models other than `AV()`
- Neumann / surface-integral boundary conditions. This is the one restriction
  that is **not** checked: whether `user_bc_neumann` returns something non-zero
  can only be known by calling it. A case with a real Neumann condition must not
  set `:ljacc`.

## Next: the shallow water equations on the sphere

`problems/ShallowWater/SWsphere` does **not** go through `rhs!`. It has its own
right-hand side (`sphere_rhs!`, `src/kernel/operators/sphere_rhs.jl`) and its own
time loop (`src/kernel/solvers/sphere_time_loop.jl`), so it needs its own JACC
entry point rather than the `:ljacc` switch.

The port is the same three kernels as here, with one change of substance: the
shell's metric is 3×2, so `a¹` and `a²` each carry a z component and the surface
divergence is `dFdx + dGdy + dHdz` with a third flux component `H`. Concretely:

- `_jacc_element_rhs_2d!` gains one term and one metric triple;
- the viscous kernel becomes Laplace-Beltrami — the same weak form, with
  `gξ = ν∇ₛq·a¹`, `gη = ν∇ₛq·a²` over three components instead of two
  (`_sphere_visc_el!` is the serial original to match);
- the CSR gather, the cache and the driver are unchanged;
- a closed shell has no boundary, so `nbnode` is zero and the boundary kernel
  never launches;
- the Lagrange projection applied after every RK stage
  (`_sphere_stage_limiter!`) and the modal filter (`sphere_filter!`) are the two
  remaining pieces that would also want porting, or the state has to come back to
  the host every stage and the transfers eat the win.

## A measured note on the `threads` backend

JACC's threads backend takes a different code path when `Threads.nthreads() == 1`
(a plain loop) than when it is greater (`Polyester.@batch`). On this test the
single-thread path is several times slower than the two-thread one — far more
than two threads can explain. Run with `-t 2` or more; a single-threaded timing
is not a meaningful serial baseline.
