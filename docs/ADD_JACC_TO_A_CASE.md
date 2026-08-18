# Adding JACC (multi-threading / GPU) to an existing case

This walks through turning an existing Jexpresso case into one that can evaluate
its right-hand side through [JACC.jl](https://github.com/JuliaORNL/JACC.jl) — one
kernel body that runs on multi-core CPU, NVIDIA, AMD, Intel or Apple hardware.

The worked example is **`problems/CompEuler/theta`** (the 2D rising thermal
bubble), which needed exactly one line changed. Background on the design lives in
`problems/JACCtestz/README.md` and in the header of
`src/kernel/operators/rhs_jacc.jl`.

> **Read step 0 before doing any of this.** For a small case the answer may well
> be "this will not make your run faster", and it is better to know that first.

---

## Step 0 — decide whether it is worth it

The JACC path replaces the **right-hand side**, and nothing else. It does not
touch mesh generation, the metric terms, the time integrator, the I/O, or the
JIT. So it can only speed up the fraction of your wall clock that is spent in
`rhs!`.

`CompEuler/theta` is a case where the answer is **no**, and it is worth seeing why
before you spend an afternoon on your own deck. Measured, 100 elements, 1681
nodes, 4 equations, one RHS evaluation:

| threads | CPU RHS | JACC RHS | JACC speedup |
|---------|---------|----------|--------------|
| 1 | 0.718 ms | 3.107 ms | **0.23×** — four times *slower* |
| 2 | 0.821 ms | 0.661 ms | 1.24× |
| 4 | 0.784 ms | 0.386 ms | 2.03× |

Now put that in context. The shipped deck is `tend = 1000`, `Δt = 0.5`, and
CarpenterKennedy2N54 has five stages — 2000 steps × 5 = 10 000 RHS calls:

| | total time in the RHS |
|---|---|
| CPU | 10 000 × 0.821 ms ≈ **8.2 s** |
| JACC, 2 threads | 10 000 × 0.661 ms ≈ **6.6 s** |

A 1.6 s difference, inside a run whose wall clock is dominated by the mesh read,
the metric terms, JIT and eleven VTK writes. Switch `:ljacc` on for this case with
`-t 2` and **the wall clock will not visibly move**. That is the expected result,
not a symptom of anything being wrong.

So, two things routinely make the answer "no":

* **The case is small.** 100 elements is not enough work per `parallel_for` to pay
  for the thread barriers. JACC starts paying when the element count is large
  enough that each kernel launch has real work in it.
* **The run is short.** A case that takes 40 s of which 35 s is compilation will
  look identical whatever you do to the RHS.

Time the RHS before you invest — not the wall clock of the whole process. The
honest check is to compare *the time-integration line* Jexpresso prints with
`:ljacc` off and on, or to call `rhs!` in a loop yourself.

And note the first row of that table: at **one** thread the JACC path is four
times slower than the serial CPU path. See step 4 — this is the single most
common way to get a disappointing result.

JACC is also the only way to reach a GPU at all, which for a large enough case is
the real reason to do this.

---

## Step 1 — check the case is covered

`jacc_check_inputs` (in `rhs_jacc.jl`) refuses, **by name and at setup**, anything
the JACC path does not implement. Your case must be:

| requirement | why |
|---|---|
| 2D (`NSD_2D`) | 1D/3D kernels are not written yet |
| `:backend => CPU()` (the default) | that key selects the separate KernelAbstractions path; the two must not both be on |
| not `:lspherical_shell` | the shell has its own RHS and time loop |
| no Laguerre elements | not implemented |
| `:lfilter => false` | not implemented |
| `:lmoist => false` | not implemented |
| `:ladapt`/`:lamr` false | the device cache is staged once, at setup |
| `:visc_model => AV()` if `:lvisc` | only constant artificial viscosity is implemented |

One restriction is **not** checked, because it cannot be: a case whose
`user_bc_neumann` returns something non-zero must not use this path. The JACC
kernels have no Neumann/surface-integral term.

`CompEuler/theta` passes: 2D, `AV()` viscosity, no filter, no moisture, no AMR.

---

## Step 2 — make sure the four device callbacks exist

The kernels call the *device* spelling of the user callbacks — the same ones the
KernelAbstractions path uses, which return values instead of writing through a
pointer:

| callback | needed when | lives in |
|---|---|---|
| `user_flux_gpu(q, qe, PhysConst, lpert)` | always | `user_flux.jl` |
| `user_source_gpu(q, qe, x, y, PhysConst, xmax, xmin, ymax, ymin, lpert)` | `:lsource` | `user_source.jl` |
| `user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)` | mesh has boundary edges | `user_bc.jl` |
| `user_primitives_gpu(u, qe, lpert)` | `:lvisc` | `user_primitives.jl` |

Check yours:

```bash
cd problems/CompEuler/theta
for f in user_flux_gpu user_source_gpu user_primitives_gpu user_bc_dirichlet_gpu; do
    printf "%-24s %s\n" "$f" "$(grep -l "function $f" *.jl)"
done
```

### `user_bc_dirichlet_gpu` takes `coords`, not `x, y, z`

The kernels hand the callback **one node's coordinates as a column view**,
`@view(coords[:, ip])` — `mesh.coords` is `(nsd, npoin)` precisely so a node costs
one cache line instead of `nsd` of them. Index it positionally:

```julia
function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, qbdy, lpert)      # 2D
    x, y = coords[1], coords[2]
    ...
```

```julia
function user_bc_dirichlet_gpu(q, qe, coords, t, nx, ny, nz, qbdy, lpert)  # 3D
    x, y, z = coords[1], coords[2], coords[3]
    ...
```

1D is `(q, qe, coords, t, lpert)`.

Getting this wrong is not a compile error and not a test failure — these callbacks
are resolved inside a kernel at run time, so a mismatched signature is a
`MethodError` the first time somebody points a GPU at that deck. `test/
test_gpu_callback_signatures.jl` checks every case in the tree for exactly this,
and is worth running after touching any `user_bc.jl`:

```bash
julia --project=. test/test_gpu_callback_signatures.jl
```

Two conventions inside that callback are worth knowing:

* `qbdy` arrives pre-filled with a sentinel. Returning it in a slot means "leave
  this component alone" — that is how `theta` leaves density and `ρθ` free at a
  wall and constrains only the momentum:
  `return T(qbdy[1]), T(u), T(v), T(qbdy[4])`.
* Whether a returned value is actually written is decided by the same
  `AlmostEqual` test the CPU path uses (Kopriva alg. 139, **ε = 1e-6, absolute**).
  Corrections smaller than that are declined, by design, to match the CPU.

---

## Step 3 — turn it on in the deck

```julia
:ljacc                => true,
```

Leave `:backend` at its `CPU()` default. That key drives the separate
KernelAbstractions path, and `rhs_jacc.jl` refuses to start if both are on.

Optionally pick the residual precision (see step 6):

```julia
:jacc_float           => :auto,     # default; Float64 where the backend has it
```

---

## Step 4 — give Julia some threads

**`:ljacc => true` does not give you threads.** It chooses which kernels run;
Julia still starts with **one** thread unless you say otherwise, and the count is
fixed at process startup — there is no way to raise it from inside a running
REPL.

```bash
julia --project=. -t 8 src/Jexpresso.jl CompEuler theta      # script form
```

```bash
julia --project=. -t auto                                    # REPL form
julia> using Jexpresso
julia> Jexpresso.run_case("CompEuler", "theta")
```

or `export JULIA_NUM_THREADS=auto` before launching. In VS Code, set
`"julia.NumThreads": "auto"`.

Running the JACC path on **one** thread is worse than not using it at all: JACC's
threads backend runs a plain loop at `nthreads() == 1` and `Polyester.@batch`
above it, and the two are not the same code. On `CompEuler/theta`'s RHS that is
3.107 ms against 0.661 ms — and against 0.821 ms for the serial CPU path it never
left. The 1→2 jump is the code-path switch, not the second thread.

Jexpresso warns if it catches you at one thread with `:ljacc` on.

---

## Step 5 — confirm it is actually running, and that it is right

The banner appears during `params_setup`:

```
 # JACC RHS ........................ backend = threads  (8 threads),  arrays = Array
 #   npoin = 1681,  nelem = 100,  ngl = 5,  neqs = 4,  viscous = true
 #   residual precision = Float64
```

If you do not see it, `:ljacc` did not reach the deck that was loaded.

Then check the answer, which is the part that matters. Run the case once with
`:ljacc => false` and once with `true`, and diff the output. The JACC RHS agrees
with the serial CPU RHS **to floating-point round-off** — not bit for bit, because
the direct stiffness summation adds the same numbers in a different order (a
gather per node instead of a scatter per element) and floating-point addition is
not associative.

Measured on `CompEuler/theta`, one RHS evaluation:

```
max|Δ| cpu vs jacc = 3.3e-16     (max|rhs| = 1.4e-01)
```

and on `ShallowWater/SoliWaveIsland`, five steps of SSPRK54 through the full
solver, HDF5 output diffed: 4.4e-16 relative on H, 7.7e-16 on Hu.

Anything larger than round-off means something is genuinely different, and the
place to look first is the boundary condition — see step 2.

---

## Step 6 — a GPU, and the Float32 question

The device is **not** a deck setting. JACC reads it from a `Preferences` value,
written once per project, which needs a restart:

```bash
julia --project=. -e 'using JACC; JACC.set_backend("cuda")'   # amdgpu|oneapi|metal|threads
```

Then run exactly the same commands. `set_backend` adds the vendor package
(CUDA.jl, AMDGPU.jl, …) to `Project.toml`; think before committing that.

**Apple GPUs have no double precision.** Metal.jl's array constructor refuses
Float64 outright and JACC declares `default_float(::MetalBackend) = Float32`, so
on an M-series Mac the residual must be single precision. `:jacc_float` controls
this:

| value | meaning |
|---|---|
| `:auto` (default) | Float64 where the backend has it, the backend's default where it does not, with a warning |
| `Float32` | opt in explicitly |
| `Float64` | demand it; error naming the backend if it has none |

`Float32` is **mixed** precision: the state, the time integrator and the MPI
assembly stay Float64 and only the residual is evaluated in single, so the error
is O(eps32) in the tendency rather than accumulating in the solution. On
SoliWaveIsland, five steps: 1.6e-07 relative on H, about 1.3 × eps32, not
compounding.

It is a **portability** feature, not a speed one — on a CPU it is slower — and
four significant digits of a small perturbation is fine for a demo and not
something to trust for a vorticity-sensitive run.

---

## Step 7 — if you had to write the `_gpu` callbacks from scratch

Rules for anything that runs inside a kernel:

* **Return values, do not allocate.** Return a tuple; never build an array.
* **No `Float64` literals if you want Metal.** Use `T = eltype(q)` and `T(0.5)`.
  Every `_gpu` callback in the tree already opens with `T = eltype(q)` — follow
  that and the same source works in both precisions.
* **No closures, no captured globals, no `try`/`catch`, no I/O, no `error()`.**
* **Branches are fine**, but remember they make the residual non-smooth: a
  threshold test that two implementations land on opposite sides of will make
  their trajectories diverge. Prefer smooth desingularisation where you can.
* **`PhysConst` follows the working precision** — the cache hands you a
  `PhysicalConst{Float32}` in Float32 mode, so do not assume Float64 fields.

Fastest way to iterate: `problems/JACCtestz/JACC_2d_swe.jl` runs the kernels
standalone in seconds with no MPI, no gmsh and no `using Jexpresso`, and
`test/test_jacc_rhs.jl` is the same checks as a testset.

---

## Checklist

```
[ ] step 0  timed the RHS and it is worth doing
[ ] step 1  case is 2D, AV() viscosity, no filter/moisture/AMR/Laguerre/shell
[ ] step 1  user_bc_neumann returns zero
[ ] step 2  the four *_gpu callbacks exist
[ ] step 2  user_bc_dirichlet_gpu takes (q, qe, coords, t, nx, ny[, nz], qbdy, lpert)
[ ] step 2  julia --project=. test/test_gpu_callback_signatures.jl passes
[ ] step 3  :ljacc => true, :backend left at CPU()
[ ] step 4  julia -t <n>
[ ] step 5  banner appears, and the answer matches the CPU run to round-off
```
