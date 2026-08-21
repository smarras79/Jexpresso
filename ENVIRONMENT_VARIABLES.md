# Environment Variables

This file lists every environment variable Jexpresso reads, plus the
third-party ones (Julia, MPI, BLAS) that affect Jexpresso runs in
ways worth documenting. A final section covers the `user_inputs.jl`
switches that change **MPI behaviour** rather than physics — they have
no environment-variable form, but they belong next to these because
they are the other place a run's parallel cost is decided.

Each entry covers: meaning, accepted values, default, the file the
variable is read in, and the equivalent `user_inputs.jl` key or CLI
flag (if any). Precedence is shown explicitly when a variable can be
overridden in multiple ways.

For all boolean-style variables below, **accepted truthy values** are
`1`, `true`, `yes`, `on` (case-insensitive). **Accepted falsy values**
are `0`, `false`, `no`, `off` (case-insensitive). Anything else is
treated as the default.

---

## Jexpresso-specific

### `JEXPRESSO_COUPLED`

Enables MPMD coupling with Alya. When set, Jexpresso splits
`MPI.COMM_WORLD` by `APPID` so it can share a world communicator with
another MPI code (typically Alya for atmosphere–CFD coupling). When
unset (the default), Jexpresso is byte-for-byte the historical
standalone code path.

- **Type:** boolean
- **Default:** unset (standalone, no coupling)
- **Read in:** `src/run.jl`
- **Example:**
  ```bash
  JEXPRESSO_COUPLED=1 mpirun -np 4 julia --project=. src/Jexpresso.jl CompEuler 3dAlya
  ```
- **Related:** `APPID` (below); `RUN-COUPLED.md` for the full
  Alya↔Jexpresso launch recipe.

### `APPID`

Application identifier used to split `MPI.COMM_WORLD` when running
under MPMD with Alya. Each peer in the MPMD launch passes a different
`APPID` so MPI's `MPI_Comm_split` can hand each app its own
sub-communicator.

- **Type:** integer
- **Default:** `2` (Jexpresso); Alya conventionally uses `1`
- **Read in:** `src/kernel/coupling/couplingStructs.jl`
- **Active only when:** `JEXPRESSO_COUPLED=1`
- **Example:** see `RUN-COUPLED.md`.

### `JEXPRESSO_ALLOC_SUMMARY`

Enables the end-of-run per-function timing and allocation summary
table (`TimerOutputs`). Off by default everywhere, CI included: it
forces one full extra RHS warm-up call so the post-JIT measurement
window is meaningful, which is a cost only a performance run should
pay.

- **Type:** boolean
- **Default:** `false` (off)
- **Read in:** `src/kernel/solvers/TimeIntegrators.jl`
  (`alloc_summary_enabled`), `src/auxiliary/timing.jl`
- **Precedence (highest first):**
  1. `JEXPRESSO_ALLOC_SUMMARY` env var
  2. `--no-alloc-summary` CLI flag in `ARGS` (the bare
     `no-alloc-summary` is accepted too)
  3. `:lalloc_summary => true` in `user_inputs.jl`
  4. Default (`false`)
- **Example:**
  ```bash
  JEXPRESSO_ALLOC_SUMMARY=1 mpirun -np 4 julia --project=. src/Jexpresso.jl CompEuler theta
  ```
  For coupled-mode mpirun pass through `mpirun -x`:
  ```bash
  mpirun -np 2 ./AlyaProxy/Alya.x : \
         -x JEXPRESSO_COUPLED=1 -x JEXPRESSO_ALLOC_SUMMARY=1 \
         -np 2 julia --project=. src/Jexpresso.jl CompEuler 3dAlya
  ```

### `JEXPRESSO_PRECOMPILE_WARMUP`

Controls whether `time_loop!` runs a one-shot RHS warm-up call (a
single fake timestep) before the real `solve(...)`. The warm-up
triggers JIT compilation of `rhs!`, `_build_rhs!`, and the rest of
the per-step kernel chain on an arbitrary initial-condition state, so
the production run starts already compiled.

Useful when launching from the command line, where every run incurs
JIT-compile cost on the first timestep; REPL users get the same
benefit by re-running inside the same session and don't need this.

- **Type:** boolean
- **Default:** `true` (on)
- **Read in:** `src/kernel/solvers/TimeIntegrators.jl`
  (`precompile_warmup_enabled`)
- **Precedence (highest first):**
  1. `JEXPRESSO_PRECOMPILE_WARMUP` env var
  2. `--no-precompile-warmup` CLI flag in `ARGS` (the bare
     `no-precompile-warmup` is accepted too)
  3. `:lprecompile_warmup => true/false` in `user_inputs.jl`
  4. Default (`true`)
- **Example — disable the warm-up:**
  ```bash
  JEXPRESSO_PRECOMPILE_WARMUP=0 mpirun -np 4 julia --project=. src/Jexpresso.jl CompEuler theta
  ```
- **Superseded by** `JEXPRESSO_PRECOMPILE_PASS` (below): when the
  pre-compilation pass is on, this warm-up is skipped entirely.
- **Note:** when `JEXPRESSO_ALLOC_SUMMARY=1`, the warm-up runs
  unconditionally regardless of this flag — the alloc summary needs
  the post-JIT measurement window to be meaningful, so the warm-up
  must run.

### `JEXPRESSO_PRECOMPILE_PASS`

Runs the simulation as a **pre-compilation pass** (exactly one timestep)
followed by the production loop, instead of one `solve(...)` call that JITs
inside its own hot loop.

This is the switch to use for a large problem on many ranks. Launched in one
go, a big case spends its first timestep compiling — the whole RHS chain, the
SciML integrator specialised on the real `CallbackSet`, the MPI halo-exchange
paths, the diagnostic/output path — and first-touching every per-rank working
array, all of that *interleaved with collective communication inside the
production time loop*. Ranks that compiled quickly sit in `MPI_Wait` waiting
for ones that did not; the step rate reported for the first hundreds of steps
is meaningless; and the long run inherits a heap full of compiler garbage.

`solve(prob, alg; kw...)` is by definition `solve!(init(prob, alg; kw...))`,
so the pass simply takes it apart:

```julia
integrator = init(prob, alg; kwargs...)
step!(integrator)      # PHASE 1 — one timestep: all JIT, all first touches
GC.gc(); MPI.Barrier() # compacted heap, all ranks aligned
solve!(integrator)     # PHASE 2 — the rest of the run, hot
```

Phase 2 continues the **same integrator object**. This is not a restart: no
second `init`, and the step-size controller, callback caches, FSAL history
and `u` all carry across untouched. The trajectory is bit-for-bit the one a
single `solve(...)` would have produced, and nothing is compiled twice — the
split costs exactly zero additional RHS evaluations
(`test/test_precompile_pass.jl` pins both down).

Phase 1's step *is* the simulation's first step: nothing is snapshotted,
restored or thrown away.

- **Type:** boolean
- **Default:** `false` (historical single-phase behaviour)
- **Read in:** `src/kernel/solvers/TimeIntegrators.jl`
  (`precompile_pass_enabled`)
- **Precedence (highest first):**
  1. `JEXPRESSO_PRECOMPILE_PASS` env var
  2. `--precompile-pass` / `--no-precompile-pass` CLI flag in `ARGS`
  3. `:lprecompile_pass => true/false` in `user_inputs.jl`
  4. Default (`false`)
- **Example — a 1536-rank LES run:**
  ```bash
  JEXPRESSO_PRECOMPILE_PASS=1 mpirun -np 1536 julia --project=. \
      src/Jexpresso.jl CompEuler LESICP2
  ```
  MPICH/Hydra propagates the environment to every rank by default; under
  OpenMPI pass it explicitly with `mpirun -x JEXPRESSO_PRECOMPILE_PASS ...`.
  `submit_Jexpresso_precompile_pass.sh` is a ready-made batch script that
  does this along with the serial package-precompilation phases that must
  precede any multi-rank launch.
- **Supersedes the throw-away warm-ups.** When the pass is on, both
  `precompile_warmup_run!` (drivers.jl) and the in-`time_loop!` integrator
  warm-up are skipped: they compile the same code but then restore the
  initial condition, so a run would pay for two extra full-mesh RHS
  evaluations whose results are discarded.
- **Works with adaptive stepping.** Because phase 2 resumes the same
  integrator rather than starting a new one, there is no controller to
  restart — `:ode_adaptive_solver => true` is fine.
- **Failure is not fatal.** If `init` or the first `step!` throws, the
  advanced history is rolled back and the run falls through to the historical
  single-phase `solve(...)`, JIT-ing on its first step as it always did. The
  integrator works on a copy of `prob.u0`, so the fallback always starts from
  a clean initial condition. The decision is collective — one rank failing
  sends every rank down the same path, since half a job in `solve!` and half
  in a fresh `solve` would hang rather than fall back.
- **What it does NOT do.** This is *runtime JIT*, which cannot be cached
  across processes. It is not a substitute for a serial `Pkg.precompile()`
  before a multi-rank launch (see `JEXPRESSO_PRECOMPILE_WORKLOAD` and
  `submit_Jexpresso_precompile_pass.sh`), nor for a PackageCompiler sysimage
  (`create_Jexpresso_sysimage.jl`), which is the only way to remove the JIT
  entirely.
- **Not to be confused with** `JEXPRESSO_PRECOMPILE_WARMUP` (throw-away
  warm-up inside a run, keeps nothing) or `JEXPRESSO_PRECOMPILE_WORKLOAD` (a
  `Pkg.precompile()`-time workload — a different stage altogether).

### `JEXPRESSO_STEP_HEARTBEAT`

Enables the per-step heartbeat callback that prints progress lines at
intervals during `solve(...)`:

```
 #   step 200   t = 4.000000   wall 00:03:22   1.004 s/step   ETA 3d 00:52:16
```

Useful when diagnostics are sparse (e.g. city2d's
`:diagnostics_at_times => 0:10:600` with `Δt = 0.004` means 2500
silent steps between user-visible writes — hard to tell from a hang).
Throttled: prints every step for the first 5, then every 100.

The wall-clock fields make the trace a performance meter as well as a
liveness check:

- **`wall`** — elapsed since the loop's first step. Under the
  pre-compilation pass (`JEXPRESSO_PRECOMPILE_PASS`) that first step *is*
  the pass, so `wall` times hot code only; JIT and first-touch allocation
  are already excluded. Comparable across rank counts, meshes and machines.
- **`s/step`** — averaged over the steps since the *previous* heartbeat
  (the last 100, once past step 5), not since `t0`, so a slowdown shows up
  in the line where it happens instead of being diluted by earlier history.
- **`ETA`** — that interval's simulated-seconds-per-wall-second carried out
  to `:tend`. A straight-line extrapolation of the last 100 steps: expect it
  to move around while the rate does.

Rank 0 does the timing and the printing, and the callback makes no MPI
calls, so it cannot deadlock a large job.

- **Type:** boolean
- **Default:** `false` (off)
- **Read in:** `src/kernel/solvers/TimeIntegrators.jl`
- **Precedence (highest first):**
  1. `JEXPRESSO_STEP_HEARTBEAT` env var
  2. `:lstep_heartbeat => true/false` in `user_inputs.jl`
  3. Default (`false`)
- **Example — enable the heartbeat for a debugging run:**
  ```bash
  JEXPRESSO_STEP_HEARTBEAT=1 mpirun -np 4 julia --project=. src/Jexpresso.jl CompEuler theta
  ```

### `JEXPRESSO_PRECOMPILE_WORKLOAD`

Opts the `@compile_workload` block at the bottom of
`src/Jexpresso.jl` into running a real (3-step) solve during package
precompilation, so `using Jexpresso` comes back with the per-step
kernel chain already compiled.

**Off by default, and deliberately so.** The workload goes through
`run.jl`, which calls `MPI.Init()`, and on a cluster login node that
is hostile: libfabric's `verbs`/`mlx5` provider aborts (RDMA queue
pair not permitted, locked memory too low) and `tcp` enumerates every
NIC and does reverse-DNS, which can stall for many minutes. Turning it
on can therefore convert a precompilation that would have succeeded
into a hang or a hard error. Without it the package still precompiles
fully — you only lose the warm-up, so the first solve of a given
problem shape pays its JIT at runtime.

- **Type:** boolean
- **Default:** `false` (off)
- **Read in:** `src/Jexpresso.jl`
- **Safe to enable on:** a compute node with a healthy fabric, or any
  machine where `MPI.Init` is cheap (a laptop).
- **Example:**
  ```bash
  JEXPRESSO_PRECOMPILE_WORKLOAD=1 julia --project=. -e 'using Pkg; Pkg.precompile()'
  ```
- **Side effect:** when opted in, and only if the user has not already
  pinned it, Jexpresso sets `FI_PROVIDER` for the duration of the
  precompile worker and restores it afterwards — see below.
- **Not to be confused with** `JEXPRESSO_PRECOMPILE_WARMUP`, which is
  about the warm-up *inside a run*; this one is about *precompilation*.

### `FI_PROVIDER`

libfabric's provider selection. Jexpresso **reads** this to check
whether the user has pinned a provider, and **sets** it only while the
precompile workload above is running, restoring the previous state (or
deleting the key) afterwards. If you have already set it, Jexpresso
leaves it alone.

- **Type:** libfabric provider name (`shm`, `tcp`, `verbs`, …)
- **Set by Jexpresso to:** `shm` on Linux, `sockets` on macOS —
  libfabric's `shm` is Linux-only, and requesting it on macOS makes
  `MPI_Init` fail with `OFI call getinfo failed (default nic=(n/a))`.
- **Read/written in:** `src/Jexpresso.jl`
- **Active only when:** `JEXPRESSO_PRECOMPILE_WORKLOAD=1`

---

## CI mode (`test/`)

These two only take effect when the run is in CI mode
(`parsed_CI_mode == "true"`, i.e. launched through `test/runtests.jl`
or with the CI-mode argument). They have no effect on an ordinary
`problems/` run.

### `JEXPRESSO_CI_OUTPUT`

CI mode overrides three output settings in the deck so the HDF5
comparison machinery can find what the run wrote:
`:outformat => "hdf5"`, `:output_dir => "none"`,
`:loverwrite_output => true`, plus a single output at `:tend` and
`:lwrite_initial => false`. Each override is announced on stdout as a
`# CI_MODE: forcing …` line.

Set this to a falsy value to keep the deck's own settings — the usual
reason is to get VTK out of a CI deck for visualisation.

- **Type:** boolean
- **Default:** `1` (overrides on)
- **Read in:** `src/run.jl`
- **Example:**
  ```bash
  JEXPRESSO_CI_OUTPUT=0 julia --project=. test/runtests.jl CompEuler/theta
  ```

### `JEXPRESSO_CI_OUTFORMAT`

Which format CI mode forces when the overrides above are active. The
HDF5 default is what `test/CI-ref/` comparison reads; the VTK smoke
test (`test/runtests.jl --vtk`, or `vtk_smoke = true` in
`test/ci_cases.jl`) sets this to `vtk` to exercise the writer that
production runs actually use, leaving every other CI convention
unchanged. `test/ci_compare.jl` sets and restores it around the smoke
test, so you rarely set it by hand.

- **Type:** string — `hdf5` or `vtk`
- **Default:** `hdf5`
- **Read in:** `src/run.jl`
- **Set in:** `test/ci_compare.jl` (`vtk_smoke_case`)
- **Note:** Jexpresso implements `write_vtk` for 2D and 3D only. A 1D
  case run with `vtk` raises a `MethodError`, which the smoke test
  catches and reports as "this case has no VTK writer".

---

## Internal cache prefetch (not user-set)

### `JEXPRESSO_PREFETCHED_MESH_CACHE` and `JEXPRESSO_PREFETCHED_SEM_CACHE`

These are *not* shell environment variables — they are Julia-level
`Ref`s set by `je_prefetch_caches!` during the coupling handshake to
ship pre-loaded mesh / SEM-preprocess caches into the `with_mpi`
block. Listed here so a `grep` for `JEXPRESSO_` finds them, but
nothing about them is settable from outside Julia.

---

## Julia runtime (third-party) — interaction notes

### `JULIA_PKG_PRECOMPILE_AUTO`

Used in `INSTALL.md` step 3a to disable Julia's automatic
precompilation during `Pkg.instantiate()` so the user can control the
precompile pass explicitly (since Jexpresso's dep tree is large and
parallel auto-precompile can deadlock on macOS file-cache locks).

- **Type:** boolean
- **Set to:** `0` for the install step
- **Example:**
  ```bash
  julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
  ```

### `JULIA_NUM_PRECOMPILE_TASKS`

Limits Julia's precompile parallelism. Set to `1` to serialise
precompilation when running multiple `mpirun` ranks against an
unpopulated precompile cache (4 ranks trying to write the same
`.ji` files at once can deadlock on macOS).

- **Type:** integer
- **Default:** automatic (Julia chooses based on CPU count)
- **Recommended:** `1` for parallel MPI launches if the cache may be
  stale; not needed after a successful single-process
  `Pkg.precompile()`.

---

## MPICH (third-party) — required on macOS for MPI.jl JLL

### `MPICH_INTERFACE_HOSTNAME`

Tells MPICH's TCP channel which hostname to bind to. Required on
macOS when not editing `/etc/hosts` — MPICH's `gethostbyname()` call
fails on the bare machine hostname because macOS only registers `.local`
mDNS names. See `INSTALL.md` step 3d.

- **Type:** hostname / IP string
- **Recommended value:** `127.0.0.1`
- **Used by:** any MPICH-based MPI (system MPICH, MPItrampoline_jll,
  Open MPI compiled with MPICH backend)
- **Example:**
  ```bash
  export MPICH_INTERFACE_HOSTNAME=127.0.0.1
  ```

---

## Thread-count guards set by Jexpresso itself

When element-learning is active (`:lelementLearning => true` in
`user_inputs.jl`) **and the model file is `.onnx`**, `EL_WorkBuffers`
in `src/kernel/elementLearningStructs.jl` sets the following just
before `ONNXRunTime.load_inference`:

| Variable                  | Set to | Why                                                  |
|---------------------------|--------|------------------------------------------------------|
| `OMP_NUM_THREADS`         | `1`    | Prevent oversubscription when MPI ranks share cores. |
| `ORT_NUM_THREADS`         | `1`    | Same, for the ONNXRunTime inference session.         |
| `OPENBLAS_NUM_THREADS`    | `1`    | Same, for the BLAS inside Gridap / LinearAlgebra.    |

These are written to `ENV[...]` rather than read from it; they
override whatever the user had configured, and unlike `FI_PROVIDER`
they are **not** restored afterwards. If you need a different value,
edit them in `elementLearningStructs.jl`. The `.jld2` model path does
not touch them, and neither does any non-element-learning run, so
there the user's shell settings stand.

---

## `user_inputs.jl` switches that change MPI behaviour

These are deck keys, not environment variables — there is no shell
form. They are documented here because, like the variables above, what
they change is a run's parallel cost rather than its physics.

### `:ldsgs_global_norms`

Scope of the two normalising scales the DynSGS (`:visc_model =>
DSGS()` / `DSGS_MHD()`) residual indicator divides by: the mean
`⟨q_i⟩` and the L∞ spread `‖q_i − ⟨q_i⟩‖`.

Marras eq. (9) and Nazarov & Hoffman eq. (3.5) write both over the
whole domain Ω. Under MPI that is a collective on the critical path of
every rank on **every RHS call** — five times per step under
`CarpenterKennedy2N54`, two or three reductions each, so 10–15
`Allreduce` per time step.

The default is rank-local and communicates nothing. These two
quantities only set the *scale* the element residual is measured
against; what the model needs from them is the order of magnitude of
the solution's variation, and a partition of a connected domain
resolves that as well as the whole domain does. μ is bounded by
`min(μ_res, μ_max)` either way, so the flow solution differs only at
the level of the usual round-off divergence.

- **Type:** boolean
- **Default:** `false` (rank-local, no communication)
- **Read in:** `src/io/mod_inputs.jl` (default),
  `src/kernel/operators/rhs.jl` (threaded to all three call sites),
  `src/kernel/physics/SGS.jl` (`_dsgs_norm_scope` block)
- **Applies to:** all four DynSGS implementations — 1D, 2D Euler-θ,
  2D total-energy, 2D GLM-MHD
- **No effect on:** serial runs, and any run whose `:visc_model` is
  not a DynSGS variant
- **Set it `true` when:**
  - μ has to be reproducible across rank counts — e.g. a regression
    test comparing a 1-rank and an N-rank run field-by-field;
  - a rank's subdomain genuinely cannot see the solution's scale (a
    partition lying entirely inside a uniform region while the
    interesting structure lives on another rank).
- **Example:**
  ```julia
  :visc_model         => DSGS(),
  :ldsgs_global_norms => true,   # paper's domain norms; costs 2-3 Allreduce/RHS
  ```
- **History:** the 2D total-energy and MHD implementations used to do
  these reductions unconditionally, and the 1D and Euler-θ ones never
  did. All four now share one switch, defaulting to rank-local. The
  2D total-energy path also reduced over a hardcoded
  `MPI.COMM_WORLD`, which would have deadlocked under Alya MPMD
  coupling (COMM_WORLD carries Alya's ranks, which never enter
  DynSGS); the global mode now uses `get_mpi_comm()` like everything
  else.

---

## Quick reference

### Environment variables

| Variable                          | Type   | Default     | Purpose                                          |
|-----------------------------------|--------|-------------|--------------------------------------------------|
| `JEXPRESSO_COUPLED`               | bool   | unset       | Enable Alya MPMD coupling                        |
| `APPID`                           | int    | `2`         | Coupling app identifier (only with COUPLED=1)    |
| `JEXPRESSO_ALLOC_SUMMARY`         | bool   | `false`     | End-of-run timing/allocation table               |
| `JEXPRESSO_PRECOMPILE_WARMUP`     | bool   | `true`      | One-step JIT warm-up before real solve           |
| `JEXPRESSO_PRECOMPILE_WORKLOAD`   | bool   | `false`     | Run a 3-step solve during package precompilation |
| `JEXPRESSO_STEP_HEARTBEAT`        | bool   | `false`     | Step/wall-clock/s-per-step prints during solve   |
| `JEXPRESSO_CI_OUTPUT`             | bool   | `1`         | CI mode forces hdf5/none/overwrite output        |
| `JEXPRESSO_CI_OUTFORMAT`          | string | `hdf5`      | Which format CI mode forces (`hdf5` \| `vtk`)    |
| `FI_PROVIDER`                     | string | (unset)     | libfabric provider; set only for the precompile workload |
| `JULIA_PKG_PRECOMPILE_AUTO`       | bool   | `1`         | Auto-precompile in `Pkg.instantiate()`           |
| `JULIA_NUM_PRECOMPILE_TASKS`      | int    | auto        | Parallelism cap for Julia's precompile pass      |
| `MPICH_INTERFACE_HOSTNAME`        | string | (system)    | macOS hostname workaround for MPICH/JLL binary   |
| `ANTHROPIC_API_KEY`               | string | (unset)     | `tools/EquationGenerator.jl` only; not read by the solver |

### `user_inputs.jl` keys documented here

| Key                    | Type | Default | Purpose                                   |
|------------------------|------|---------|-------------------------------------------|
| `:ldsgs_global_norms`  | bool | `false` | DynSGS norms: rank-local vs domain-global |
