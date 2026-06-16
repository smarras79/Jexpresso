# Environment Variables

This file lists every environment variable Jexpresso reads, plus the
third-party ones (Julia, MPI, BLAS) that affect Jexpresso runs in
ways worth documenting.

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
table (`TimerOutputs`). On by default in CI, off by default in normal
runs because it adds one full extra RHS warm-up call to make the
post-JIT measurement window meaningful.

- **Type:** boolean
- **Default:** `false` (off)
- **Read in:** `src/kernel/solvers/TimeIntegrators.jl`,
  `src/auxiliary/timing.jl`
- **Precedence (highest first):**
  1. `JEXPRESSO_ALLOC_SUMMARY` env var
  2. `--no-alloc-summary` CLI flag in `ARGS`
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
- **Precedence (highest first):**
  1. `JEXPRESSO_PRECOMPILE_WARMUP` env var
  2. `--no-precompile-warmup` CLI flag in `ARGS`
  3. `:lprecompile_warmup => true/false` in `user_inputs.jl`
  4. Default (`true`)
- **Example — disable the warm-up:**
  ```bash
  JEXPRESSO_PRECOMPILE_WARMUP=0 mpirun -np 4 julia --project=. src/Jexpresso.jl CompEuler theta
  ```
- **Note:** when `JEXPRESSO_ALLOC_SUMMARY=1`, the warm-up runs
  unconditionally regardless of this flag — the alloc summary needs
  the post-JIT measurement window to be meaningful, so the warm-up
  must run.

### `JEXPRESSO_STEP_HEARTBEAT`

Enables the per-step heartbeat callback that prints
`#   step N   t = X.X` lines at intervals during `solve(...)`.
Useful when diagnostics are sparse (e.g. city2d's
`:diagnostics_at_times => 0:10:600` with `Δt = 0.004` means 2500
silent steps between user-visible writes — hard to tell from a hang).
Throttled: prints every step for the first 5, then every 100.

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
`user_inputs.jl`), Jexpresso programmatically sets the following at
load time inside `src/kernel/elementLearningStructs.jl`:

| Variable                  | Set to | Why                                                  |
|---------------------------|--------|------------------------------------------------------|
| `OMP_NUM_THREADS`         | `1`    | Prevent oversubscription when MPI ranks share cores. |
| `ORT_NUM_THREADS`         | `1`    | Same, for the ONNXRunTime inference session.         |
| `OPENBLAS_NUM_THREADS`    | `1`    | Same, for the BLAS inside Gridap / LinearAlgebra.    |

These are written to `ENV[...]` rather than read from it; they
override whatever the user had configured. If you need a different
value, edit the constants in `elementLearningStructs.jl`. For
non-element-learning runs Jexpresso does not touch these variables,
so the user's shell settings stand.

---

## Quick reference

| Variable                       | Type   | Default     | Purpose                                          |
|--------------------------------|--------|-------------|--------------------------------------------------|
| `JEXPRESSO_COUPLED`            | bool   | unset       | Enable Alya MPMD coupling                        |
| `APPID`                        | int    | `2`         | Coupling app identifier (only with COUPLED=1)    |
| `JEXPRESSO_ALLOC_SUMMARY`      | bool   | `false`     | End-of-run timing/allocation table               |
| `JEXPRESSO_PRECOMPILE_WARMUP`  | bool   | `true`      | One-step JIT warm-up before real solve           |
| `JEXPRESSO_STEP_HEARTBEAT`     | bool   | `false`     | Per-step progress prints during solve            |
| `JULIA_PKG_PRECOMPILE_AUTO`    | bool   | `1`         | Auto-precompile in `Pkg.instantiate()`           |
| `JULIA_NUM_PRECOMPILE_TASKS`   | int    | auto        | Parallelism cap for Julia's precompile pass      |
| `MPICH_INTERFACE_HOSTNAME`     | string | (system)    | macOS hostname workaround for MPICH/JLL binary   |
