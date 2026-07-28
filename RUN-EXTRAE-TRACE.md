# Tracing Jexpresso with Extrae and visualizing it in Paraver

This guide walks through, end to end:

1. **running a test** that exercises the Extrae instrumentation,
2. **extracting an MPI trace** of a real Jexpresso case, and
3. **visualizing** that trace in Paraver.

Jexpresso uses [**Extrae.jl**](https://github.com/bsc-quantic/Extrae.jl), the
Julia bindings to the [Extrae](https://tools.bsc.es/extrae) HPC profiler from
the Barcelona Supercomputing Center. Extrae intercepts every MPI call and
writes a [Paraver](https://tools.bsc.es/paraver) trace, so you can see — per
rank, over time — how much of the run is real computation and how much is
`MPI_Waitall` / `MPI_Allreduce` / halo exchange.

> Reference: S. Sánchez-Ramírez and M. Giordano, *"Extrae.jl: Julia bindings
> for the Extrae HPC Profiler"*, Proceedings of JuliaCon
> (arXiv:2504.12087v1, 2025).

For the design of the in-tree instrumentation and the standalone examples, see
[`tools/Extrae/README.md`](tools/Extrae/README.md). This document is the
task-oriented "how do I get a trace" companion to it.

---

## ⚠️ Prerequisite: Linux

The Extrae native library (`Extrae_jll`) is built **for Linux only**
(`x86_64`, `aarch64`, `powerpc64le`, glibc). **There is no macOS build**, so
you cannot produce a `.prv` trace on a Mac.

That is not a blocker for development: all instrumentation is written so that
it **degrades to no-ops** when Extrae is unavailable. The same instrumented
code runs to completion on macOS — it just produces no trace. So you can
develop and test on a laptop, then run the identical code on a Linux cluster to
capture the trace.

Two things are always off unless you ask for them:

| Switch | Effect |
|--------|--------|
| `JEXPRESSO_EXTRAE=1` | Turns the **in-solver** instrumentation on. Unset (default) ⇒ every instrumentation call is a single `Ref{Bool}` read and Extrae is never even loaded. |
| `LD_PRELOAD=…/libmpitrace.so` | Turns **automatic MPI interception** on. Without it you get the custom Jexpresso regions but no MPI calls. |

The launcher scripts below set both for you.

---

## Step 1 — Run the test (works anywhere, including macOS)

This is the portable smoke test. It verifies the shim API never throws, and
that both the serial and the MPI example run to completion:

```bash
# from the Jexpresso project root
julia --project=. tools/Extrae/runtests.jl
```

It runs three test sets: the shim API, the serial `axpy` example, and the MPI
example on 4 ranks. Override the rank count with:

```bash
JEXPRESSO_EXTRAE_NRANKS=8 julia --project=. tools/Extrae/runtests.jl
```

On macOS you will see

```
Ran in no-op shim mode (Extrae native library not available on this platform).
```

which is expected — the test still **passes**, confirming the instrumented code
paths are correct. On Linux with Extrae installed, the same runs additionally
produce a real trace.

The pieces can also be run individually:

```bash
julia --project=. tools/Extrae/extrae_axpy.jl     # serial example
./tools/Extrae/run_extrae_example.sh 4            # MPI example, 4 ranks
```

---

## Step 2 — Install Extrae (Linux, once per project)

Extrae is deliberately **not** a `Project.toml` dependency — it is a
Linux-only binary and the instrumentation loads it lazily at runtime, so the
project stays installable on macOS. Add it explicitly:

```bash
julia --project=. -e 'import Pkg; Pkg.add("Extrae"); Pkg.precompile()'
```

Do this **serially, before any MPI launch**. If several ranks race to
precompile Extrae under an active MPI session, the run can hang. (The
launchers in Step 3 do this for you.)

---

## Step 3 — Extract the trace

### 3a. The easy path: `run_jexpresso_traced.sh`

One command traces a real case. It derives the Extrae paths from `Extrae_jll`,
precompiles Extrae serially, sets `JEXPRESSO_EXTRAE=1`, and attaches the
preload to the rank processes:

```bash
./tools/Extrae/run_jexpresso_traced.sh 4 CompEuler theta
#                                      │ │         └ case dir under problems/CompEuler/
#                                      │ └ equations dir under problems/
#                                      └ number of MPI ranks
```

Defaults are `4 CompEuler theta`. Another example:

```bash
./tools/Extrae/run_jexpresso_traced.sh 16 CompEuler 3d
```

On success rank 0 prints

```
Jexpresso: Extrae tracing ACTIVE (JEXPRESSO_EXTRAE set; already_initialised=true)
```

and the script lists the trace files at the end.

### 3b. Overriding the Extrae installation

Export any of these before calling the launcher and it will use yours instead
of the `Extrae_jll` artifact — e.g. to use a system `module load extrae`:

```bash
module load extrae
export EXTRAE_LIB=$EXTRAE_HOME/lib/libmpitrace.so    # the C lib, see warning below
export EXTRAE_LIBPATH=$EXTRAE_HOME/lib               # so libunwind/PAPI resolve
export EXTRAE_CONFIG_FILE=$PWD/tools/Extrae/extrae.xml
./tools/Extrae/run_jexpresso_traced.sh 16 CompEuler theta
```

> **Use `libmpitrace.so`, not `libmpitracef.so`.** MPI.jl calls the **C** MPI
> ABI; the `f` variant is the Fortran interception library and will not trace
> MPI.jl's calls. The launchers warn if you point at the Fortran one.

### 3c. Launching by hand

The preload must be attached to the **rank processes**, not to the `julia`
binary — Extrae has to load *after* Julia but *before* the MPI library.
Preloading onto `julia` itself causes libstdc++ version clashes. Build the
`env` prefix in **bash**, not inside the Julia `-e` string (Julia would try to
interpolate `$LD_LIBRARY_PATH` and fail with `UndefVarError`):

```bash
export JEXPRESSO_EXTRAE=1
export EXTRAE_CONFIG_FILE=$PWD/tools/Extrae/extrae.xml
ART=$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.artifact_dir)')
export EXTRAE_LIB=$ART/lib/libmpitrace.so
export EXTRAE_LIBPATH=$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.LIBPATH[])')

PREFIX="env -u OMP_NUM_THREADS LD_LIBRARY_PATH=$EXTRAE_LIBPATH:${LD_LIBRARY_PATH:-} LD_PRELOAD=$EXTRAE_LIB"

julia --project=. -e "
  using MPI
  run(\`\$(mpiexec()) -n 4 $PREFIX \$(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler theta\`)"
```

### 3d. As a SLURM batch job

```bash
sbatch tools/Extrae/submit_extrae.sh
```

The script requests one CPU and ~4 GB per rank, precompiles once, sets up the
preload from `Extrae_jll`, and launches. **Edit the `#SBATCH` header** for your
site: `--account`, `--partition`, `--qos`, `--ntasks` (= rank count), and
`--mem-per-cpu`. As shipped it profiles the standalone MPI *example*
(`extrae_mpi_jexpresso_pattern.jl`) on 32 ranks; to profile a real case
instead, replace the `run_extrae_example.sh` line at the bottom with a
`run_jexpresso_traced.sh <NRANKS> <EQ> <CASE>` call. Progress goes to
`jexp-extrae.<jobid>.out` / `.err`.

---

## Step 4 — Check what you got

A successful run leaves **three files** in the project root that together form
one Paraver trace. Keep them together — they are matched by basename:

```
jexpresso-extrae.prv     # the trace data
jexpresso-extrae.pcf     # event/value labels (so you see "time_loop", not "4")
jexpresso-extrae.row     # the rank/thread layout
```

The merge is done automatically by the `<merge>` block in
[`tools/Extrae/extrae.xml`](tools/Extrae/extrae.xml). If you disable it, merge
the per-process `.mpit` files by hand:

```bash
mpi2prv -f TRACE.mpits -o jexpresso-extrae.prv
```

### What is in the trace

Besides every intercepted MPI call, the trace carries one custom event type,
**`Jexpresso phase`** (numeric type `6700000`), whose value says which part of
the solver is executing. Currently emitted:

| Phase label | What it covers | Instrumented in |
|-------------|----------------|-----------------|
| `time_loop` | the whole main workload | `problems/drivers.jl` |
| `rhs` | every RHS evaluation (each RK stage) | `src/kernel/operators/rhs.jl` |
| `rhs_comm` | the RHS inter-rank assembly (`DSS_global_RHS!` → `assemble_mpi!`) *inside* each `rhs` region | `src/kernel/operators/rhs.jl` |
| `coupling_setup` | one-time Jexpresso↔Alya handshake / data receive | `src/run.jl` |
| `coupling_interp` | per-step interpolation of the solution to Alya points | `src/kernel/coupling/couplingStructs.jl` |
| `coupling_comm` | per-step MPI send of the packed data to Alya | `src/kernel/coupling/couplingStructs.jl` |

Reading **compute vs. communication inside the RHS**: within an `rhs` region
the phase event reads `rhs` during the volume/flux compute and switches to
`rhs_comm` for the MPI assembly, then back. The assembly is marked with
`Profiling.mark`, which only changes the event value — it does **not** open a
nested user-function region, so region nesting stays flat.

The `coupling_*` phases only appear in a **coupled** run (see Step 6). In a
standalone run they are simply absent.

> The labels `idle`, `sem_setup`, `initialize`, `params_setup` and
> `halo_exchange` are already registered in the `.pcf` (so Paraver will show
> names for them) but **no call site emits them yet** — the instrumentation is
> being layered on incrementally. Do not expect them in a timeline today.

---

## Step 5 — Visualize in Paraver

### 5a. Get Paraver

`wxparaver` is BSC's free GUI viewer. Either load a cluster module:

```bash
module avail paraver && module load paraver
```

…or install it on your laptop from <https://tools.bsc.es/paraver> and copy the
three files across (note the glob — you need all three):

```bash
scp 'mycluster:/path/to/Jexpresso/jexpresso-extrae.*' .
```

### 5b. Open the trace

```bash
wxparaver jexpresso-extrae.prv
```

The `.pcf` and `.row` are picked up automatically from the shared basename. If
your timeline shows bare integers instead of names, the `.pcf` is missing.

### 5c. The views worth loading first

Paraver ships ready-made **configuration files** (`.cfg`) that build the
standard views in one click — `File ▸ Load Configuration…`, then from the
bundled `cfgs/` directory:

| Configuration | What it shows |
|---------------|---------------|
| `mpi/views/MPI_call.cfg` | each rank's timeline coloured by MPI routine — the `MPI_Waitall` / `MPI_Allreduce` blocks |
| `General/views/user_functions.cfg` | the Jexpresso user regions (`time_loop`, `rhs`, `coupling_*`) |
| `mpi/analysis/2dp_MPI_call_profile.cfg` | a 2D table: % of time each rank spends in each MPI call |

To see the `Jexpresso phase` event itself, right-click a timeline ▸
*View ▸ Event Flags*, or open a 2D **Analyzer** histogram on the
`Jexpresso phase` event type — that is what gives you per-phase time
per rank, including the `rhs` vs `rhs_comm` split.

### 5d. Headless / CLI

If you only have a terminal, `prv2dim` and `paramedir` (both shipped with
Paraver) compute the same profiles without the GUI and dump them to CSV.

---

## Step 6 — Tracing a coupled (Jexpresso ↔ Alya) run

`run_jexpresso_traced.sh` launches a **single-program** job, so it cannot start
a coupled MPMD run. Set the coupling up as described in
[RUN-COUPLED.md](RUN-COUPLED.md) — in particular, **both codes must use the
same MPI library** — and then add the Extrae `env` prefix to the *Jexpresso
side* of the MPMD command line:

```bash
export JEXPRESSO_COUPLED=1
export JEXPRESSO_EXTRAE=1
export EXTRAE_CONFIG_FILE=$PWD/tools/Extrae/extrae.xml
ART=$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.artifact_dir)')
export EXTRAE_LIB=$ART/lib/libmpitrace.so
export EXTRAE_LIBPATH=$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.LIBPATH[])')

mpirun -np 2 ./AlyaProxy/Alya.x \
     : -np 2 env LD_LIBRARY_PATH="$EXTRAE_LIBPATH:${LD_LIBRARY_PATH:-}" \
                 LD_PRELOAD="$EXTRAE_LIB" \
       julia --project=. ./src/Jexpresso.jl CompEuler 3dAlya
```

In the timeline you then see, for every coupled step, a `coupling_interp` block
immediately followed by a `coupling_comm` block, which is what lets you
quantify interpolation vs. communication cost of the coupling.

Two caveats specific to the coupled case:

- **Only the preloaded ranks are traced.** The prefix above attaches Extrae to
  the Jexpresso ranks only, so the Alya ranks contribute no events. The
  Jexpresso-side MPI calls that talk to Alya are still recorded, so the coupling
  cost is visible from Jexpresso's point of view; you just do not see what Alya
  was doing on the other side.
- **This MPMD + Extrae combination is less well exercised** than the
  single-program path. Verify you actually get a `.prv` on a short case before
  spending a large allocation on it.

---

## Troubleshooting

| Symptom | Cause and fix |
|---------|---------------|
| `julia: error while loading shared libraries: libunwind.so.8: cannot open shared object file` | `EXTRAE_LIBPATH` not set, so the preload's own dependencies aren't on the loader path. Set it (Step 3b/3c). |
| MPI calls missing from the trace, only user regions present | `EXTRAE_LIB` points at `libmpitracef.so` (Fortran). Use the C library `libmpitrace.so`. |
| Run **hangs right after the Extrae banner** | Usually ranks racing to precompile Extrae under an active MPI session. Precompile serially first (Step 2). The double-`Extrae_init` deadlock is already handled — `Profiling.init` skips initialisation when the preload has already done it. |
| `UndefVarError` mentioning `LD_LIBRARY_PATH` | You built the `env …` prefix *inside* the Julia `-e` string, so Julia tried to interpolate it. Build it in bash (Step 3c). |
| `oom_kill` / `Out Of Memory` / `EXIT CODE: 9` | Each rank is a full Julia process plus an Extrae buffer — budget **~1 GB per rank**. Request the memory (`salloc --mem=0`, or `--mem-per-cpu=2G`) or use fewer ranks while validating. The buffer in `extrae.xml` is already small (100 000 events); a smaller buffer never loses data, it just flushes more often. |
| No `.prv` produced at all | Check the `<merge>` block in `extrae.xml` is `enabled="yes"`, or merge manually with `mpi2prv -f TRACE.mpits -o jexpresso-extrae.prv`. |
| Timeline shows numbers instead of phase names | The `.pcf` is missing or has a different basename than the `.prv`. Copy all three files together. |
| No `Jexpresso: Extrae tracing ACTIVE` message | `JEXPRESSO_EXTRAE` is unset, or Extrae isn't installed/loadable. Set the variable and run Step 2. |
| Need to see *where* startup stalls | `export JEXPRESSO_EXTRAE_DEBUG=1` for per-rank, flushed diagnostics through each Extrae start-up step. |
| Extrae warns `OMP_NUM_THREADS is set but OpenMP is not supported!` | Harmless. The launchers already `env -u OMP_NUM_THREADS` per rank; Julia threading uses `JULIA_NUM_THREADS`. |

---

## Related documentation

- [`tools/Extrae/README.md`](tools/Extrae/README.md) — instrumentation design, the standalone examples, how the phases map onto Jexpresso
- [`src/kernel/infrastructure/Profiling.jl`](src/kernel/infrastructure/Profiling.jl) — the opt-in `Profiling` submodule (phase codes, `region`/`mark` API)
- [`tools/Extrae/extrae.xml`](tools/Extrae/extrae.xml) — Extrae runtime configuration (buffer size, merge, counters)
- [RUN-COUPLED.md](RUN-COUPLED.md) — setting up and launching a coupled Jexpresso ↔ Alya run
- [ENVIRONMENT_VARIABLES.md](ENVIRONMENT_VARIABLES.md) — all environment variables Jexpresso reads
- [INSTALL.md](INSTALL.md) — MPI setup (§5, "Running in parallel with MPI")
- Extrae user guide: <https://tools.bsc.es/doc/html/extrae/> · Paraver: <https://tools.bsc.es/paraver>
