# Profiling Jexpresso with Extrae.jl

This directory contains a **working, self-contained example** of how to
instrument Jexpresso-style code with
[**Extrae.jl**](https://github.com/bsc-quantic/Extrae.jl) — the Julia
bindings to the [Extrae](https://tools.bsc.es/extrae) HPC profiler from the
Barcelona Supercomputing Center.

Extrae produces [Paraver](https://tools.bsc.es/paraver) traces, which let you
see, per MPI rank and over time, where your program spends its time
(computation vs. `MPI_Waitall` vs. `MPI_Allreduce`, etc.). This is exactly the
analysis carried out in the reference paper:

> S. Sánchez-Ramírez and M. Giordano, *“Extrae.jl: Julia bindings for the
> Extrae HPC Profiler”*, Proceedings of JuliaCon (arXiv:2504.12087v1, 2025).

---

## ⚠️ Important: macOS vs. Linux

The Extrae native library (`Extrae_jll`) is built **only for Linux**
(`x86_64`, `aarch64`, `powerpc64le`, glibc). **There is no macOS build.**

So on your MacBook Air you **cannot generate a real Paraver trace**. To keep
the example usable everywhere, all the instrumentation goes through a small
shim (`ExtraeShim.jl`) that:

- uses the **real Extrae** when it is installed and the native library
  initialises (Linux / HPC), and
- **degrades to no-ops** otherwise (macOS / Windows), so the *exact same
  instrumented code* still runs to completion — it just produces no trace.

This means you can **develop and run the example on your Mac**, and then run
the **byte-for-byte identical code on a Linux cluster** to get the real
`.prv` trace.

---

## Files

| File | Purpose |
|------|---------|
| `ExtraeShim.jl` | Portable wrapper around Extrae.jl (real on Linux, no-op elsewhere). |
| `extrae_axpy.jl` | Serial example — reproduces Listings 1 & 2 of the paper (`@user_function` + custom events). |
| `extrae_mpi_jexpresso_pattern.jl` | MPI example — a halo-exchange + `Allreduce` stencil solver mirroring `src/kernel/mpi/mpi_communications.jl`. |
| `runtests.jl` | The test to run on your MacBook Air. |
| `run_extrae_example.sh` | Convenience launcher for the MPI example. |
| `extrae.xml` | Extrae runtime configuration (used on Linux only). |

---

## Quick start (works on your MacBook Air)

From the **Jexpresso project root**:

```bash
# 1. Run the whole test (serial + MPI on 4 ranks — one per core):
julia --project=. tools/Extrae/runtests.jl

# 2. Or run the pieces individually:
julia --project=. tools/Extrae/extrae_axpy.jl

./tools/Extrae/run_extrae_example.sh 4        # MPI example, 4 ranks
```

On macOS you will see messages such as:

```
Ran in no-op shim mode (Extrae native library not available on this platform).
```

That is expected — it confirms the instrumentation is wired correctly even
though no trace is produced. The test **passes** on macOS.

> The MPI example uses MPI.jl's bundled `mpiexec`, so it honours whatever your
> `MPIPreferences` point at. On macOS the recommended route is **MPICH_jll**
> (see the repo's `INSTALL.md`, "Route C"). 4 ranks fit one-per-core on a
> 4-core MacBook Air; override the rank count with
> `JEXPRESSO_EXTRAE_NRANKS=N julia --project=. tools/Extrae/runtests.jl`.

---

## Getting a REAL Paraver trace (on Linux / an HPC cluster)

1. **Install Extrae.jl** into the project (Linux only):

   ```bash
   julia --project=. -e 'import Pkg; Pkg.add("Extrae")'
   ```

2. **Locate the Extrae MPI tracing library.** Use the **C** library
   `libmpitrace.so` — *not* the Fortran `libmpitracef.so`, since MPI.jl calls
   the C MPI ABI. Also capture the dependency search path (`libunwind`, PAPI,
   …) so the dynamic loader can start `julia` with the preload attached:

   ```bash
   # (A) system module — preferred on HPC (deps resolve via system libs):
   module load extrae
   export EXTRAE_LIB=$EXTRAE_HOME/lib/libmpitrace.so
   export EXTRAE_LIBPATH=$EXTRAE_HOME/lib

   # (B) or the Julia artifact (note the exact C lib name, and use the jll's
   #     own LIBPATH so libunwind.so.8 etc. are found):
   ART=$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.artifact_dir)')
   export EXTRAE_LIB=$ART/lib/libmpitrace.so
   export EXTRAE_LIBPATH=$(julia --project=. -e 'using Extrae_jll; print(Extrae_jll.LIBPATH[])')
   ```

   > If you skip `EXTRAE_LIBPATH` with the artifact route you'll see
   > `julia: error while loading shared libraries: libunwind.so.8: cannot open
   > shared object file` — that just means the preload's dependencies aren't on
   > the loader path yet.

3. **Let Extrae instrument MPI automatically.** Extrae must be loaded *after*
   Julia but *before* the MPI library. The paper (§3.1) suggests MPI.jl's
   `preloads` preference — but that keyword only exists in newer
   `MPIPreferences`, and Jexpresso pins `MPIPreferences = "=0.1.11"`, which
   does **not** support it. The portable equivalent is to attach the preload
   to the **rank processes only**, via `env`, so it loads after Julia but
   before MPI (do *not* `LD_PRELOAD` the `julia` binary itself — that loads
   Extrae before Julia and triggers libstdc++ version clashes). The launcher
   does this for you when `EXTRAE_LIB` is set:

   ```bash
   export EXTRAE_CONFIG_FILE=$PWD/tools/Extrae/extrae.xml
   ./tools/Extrae/run_extrae_example.sh 32          # 32 ranks on a 32-core node
   ```

   Equivalent explicit command:

   ```bash
   julia --project=. -e '
     using MPI
     run(`$(mpiexec()) -n 32 env LD_PRELOAD='"$EXTRAE_LIB"' \
         $(Base.julia_cmd()) --project=. tools/Extrae/extrae_mpi_jexpresso_pattern.jl`)'
   ```

4. **Merge** the per-process trace files into a single Paraver trace (this is
   done automatically by the `<merge>` block in `extrae.xml`; if you disable
   it, merge manually):

   ```bash
   mpi2prv -f TRACE.mpits -o jexpresso-extrae.prv
   ```

5. **Open `jexpresso-extrae.prv` in Paraver** and inspect the timeline. As in
   the paper, you will be able to see the compute (`compute RHS`) regions, the
   `halo exchange` communication, and the `allreduce` collective, and quantify
   how much time each MPI rank spends in each.

### Memory / OOM on a single node

Each rank is a full Julia process (base runtime + MPI + Extrae ≈ 0.5–1 GB)
**plus** Extrae's per-process trace buffer. Running many ranks on one node can
exceed the memory granted to an interactive allocation (which often defaults
to only a few GB total, not the whole node), giving an `oom_kill` /
`Out Of Memory` / `EXIT CODE: 9` failure.

If you hit OOM:

- **Give the allocation the node's memory.** In your interactive session
  request enough RAM, e.g. `salloc --nodes=1 --ntasks=32 --mem=0` (`--mem=0`
  asks for all memory on the node) or `--mem-per-cpu=2G`.
- **Use fewer ranks** while validating, then scale up:
  `./tools/Extrae/run_extrae_example.sh 8`.
- The trace buffer in `extrae.xml` is already set small (`<buffer><size>` =
  100000 events ≈ a few MB/rank); keep it modest for many-rank single-node
  runs.

A practical rule of thumb: budget ~1 GB per rank, so 32 ranks ⇒ ask for
~32 GB.

---

## Submitting as a SLURM batch job

Instead of hand-setting an interactive allocation every time, use the ready
made job script, which requests the right `--ntasks`/`--mem`, precompiles the
project once, sets up the Extrae preload from `Extrae_jll`, and launches the
example:

```bash
sbatch tools/Extrae/submit_extrae.sh
```

Edit the `#SBATCH` header to change the rank count (`--ntasks`), memory
(`--mem-per-cpu`), and your `--account`. The trace files appear in the submit
directory; watch progress in `jexp-extrae.<jobid>.out` / `.err`.

---

## Visualizing the trace in Paraver

A run produces three files that together form one Paraver trace — keep them
together (same basename):

```
jexpresso-extrae.prv     # the trace data
jexpresso-extrae.pcf     # event/value labels (so you see "compute RHS", etc.)
jexpresso-extrae.row     # the rank/thread layout
```

1. **Get Paraver (`wxparaver`).** It is BSC's free GUI viewer. Either load a
   cluster module if available (`module avail paraver` → `module load paraver`)
   or download it for your laptop from <https://tools.bsc.es/paraver> and copy
   the three files over:

   ```bash
   scp 'wulver:/project/smarras/smarras/Jexpresso/jexpresso-extrae.*' .
   ```

2. **Open the trace:**

   ```bash
   wxparaver jexpresso-extrae.prv
   ```

   The `.pcf`/`.row` are picked up automatically because they share the
   basename.

3. **Look at the timeline.** Paraver ships ready-made *configuration files*
   (`.cfg`) that build the standard views in one click:
   - `File ▸ Load Configuration…` → from the bundled `cfgs/` pick
     **`mpi/views/MPI_call.cfg`** to colour each rank's timeline by MPI routine
     (you'll see `MPI_Waitall` / `MPI_Allreduce` blocks, like Fig. 2 in the
     paper), and **`General/views/user_functions.cfg`** to see the
     `compute RHS` / `halo exchange` / `allreduce` user regions you annotated.
   - The custom events you `emit` (Solver phase, Iteration, residual) show up
     under `<event-type>` — right-click a timeline ▸ *View ▸ Event Flags*, or
     open a 2D *Analyzer* histogram on the "Solver phase" event type.

4. **Quantify with the Analyzer (2D tables).** `File ▸ Load Configuration…` →
   **`mpi/analysis/2dp_MPI_call_profile.cfg`** gives the % of time each rank
   spends in each MPI call — the table behind the paper's conclusion that
   `MPI_Waitall`/`MPI_Allreduce` dominate.

Tip: if you only have a terminal, `prv2dim` and the `paramedir` CLI (shipped
with Paraver) can compute the same profiles headless and dump them to CSV.

---

## How the instrumentation maps to Jexpresso

The MPI example is intentionally shaped like Jexpresso's real solver loop:

| Example phase | Extrae annotation | Jexpresso analogue |
|---------------|-------------------|--------------------|
| Local stencil update | `@user_function` region "compute RHS" | RHS / element-matrix kernels |
| `Isend`/`Irecv!`/`Waitall` | `@user_function` region "halo exchange" | `src/kernel/mpi/mpi_communications.jl` assembler exchange |
| `Allreduce` | `@user_function` region "allreduce" | global residual / norm reductions |
| Iteration counter, residual | `emit(...)` events | time-step / convergence diagnostics |

To instrument the **real** Jexpresso solver, follow the same recipe: `include`
`ExtraeShim.jl`, call `ExtraeShim.init()` once after `MPI.Init()`, wrap the hot
regions of the time loop in `@user_function`, `emit` the step number, and call
`ExtraeShim.finish()` at the end. Because the shim is a no-op off Linux, this
can live in the codebase without breaking laptop/macOS runs.
