# Jexpresso FAQ

Frequently asked questions and common errors, grouped by topic. For full setup
details see [INSTALL.md](INSTALL.md); for coupled Jexpresso+Alya runs see
[RUN-COUPLED.md](RUN-COUPLED.md).

## Table of contents

- [Run](#run)
- [Installation](#installation)

---

## Run

### `prterun` cannot find the `julia` executable

**Q.** If you run the code with command

```bash
mpirun -n 4 julia --project=. src/Jexpresso.jl CompEuler theta
```

and receive the error

```
prterun was unable to find the specified executable file, and therefore did
not launch the job.  This error was first reported for process rank
2; it may have occurred for other processes as well.

NOTE: A common cause for this error is misspelling a prterun command
   line parameter option (remember that prterun interprets the first
   unrecognized command line token as the executable).

Executable: julia
```

**A.** You simply need to replace `julia` in your command line call with its
full path; for example:

```bash
mpirun -n 4 /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia --project=. src/Jexpresso.jl CompEuler theta
```

This usually happens with OpenMPI's `prterun` launcher when `julia` is provided
through a shim (e.g. `juliaup`) that the launcher cannot resolve on the worker
processes. To find the real binary automatically:

```bash
JULIA_BIN=$(julia -e 'print(joinpath(Sys.BINDIR, "julia"))')
echo "$JULIA_BIN"   # e.g. ~/.julia/juliaup/julia-1.11.x+.../bin/julia
mpirun -n 4 "$JULIA_BIN" --project=. src/Jexpresso.jl CompEuler theta
```

### A parallel run hangs forever at `MPI.Init` / first collective

**A.** The most common cause is that Julia's `MPI.jl` is bound to a different MPI
than the launcher (`mpiexec`/`mpirun`) you used — or, in a coupled run, a
different MPI than Alya was compiled with. They must be the same implementation,
version, and ABI. Check what `MPI.jl` is bound to:

```bash
julia --project=. -e '
  using MPIPreferences; println("binary = ", MPIPreferences.binary)
  using MPI;            println(MPI.identify_implementation())'
```

Then make sure you launch with the matching `mpiexec`. See
[INSTALL.md, Section 5](INSTALL.md#5-running-in-parallel-with-mpi) for standalone
runs and [RUN-COUPLED.md](RUN-COUPLED.md) for coupled runs.

### `MPI_Init` fails with `gethostbyname failed` on macOS

**A.** This affects any MPICH-based MPI on macOS (including the bundled
`MPItrampoline_jll`). Register your hostname in `/etc/hosts`:

```bash
echo "127.0.0.1   $(hostname -s)" | sudo tee -a /etc/hosts
echo "127.0.0.1   $(hostname)"    | sudo tee -a /etc/hosts
```

Full explanation:
[INSTALL.md, Section 5.5](INSTALL.md#55-macos-hostname-fix-mpich-and-mpich_jll-only).

### Precompiling on an InfiniBand cluster fails with `Failed to modify UD QP to INIT on mlx5_0`

**A.** When `using Jexpresso` precompiles, its `@compile_workload` runs one serial
driver pass that calls `MPI.Init()`. On an InfiniBand cluster the default
libfabric (OFI) provider is `verbs;ofi_rxm` (mlx5), which allocates an RDMA queue
pair at init time. On a **login node** (verbs disabled) or when the locked-memory
limit is low (`ulimit -l`), that allocation is refused and precompilation aborts:

```
n0096:rank0.julia.bin: Failed to modify UD QP to INIT on mlx5_0: Operation not permitted
MPIDI_OFI_init_local ... create_vni_context: Cannot allocate memory
ERROR: The following 1 direct dependency failed to precompile: Jexpresso
```

Jexpresso now steers the precompile worker onto the shared-memory `shm` provider
automatically (`src/Jexpresso.jl`, `@setup_workload`), so a plain
`julia --project=. -e 'using Pkg; Pkg.precompile()'` succeeds out of the box on a
login node. `shm` is intra-node only and does no network-interface probing, so it
avoids both the verbs/mlx5 failure above *and* the multi-minute stall the `tcp`
provider can hit while enumerating every NIC (IPoIB, bonded, …) during
`MPI.Init`. The override is scoped to precompilation only and is skipped if you
have already pinned `FI_PROVIDER` yourself — so your compute-node runs keep using
verbs/mlx5.

If you still hit this (e.g. with a customized `FI_PROVIDER`), force a benign
provider just for the precompile step, or precompile inside a compute-node
allocation where verbs and locked memory are permitted:

```bash
# Login node: sidestep the RDMA fabric for precompilation (shm = no network)
FI_PROVIDER=shm julia --project=. -e 'using Pkg; Pkg.precompile()'

# tcp also works but probes every NIC and can stall for minutes on some
# login nodes; prefer shm for the single-process precompile.

# Or precompile on a compute node where verbs/locked memory are allowed
salloc -N1 -n1 -t 0:30:00      # your partition/account flags
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

Quick way to find a provider that initializes fast on your node:

```bash
for p in shm tcp; do
  echo -n "$p: "
  timeout 30 env FI_PROVIDER=$p julia --project=. -e \
    'using MPI; MPI.Init(); println("ok, size=", MPI.Comm_size(MPI.COMM_WORLD)); MPI.Finalize()' \
    || echo "HANG/FAIL"
done
```

### A run is "stuck" for ~30–60 s before the time loop advances

**A.** That is the one-time JIT compilation cost (`sem_setup`, the SciML
integrator, the VTK writer, etc.), not a hang. For serial development, use the
interactive REPL workflow so you pay it only once per session; see
[INSTALL.md, Section 6](INSTALL.md#6-daily-workflow--interactive-repl-for-fast-iteration).
Every fresh `mpirun` is its own cold start.

---

## Installation

### Which MPI route should I use — OpenMPI, MPICH, or the bundled `MPICH_jll`?

**A.** All three are supported for standalone Jexpresso; you pick one and tell
`MPIPreferences` about it:

- **OpenMPI / MPICH (system binary):** install via your package manager and run
  `MPIPreferences.use_system_binary()`. Best on Linux/HPC.
- **`MPICH_jll` (native, bundled with MPI.jl):** nothing to install, run
  `MPIPreferences.use_jll_binary()`. Best on laptops, or when system OpenMPI
  deadlocks on macOS.

The full comparison and step-by-step setup are in
[INSTALL.md, Section 5](INSTALL.md#5-running-in-parallel-with-mpi). **Note:** for
**coupled** runs you must use a *system* MPI for both codes — see
[RUN-COUPLED.md](RUN-COUPLED.md).

### OpenMPI 5 + macOS deadlocks on `MPI.Init`

**A.** A known sharp edge where the PMIx handshake never completes. Switch to the
bundled JLL MPI:

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
```

See [INSTALL.md, Section 5.7](INSTALL.md#57-troubleshooting).

### I switched MPI and now things behave strangely / won't bind

**A.** Remove the recorded preference and reconfigure from scratch:

```bash
rm -f LocalPreferences.toml
```

Then redo the configure → rebuild → precompile steps in
[INSTALL.md, Section 5.2–5.3](INSTALL.md#52-point-mpipreferences-at-your-mpi).

### How do I verify which MPI `MPI.jl` is actually using?

**A.**

```bash
julia --project=. -e '
  using MPIPreferences; println("binary = ", MPIPreferences.binary)
  using MPI;            println(MPI.identify_implementation())'
```

`binary = "system"` means a system MPI (the implementation/version follows);
`binary = "MPItrampoline_jll"` means the bundled native MPI.

### Precompilation takes a long time / packages fail to resolve

**A.** The first `Pkg.precompile()` is expected to be slow as Julia compiles all
dependencies. If packages fail to resolve, make sure you instantiated first
(`Pkg.instantiate()`) and that you are using a supported Julia version
(1.11.9 recommended). The pinned compatible package versions are listed at the
bottom of [INSTALL.md](INSTALL.md#notes-on-package-list).
