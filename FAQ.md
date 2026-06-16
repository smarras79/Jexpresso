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

### Precompiling on an InfiniBand cluster hangs or fails at `MPI.Init` (`Failed to modify UD QP to INIT on mlx5_0`)

**A.** Historically Jexpresso's `@compile_workload` ran one serial driver pass that
called `MPI.Init()` to bake in the integrator/RHS specializations. On an
InfiniBand cluster that init brings up the libfabric (OFI) fabric, which is
hostile to precompilation on a **login node** in several ways:

- the `verbs;ofi_rxm` (mlx5) default allocates an RDMA queue pair and aborts when
  verbs are disabled or locked memory is low (`ulimit -l`):

  ```
  n0096:rank0.julia.bin: Failed to modify UD QP to INIT on mlx5_0: Operation not permitted
  MPIDI_OFI_init_local ... create_vni_context: Cannot allocate memory
  ERROR: The following 1 direct dependency failed to precompile: Jexpresso
  ```

- the `tcp` provider enumerates every NIC (IPoIB, bonded, …) and does reverse-DNS
  during `MPI.Init`, which can **stall for many minutes**;
- even `shm` can stall during fabric bring-up on some builds.

Jexpresso handles both problems so the warm-up still happens (and your first
`run_case` stays fast) without precompilation hanging:

1. **Lean workload.** The precompile pass no longer runs a full 2000-step sod1d
   simulation — `drivers.jl` caps it to 3 timesteps while `jl_generating_output`
   is set. That still triggers every RHS / integrator / callback specialization,
   so the cache is fully warmed, but precompilation finishes in seconds-to-a-minute
   instead of many minutes. (If you saw it "hang at `Progress 0/1`" for 7+ minutes,
   that was the full simulation running, not a true hang — the lean pass fixes it.)
2. **shm fabric.** For the single-process precompile worker, libfabric is steered
   onto the shared-memory `shm` provider (unless you've pinned `FI_PROVIDER`),
   which needs no network and avoids both the verbs/mlx5 abort and the tcp NIC
   stall. It's restored afterwards, so compute-node runs are unaffected.

So a plain `julia --project=. -e 'using Pkg; Pkg.precompile()'` should now complete
on a login node **and** keep cold runs fast.

If `MPI.Init` still can't come up at all on your node (no working fabric
provider), disable the workload — the package still precompiles fully, you just
lose the warm-up so the first run pays more JIT:

```bash
JEXPRESSO_PRECOMPILE_WORKLOAD=0 julia --project=. -e 'using Pkg; Pkg.precompile()'
```

Or precompile inside a compute-node allocation where the fabric is healthy:

```bash
salloc -N1 -n1 -t 0:30:00                 # your partition/account flags
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

To check which provider initializes fast on your node (the default is `shm`):

```bash
for p in shm tcp; do
  echo -n "$p: "
  timeout 30 env FI_PROVIDER=$p julia --project=. -e \
    'using MPI; MPI.Init(); println("ok, size=", MPI.Comm_size(MPI.COMM_WORLD)); MPI.Finalize()' \
    || echo "HANG/FAIL"
done
# if shm is slow but another provider is fast, pin it for precompile:
FI_PROVIDER=<fast one> julia --project=. -e 'using Pkg; Pkg.precompile()'
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
