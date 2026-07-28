# Running Jexpresso on a Slurm cluster

`submit_jexpresso_parallel.sh` is a generic, site-agnostic batch script for
running Jexpresso across many nodes. The site-specific bits (modules, account,
partition, QoS) are marked `CHANGE_ME` — edit those, and nothing else, to get
started.

```bash
sbatch --export=ALL,EQS=CompEuler,CASE=theta \
       auxiliary/slurm/submit_jexpresso_parallel.sh
```

## Why the script has two phases

The one thing that reliably ruins a large Julia MPI job is letting every rank
precompile. With a cold depot, all *N* ranks discover the same missing
`.ji`/`.so` files at the same instant, all of them start compiling, and all of
them write into the same shared depot. At 512 ranks that is 512 compilers on
one filesystem: minutes to hours of wall time, lock contention, corrupt cache
files, sometimes a job that never reaches the time loop.

So the script does this instead:

| Phase | Ranks | What happens |
|-------|-------|--------------|
| 1 | **1 process, 1 node** | `Pkg.instantiate()`, `Pkg.precompile()`, optionally the sysimage build. All compilation happens here, exactly once. |
| 2 | all `$SLURM_NTASKS` | The actual run, launched with precompilation **switched off** at the Julia level. |

Phase 2 is locked down three ways:

```
JULIA_PKG_PRECOMPILE_AUTO=0    # Pkg never starts a precompile pass
--compiled-modules=existing    # load existing caches, never create them
--pkgimages=existing           # ditto for the native-code images
```

A rank that still misses a cache falls back to compiling privately in its own
memory. It can no longer write to the shared depot, and it can no longer race
the other ranks for the same file.

Phase 1 does leave the rest of the allocation idle for its duration. That is
the intended trade — it is far cheaper than the alternative. If your
precompile is slow, warm the depot once from a login node or a small
interactive allocation:

```bash
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

and then submit with `PRECOMPILE=0`.

## Knobs

All of these are overridable at submit time via
`sbatch --export=ALL,NAME=value,...`:

| Variable | Default | Meaning |
|----------|---------|---------|
| `EQS` / `CASE` | `CompEuler` / `theta` | the case, i.e. `problems/$EQS/$CASE` |
| `JEXPRESSO_ROOT` | `$SLURM_SUBMIT_DIR` | repo root |
| `PRECOMPILE` | `1` | run phase 1; set `0` if the depot is already warm |
| `PRECOMPILE_TASKS` | `8` | Julia's internal precompile workers (clamped to the node's core count; set `1` for a strictly single-core pass) |
| `BUILD_SYSIMAGE` | `0` | build `jexpresso.so` during phase 1 |
| `SYSIMAGE` | `jexpresso.so` | sysimage to use if present |
| `LAUNCHER` | `auto` | `srun`, `mpirun`, or auto-detect from `MPIPreferences.binary` |
| `SRUN_MPI` | site default | value for `srun --mpi=` (`pmix` for Open MPI, `pmi2` for MPICH/Intel) |
| `TAG_OUTPUT` | `1` | prefix each output line with its rank |
| `JULIA` | `julia` | which Julia to use |

## Things that bite on a big machine

**MPI.jl must be bound to the system MPI.** The script prints
`MPIPreferences.binary` and warns if it is anything else. A bundled JLL binary
cannot be launched by `srun`'s PMI, and it is not using the cluster's
high-speed fabric either — you would be running an InfiniBand machine over
TCP. Fix it with `MPIPreferences.use_system_binary()` (INSTALL.md §5.2), and
make sure the MPI module loaded in the script is the same one you bound
against. A mismatch does not error: it hangs inside `MPI_Init`.

**Do not stage the depot to node-local scratch.** It is a tempting way to cut
filesystem load, and it silently defeats the whole design: Julia's cache files
record the absolute paths of the sources they were built from, so a depot at a
different path is a depot with no valid caches, and every rank recompiles.
Keep `JULIA_DEPOT_PATH` on a shared filesystem all compute nodes can see.

**Heterogeneous partitions need `JULIA_CPU_TARGET`.** Caches built on one
microarchitecture are rejected on another, per rank, at load time. If your
nodes are not identical, uncomment and set the multiversioned target in the
script.

**One core per rank.** Jexpresso is essentially pure MPI, so the script pins
`JULIA_NUM_THREADS=OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=1`. Without that, each
of the 64 ranks on a node starts 64 BLAS threads and they fight over the same
cores.

**A sysimage is the next big win.** Without one, every rank JIT-compiles the
solver on startup — tens of seconds, on every rank, on every launch.
`BUILD_SYSIMAGE=1` builds it in phase 1; it needs `using Revise` commented out
in `src/Jexpresso.jl` (INSTALL.md §6.1). Rebuild it after changing
`Project.toml` or `src/`, not after editing `user_inputs.jl`.

**The mesh/SEM cache is per-rank-count.** `<case>/.jexpresso_cache/` files are
keyed by rank and by the total number of parts, so a serial warm-up run does
not populate them for a parallel one, and changing the rank count rebuilds
them. The first run at a given rank count pays that cost.

## Site-specific examples

`auxiliary/wulver/` holds the older, hardcoded scripts for NJIT's Wulver
cluster. They are kept for reference; prefer this one for anything new.
