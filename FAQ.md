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

### AMR (`theta_amr` and other `:lamr`/`:linitial_refine` cases) fails with `AssertionError: A check failed` in `OctreeDistributedDiscreteModel`

**Q.** Loading a mesh into an AMR case fails with something like:

```
ERROR: LoadError: AssertionError: A check failed
Stacktrace:
  [1] macro expansion
    @ ~/.julia/packages/Gridap/.../src/Helpers/Macros.jl:61 [inlined]
  [2] GridapP4est.OctreeDistributedDiscreteModel(Dc::Int64, Dp::Int64, ...)
    @ GridapP4est ~/.../GridapP4est.jl/src/OctreeDistributedDiscreteModels.jl:325
```

**A.** This is a mesh-file issue, not a code bug: `GridapP4est`/p4est is not
manifold-aware — it requires the cell dimension (`Dc`) and the point/embedding
dimension (`Dp`) of the coarse mesh to be equal. Your `.msh` file has
`Dp != Dc` (e.g. a mesh that is logically 2D but whose vertices carry a
nonzero z-coordinate, so Gridap's GMSH reader reports `Dp=3`).

The same mesh loads without any problem for a **non-AMR** case (e.g. `theta`)
because plain `GridapDistributed` happily assembles on `Dc=2, Dp=3` embedded
meshes — the mismatch only breaks the octree/AMR path.

To fix it at the source, check whether your mesh's z-coordinate is truly zero
everywhere:

```bash
julia --project=. -e '
  using GridapGmsh
  m = GmshDiscreteModel("path/to/your.msh")
  println(typeof(m))   # look for e.g. UnstructuredDiscreteModel{2,2,...}
                        # if you see {2,3,...} instead, Dp != Dc
'
```

If it prints `{Dc,Dp,...}` with `Dp > Dc`, regenerate the mesh (or edit the
`.geo`) so every `Point(...)` uses `0` for the extra coordinate rather than a
nonzero constant — a common cause is a copy-paste bug where a mesh-spacing
variable (e.g. `gridsize`) ends up in the z slot instead of `0`.

### A package fails to precompile with a missing file inside a Julia artifact

**Q.** `Pkg.instantiate()`/`Pkg.precompile()` (or a run that triggers a lazy
precompile) fails with a missing file somewhere under `~/.julia/artifacts/`.
Two examples seen in the wild:

```
ERROR: LoadError: could not load library ".../artifacts/.../lib/libp4est.4.dylib"
dlopen(.../libp4est.4.dylib, 0x0001): Library not loaded: @rpath/libjansson.4.dylib
  Referenced from: <...> .../libp4est.4.dylib
  Reason: tried: '.../lib/./libjansson.4.dylib' (no such file), ...
```

```
ERROR: LoadError: SystemError: opening file
".../artifacts/.../lib/gmsh.jl": No such file or directory
```

The first is `P4est_jll` missing its `Jansson_jll` dependency's dylib; the
second is `Gmsh_jll` (pulled in by `GridapGmsh`) missing its own bundled
`gmsh.jl`. Different packages, same root cause — see below.

**A.** This is a corrupted/incomplete Julia artifact on disk, not a code bug.
Artifacts are content-addressed directories under `~/.julia/artifacts/<hash>/`
that Pkg downloads once and then trusts are complete forever — if the
download was interrupted, partial, or a file inside got deleted/modified
after the fact, Pkg has no way to notice and will keep handing out the
broken directory. On macOS specifically, Gatekeeper quarantining/stripping a
freshly-downloaded artifact is a common trigger.

Note Pkg.jl's live-updating progress bar can also *hide* which package
actually failed and why — if you only see a generic "Failed to precompile
X" with no real error above it, re-run that one package in isolation to see
past the redraw, e.g. `julia --project=. -e 'using Pkg; Pkg.precompile("GridapGmsh")'`
or just `julia --project=. -e 'using GridapGmsh'`.

Fix, step by step:

1. Quit any running Julia session first.
2. Clear the artifact cache so Pkg redownloads everything cleanly (this
   affects every Julia project on the machine, but is the most reliable
   fix — artifacts are content-addressed and just redownload on demand):
   ```bash
   rm -rf ~/.julia/artifacts
   ```
3. From your Jexpresso project directory, reinstantiate:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   ```
4. Rebuild the packages that link against native libraries:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.build()'
   ```
5. Precompile:
   ```bash
   julia --project=. -e 'using Pkg; Pkg.precompile()'
   ```
6. Re-run whatever case failed originally.

If the exact same error comes back immediately after step 2–3 (i.e. the
redownloaded artifact is *still* missing the file), macOS Gatekeeper may be
quarantining it. Check and clear the quarantine flag, then repeat steps 3–5:
```bash
xattr -lr ~/.julia/artifacts | grep -i quarantine   # check first
xattr -dr com.apple.quarantine ~/.julia/artifacts   # clear if present
```

If you're on Apple Silicon and running an AMR case, also make sure you're on
the patched `GridapP4est` fork, not the registered `=0.3.11` release — see
[INSTALL.md, Section 7](INSTALL.md#7-amr-on-macos-apple-silicon-use-the-patched-gridapp4est-fork).
That's a separate issue (missing ARM64 `@cfunction`/struct-stride support)
from this artifact-corruption one, but both surface on the same machines.

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
