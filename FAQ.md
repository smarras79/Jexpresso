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

**A.** That `@check` is a **dimensional guard**, not a corruption symptom: line
325 asserts that the coarse model's point/embedding dimension `Dp` equals its
cell dimension `Dc`. p4est is not manifold-aware and requires `Dp == Dc`. The
mismatch comes from GridapGmsh: its `GmshDiscreteModel` reports `Dp=3` for **any**
`.msh` whose nodes carry a non-zero z-coordinate — see `_setup_point_dim`, which
returns 3 the moment one node has `z !≈ 0`, with no keyword to override it. So a
logically-2D mesh with a stray non-zero z reads back as `{Dc=2, Dp=3}`: it loads
fine for a **non-AMR** case (e.g. `theta` — plain `GridapDistributed` assembles
happily on a `Dc=2, Dp=3` embedded model) but every AMR run aborts here. Because
the forest is 2D (the p4est log shows `121 = 11×11` nodes), the *cell* check
passes and it is always the *point-dim* check that fires.

**This is now handled in code.** `mod_mesh_read_gmsh!` projects the coarse Gmsh
model down to `Dp=Dc` (`_flatten_model_to_cell_dim` in `src/kernel/mesh/mesh.jl`)
before building the octree model, dropping the trailing coordinate component.
Face/boundary numbering is derived from cell connectivity, not coordinates, so
boundary tags are preserved; the flatten is a no-op when `Dp == Dc`. AMR cases
on 2D-embedded-in-3D meshes therefore run without touching the `.msh`. If you
still hit this after pulling, confirm you're running the patched
`mod_mesh_read_gmsh!` and check what dimension your mesh reports:

```bash
julia --project=. -e '
  using GridapGmsh
  m = GmshDiscreteModel("path/to/your.msh")
  println(typeof(m))   # {2,3,...} means Dp=3 (the flatten handles it);
                        # {2,2,...} means Dp already equals Dc
'
```

**Separately, on macOS you still need the patched `GridapP4est` fork** — not for
*this* assertion (it fires before any p4est callback runs) but for the actual
refinement steps that follow it: the registered `0.3.11` lacks ARM64
`@cfunction` support and has a p4est-iterator struct-stride mismatch (ARM64 on
Julia ≥ 1.11, x86_64 on Julia ≥ 1.12). The fork is now **pinned automatically**
via a `[sources]` block in `Project.toml`, so a fresh `Pkg.instantiate()`
resolves it. Verify with:

```bash
julia --project=. -e 'using Pkg; Pkg.status("GridapP4est")'
# expected: "...#arm64-cfunction-fix"
```

See [INSTALL.md, Section 7](INSTALL.md#7-amr-on-macos-apple-silicon-the-patched-gridapp4est-fork)
for details.

The in-code flatten means you no longer *have* to fix the mesh, but if you'd
rather remove the stray z at the source, regenerate it (or edit the `.geo`) so
every `Point(...)` uses `0` for the extra coordinate rather than a nonzero
constant — a common cause is a copy-paste bug where a mesh-spacing variable
(e.g. `gridsize`) ends up in the z slot instead of `0`.

### AMR fails with `could not load library "libp4est.4.dylib"` / `Library not loaded: @rpath/libmpi.12.dylib`

**Q.** An AMR case (`theta_amr`, or anything with `:lamr`/`:linitial_refine`)
aborts immediately after reading the mesh:

```
Info    : Done reading './meshes/gmsh_grids/hexa_TFI_10x10.msh'
ERROR: LoadError: could not load library ".../artifacts/904c551e.../lib/libp4est.4.dylib"
dlopen(...): Library not loaded: @rpath/libmpi.12.dylib
  Reason: tried: '.../lib/./libmpi.12.dylib' (no such file),
          '/Applications/Julia-1.11.app/Contents/Resources/julia/lib/libmpi.12.dylib' (no such file),
          '/usr/local/lib/libmpi.12.dylib' (no such file),
          '/usr/lib/libmpi.12.dylib' (no such file, not in dyld cache)
Stacktrace:
  [1] p4est_connectivity_new
    @ ~/.julia/packages/P4est_wrapper/.../src/bindings/p4est_api.jl:361 [inlined]
  [2] setup_pXest_connectivity_from_geometry(...)
    @ GridapP4est ...
```

Non-AMR cases run fine on the same machine, and often the *same* case runs on a
different machine.

**A.** This is an **MPI binding mismatch**, not a p4est or mesh problem. Read the
soname: `libmpi.12` is the **MPICH** ABI (OpenMPI is `libmpi.40`). `P4est_jll` is
built against `MPICH_jll`, and the `MPICH_jll` artifact directory only lands on
the loader's `@rpath` when `MPICH_jll` is *actually loaded* — which happens only
when `MPI.jl` is bound to the JLL binary. If `MPIPreferences` points at a
**system** MPI, `MPICH_jll` is never loaded, nothing supplies
`libmpi.12.dylib`, and `dlopen` fails exactly as above. The dyld search list in
the error is the giveaway: artifact dir, Julia's own lib dirs, `/usr/local/lib`,
`/usr/lib` — no `MPICH_jll` artifact, no Homebrew prefix.

Two things follow from this, and they explain the usual symptoms:

- **Only AMR breaks.** `GridapP4est`/`P4est_wrapper` are the only packages that
  dlopen `libp4est`. Every non-AMR case avoids p4est entirely and runs happily
  on a system-MPI binding.
- **It differs between your own machines.** `LocalPreferences.toml` — the file
  that records the binding — is in `.gitignore`, so it is per-machine and never
  travels with a clone. One laptop on `use_jll_binary()` and one desktop on
  `use_system_binary()` is the common way to see this.

First confirm the binding ([INSTALL.md §5.4](INSTALL.md#54-verify-the-binding)):

```bash
julia --project=. -e '
  using MPIPreferences; println("binary  = ", MPIPreferences.binary)
  using MPI;            println("impl    = ", MPI.identify_implementation())
                        println("libmpi  = ", MPI.API.libmpi)'
```

If `binary` is anything other than `MPICH_jll`, that is the cause. Rebind with
[INSTALL.md §5.2 Route C](INSTALL.md#route-c--mpich_jll-native) and rebuild —
**the cache clearing is not optional**, or you keep the old library under the
new preference:

```bash
rm -f LocalPreferences.toml
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

If you must keep a system MPI on that machine (e.g. for coupled Alya runs, which
need a system `mpif90` anyway), it has to be **MPICH**, not OpenMPI: Homebrew
`mpich` provides `libmpi.12.dylib`, `open-mpi` provides `libmpi.40.dylib` and can
never satisfy `P4est_jll`. Even with MPICH you will usually also need its lib
directory on the loader path, since Homebrew's prefix is not in the list dyld
searched above:

```bash
export DYLD_FALLBACK_LIBRARY_PATH="$(brew --prefix mpich)/lib:${DYLD_FALLBACK_LIBRARY_PATH:-}"
```

`bash tools/check_mpi_setup.sh` walks the whole binding end to end and reports
which implementation each layer resolved to. See also
[Which MPI route should I use](#which-mpi-route-should-i-use--openmpi-mpich-or-the-bundled-mpich_jll)
and [I switched MPI and now things behave strangely](#i-switched-mpi-and-now-things-behave-strangely--wont-bind).

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

If you're on Apple Silicon and running an AMR case, the patched `GridapP4est`
fork is already pinned for you via a `[sources]` block in `Project.toml` — see
[INSTALL.md, Section 7](INSTALL.md#7-amr-on-macos-apple-silicon-the-patched-gridapp4est-fork).
Its missing ARM64 `@cfunction`/struct-stride support is a separate issue from
this artifact-corruption one, but both surface on the same machines, so verify
`Pkg.status("GridapP4est")` shows the `#arm64-cfunction-fix` fork after
reinstantiating.

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
