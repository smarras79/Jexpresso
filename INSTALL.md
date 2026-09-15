# Installing Jexpresso

This guide takes you from a fresh clone to a running **Jexpresso**, first in
serial and then in parallel with MPI.

> **Hit an error?** See [FAQ.md](FAQ.md) for fixes to common installation and
> run problems.

**Which sections do you need?**

| I want to… | Read |
|---|---|
| Run serially on a laptop | §1–§4 |
| Run in parallel with MPI | §1–§4, then §5 |
| Run AMR cases (`theta_amr`, …) on macOS | §1–§4, §5, then §7 |
| Iterate quickly on code | §6 |

---

## Prerequisites

- **Julia 1.11.9** (recommended; 1.12.6 also works, but we stay on 1.11 for now).
- **Git** with SSH access to GitHub ([setup guide](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)).

## 1. Download the repositories

```bash
git clone git@github.com:smarras79/Jexpresso.git
git clone git@github.com:smarras79/JexpressoMeshes.git
```

`JexpressoMeshes` holds the sample meshes used by the bundled test cases.

## 2. Link the sample meshes

```bash
cd Jexpresso
ln -s ../JexpressoMeshes/meshes .
```

## 3. Build and precompile (serial)

From inside `Jexpresso/`:

**3a. Instantiate the dependencies** (precompilation deferred so we control it
explicitly):

```bash
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
```

> **Going parallel?** Stop here and do [§5](#5-running-in-parallel-with-mpi)
> *before* precompiling. Binding MPI rebuilds `MPI` and the C libraries that
> link against it, which forces a recompile anyway — doing MPI first saves you
> a full precompilation pass.

**3b. Precompile:**

```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

The first pass takes a while.

## 4. Test the installation

```bash
julia --project=.
```

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "sod1d")
```

If the case runs to completion, your serial Jexpresso is ready. 🎉

---

## 5. Running in parallel with MPI

`MPI.jl` does not contain an MPI implementation; it binds to one at build time,
and `MPIPreferences` records which. The workflow is the same for every route:

> **install MPI → point `MPIPreferences` at it → rebuild `MPI` and the
> MPI-linked C libraries → precompile → launch with the matching `mpiexec`.**

You may have several MPIs installed at once and rebind whenever you like
(§5.2, case 2). What is never allowed is mixing them **within one run** — the
`libmpi` MPI.jl loads, the `mpiexec` that launches the job, the C libraries
compiled against MPI (§5.3), and, for coupled runs, the MPI `Alya.x` was linked
with must all come from the same installation.

| Route | What provides MPI | When to choose it |
|-------|-------------------|-------------------|
| **A. OpenMPI** (system) | An OpenMPI you install | Linux clusters where OpenMPI is the site default |
| **B. MPICH** (system) | An MPICH you install | You prefer MPICH, or the cluster ships it |
| **C. MPICH_jll** (bundled) | Shipped inside Julia's package environment | **Laptops/desktops (recommended on macOS)**, no admin rights, or OpenMPI deadlocks on macOS |

### 5.1 Install an MPI implementation

**Route A — OpenMPI**

```bash
sudo apt install libopenmpi-dev openmpi-bin   # Ubuntu/Debian
brew install open-mpi                         # macOS
mpiexec --version                             # verify
```

**Route B — MPICH**

```bash
sudo apt install mpich libmpich-dev           # Ubuntu/Debian
brew install mpich                            # macOS
mpiexec --version                             # verify
```

**Route C — MPICH_jll**

Nothing to install. MPI.jl ships MPItrampoline (MPICH-based) together with its
own `mpiexec`, which is the launcher you must use (§5.6). It is unaffected by
any system MPI on the machine.

> Keeping OpenMPI *and* MPICH installed side by side is fine and often useful
> (OpenMPI 5 on macOS is prone to hanging in `MPI_Init`). Just run both install
> commands; §5.2 case 2 shows how to select between them.

### 5.2 Point `MPIPreferences` at your MPI

Run **one** of these. The choice is recorded in `LocalPreferences.toml` in the
project root.

> **Two things get selected, not one.** `MPIPreferences` records which
> **library** Julia loads; your shell decides which **launcher** starts the
> job. A run whose launcher and library come from different MPIs does not
> error — it hangs in `MPI_Init`, or gives every rank its own world of size 1.
> Whenever you bind one, pin the other to match.

**Route C — MPICH_jll (simplest)**

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'
```

**Route A or B, case 1 — only one system MPI installed**

Standard system paths (`/usr/bin`, `/usr/local/bin`):

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
```

Elsewhere (`/opt/...`, Homebrew on Apple Silicon), pass the `lib` directory:

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/PATH/TO/MPI/lib"])'
```

For a single Homebrew MPI on Apple Silicon that is `/opt/homebrew/lib` —
**but not if both OpenMPI and MPICH are installed**; use case 2.

**Route A or B, case 2 — both OpenMPI and MPICH installed**

The bare `use_system_binary()` binds whichever `libmpi` the loader finds first,
which need not match the `mpiexec` on your `PATH`. Name the installation
explicitly.

*Step 1 — find each prefix.*

| Platform | OpenMPI prefix | MPICH prefix |
|---|---|---|
| macOS, Homebrew (Apple Silicon) | `$(brew --prefix open-mpi)` → `/opt/homebrew/opt/open-mpi` | `$(brew --prefix mpich)` → `/opt/homebrew/opt/mpich` |
| macOS, Homebrew (Intel) | `/usr/local/opt/open-mpi` | `/usr/local/opt/mpich` |
| Ubuntu/Debian | `/usr/lib/x86_64-linux-gnu/openmpi` | `/usr/lib/x86_64-linux-gnu/mpich` |
| HPC modules | `$MPI_HOME` after `module load openmpi` | `$MPI_HOME` after `module load mpich` |

*Step 2 — bind library **and** launcher.* `extra_paths` is searched before the
default loader paths, and passing `mpiexec` records the matching launcher so the
two halves cannot drift apart.

> ⚠️ **Paste one block at a time and never put a trailing `# comment` on the
> `PREFIX=` line.** Interactive `zsh` (the macOS default) does not treat `#`
> as a comment unless `setopt interactivecomments` is set; it runs the command
> `#` with the variable set only for that command, so the variable is never
> set and the `julia` call fails with `extra_paths = ["/lib"]`.

To use OpenMPI:

```bash
OMPI_PREFIX=$(brew --prefix open-mpi)
echo "prefix = $OMPI_PREFIX"
ls "$OMPI_PREFIX"/lib/libmpi.*
```

```bash
julia --project=. -e "using MPIPreferences; MPIPreferences.use_system_binary(
        extra_paths = [\"$OMPI_PREFIX/lib\"],
        mpiexec     = \"$OMPI_PREFIX/bin/mpiexec\")"
```

To use MPICH:

```bash
MPICH_PREFIX=$(brew --prefix mpich)
echo "prefix = $MPICH_PREFIX"
ls "$MPICH_PREFIX"/lib/libmpi.*
```

```bash
julia --project=. -e "using MPIPreferences; MPIPreferences.use_system_binary(
        extra_paths = [\"$MPICH_PREFIX/lib\"],
        mpiexec     = \"$MPICH_PREFIX/bin/mpiexec\")"
```

On Linux set the prefix by hand (e.g.
`MPICH_PREFIX=/usr/lib/x86_64-linux-gnu/mpich`) and run the same `julia` line.

> `MPI library could not be found ... extra directories: ["/lib"]` means the
> prefix variable was empty — either the `zsh` comment trap above, or
> `brew --prefix <formula>` failed because the formula is not installed.

*Step 3 — put the matching launcher first on `PATH`.* Define both switches in
`~/.zshrc` / `~/.bashrc` and call the one you need before launching:

```bash
# macOS / Homebrew
use-openmpi() { export PATH="$(brew --prefix open-mpi)/bin:$PATH"; hash -r; }
use-mpich()   { export PATH="$(brew --prefix mpich)/bin:$PATH";    hash -r; }

# Linux — substitute the prefixes from step 1
use-openmpi() { export PATH="/usr/lib/x86_64-linux-gnu/openmpi/bin:$PATH"; hash -r; }
use-mpich()   { export PATH="/usr/lib/x86_64-linux-gnu/mpich/bin:$PATH";   hash -r; }
```

Before every parallel run, confirm all three agree: `which mpiexec mpirun mpif90`.

> **Homebrew:** `open-mpi` and `mpich` both provide `mpicc`/`mpiexec`/`libmpi`,
> so only one is *linked* into `/opt/homebrew/{bin,lib}` at a time. Both stay
> usable under their own `opt` prefixes, which is why the commands above never
> use the shared `/opt/homebrew/lib`. Rebinding MPI.jl does **not** require
> `brew link`.
>
> **Ubuntu/Debian:** `update-alternatives --config mpi` only changes what the
> bare `mpicc`/`mpirun` names resolve to; it has no effect on MPI.jl.

### 5.3 Rebuild `MPI` **and the MPI-linked C libraries**, then precompile

Two things in the dependency tree are compiled against a specific `libmpi`:

- `MPI` itself.
- **`P4est_wrapper`** — a small C shim that `GridapP4est` compiles on your
  machine (its `deps/build.jl` links against the `libp4est` and `libmpi`
  present *at build time*). `Pkg.instantiate()` builds it once, in §3a, against
  whatever MPI was configured then — and never rebuilds it on its own. If you
  bind a different MPI afterwards and skip this step, the shim stays linked to
  the old `libmpi`, and the first AMR call `dlopen`s two MPIs into one process.
  On macOS that shows up as a **`Segmentation fault: 11` in
  `_platform_memmove` right after the gmsh mesh is read**, with a backtrace
  that is nothing but `dyld` frames.

So after **every** `MPIPreferences` change, rebuild all of them, in this order:

```bash
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.build("P4est_wrapper"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.build("GridapP4est"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

`Pkg.build("P4est_wrapper")` prints the path of its `build.log`; the log
records which `libmpi` was linked, so you can confirm it matches §5.4.

### 5.4 Verify the binding

```bash
julia --project=. -e '
  using MPIPreferences; println("binary  = ", MPIPreferences.binary)
  using MPI;            println("impl    = ", MPI.identify_implementation())
                        println("libmpi  = ", MPI.API.libmpi)
                        println("mpiexec = ", MPI.mpiexec())'
```

- **Route A:** an Open MPI implementation, `binary = "system"`.
- **Route B:** MPICH, `binary = "system"`.
- **Route C:** `binary = "MPItrampoline_jll"`.

With both OpenMPI and MPICH installed, also check the paths: `libmpi` must sit
under the prefix chosen in §5.2 (not the ambiguous `/opt/homebrew/lib`), and
`mpiexec` must come from the same prefix **and** match `which mpiexec` in the
shell you launch from.

To also confirm `P4est_wrapper` links the same MPI (macOS shown; use `ldd` on
Linux):

```bash
otool -L $(find ~/.julia/scratchspaces -name 'libp4est_wrapper*.dylib') | grep mpi
```

`bash tools/check_mpi_setup.sh` runs all of these checks, plus the launcher,
and times out instead of hanging.

### 5.5 macOS hostname fix (MPICH and MPICH_jll)

Required for any MPICH-based MPI on macOS (Routes B and C), and harmless — and
recommended — on OpenMPI, whose PRRTE launcher also resolves the hostname and
hangs silently when it cannot. Without it `MPI_Init` fails with
`gethostbyname failed, <hostname> (errno 0)`.

```bash
echo "127.0.0.1   $(hostname -s)" | sudo tee -a /etc/hosts
echo "127.0.0.1   $(hostname)"    | sudo tee -a /etc/hosts
ping -c 1 $(hostname -s)     # should answer from 127.0.0.1
```

No `sudo`? Export `MPICH_INTERFACE_HOSTNAME=127.0.0.1` in every shell instead
(or add it to `~/.zshrc`).

### 5.6 Launch a parallel run

Use the launcher that matches your route. Mixing a system `mpiexec` with the
JLL library (or vice versa) is the most common "it won't start" failure.

**Route A or B — system MPI**

```bash
mpiexec -n 4 julia --project=. src/Jexpresso.jl CompEuler 3d
```

With several MPIs installed, use absolute paths to both launcher and `julia`:

```bash
/opt/homebrew/opt/mpich/bin/mpiexec -n 4 \
  /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia \
  --project=. src/Jexpresso.jl CompEuler theta
```

**Route C — MPICH_jll**

Do **not** use a system `mpirun`. Launch with the `mpiexec` MPI.jl ships:

```bash
julia --project=. -e '
  using MPI
  run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler city2d`)'
```

For daily use, [`jexp_mpich.sh`](jexp_mpich.sh) in the repo root wraps this
(edit the `JULIA=` path at its top):

```bash
./jexp_mpich.sh 4 CompEuler city2d
```

### 5.7 Troubleshooting

- **Hang in `MPI.Init` with OpenMPI 5 on macOS (Route A).** A known PMIx
  handshake deadlock. Switch to Route C: redo §5.2–§5.6 with the JLL binary.

- **Segfault in `_platform_memmove` right after `Done reading … .msh` (AMR
  cases, macOS).** `P4est_wrapper` was built against a different MPI than the
  one MPI.jl now loads. Rebuild per §5.3.

- **MPICH aborts in `MPI_Finalize` with an `OFI` / `nic=utunN` error.**
  MPICH's libfabric netmod picked a VPN/tunnel device as its NIC. Pin the
  provider to the machine's real default interface:
  ```bash
  route get default | awk '/interface:/{print $2}'   # usually en0
  export FI_PROVIDER=tcp
  export FI_TCP_IFACE=en0
  ```
  Do not use `lo0` (the `tcp` provider often cannot open an endpoint on it,
  giving `ep_enable failed … Bad file descriptor`). If `tcp` fails,
  `export FI_PROVIDER=sockets` needs no interface at all.
  `bash tools/check_mpi_setup.sh` tries the combinations for you.

  To make it permanent, put the winning pair in `~/.zshrc`, guarded so it does
  not leak onto a cluster with a real fabric:
  ```bash
  if [[ "$(hostname -s)" == "my-laptop" ]]; then
    export FI_PROVIDER=tcp
    export FI_TCP_IFACE=en0
  fi
  ```
  Notes: `FI_*` affects every libfabric user on the machine (OpenMPI's OFI
  components included); `FI_TCP_IFACE` must exist at `MPI_Init`, so Wi-Fi off
  or a dock/VPN renumbering can bring the failure back — `FI_PROVIDER=sockets`
  with no interface pinned is more robust if that bites.

  > **Never `FI_PROVIDER=shm` on macOS** — it is Linux-only. Asking for it
  > leaves libfabric with no provider and `MPI_Init` fails with
  > `getinfo failed (default nic=(n/a) …)`. `nic=(n/a)` means "no provider
  > matched"; check `echo $FI_PROVIDER` and switch to `tcp`/`sockets` or
  > `unset` it. `fi_info -l` lists the providers your build has.

- **Stale binding / library conflicts.** `rm -f LocalPreferences.toml`, then
  reconfigure from §5.2 and rebuild per §5.3.

- **Wrong launcher picked up.** `which mpiexec mpirun`; use absolute paths
  (§5.6) if needed.

- **Version mismatch (system MPI).** `mpicc --version` and `mpif90 --version`
  must agree with the runtime.

#### Switching MPI implementations later

Rebinding is the `MPIPreferences` call **plus** a clean rebuild of everything
linked against MPI; otherwise you get the old library under the new preference:

```bash
rm -f LocalPreferences.toml
# → re-run the §5.2 block for the MPI you want
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true); Pkg.build("P4est_wrapper"; verbose=true); Pkg.build("GridapP4est"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'
bash tools/check_mpi_setup.sh
```

> **Coupled runs need one more thing to match.** `AlyaProxy/Alya.x` has its own
> MPI baked in at link time. After switching, rebuild it with the matching
> `mpif90`: `cd AlyaProxy && MPIF90=<prefix>/bin/mpif90 bash compilef90.sh`.
> See [RUN-COUPLED.md](RUN-COUPLED.md) §3 and §5.

### 5.8 Running on a Slurm cluster

```bash
sbatch --export=ALL,EQS=CompEuler,CASE=theta \
       auxiliary/slurm/submit_jexpresso_parallel.sh
```

Edit the `CHANGE_ME` lines in its `#SBATCH` block and the `module load` lines;
the rest is site-agnostic. The script runs a **serial precompile phase**
(`Pkg.instantiate()` + `Pkg.precompile()` on one process) before the
**parallel run phase** (launched with `JULIA_PKG_PRECOMPILE_AUTO=0`,
`--compiled-modules=existing`, `--pkgimages=existing`), so hundreds of ranks
do not all compile into the same depot at once. Full details and
cluster-specific pitfalls (system-MPI binding, `JULIA_CPU_TARGET`, why not to
stage the depot to node-local disk): [`auxiliary/slurm/README.md`](auxiliary/slurm/README.md).

---

## 6. Daily workflow — fast iteration

Every fresh Julia process pays a one-time JIT cost (`sem_setup`, the `with_mpi`
closure, the SciML integrator, the VTK writer…), typically 30–60 s of silent
wall time after `# Read inputs dict ... DONE`.

### Serial: stay in one REPL

```bash
julia --project=.
```

```julia
julia> using Jexpresso
julia> Jexpresso.run_case("CompEuler", "theta")   # slow the first time
julia> Jexpresso.run_case("CompEuler", "theta")   # near-instant after edits
```

Use `Revise.jl` so edits to source files outside `user_inputs.jl` are picked up
without restarting.

### Parallel from the same REPL

```julia
julia> using MPI
julia> run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler city2d`)
```

Each `run(...)` spawns fresh ranks, so every parallel launch is a cold start —
the REPL does not help here. The fix is a sysimage (§6.1).

### 6.1 Fastest parallel execution on macOS with MPICH

On a Mac the cost of a parallel run is dominated by Julia startup, not MPI.

**One-time shell setup** (`~/.zshrc`; see §5.7 for the libfabric pair):

```bash
export FI_PROVIDER=tcp
export FI_TCP_IFACE=en0
export PATH="$(brew --prefix mpich)/bin:$PATH"
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

The thread settings matter: Jexpresso is essentially pure MPI, and without
them *N* ranks each start *N* BLAS threads and fight over the same cores.

**One-time sysimage** — the single biggest win:

```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
julia --project=. create_Jexpresso_sysimage.jl        # writes jexpresso.so (minutes, once)
```

```bash
mpiexec -np 4 -prepend-rank \
    julia --project=. --sysimage jexpresso.so --startup-file=no \
    ./src/Jexpresso.jl CompEuler theta
```

Rebuild it after changing `Project.toml` or `src/` (not after editing
`user_inputs.jl`). `run_coupled.sh` picks it up automatically;
`REBUILD_SYSIMAGE=1 ./run_coupled.sh` regenerates it. The mesh/SEM caches
under `<case>/.jexpresso_cache/` are the second-biggest win and are automatic —
leave `:luse_mesh_cache` at its default.

**How many ranks.** Only the performance cores are worth a rank:

```bash
sysctl -n hw.perflevel0.physicalcpu   # your rank ceiling: 4 on an M1/M2 Air, 8–12 on Pro/Max
```

Oversubscribing makes runs slower. Rank pinning (`-bind-to`, `-map-by`) does
nothing on macOS; don't pass it.

**Launch flags.** `-prepend-rank` tags each rank's output (otherwise stdout is
interleaved and block-buffered and looks like a hang); `--startup-file=no`
skips `startup.jl` on every rank; `-genv VAR VAL` sets a variable for all ranks.

**Coupled runs.** Split the P-cores, e.g. 2 + 2. Alya needs **at least two
ranks** (its rank 0 is a master with no grid points); beyond that, spend extra
ranks on Jexpresso:

```bash
export JEXPRESSO_COUPLED=1
mpiexec -np 2 ./AlyaProxy/Alya.x \
      : -np 2 julia --project=. --sysimage jexpresso.so --startup-file=no \
        ./src/Jexpresso.jl CompEuler 3dAlya
```

or `./run_coupled.sh 2 2`, which applies the sysimage and preflight checks.

**Things that cost you time:** a VPN connected during a run (renumbers the
interface `FI_TCP_IFACE` is pinned to); the repo on iCloud Drive or a network
share (every rank reads mesh and cache files — keep it on the internal SSD);
frequent output (`:diagnostics_at_times` and VTK writes are synchronous);
rebuilding the sysimage more often than needed.

---

## 7. AMR on macOS (Apple Silicon): the patched GridapP4est fork

The registered `GridapP4est 0.3.11` does not support AMR on macOS:
`@cfunction` closures aren't supported on ARM64, and there is a Julia/C struct
stride mismatch for the p4est iterator structs (ARM64 with Julia ≥ 1.11, x86_64
with Julia ≥ 1.12). Both surface in the refinement callbacks that run *after*
the coarse octree model is built, so every `:lamr`/`:lpreadapt`/
`:linitial_refine` case (see [docs/amr_setup.md](docs/amr_setup.md)) needs the
patched fork on macOS:
<https://github.com/Hwang1229/GridapP4est.jl/tree/arm64-cfunction-fix>.

**You do not install the fork by hand.** `Project.toml` pins it in a
`[sources]` block, so `Pkg.instantiate()` resolves `GridapP4est` to the fork on
every platform. Verify:

```bash
julia --project=. -e 'using Pkg; Pkg.status("GridapP4est")'
# expect: https://github.com/Hwang1229/GridapP4est.jl#arm64-cfunction-fix
```

To work from an editable local clone instead (overrides the pin on your machine
only):

```bash
git clone -b arm64-cfunction-fix git@github.com:Hwang1229/GridapP4est.jl.git ~/GridapP4est.jl
julia --project=. -e 'using Pkg; Pkg.develop(path=expanduser("~/GridapP4est.jl"))'
```

> The `AssertionError: A check failed` at
> `OctreeDistributedDiscreteModels.jl:325` is a *different*, platform-independent
> issue (a `Dp != Dc` coarse mesh), handled in code by
> `_flatten_model_to_cell_dim` in `src/kernel/mesh/mesh.jl` — see the
> [FAQ entry](FAQ.md#amr-theta_amr-and-other-lamrlinitial_refine-cases-fails-with-assertionerror-a-check-failed-in-octreedistributeddiscretemodel).
>
> `dlopen` failures for `libp4est`/`libjansson`, or a `SystemError: opening
> file ".../artifacts/..."` during precompile, are an unrelated
> artifact-corruption issue — see
> [FAQ: missing file inside a Julia artifact](FAQ.md#a-package-fails-to-precompile-with-a-missing-file-inside-a-julia-artifact).

### 7.1 AMR quick-start on macOS (clean build → `theta_amr`)

```bash
# 0. Prereqs
julia --version                              # 1.11.9

# 1. Clone side by side
git clone git@github.com:smarras79/Jexpresso.git
git clone git@github.com:smarras79/JexpressoMeshes.git
cd Jexpresso
git checkout sm/newmaster

# 2. Link the sample meshes (provides hexa_TFI_10x10.msh)
ln -s ../JexpressoMeshes/meshes .
ls -l meshes/gmsh_grids/hexa_TFI_10x10.msh   # must resolve

# 3. Instantiate (clones the fork via [sources]; needs network). No precompile yet.
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'

# 4. Confirm the fork took effect
julia --project=. -e 'using Pkg; Pkg.status("GridapP4est")'
#    expect: ...GridapP4est.jl#arm64-cfunction-fix

# 5. Bind MPI — Route C (bundled MPICH_jll) is the reliable choice on macOS
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'

# 6. macOS hostname fix (§5.5)
echo "127.0.0.1   $(hostname -s)" | sudo tee -a /etc/hosts
echo "127.0.0.1   $(hostname)"    | sudo tee -a /etc/hosts

# 7. Rebuild MPI AND the MPI-linked C shim (§5.3), then precompile
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.build("P4est_wrapper"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.build("GridapP4est"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'

# 8. Sanity-check a non-AMR case first
julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler","theta")'

# 9. Run the AMR case
julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler","theta_amr")'
```

**Updating an existing clone** (e.g. you hand-added the fork earlier): force a
clean re-resolve so the `[sources]` pin and new mesh code are picked up.

```bash
cd Jexpresso
git pull
rm -f Manifest.toml            # gitignored; forces a from-scratch resolve
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true); Pkg.build("P4est_wrapper"; verbose=true); Pkg.build("GridapP4est"; verbose=true); Pkg.precompile()'
```

**If `theta_amr` still misbehaves**, check in this order:

```bash
# a) Segfault in _platform_memmove right after "Done reading ... .msh"?
#    → P4est_wrapper is linked to a different MPI than MPI.jl loads. Redo step 7
#      and compare (§5.4):
otool -L $(find ~/.julia/scratchspaces -name 'libp4est_wrapper*.dylib') | grep mpi
julia --project=. -e 'using MPI; println(MPI.API.libmpi)'

# b) Is the fork really in use? (step 4) — must show #arm64-cfunction-fix

# c) Embedding dimension of the mesh (Dp=3 is flattened at runtime by _flatten_model_to_cell_dim)
julia --project=. -e 'using GridapGmsh; println(typeof(GmshDiscreteModel("./meshes/gmsh_grids/hexa_TFI_10x10.msh")))'
#    {2,3,...} = Dp=3 (flattened);  {2,2,...} = already flat
```

Other macOS stumbles (`dlopen`/`libjansson`, `gethostbyname failed`, hangs at
`MPI.Init`) are in [FAQ.md](FAQ.md#run).

---

## Adding your own test cases

See [ADD_A_NEW_TEST.md](ADD_A_NEW_TEST.md).

## Notes on package versions

Jexpresso pins several packages whose latest versions are incompatible. This
list is an informational snapshot; the authoritative versions are in
[`Project.toml`](Project.toml). `GridapP4est` in particular is not plain
`=0.3.11` but the patched fork via `[sources]` (§7).

```
[compat]
BenchmarkTools = "1.8.0"
CSV = "0.10.16"
Crayons = "=4.1.1"
Gridap = "=0.18.12"
GridapDistributed = "=0.4.7"
GridapGmsh = "=0.7.2"
GridapP4est = "=0.3.11"
JACC = "1.0.0"
JLD2 = "0.5.15"
KrylovPreconditioners = "0.3.5"
LinearOperators = "2.11.0"
MPI = "=0.20.22"
MPIPreferences = "=0.1.11"
ONNXRunTime = "1.3.1"
PProf = "3.2.0"
Preferences = "1.5.2"
PrettyTables = "=2.4.0"
Profile = "1.11.0"
QuadGK = "2.11.2"
Roots = "2.2.13"
SciMLBase = "2.148.0"
Serialization = "1.11.0"
Thermodynamics = "=0.12.7"
TimerOutputs = "0.5.29"
TrixiBase = "0.1.8"
UUIDs = "1.11.0"
UnicodePlots = "=3.7.2"
```
