# Installing Jexpresso

This guide walks you through downloading, building, and testing **Jexpresso**.

> **Hit an error?** See [FAQ.md](FAQ.md) for fixes to common installation and run problems.

## Prerequisites

- **Julia 1.11.9** is the recommended version.
  > Julia 1.12.6 also works, but for now we prefer to stay on the earlier release.
- **Git** with SSH access to GitHub configured ([setup guide](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)).

## 1. Download the repositories

Clone Jexpresso and the companion mesh repository:

```julia
git clone git@github.com:smarras79/Jexpresso.git
git clone git@github.com:smarras79/JexpressoMeshes.git
```

> **Note:** `JexpressoMeshes` contains sample meshes used to run the existing
> tests without having to build them from scratch.

## 2. Link the sample meshes

From inside the `Jexpresso` directory, create a symbolic link to the meshes:

```bash
cd Jexpresso
ln -s ../JexpressoMeshes/meshes .
```

## 3. Build and precompile (serial)

These steps give you a working **serial** Jexpresso. If you want to run in
parallel, do these steps first and then continue with
[Section 5 — Running in parallel with MPI](#5-running-in-parallel-with-mpi).

If you are not already inside the project directory, move into it first:

```bash
cd PATH/TO/Jexpresso
```

**3a. Instantiate the dependencies** (with automatic precompilation disabled so
it can be controlled explicitly below):

```bash
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
```

**3b. Precompile everything:**

```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

This step may take a while the first time as Julia compiles all dependencies.

> **Going parallel?** Do **not** precompile yet if you already know you will
> run with MPI — configuring MPI first (Section 5) rebuilds the `MPI` package,
> which triggers a recompile anyway. Either order works, but configuring MPI
> first saves you one precompilation pass.

## 4. Test the installation

Once compilation has finished, verify everything works by running one of the
bundled cases.

Start a Julia session scoped to the project:

```bash
julia --project=.
```

Then, at the Julia prompt, run the `CompEuler` / `sod1d` test:

```julia
using Jexpresso
Jexpresso.run_case("CompEuler", "sod1d")
```

If the test runs to completion, your serial Jexpresso installation is ready. 🎉

---

## 5. Running in parallel with MPI

Everything you need to install, configure, and run Jexpresso in parallel lives
in this section. **Skip it entirely if you only run serially.**

Julia's [`MPI.jl`](https://juliaparallel.org/MPI.jl/stable/) does not contain an
MPI implementation itself — it binds to one at build time. You tell
`MPIPreferences` which of the three routes below to use.

You may have **several MPIs installed at once** — MPI.jl is bound to exactly
one of them *at a time*, and you can rebind whenever you like (§5.2, case 2).
What is never allowed is mixing them *within a single run*: the library MPI.jl
loads, the `mpiexec` that launches the job, and — for coupled runs — the MPI
`Alya.x` was linked against must all be the same installation.

| Route | What provides MPI | When to choose it |
|-------|-------------------|-------------------|
| **A. OpenMPI** (system binary) | An OpenMPI you install on the machine | Linux clusters / HPC where OpenMPI is the site default |
| **B. MPICH** (system binary) | An MPICH you install on the machine | You prefer MPICH, or your cluster ships MPICH |
| **C. MPICH_jll** (native, bundled) | MPI shipped *inside* Julia's package env — nothing to install | Laptops/desktops, no admin rights, or OpenMPI deadlocks on macOS (see [Troubleshooting](#57-troubleshooting)) |

The three routes share the same workflow: **install MPI → point `MPIPreferences`
at it → rebuild `MPI` → launch with the matching `mpiexec`.** Only the details
differ, and they are spelled out per route below.

### 5.1 Install an MPI implementation

#### Route A — OpenMPI (system binary)

```bash
# Ubuntu/Debian
sudo apt install libopenmpi-dev openmpi-bin

# macOS (Homebrew)
brew install open-mpi

# Verify
mpiexec --version
```

#### Route B — MPICH (system binary)

```bash
# Ubuntu/Debian
sudo apt install mpich libmpich-dev

# macOS (Homebrew)
brew install mpich

# Verify
mpiexec --version
```

#### Route C — MPICH_jll (native, bundled with MPI.jl)

**Nothing to install.** MPI ships *with* Julia's package environment as a JLL
(MPItrampoline, which defaults to MPICH). This is the most reliable route on a
laptop and the recommended fallback when a system OpenMPI deadlocks on macOS.
Proceed straight to the configuration step below. It is unaffected by whatever
system MPIs you have — it brings its own library *and* its own `mpiexec`
(`julia --project=. -e 'using MPI; println(mpiexec())'`), which is the one you
must launch with.

#### Installing both OpenMPI and MPICH side by side

You do not have to choose once and for all — keeping both installed is
supported, and useful (OpenMPI 5 on macOS is prone to hanging in `MPI_Init`,
so having MPICH ready to switch to saves a reinstall). Just run both install
commands above. The only rule is that **each individual run must use one
implementation end to end**: library, launcher, and — for coupled runs —
`Alya.x`. Section 5.2, case 2 below shows how to select between them.

### 5.2 Point `MPIPreferences` at your MPI

Run **one** of these, matching the route you chose. The setting is recorded in
`LocalPreferences.toml` in the project root.

> **Two things get selected, not one.** `MPIPreferences` records which
> **library** (`libmpi`) Julia loads. Your shell separately decides which
> **launcher** (`mpiexec`/`mpirun`) starts the job. They are independent, and a
> job whose launcher and library come from different MPIs does not fail with an
> error — it hangs in `MPI_Init`, or gives every rank its own world of size 1.
> Whenever you bind one, pin the other to match.

#### Route A or B, case 1 — only ONE system MPI installed

If MPI is installed in standard system paths (`/usr/bin`, `/usr/local/bin`):

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
```

If it lives somewhere else (`/opt/...`, Homebrew on Apple Silicon), pass its
`lib` directory:

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/PATH/TO/MPI/lib"])'
```

> For a single Homebrew MPI on Apple Silicon that directory is `/opt/homebrew/lib`.
> **Do not use `/opt/homebrew/lib` if you have both OpenMPI and MPICH installed** —
> see case 2.

#### Route A or B, case 2 — BOTH OpenMPI and MPICH installed

This is fine and common; the two just have to be kept apart. The bare call
`use_system_binary()` searches the default loader paths and binds whichever
`libmpi` it finds first — which is not necessarily the one whose `mpiexec` is
first on your `PATH`. Name the install explicitly instead.

**Step 1 — find each installation's prefix.**

On macOS/Homebrew (these print stable symlinks into the versioned Cellar dirs,
e.g. `/opt/homebrew/opt/open-mpi`):

```bash
brew --prefix open-mpi
brew --prefix mpich
```

On Ubuntu/Debian, both packages coexist under per-implementation subdirs:

```bash
ls -d /usr/lib/*/openmpi /usr/lib/*/mpich 2>/dev/null
```

Anywhere else, find the libraries themselves:

```bash
ls -d /opt/*/lib/libmpi.* /usr/local/*/lib/libmpi.* 2>/dev/null
```

Typical layouts:

| Platform | OpenMPI prefix | MPICH prefix |
|---|---|---|
| macOS, Homebrew (Apple Silicon) | `/opt/homebrew/opt/open-mpi` | `/opt/homebrew/opt/mpich` |
| macOS, Homebrew (Intel) | `/usr/local/opt/open-mpi` | `/usr/local/opt/mpich` |
| Ubuntu/Debian | `/usr/lib/x86_64-linux-gnu/openmpi` | `/usr/lib/x86_64-linux-gnu/mpich` |
| HPC modules | `$MPI_HOME` after `module load openmpi` | `$MPI_HOME` after `module load mpich` |

**Step 2 — bind MPI.jl to the one you want, library *and* launcher.**

`extra_paths` is searched *before* the default loader paths, so it decides the
tie. Passing `mpiexec` at the same time records the matching launcher in
`LocalPreferences.toml`, which is what keeps the two halves from drifting apart:

> ⚠️ **Paste these one block at a time, and never add a trailing `# comment`
> to the `MPI_PREFIX=` line.** Interactive `zsh` (the macOS default) does not
> treat `#` as a comment unless you `setopt interactivecomments`. It parses
> `MPI_PREFIX=$(...)  # note` as the command `#` with `MPI_PREFIX` set only for
> *that* command — so you get `zsh: command not found: #` **and the variable is
> never set**. The `use_system_binary` call then runs with an empty prefix and
> fails with `extra_paths = ["/lib"]`. The blocks below are comment-free for
> exactly this reason.

**To use OpenMPI**, set the prefix and check it is non-empty:

```bash
OMPI_PREFIX=$(brew --prefix open-mpi)
echo "prefix = $OMPI_PREFIX"
ls "$OMPI_PREFIX"/lib/libmpi.*
```

Both commands must print something. Then bind:

```bash
julia --project=. -e "using MPIPreferences; MPIPreferences.use_system_binary(
        extra_paths = [\"$OMPI_PREFIX/lib\"],
        mpiexec     = \"$OMPI_PREFIX/bin/mpiexec\")"
```

**To use MPICH**, the same in its prefix:

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

On Linux, set the prefix by hand instead of with `brew --prefix` — e.g.
`OMPI_PREFIX=/usr/lib/x86_64-linux-gnu/openmpi` or
`MPICH_PREFIX=/usr/lib/x86_64-linux-gnu/mpich` — and run the identical `julia`
command.

> **`MPI library could not be found ... in the following extra directories:
> ["/lib"]`** means the prefix variable was empty when the `julia` line ran —
> the `extra_paths` entry is literally `"" * "/lib"`. Re-run the `echo` above:
> if it prints `prefix = ` with nothing after it, either you hit the `zsh`
> comment trap described above, or `brew --prefix <formula>` failed because
> that formula is not installed.

Then continue with [5.3](#53-rebuild-mpi-and-precompile) — the rebuild is
**mandatory** after every rebinding.

**Step 3 — put the matching launcher first on `PATH` for the same shell.**

`MPIPreferences` cannot change your `PATH`; you must. Define both switches
once in `~/.zshrc` / `~/.bashrc` and call the one you need before launching:

```bash
# macOS / Homebrew
use-openmpi() { export PATH="$(brew --prefix open-mpi)/bin:$PATH"; hash -r; }
use-mpich()   { export PATH="$(brew --prefix mpich)/bin:$PATH";    hash -r; }

# Linux — substitute the prefixes from step 1
use-openmpi() { export PATH="/usr/lib/x86_64-linux-gnu/openmpi/bin:$PATH"; hash -r; }
use-mpich()   { export PATH="/usr/lib/x86_64-linux-gnu/mpich/bin:$PATH";   hash -r; }
```

Confirm before every parallel run that all three agree:

```bash
which mpiexec mpirun mpif90
```

> **Homebrew caveat.** The `open-mpi` and `mpich` formulas both provide
> `mpicc`, `mpiexec`, `libmpi.dylib`, so Homebrew will only *link* one of them
> into `/opt/homebrew/{bin,lib}` at a time. Both stay fully usable under their
> own `opt`/`Cellar` prefixes — which is exactly why the commands above use
> `brew --prefix <formula>` and never the shared `/opt/homebrew/lib`. If you
> want to change which one owns the shared symlinks:
> ```bash
> brew unlink mpich && brew link --overwrite open-mpi
> ```
> Rebinding MPI.jl does **not** require this — the prefix paths work whether or
> not the formula is linked.

> **Ubuntu/Debian caveat.** `update-alternatives --config mpi` (and
> `mpi-default-bin`) switches which `mpicc`/`mpirun` the bare names resolve to,
> system-wide. It has no effect on what MPI.jl is bound to. Use the explicit
> prefixes above and treat the alternative purely as a `PATH` convenience.

#### Switching from one to the other later

Rebinding is not just the `MPIPreferences` call — the stale build and
precompile caches must go too, or you get the old library with the new
preference:

```bash
rm -f LocalPreferences.toml
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

Then verify with [5.4](#54-verify-the-binding), and check the whole setup end
to end with:

```bash
bash tools/check_mpi_setup.sh
```

> **Coupled runs need a third thing to match.** `AlyaProxy/Alya.x` is a
> separately compiled Fortran binary with its own MPI baked in at link time.
> After switching implementations you must rebuild it with the matching
> `mpif90` — `cd AlyaProxy && MPIF90=<prefix>/bin/mpif90 bash compilef90.sh`.
> See [RUN-COUPLED.md](RUN-COUPLED.md) §3 and §5.

#### Route C — MPICH_jll (native)

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'
```

### 5.3 Rebuild `MPI` and precompile

`MPI` must be rebuilt against whatever you just selected, then everything
precompiled:

```bash
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

### 5.4 Verify the binding

```bash
julia --project=. -e '
  using MPIPreferences; println("binary  = ", MPIPreferences.binary)
  using MPI;            println("impl    = ", MPI.identify_implementation())
                        println("libmpi  = ", MPI.API.libmpi)
                        println("mpiexec = ", MPI.mpiexec())'
```

- **Route A** should report an Open MPI implementation and `binary = "system"`.
- **Route B** should report MPICH and `binary = "system"`.
- **Route C** should print `binary = "MPItrampoline_jll"` (MPICH-based).

With **both** OpenMPI and MPICH installed, the implementation name alone is not
enough — check the two paths as well:

- `libmpi` must sit under the prefix you chose in 5.2 (e.g.
  `/opt/homebrew/opt/mpich/lib/...`, not `/opt/homebrew/lib/...`, which is
  ambiguous because it belongs to whichever formula is currently linked).
- `mpiexec` must come from that **same** prefix, and must match what
  `which mpiexec` reports in the shell you launch from. If those disagree, your
  runs will hang in `MPI_Init` rather than report an error.

`bash tools/check_mpi_setup.sh` checks all of this for you, plus the launcher
itself, and times out instead of hanging.

### 5.5 macOS hostname fix (MPICH and MPICH_jll only)

Required for **any MPICH-based MPI on macOS** — i.e. Route B on macOS and
Route C on macOS.

> It is worth applying on **OpenMPI too**. Open MPI 5's PRRTE launcher also
> resolves the local hostname when it brings up its out-of-band channel, and an
> unresolvable hostname makes `mpirun` hang before your program starts rather
> than print the MPICH error below. The fix is harmless either way.

MPICH's TCP channel resolves the machine hostname via the C `gethostbyname()`
call, which on macOS only returns mDNS names (`*.local`) and fails on the bare
hostname, producing `MPI_Init` errors like:

```
GetSockInterfaceAddr ... gethostbyname failed, <your-hostname> (errno 0)
```

Fix once, permanently:

```bash
echo "127.0.0.1   $(hostname -s)" | sudo tee -a /etc/hosts
echo "127.0.0.1   $(hostname)"    | sudo tee -a /etc/hosts
```

Verify:

```bash
ping -c 1 $(hostname -s)     # should respond from 127.0.0.1
```

If you cannot `sudo`, export `MPICH_INTERFACE_HOSTNAME=127.0.0.1` in every shell
session before launching MPI jobs (or add it to `~/.zshrc`).

### 5.6 Launch a parallel run

> **Use the launcher that matches your route.** A system `mpiexec`/`mpirun`
> belongs to the system MPI; the JLL provides its own `mpiexec` through MPI.jl.
> Mixing them is the most common cause of "it won't start" failures.

#### Route A or B — system MPI (OpenMPI / MPICH)

Use the system launcher directly:

```bash
mpiexec -n <NPROCS> julia --project=. src/Jexpresso.jl <EQUATIONS> <CASE_NAME>
```

For example, 4 ranks of the 3D Euler case:

```bash
mpiexec -n 4 julia --project=. src/Jexpresso.jl CompEuler 3d
```

If `mpiexec` is not on your `PATH`, or you have several MPIs installed, use
absolute paths to both the launcher and `julia`:

```bash
/opt/homebrew/Cellar/open-mpi/5.0.6/bin/mpirun -n 4 \
  /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia \
  --project=. src/Jexpresso.jl CompEuler theta
```

#### Route C — MPICH_jll (native)

Do **not** use the system `mpirun`. Launch with the `mpiexec` that MPI.jl ships,
which you reach from inside Julia:

```bash
julia --project=. -e '
  using MPI
  run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler city2d`)'
```

For daily use, wrap the launcher in a small shell script (a ready-made copy
ships as [`jexp_mpich.sh`](jexp_mpich.sh) in the repo root):

```bash
#!/bin/bash

# USER: change to your julia path: ---------------------------------------
JULIA=/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia
# END USER ---------------------------------------------------------------

jexp_mpich() {
    $JULIA --project=. -e "
      using MPI
      run(\`\$(mpiexec()) -n $1 \$(Base.julia_cmd()) --project=. src/Jexpresso.jl $2 $3\`)"
}

jexp_mpich "$@"
```

Usage:

```bash
./jexp_mpich.sh 4 CompEuler city2d
```

### 5.7 Troubleshooting

- **OpenMPI 5 + macOS deadlock on `MPI.Init` (Route A).** A known sharp edge
  where the PMIx handshake never completes and the run hangs forever. The fix is
  to switch to **Route C (MPICH_jll)** — redo Sections 5.2–5.6 with the JLL
  binary. This needs no system MPI at all.
- **MPICH aborts in `MPI_Finalize` with an `OFI` / `nic=utunN` error.** For
  example:
  ```
  Fatal error in internal_Finalize: Other MPI error
  MPIDI_OFI_mpi_finalize_hook(861): flush_send(771):
  OFI call tsenddata failed (default nic=utun7: No such file or directory
  ```
  MPICH's OFI (libfabric) netmod enumerates the machine's network interfaces
  and picked a **VPN / tunnel device** (`utun0`, `utun7`, …) as its NIC. That
  interface disappears or refuses the send, and the job dies — often at
  finalize, after the science already ran, and on a laptop where no real
  network fabric is needed at all.

  Pin the provider, and pin it to the machine's **real** default interface —
  which you can read off the routing table rather than guess:
  ```bash
  route get default | awk '/interface:/{print $2}'
  ```
  Then, with that name (`en0` on most Macs):
  ```bash
  export FI_PROVIDER=tcp
  export FI_TCP_IFACE=en0
  ```
  Do not reach for `lo0` first. Loopback looks like the obvious choice for a
  single-machine run, but libfabric's `tcp` provider often cannot open an
  endpoint on it, which trades the original error for this one:
  ```
  MPIDI_OFI_create_vci_context(1054): OFI call ep_enable failed
  (default nic=tcp: Bad file descriptor)
  ```
  If neither works, the older `sockets` provider is the reliable fallback and
  needs no interface at all:
  ```bash
  export FI_PROVIDER=sockets
  ```
  `bash tools/check_mpi_setup.sh` tries all of these combinations and tells you
  which one works on your machine. Disconnecting the VPN also makes the
  original symptom go away, which is a quick way to confirm the diagnosis.

  **Making it permanent.** Put the winning pair — and only that pair, since
  they interact — in `~/.zshrc`. Two things to know first:

  * `FI_*` is read by *any* libfabric user, so this affects every MPI job on
    the machine, OpenMPI's OFI components included. On a laptop that is
    normally what you want.
  * `FI_TCP_IFACE` names an interface that has to exist at `MPI_Init`. On a
    laptop `en0` is usually Wi-Fi, so turning Wi-Fi off, or a dock/VPN change
    that renumbers interfaces, can bring the failure back. `FI_PROVIDER=sockets`
    with no interface pinned is the more robust choice if that bites.

  If you sync dotfiles to a cluster, guard it — these settings are wrong on a
  machine with a real fabric (InfiniBand, Slingshot), where libfabric should
  pick `verbs`, `psm2` or `cxi` for itself:
  ```bash
  if [[ "$(hostname -s)" == "my-laptop" ]]; then
    export FI_PROVIDER=tcp
    export FI_TCP_IFACE=en0
  fi
  ```

  > **Do not use `FI_PROVIDER=shm` on macOS.** The libfabric `shm` provider is
  > **Linux-only**. Asking for it on macOS leaves libfabric with no provider to
  > return, and `MPI_Init` then fails outright with a *different* error that is
  > easy to mistake for a broken install:
  > ```
  > MPIDI_OFI_find_provider(84): find_provider(126):
  > OFI call getinfo failed (default nic=(n/a): No message available on STREAM)
  > ```
  > `nic=(n/a)` means "no provider matched your filter". If you see it, check
  > whether `FI_PROVIDER` is set to something unavailable — `echo $FI_PROVIDER`
  > — and switch to `tcp` or `sockets`, or `unset FI_PROVIDER` to let
  > libfabric choose again. `fi_info -l` lists the providers your build
  > actually has.
- **Library conflicts / stale binding.** Remove the recorded preference and
  reconfigure from Section 5.2:
  ```bash
  rm -f LocalPreferences.toml
  ```
- **Path issues (system MPI).** Confirm which launcher you are actually calling:
  ```bash
  which mpiexec
  which mpirun
  ```
  Use absolute paths (see Section 5.6) if the wrong one is picked up.
- **Version mismatches (system MPI).** Make sure the compiler wrappers and the
  runtime agree:
  ```bash
  mpicc --version
  mpif90 --version
  ```

### 5.8 Running on a Slurm cluster

For batch runs on many nodes, use the ready-made job script:

```bash
sbatch --export=ALL,EQS=CompEuler,CASE=theta \
       auxiliary/slurm/submit_jexpresso_parallel.sh
```

Edit the `CHANGE_ME` lines in its `#SBATCH` block (account, partition, QoS)
and the `module load` lines at the top; the rest is site-agnostic.

The script splits the job into a **serial precompile phase** — one process on
one node running `Pkg.instantiate()` + `Pkg.precompile()` — and a **parallel
run phase** launched with precompilation disabled
(`JULIA_PKG_PRECOMPILE_AUTO=0`, `--compiled-modules=existing`,
`--pkgimages=existing`). That split is the whole point: with a cold depot,
every rank otherwise discovers the same missing cache files at the same
moment and they all compile into the same shared depot at once, which on a
few hundred ranks costs far more wall time than the idle nodes during the
serial phase.

See [`auxiliary/slurm/README.md`](auxiliary/slurm/README.md) for the full set
of knobs and the cluster-specific pitfalls (system-MPI binding,
`JULIA_CPU_TARGET` on heterogeneous partitions, why you must *not* stage the
depot to node-local disk).

## 6. Daily workflow — interactive REPL for fast iteration

Every Julia process pays a one-time JIT compilation cost on first use of
`sem_setup`, the `with_mpi` closure, the SciML integrator, the VTK
writer, etc. Cold starts (a fresh `mpirun`, a fresh `julia src/Jexpresso.jl ...`)
re-pay this cost every time — typically ~30–60 s of silent wall time
between `# Read inputs dict ... DONE` and the time loop visibly
advancing.

**The cheapest way to escape that is the REPL workflow**: launch Julia
once, run the same case (or different cases) repeatedly inside the same
process. The first invocation in a session is slow; every subsequent
invocation is essentially instant for the JIT-related work — only your
actual integration time remains.

### Single-rank (serial) workflow

This is the recommended development workflow:

```bash
julia --project=.
```

Then at the prompt, run a case:

```julia
julia> using Jexpresso
julia> Jexpresso.run_case("CompEuler", "theta")
# ... lots of JIT on the first run; the simulation completes ...
```

Edit your code, then re-run in the SAME session:

```julia
julia> Jexpresso.run_case("CompEuler", "theta")
# ... starts almost immediately; only your edits get JIT-compiled ...
```

For source files outside `user_inputs.jl` etc., use `Revise.jl` so edits
are picked up without restarting the REPL.

### MPI runs from the same Julia process

You can also drive parallel `mpiexec` runs from within an interactive
Julia session — Julia's `run(...)` keeps `mpiexec` happy (this works for
all three routes; with Route C be sure `mpiexec()` comes from `using MPI`)
and you can re-run as many times as you like without restart:

```julia
julia> using MPI
julia> JULIA = Base.julia_cmd()
julia> run(`$(mpiexec()) -n 4 $JULIA --project=. src/Jexpresso.jl CompEuler city2d`)
```

Each `run(...)` invocation still spawns FRESH MPI ranks, so each
parallel run pays the cold-JIT cost again. The REPL-resident Julia
process does not save you here — the JIT cost lives in the child
processes that mpiexec spawns.

In short: **interactive REPL eliminates cold-start cost for serial
development**, but every `mpirun` is its own cold start. The next-tier
win for parallel cold starts is `PackageCompiler.create_sysimage`,
covered in [Section 6.1](#61-fastest-parallel-execution-on-macos-with-mpich).

### 6.1 Fastest parallel execution on macOS with MPICH

On a Mac the wall-clock cost of a parallel Jexpresso run is dominated by Julia
startup, not by MPI. Address that first; the MPI-side tuning is small by
comparison and mostly consists of not getting in your own way.

#### One-time: shell environment

Put these in `~/.zshrc` (see [§5.7](#57-troubleshooting) for why the libfabric
pair is needed, and how to find your own interface name):

```bash
export FI_PROVIDER=tcp
export FI_TCP_IFACE=en0
export PATH="$(brew --prefix mpich)/bin:$PATH"
```

#### One-time: build a sysimage

This is the single biggest win. Without it every rank of every `mpirun` JIT-compiles
Jexpresso from scratch — tens of seconds per launch, paid again on the next run.
A sysimage bakes that into a `.so` the ranks load directly:

```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
julia --project=. create_Jexpresso_sysimage.jl
```

That writes `jexpresso.so` (several minutes, once). Use it with `--sysimage`:

```bash
mpiexec -np 4 julia --project=. --sysimage jexpresso.so --startup-file=no \
    ./src/Jexpresso.jl CompEuler theta
```

Rebuild it after changing `Project.toml` or `src/`. `run_coupled.sh` picks up
`jexpresso.so` automatically when it is present, and `REBUILD_SYSIMAGE=1
./run_coupled.sh` regenerates it.

The second-biggest win is already automatic: the mesh and SEM caches under
`<case>/.jexpresso_cache/`. The first run of a case builds them, later runs of
the same case load them. Leave `:luse_mesh_cache` at its default.

#### How many ranks

Apple Silicon has performance and efficiency cores; only the P-cores are worth
giving MPI ranks. Count them:

```bash
sysctl -n hw.perflevel0.physicalcpu
sysctl -n hw.perflevel1.physicalcpu
```

The first number is your rank ceiling — 4 on an M1/M2 Air, 8–12 on a Pro/Max.
Oversubscribing past it makes runs slower, not faster, and MPICH will busy-wait
while it happens.

Jexpresso is essentially pure MPI (one `@threads` loop, in matrix assembly), so
give the ranks the cores and keep threads out of the way:

```bash
export JULIA_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
```

Without this, *N* ranks each start *N* BLAS threads and fight over the same
cores.

> **Rank pinning does not work on macOS.** Hydra's `-bind-to core` /
> `-map-by` rely on affinity APIs Darwin does not expose. Do not bother
> passing them; the scheduler places threads for you.

#### Launch flags worth using

```bash
mpiexec -np 4 -prepend-rank \
    julia --project=. --sysimage jexpresso.so --startup-file=no \
    ./src/Jexpresso.jl CompEuler theta
```

* `-prepend-rank` — Hydra's per-rank output tagging (MPICH's equivalent of
  OpenMPI's `--output tag`). Without it, a parallel run's stdout is
  interleaved and block-buffered, and looks like a hang.
* `--startup-file=no` — skips `~/.julia/config/startup.jl` on every rank.
* `-genv VAR VAL` — set an environment variable for all ranks. Hydra forwards
  your whole environment by default, so exported variables already arrive;
  this is for one-off overrides.

#### For coupled runs

Split the P-cores between the two codes rather than exceeding them — on a
4-P-core machine that means 2 + 2. **Alya needs at least two ranks**: its rank
0 is a master that owns no grid points, so a single Alya rank assigns zero
points and the codes do not couple at all. Beyond that minimum, Alya does very
little work, so spend extra ranks on Jexpresso:

```bash
export JEXPRESSO_COUPLED=1
mpiexec -np 2 ./AlyaProxy/Alya.x \
      : -np 2 julia --project=. --sysimage jexpresso.so --startup-file=no \
        ./src/Jexpresso.jl CompEuler 3dAlya
```

or just `./run_coupled.sh 2 2`, which applies the sysimage and the preflight
checks for you.

#### Things that cost you time

* **A VPN connected during a run.** It renumbers interfaces, which is what
  `FI_TCP_IFACE` is pinned against; a reconnect mid-run can kill the job.
* **The repo on iCloud Drive or a network share.** Every rank reads the mesh
  and cache files. Keep the working copy on the internal SSD.
* **Frequent output.** `:diagnostics_at_times` and the VTK writes are
  synchronous. Widen the interval while you are benchmarking.
* **Rebuilding the sysimage more often than needed.** It is only invalidated
  by dependency or source changes, not by editing `user_inputs.jl`.

## 7. AMR on macOS (Apple Silicon): the patched GridapP4est fork

The registered `GridapP4est` (`0.3.11`, see the package list at the bottom of
this file) does not support adaptive mesh refinement on macOS — `@cfunction`
closures aren't supported on ARM64, and there's a Julia/C struct stride
mismatch for the p4est iterator structs on both ARM64 (Julia ≥ 1.11) and
x86_64 (Julia ≥ 1.12). These surface during the refinement callbacks that run
*after* the coarse octree model is built, so every
`:lamr`/`:lpreadapt`/`:linitial_refine` case (see
[docs/amr_setup.md](docs/amr_setup.md)) needs the patched fork on macOS.

Both bugs are fixed on a branch of a fork:
<https://github.com/Hwang1229/GridapP4est.jl/tree/arm64-cfunction-fix>.

> **Note:** the `AssertionError: A check failed` at
> `OctreeDistributedDiscreteModels.jl:325` is a *different*, platform-independent
> issue — a `Dp != Dc` coarse mesh — and is handled in code by
> `_flatten_model_to_cell_dim` in `src/kernel/mesh/mesh.jl`; the fork does not
> address it. See the
> [FAQ entry](FAQ.md#amr-theta_amr-and-other-lamrlinitial_refine-cases-fails-with-assertionerror-a-check-failed-in-octreedistributeddiscretemodel).

**You no longer need to install this fork by hand.** `Project.toml` pins it in a
`[sources]` block (same package UUID and version `0.3.11`, plus the fixes), so a
plain `Pkg.instantiate()` already resolves `GridapP4est` to the fork on every
platform. Verify it took effect:

```bash
julia --project=. -e 'using Pkg; Pkg.status("GridapP4est")'
# should show "https://github.com/Hwang1229/GridapP4est.jl#arm64-cfunction-fix"
```

If you'd rather work from an editable local clone (e.g. to make further fixes
yourself), the `Pkg.develop` below overrides the `[sources]` pin for your
machine only:

```bash
git clone -b arm64-cfunction-fix git@github.com:Hwang1229/GridapP4est.jl.git ~/GridapP4est.jl
julia --project=. -e 'using Pkg; Pkg.develop(path=expanduser("~/GridapP4est.jl"))'
```

If you hit a missing-file error (e.g. `libp4est`/`libjansson` failing to
`dlopen`, or `Pkg.instantiate()` failing to precompile a package with a
`SystemError: opening file ".../artifacts/.../..."`) right after switching,
that's a separate, unrelated artifact-corruption issue — see the FAQ entry
["A package fails to precompile with a missing file inside a Julia artifact"](FAQ.md#a-package-fails-to-precompile-with-a-missing-file-inside-a-julia-artifact).

### 7.1 AMR quick-start on macOS (clean build → `theta_amr`)

A copy-pasteable path from nothing to a running AMR case on macOS. It reuses the
detailed steps above ([§2](#2-link-the-sample-meshes) for the mesh link,
[§5.2](#52-point-mpipreferences-at-your-mpi) / [§5.5](#55-macos-hostname-fix-mpich-and-mpich_jll-only)
for MPI); the notes here call out only what AMR adds.

```bash
# 0. Prereqs: Julia 1.11.9, git+SSH to GitHub.
julia --version                              # 1.11.9

# 1. Clone Jexpresso + the mesh repo side by side, get on this branch
git clone git@github.com:smarras79/Jexpresso.git
git clone git@github.com:smarras79/JexpressoMeshes.git
cd Jexpresso
git checkout claude/theta-amr-osx-assertion-bi8ubv

# 2. Link the sample meshes (provides hexa_TFI_10x10.msh — see Section 2)
ln -s ../JexpressoMeshes/meshes .
ls -l meshes/gmsh_grids/hexa_TFI_10x10.msh   # must resolve

# 3. Instantiate. The [sources] pin clones the patched GridapP4est fork here,
#    so this step needs network access. Precompile is deferred until after MPI.
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'

# 4. Confirm the fork took effect (this is what makes the refinement steps work)
julia --project=. -e 'using Pkg; Pkg.status("GridapP4est")'
#    expect: ...GridapP4est.jl#arm64-cfunction-fix

# 5. Configure MPI — Route C (bundled MPICH_jll) is the reliable choice on macOS
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'

# 6. macOS hostname fix (required for any MPICH-based MPI, incl. the JLL; Section 5.5)
echo "127.0.0.1   $(hostname -s)" | sudo tee -a /etc/hosts
echo "127.0.0.1   $(hostname)"    | sudo tee -a /etc/hosts

# 7. Rebuild MPI against that choice, then precompile everything
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
julia --project=. -e 'using Pkg; Pkg.precompile()'

# 8. Sanity-check a non-AMR case first (fast, proves the base build)
julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler","theta")'

# 9. Run the AMR case
julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler","theta_amr")'
```

**Updating an existing clone** (you previously hand-added the fork with
`Pkg.add`/`Pkg.develop`): force a clean re-resolve so the committed `[sources]`
pin and the new mesh code are picked up.

```bash
cd Jexpresso
git pull
rm -f Manifest.toml            # gitignored; forces a from-scratch resolve
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true); Pkg.precompile()'
```

**Two AMR-specific checks** if `theta_amr` still misbehaves:

```bash
# a) What embedding dimension does the mesh report? Dp=3 is why the old
#    AssertionError fired; _flatten_model_to_cell_dim now collapses it to Dp=Dc.
julia --project=. -e 'using GridapGmsh; println(typeof(GmshDiscreteModel("./meshes/gmsh_grids/hexa_TFI_10x10.msh")))'
#    {2,3,...} = Dp=3 (flattened at runtime);  {2,2,...} = already flat

# b) Is the fork really in use? (step 4) — must show #arm64-cfunction-fix
```

Common macOS stumbles are covered in the FAQ: `dlopen`/`libjansson`/missing
artifacts on instantiate, `MPI_Init … gethostbyname failed` (the hostname fix in
step 6), and hangs at `MPI.Init` (stay on Route C, step 5). See
[FAQ.md](FAQ.md#run).

# To run other tests that are already in Jexpresso or to add your own new problem,
see [ADD_A_NEW_TEST.md](ADD_A_NEW_TEST.md)

# NOTES ON PACKAGE LIST:


Jexpresso uses a few packages whose latest version may be incompatible. Please, enfornce the installation of the following versions:

> **Note:** this list is an informational snapshot — the authoritative versions
> live in [`Project.toml`](Project.toml). In particular `GridapP4est` is no
> longer plain `=0.3.11`: it is pinned to the patched fork via a `[sources]`
> block (see [Section 7](#7-amr-on-macos-apple-silicon-the-patched-gridapp4est-fork)),
> so a `Pkg.instantiate()` resolves the fork automatically.

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
</content>
</invoke>