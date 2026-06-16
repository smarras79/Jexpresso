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
MPI implementation itself — it binds to one at build time. You pick exactly
**one** of the three routes below and tell `MPIPreferences` which one to use.

| Route | What provides MPI | When to choose it |
|-------|-------------------|-------------------|
| **A. OpenMPI** (system binary) | An OpenMPI you install on the machine | Linux clusters / HPC where OpenMPI is the site default |
| **B. MPICH** (system binary) | An MPICH you install on the machine | You prefer MPICH, or your cluster ships MPICH |
| **C. MPICH_jll** (native, bundled) | MPI shipped *inside* Julia's package env — nothing to install | Laptops/desktops, no admin rights, or OpenMPI deadlocks on macOS (see [Troubleshooting](#56-troubleshooting)) |

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
Proceed straight to the configuration step below.

### 5.2 Point `MPIPreferences` at your MPI

Run **one** of these, matching the route you chose. The setting is recorded in
`LocalPreferences.toml` in the project root.

#### Route A or B — system binary (OpenMPI / MPICH)

If MPI is installed in standard system paths (`/usr/bin`, `/usr/local/bin`):

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary()'
```

If MPI lives in a non-standard location (multiple installs, `/opt/...`,
Homebrew on Apple Silicon), pass its `lib` directory explicitly:

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/PATH/TO/MPILIB/lib"])'
```

> For Homebrew-installed OpenMPI/MPICH on Apple Silicon the path is usually
> `/opt/homebrew/lib`:
>
> ```bash
> julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/opt/homebrew/lib"])'
> ```

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
  using MPIPreferences; println("binary = ", MPIPreferences.binary)
  using MPI;            println(MPI.identify_implementation())'
```

- **Route A** should report an Open MPI implementation and `binary = "system"`.
- **Route B** should report MPICH and `binary = "system"`.
- **Route C** should print `binary = "MPItrampoline_jll"` (MPICH-based).

### 5.5 macOS hostname fix (MPICH and MPICH_jll only)

Required for **any MPICH-based MPI on macOS** — i.e. Route B on macOS and
Route C on macOS. Not needed for OpenMPI (Route A).

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
win for parallel cold starts is `PackageCompiler.create_sysimage` — a
larger one-time investment that ships a pre-compiled `.dylib` and
makes every cold `mpiexec` rank skip JIT too — but that is out of
scope for this guide.

# To run other tests that are already in Jexpresso or to add your own new problem,
see [ADD_A_NEW_TEST.md](ADD_A_NEW_TEST.md)

# NOTES ON PACKAGE LIST:


Jexpresso uses a few packages whose latest version may be incompatible. Please, enfornce the installation of the following versions:

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
