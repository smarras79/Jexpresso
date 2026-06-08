# Installing Jexpresso

This guide walks you through downloading, building, and testing **Jexpresso**.

## Prerequisites

- **Julia 1.11.9** is the recommended version.
  > Julia 1.12.6 also works, but for now we prefer to stay on the earlier release.
- **Git** with SSH access to GitHub configured ([setup guide](https://docs.github.com/en/authentication/connecting-to-github-with-ssh)).

## 1. Download the repositories

Clone Jexpresso and the companion mesh repository:

```bash
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

## 3. Build and precompile

If you are not already inside the project directory, move into it first:

```bash
cd PATH/TO/Jexpresso
```

Then build the project in three steps.

**3a. Instantiate the dependencies** (with automatic precompilation disabled so
it can be controlled explicitly below):

```bash
julia --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'
```

**3b. Point MPI at your system binary**, replacing the path with the location of
your MPI library: You only need this step if you are planning to run Jexpresso in parallel and must 
have some version of MPI installed first:

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/PATH/TO/MPILIB/lib"])'
```

> For example, if you use OpenMPI installed with Homebrew, the path is likely
> `/opt/homebrew/lib`.

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/opt/homebrew/lib"])'
```

**3c. Alternative: use MPI.jl's bundled MPI (MPICH-based JLL).**

Use this route if step 3b deadlocks on `MPI.Init` (a known sharp edge
with Open MPI 5 + macOS + MPI.jl, where the PMIx handshake never
completes), or if you don't want to install a system MPI at all. With
this route MPI ships *with* Julia's package environment — no system
MPI is needed.

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_jll_binary()'
julia --project=. -e 'using Pkg; Pkg.build("MPI"; verbose=true)'
```

Verify the bind:

```bash
julia --project=. -e '
  using MPIPreferences; println("binary = ", MPIPreferences.binary)
  using MPI;            println(MPI.identify_implementation())'
```

`binary` should now print `"MPItrampoline_jll"` (which defaults to
MPICH on macOS).

> **Important — launch with the bundled `mpiexec`, NOT system `mpirun`.**
> After switching to the JLL binary, the system `mpirun` will not work
> because it belongs to a different MPI. Use the launcher MPI.jl ships:

```bash
julia --project=. -e '
  using MPI
  run(`$(mpiexec()) -n 4 $(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler city2d`)'
```

For daily use you can wrap the launcher in a small shell file `jexp_mpich.sh`:

```bash
#!/bin/bash

# USER: change to your julia path: ---------------------------------------
JULIA=/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia
# END USER ---------------------------------------------------------------

#!/bin/bash
jexp_mpich() {
    local msg="I am starting Julia; please be patient. Jexpresso hasn't started yet!"
    local border
    border=$(printf '═%.0s' $(seq 1 $(( ${#msg} + 2 ))))
    printf '\033[1;31m╔%s╗\n║ %s ║\n╚%s╝\033[0m\n' "$border" "$msg" "$border"

    # Launch Julia in the background so we can animate while it starts
    $JULIA --project=. -e "
      using MPI
      run(\`\$(mpiexec()) -n $1 \$(Base.julia_cmd()) --project=. src/Jexpresso.jl $2 $3\`)" &
    local pid=$!

    # Ctrl-C should kill Julia and restore the cursor
    trap 'kill "$pid" 2>/dev/null; tput cnorm 2>/dev/null; printf "\r\033[K"; trap - INT; return 130' INT

    local frames='⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'
    local n=${#frames}
    local i=0
    tput civis 2>/dev/null                  # hide cursor
    while kill -0 "$pid" 2>/dev/null; do
        printf '\r\033[1;31m%s waiting…\033[0m\033[K' "${frames:i++%n:1}"
        sleep 0.1
    done
    tput cnorm 2>/dev/null                   # restore cursor
    printf '\r\033[K'                        # erase the spinner line
    trap - INT
    wait "$pid"                              # propagate Julia's exit status
}
jexp_mpich "$@"
}
# Usage:
jexp_mpich 4 CompEuler city2d
```

**3d. macOS-specific: register your hostname in `/etc/hosts`.**

Required for *any* MPICH-based MPI on macOS (so: required if you took
step 3c, optional otherwise). MPICH's TCP channel resolves the machine
hostname via the C `gethostbyname()` call, which on macOS only returns
mDNS names (`*.local`) and fails on the bare hostname, producing
`MPI_Init` errors of the form:

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

If you cannot sudo, the equivalent env-var workaround is to export
`MPICH_INTERFACE_HOSTNAME=127.0.0.1` in every shell session before
launching MPI jobs (or add it to `~/.zshrc`).

**3e. Precompile everything:**

```bash
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

This last step may take a while the first time as Julia compiles all
dependencies.

## 4. Test the installation

Once compilation has finished, verify everything works by running one of the
bundled cases.

Start a Julia session scoped to the project:

```bash
julia --project=.
```

Then, at the Julia prompt, run the `CompEuler` / `sod1d` test:

```julia
push!(empty!(ARGS), "CompEuler", "sod1d");
include("./src/Jexpresso.jl")
```

If the test runs to completion, your Jexpresso installation is ready to go. 🎉

## 5. Daily workflow — interactive REPL for fast iteration

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
julia> push!(empty!(ARGS), "CompEuler", "theta");
julia> include("./src/Jexpresso.jl")
# ... lots of JIT on the first run; the simulation completes ...
```

Edit your code, then re-run in the SAME session:

```julia
julia> include("./src/Jexpresso.jl")
# ... starts almost immediately; only your edits get JIT-compiled ...
```

For source files outside `user_inputs.jl` etc., use `Revise.jl` so edits
are picked up without restarting the REPL.

### MPI runs from the same Julia process

You can also drive parallel `mpiexec` runs from within an interactive
Julia session — Julia's `run(...)` keeps the bundled `mpiexec` happy
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
