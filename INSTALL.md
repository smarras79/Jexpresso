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

**3c. Precompile everything:**

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
