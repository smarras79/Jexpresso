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
your MPI library:

```bash
julia --project=. -e 'using MPIPreferences; MPIPreferences.use_system_binary(extra_paths=["/PATH/TO/MPILIB/lib"])'
```

> For example, if you use OpenMPI installed with Homebrew, the path is likely
> `/opt/homebrew/lib`.

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