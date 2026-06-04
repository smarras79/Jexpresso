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

Then instantiate the project and precompile its dependencies:

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```

This step may take a while the first time as Julia downloads and compiles all
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
