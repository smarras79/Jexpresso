# Building and Using a Julia Sysimage for Jexpresso

A sysimage bakes precompiled code into a `.so` file, eliminating Julia's startup latency on subsequent runs.

## Prerequisites

Before building the sysimage, comment out the following lines in [src/Jexpresso.jl](../src/Jexpresso.jl) (they are incompatible with sysimage builds):

```julia
# using Revise
# using SnoopCompile
```

## Building the Sysimage

Start a fresh Julia REPL from the project root:

```bash
julia --project=.
```

Then inside the REPL, run:

```julia
julia> include("create_Jexpresso_sysimage.jl")
```

This script:
1. Runs `precompile_jexpresso.jl` to execute a real Jexpresso run, loading all packages into the session
2. Detects all loaded project packages via `Base.loaded_modules`
3. Writes the package list to `sysimage_packages.txt`
4. Calls `PackageCompiler.create_sysimage()` to produce `jexpresso.so`

## When to Rebuild the Sysimage

**Rebuild `jexpresso.so` whenever you run `Pkg.instantiate()`** (i.e., after adding, removing, or updating packages). The sysimage encodes the exact compiled state of the packages at build time, so any change to the package environment will make the existing sysimage stale.

## Running Jexpresso with the Sysimage

Pass `--sysimage jexpresso.so` to the Julia invocation. For MPI runs:

```bash
mpirun -np <N> julia --project=. --sysimage jexpresso.so -e \
  'push!(empty!(ARGS), "CompEuler", "3d"); include("src/Jexpresso.jl")'
```

## Relevant Files

| File | Purpose |
|------|---------|
| `create_Jexpresso_sysimage.jl` | Main script to build the sysimage |
| `precompile_jexpresso.jl` | Precompile trace -- runs `CompEuler/3d` to warm up all code paths |
| `sysimage_packages.txt` | Auto-generated list of packages baked into the sysimage |
| `jexpresso.so` | The compiled sysimage output |
