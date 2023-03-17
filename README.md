# JEXPRESSO
A research and educational software for the numerical solution of 1D, 2D, and 3D PDEs using spectral and spectral element methods on CPUs and GPUs. DISCLAIMER: this is WIP.

If you are interested in contributing, please get in touch.

# Some notes on using JEXPRESSO

To install and run the code assume Julia
version 1.7.2 or higher (tested up to 1.8.5)

The [MPI.jl][0] package that is used assumes that you have a working MPI installation

## Setup with CPUs

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```
followed by the following:

Push problem name to ARGS
You need to do this only when you run a new problem
```bash
julia> push!(empty!(ARGS), PROBLEM_NAME::String);
julia> include(./src/run.jl)
```

PROBLEM_NAME is the name of your problem directory as $JEXPRESSO/src/problems/problem_name
Ex. If you run the Advection Diffusion problem in $JEXPRESSO/src/problems/AdvDiff
```bash
julia> push!(empty!(ARGS), "AdvDiff");
julia> include(./src/run.jl)
```

Currently available problem names:

* AdvDiff
* Elliptic


## Plotting
For plotting we rely on `PlotlyJS`. If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)

