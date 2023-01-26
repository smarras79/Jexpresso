# JEXPRESSO
A research and educational software for the numerical solution of 1D, 2D, and 3D PDEs using spectral and spectral element methods on CPUs and GPUs. DISCLAIMER: this is WIP.

If you are interested in contributing, please get in touch.

# Some notes on using JEXPRESSO

To install and run the code assume Julia
version 1.7.2.

The [MPI.jl][0] package that is used assumes that you have a working MPI installation

## Setup with CPUs

```bash
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```
You can test that things were installed properly with
```bash
julia --project=. $JEXPRESSO_HOME/src/run.jl PROBLEM_NAME
```

`$JEXPRESSO_HOME` is the path to the base JEXPRESSO directory on your computer (you can export it in your .bashrc or simply replace its value with the explicit name of the path)

`PROBLEM_NAME` must be the same as the problem directory in `$JEXPRESSO_HOME/src/problems/PROBLEM_NAME`
Currently available problem names:

* AdvDiff
* LinearCLaw



## Plotting
For plotting we rely on `PlotlyJS`. If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)

