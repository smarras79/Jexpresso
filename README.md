# JEXPRESSO
A research and educational software for the numerical solution of 1D, 2D, and 3D PDEs using spectral and spectral element methods on CPUs and GPUs. DISCLAIMER: this is WIP.

If you are interested in contributing, please get in touch.

# Some notes on using JEXPRESSO

To install and run the code assume Julia
version 1.7.2 or higher (tested up to 1.8.5)

## Setup with CPUs

```bash
>> cd $JEXPRESSO_HOME
```

```bash
>> julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```


```bash
julia > include("/src/run.jl")
```

`$JEXPRESSO_HOME` is the path to the base JEXPRESSO directory on your computer (you can export it in your .bashrc or simply replace its value with the explicit name of the path)

Notice: command line arguments have been de-activated while we add several solvers beyond the AdvDiff equation
`PROBLEM_NAME` must be the same as the problem directory in `$JEXPRESSO_HOME/src/problems/PROBLEM_NAME`
Currently available problem names:

* AdvDiff

## Plotting
For plotting we rely on `PlotlyJS`. If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)

