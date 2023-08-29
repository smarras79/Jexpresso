# JEXPRESSO

IMPORTANT NOTICE: this master branch is NOT the latest. We are working on the performant and optimized version of the code.
You can try  CompEuler/theta with branch ```sm/efficientMemoryAccess```

A research and educational software for the numerical solution of 1D, 2D, and 3D PDEs using spectral and spectral element methods on CPUs and GPUs. DISCLAIMER: this is WIP.

If you are interested in contributing, please get in touch.

# Some notes on using JEXPRESSO

To install and run the code assume Julia
version 1.7.2 or higher (tested up to 1.8.5)

## Setup with CPUs

```bash
>> cd $JEXPRESSO_HOME
>> julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```
followed by the following:

Push equations name to ARGS
You need to do this only when you run a new equations
```bash
julia> push!(empty!(ARGS), EQUATIONS::String, EQUATIONS_CASE_NAME::String);
julia> include("./src/Jexpresso.jl")
```

* EQUATIONS is the name of your equations directory as $JEXPRESSO/src/equations/equations
* EQUATIONS_CASE_NAME is the name of the subdirectory containing the specific setup that you want to run: 

The path would look like 
```$JEXPRESSO/src/equations/EQUATIONS/EQUATIONS_CASE_NAME```

For example, if you wanted to run `AdvDiff` with the setup defined inside the case directory `case1`, then you would do the following:
```bash
julia> push!(empty!(ARGS), "AdvDiff", "case1");
julia> include("./src/Jexpresso.jl")
```

For ready to run tests, there are the currently available equations names:

* AdvDiff
* Elliptic
* CompEuler (option with total energy and theta formulation)

The code is designed to create any system of conservsation laws. See CompEuler/case1 to see an example of each file.
Details will be given in the documentation (still WIP). Write us if you need help.

More are already implemented but currently only in individual branches. They will be added to master after proper testing.

## Plotting
For plotting we rely on [Makie](https://github.com/MakieOrg/Makie.jl). If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

For non-periodic 2D tests, the output can also be written to VTK files by setting the value "vtk" for the usier_input key :outformat

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)
