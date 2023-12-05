# Documentation

```@meta
CurrentModule = Jexpresso
DocTestSetup = quote
    using Jexpresso
end
```

```@autodocs
Modules = [Jexpresso]
```

# Jexpresso.jl

Documentation of `Jexpresso.jl`.

!!! note

     This documentation is and will always be WIP!
     
     A research software for the numerical solution of a system of an arbitrary number of conservation laws using continuous spectral elements. DISCLAIMER: this is WIP and only 2D is being maintained until parallelization is complete.

     Jexpresso uses arbitrarily high-order (3rd and above) continuous spectral elements to solve

$$\frac{\partial \bf q}{\partial t} + \sum_{i=1}^{nd}\nabla\cdot{{\bf F}_i({\bf q})} = \mu\nabla^2**q** + {\bf S}({\bf q}) + ~{\rm b.c.}$$

where the vectors **q**, **F**, and **S** are problem-dependent as shown below,
and are taken to be zero vectors of the appropriate size when not explicitly stated otherwise.

     In order, we provide tests and results for the following equations:
1. 1D wave equation:
   
$$**q**=\begin{bmatrix}
u \\
v
\end{bmatrix} **F**=\begin{bmatrix}
v\\
u
\end{bmatrix}$$

2: 1D shallow water:

$$**q**=\begin{bmatrix}
h \\
u
\end{bmatrix} **F**=\begin{bmatrix}
Uh + Hu\\
gh + Uu
\end{bmatrix}$$
   
3. 2D Helmholtz:
   
$$**S**=\begin{bmatrix}
\alpha^2 u + f(x,z)
\end{bmatrix}, \mu\nabla^2**q**=\mu\begin{bmatrix}
u_{xx} + u_{zz}
\end{bmatrix}.$$

4. 2D scalar advection-diffusion:

$$**S**=\begin{bmatrix}
q\\
\end{bmatrix} **F**=\begin{bmatrix}
qu\\
\end{bmatrix} **F**=\begin{bmatrix}
qv\\
\end{bmatrix}, \mu\nabla^2**q**=\mu\begin{bmatrix}
q_{xx} + q_{zz}
\end{bmatrix},$$
where $H$ and $U$ are a reference height and velocity, respectively.



for a constant value of $\mu$ which is case-dependent.

If you are interested in contributing, please get in touch:
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)

# Some notes on using JEXPRESSO

To install and run the code assume Julia 1.9.3

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

For example, if you wanted to run `CompEuler` with the setup defined inside the case directory `theta`, then you would do the following:
```bash
julia> push!(empty!(ARGS), "CompEuler", "theta");
julia> include("./src/Jexpresso.jl")
```

For ready to run tests, there are the currently available equations names:

* CompEuler (option with total energy and theta formulation)

The code is designed to create any system of conservsation laws. See CompEuler/case1 to see an example of each file.
Details will be given in the documentation (still WIP). Write us if you need help.

More are already implemented but currently only in individual branches. They will be added to master after proper testing.

## Plotting
For plotting we rely on [Makie](https://github.com/MakieOrg/Makie.jl). If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

For non-periodic 2D tests, the output can also be written to VTK files by setting the value "vtk" for the usier_input key :outformat


## Manual

```@contents
Pages = ["Jexpresso.md"]
```
