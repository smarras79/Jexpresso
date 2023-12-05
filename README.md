| **Documentation** |
|:------------ |
| [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://smarras79.github.io/Jexpresso/dev/)

# JEXPRESSO
A research software for the numerical solution of a system of arbitrary conservation laws using continuous spectral elements. DISCLAIMER: this is WIP. Its GPU expansion is also under development. 

Jexpresso uses arbitrarily high-order (3rd and above) continuous spectral elements to solve

$$\frac{\partial \bf q}{\partial t} + \sum_{i=1}^{nd}\nabla\cdot{{\bf F}_i({\bf q})} = \mu\nabla^2**V** + {\bf S}({\bf q}) + ~{\rm b.c.}$$

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
\end{bmatrix}, \mu\nabla^2**V**=\mu\begin{bmatrix}
u_{xx} + u_{zz}
\end{bmatrix}.$$


for a constant value of $\mu$ which is case-dependent.

NOTICE: PLEASE CONTACT ME IF YOU ARE INTERESTED IN TESTING THIS WIP. 
I WILL POINT YOU TO THE MOST EFFICIENT, but less general BRANCH OF THE CODE!

A research software for the numerical solution of conservation laws using spectral element methods. DISCLAIMER: this is WIP and only 2D is being maintained until parallelization is complete.

If you are interested in contributing, please get in touch.

# Some notes on using JEXPRESSO

To install and run the code assume Julia 1.9.3

## Setup with CPUs

```bash
>> cd $JEXPRESSO_HOME
>> julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```
followed by the following:

Push problem name to ARGS
You need to do this only when you run a new problem
```bash
julia> push!(empty!(ARGS), EQUATIONS::String, EQUATIONS_CASE_NAME::String);
julia> include("./src/Jexpresso.jl")
```

* PROBLEM_NAME is the name of your problem directory as $JEXPRESSO/problems/equations/problem_name
* PROBLEM_CASE_NAME is the name of the subdirectory containing the specific setup that you want to run: 

The path would look like 
```$JEXPRESSO/problems/equations/PROBLEM_NAME/PROBLEM_CASE_NAME```

Example 1: to solve the 2D Euler equations with buyoancy and two passive tracers defined in `problems/equations/CompEuler/thetaTracers` you would do the following:
```bash
julia> push!(empty!(ARGS), "CompEuler", "thetaTracers");
julia> include("./src/Jexpresso.jl")
```

<img src="assets/thetaTracersMesh.png"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />


Example 2: to solve the 2D Euler equations leading to a density current defined in `problems/equations/CompEuler/dc` you would do the following:
```bash
julia> push!(empty!(ARGS), "CompEuler", "dc");
julia> include("./src/Jexpresso.jl")
```

<img src="assets/dc.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Example 3: to solve the 1D wave equation  defined in `problems/equations/CompEuler/wave1d` you would do the following:
```bash
julia> push!(empty!(ARGS), "CompEuler", "wave1d");
julia> include("./src/Jexpresso.jl")
```

<img src="assets/wave1d-v.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />



For ready to run tests, there are the currently available equations names:

* CompEuler (option with total energy and theta formulation)

The code is designed to create any system of conservsation laws. See CompEuler/case1 to see an example of each file.
Details will be given in the documentation (still WIP). Write us if you need help.

More are already implemented but currently only in individual branches. They will be added to master after proper testing.

## Plotting
Files can be written to VTK (recommended) or png. For the png plots, we use [Makie](https://github.com/MakieOrg/Makie.jl). If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

For non-periodic 2D tests, the output can also be written to VTK files by setting the value "vtk" for the usier_input key :outformat

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)



$
\begin{matrix}
\rho  \\
\rho u \\
\rho v \\
\rho \theta 
\end{matrix}_t
+
\begin{matrix}
\rho u  \\
\rho uu + p\\
\rho vu \\
\rho \theta u
\end{matrix}_x
+
\begin{matrix}
\rho v  \\
\rho uv \\
\rho vv + p\\
\rho \theta v
\end{matrix}_y
= 
\begin{matrix}
0  \\
0 \\
\rho g \\
0
\end{matrix}
$
