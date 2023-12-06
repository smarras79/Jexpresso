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

# Equations:
Jexpresso uses arbitrarily high-order (3rd and above) continuous spectral elements to approximate in space the equations below:

$$\frac{\partial \bf q}{\partial t} + \sum_{i=1}^{nd}\nabla\cdot{{\bf F}_i({\bf q})} = \mu\nabla^2{\bf q} + {\bf S}({\bf q}) + ~{\rm b.c.}$$

where the vectors ${\bf q}$, ${\bf F}$, and ${\bf S}$ are problem-dependent as shown below,
and are taken to be zero vectors of the appropriate size when not explicitly stated otherwise.

The Julia package [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) is used for time discretization and stepping.

In order, we provide tests and results for the following equations:


1. 1D wave equation:
   
$${\bf q}=\begin{bmatrix}
u \\
v
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
v\\
u
\end{bmatrix}$$

2: 1D shallow water:

$${\bf q}=\begin{bmatrix}
h \\
u
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
Uh + Hu\\
gh + Uu
\end{bmatrix},$$

where $H$ and $U$ are a reference height and velocity, respectively.

3. 2D Helmholtz:
   
$${\bf S}=\begin{bmatrix}
\alpha^2 u + f(x,z)
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
u_{xx} + u_{zz}
\end{bmatrix},$$

for a constant value of $\alpha$ and $\mu$, which are case-dependent.

4. 2D scalar advection-diffusion:

$${\bf q}=\begin{bmatrix}
q\\
\end{bmatrix}\quad {\bf F}_1=\begin{bmatrix}
qu\\
\end{bmatrix}\quad {\bf F}_2=\begin{bmatrix}
qv\\
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
q_{xx} + q_{zz}
\end{bmatrix},$$

5. 2D Euler equations of compressible flows with gravity and N passive chemicals ``c_i, \forall i=1,...,N`` 

$${\bf q}=\begin{bmatrix}
\rho \\
\rho u\\
\rho v\\
\rho \theta\\
\rho c1\\
...\\
\rho cN
\end{bmatrix}\quad {\bf F}_1=\begin{bmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
\rho u \theta\\
\rho u c1\\
...\\
\rho u cN
\end{bmatrix}\quad {\bf F}_2=\begin{bmatrix}
\rho v\\
\rho v u\\
\rho v^2 + p\\
\rho v \theta\\
\rho v c1\\
...\\
\rho v cN
\end{bmatrix}\quad {\bf S}=\begin{bmatrix}
0\\
0\\
-\rho g\\
0\\
0\\
...\\
0
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
0\\
u_{xx} + u_{zz}\\
v_{xx} + v_{zz}\\
\theta_{xx} + \theta_{zz}\\
c1_{xx} + c1_{zz}\\
...\\
cN_{xx} + cN_{zz}
\end{bmatrix}.$$

The equation of state for a perfect gas is used to close the system.

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
