# Jexpresso.jl

Documentation of `Jexpresso.jl`.

!!! note

     This documentation is and will always be WIP!
          
## Introduction

Jexpresso is a CPU/GPU research software for the numerical solution of a system of arbitrary conservation laws in 1D, 2D, 3D using continuous spectral elements (SEM). Nevertheless, the code is built so that any other numerical method can be added. For example, the Jexpresso already contains a 1D finite difference implementation as well.

Jexpresso is written in the [Julia programming language](https://julialang.org/) and was thought to be modular and allow any user to add any equations in any dimensions without knowing anything about numerical methods. 

## Do I need to know Julia to use Jexpresso?
Yes and no. It depends how much you are interested in adding your own equation set in the code rather than using it as a black box. 

The following are useful resources about Julia:
* Julia webpage [docs.julialang.org](https://docs.julialang.org/)
* Official list of learning resources [julialang.org/learning](https://julialang.org/learning/)

```@contents
Pages = [
  "Jexpresso.md",
  ]
```

## Equations:
Jexpresso uses arbitrarily high-order (3rd and above) **continuous spectral elements** to solve

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
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
qu\\
\end{bmatrix}\quad {\bf F}=\begin{bmatrix}
qv\\
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
q_{xx} + q_{zz}
\end{bmatrix},$$

5. 2D Euler equations of compressible flows with gravity and N passive chemicals $c_i, \forall i=1,...,N$ 

$${\bf q}=\begin{bmatrix}
\rho \\
\rho u\\
\rho v\\
\rho \theta\\
\rho c1\\
...\\
\rho cN
\end{bmatrix}\quad {\bf F}1=\begin{bmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
\rho u \theta\\
\rho u c1\\
...\\
\rho u cN
\end{bmatrix}\quad {\bf F}2=\begin{bmatrix}
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


6. 3D Euler equations of compressible flows with gravity

$${\bf q}=\begin{bmatrix}
\rho \\
\rho u\\
\rho v\\
\rho w\\
\rho \theta\\
\end{bmatrix}\quad {\bf F}1=\begin{bmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
\rho u w\\
\rho u \theta\\
\end{bmatrix}\quad {\bf F}2=\begin{bmatrix}
\rho v\\
\rho v u\\
\rho v^2 + p\\
\rho v w\\
\rho v \theta\\
\end{bmatrix}\quad {\bf F}3=\begin{bmatrix}
\rho w\\
\rho w u\\
\rho w v\\
\rho w^2 + p\\
\rho w \theta\\
\end{bmatrix}\quad {\bf S}=\begin{bmatrix}
0\\
0\\
0\\
-\rho g\\
0\\
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
0\\
u_{xx} + u_{yy} + u_{zz}\\
v_{xx} + v_{yy} + v_{zz}\\
w_{xx} + w_{yy} + w_{zz}\\
\theta_{xx} + \theta_{yy} + \theta_{zz}\\
\end{bmatrix}.$$


If you are interested in contributing, please get in touch:
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)


# Some notes on using JEXPRESSO

To install and run the code assume Julia 1.10.0

## Setup with CPUs

```bash
cd $JEXPRESSO_HOME
julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
```
followed by the following:

Push problem name to ARGS
You need to do this only when you run a new problem
```bash
push!(empty!(ARGS), EQUATIONS::String, EQUATIONS_CASE_NAME::String);
include("./src/Jexpresso.jl")
```

* PROBLEM_NAME is the name of your problem directory as $JEXPRESSO/problems/equations/problem_name
* PROBLEM_CASE_NAME is the name of the subdirectory containing the specific setup that you want to run: 

The path would look like 
```$JEXPRESSO/problems/equations/PROBLEM_NAME/PROBLEM_CASE_NAME```

## Tutorials

The following tutorials will introduce you to the functionality of
Jexpresso.jl.

```@contents
Pages = [
    "features/performance.md",
    "tutorials/user_inputs.md",
    "tutorials/theta.md",
    ]
Depth = 2
```

Example 1: to solve the 2D Euler equations with buyoancy and two passive tracers defined in `problems/equations/CompEuler/thetaTracers` you would do the following:
```bash
push!(empty!(ARGS), "CompEuler", "thetaTracers");
include("./src/Jexpresso.jl")
```

<img src="assets/thetaTracersMesh.png"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />


Example 2: to solve the 2D Euler equations leading to a density current defined in `problems/equations/CompEuler/dc` you would do the following:
```bash
push!(empty!(ARGS), "CompEuler", "dc");
include("./src/Jexpresso.jl")
```

<img src="assets/dc.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Example 3: to solve the 1D wave equation  defined in `problems/equations/CompEuler/wave1d` you would do the following:
```bash
push!(empty!(ARGS), "CompEuler", "wave1d");
include("./src/Jexpresso.jl")
```

<img src="assets/wave1d-v.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />



For ready to run tests, there are the currently available equations names:

* CompEuler (option with total energy and theta formulation)

The code is designed to create any system of conservsation laws. See CompEuler/case1 to see an example of each file.
Details will be given in the documentation (still WIP). Write us if you need help.

More are already implemented but currently only in individual branches. They will be added to master after proper testing.

## Laguerre semi-infinite element test suite
This section contains instructions to run all of the test cases presented in

```
@article{tissaoui2024,
  doi = {},
  url = {},
  year = {2020},
  volume = {},
  number = {},
  pages = {},
  author = {Yassine Tissaoui and James F. Kelly and Simone Marras}
  title = {Efficient Spectral Element Method for the Euler Equations on Unbounded Domains in Multiple Dimensions},
  journal = {arXiv},
}
```

Test 1: 1D wave equation with Laguerre semi-infinite element absorbing layers

The problem is defined in [`problems/CompEuler/wave1d_lag`](https://github.com/smarras79/Jexpresso/tree/yt/Laguerre_test_suite/problems/equations/CompEuler/wave1d_lag) and by default output will be written to `output/CompEuler/wave1d_lag`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "CompEuler", "wave1d_lag");
include("./src/Jexpresso.jl")
```

<img src="assets/wave_v_4.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 2: 1D wave train for linearized shallow water equations

The problem is defined in [`problems/equations/AdvDiff/Wave_Train`](https://github.com/smarras79/Jexpresso/tree/yt/Laguerre_test_suite/problems/equations/AdvDiff/Wave_Train) and by default output will be written to `output/AdvDiff/Wave_Train`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "AdvDiff", "Wave_Train");
include("./src/Jexpresso.jl")
```

<img src="assets/Wave_Train_final.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

A second version of this tests generate images with the solutions at different times overlapped.

This version is defined in [`problems/equations/AdvDiff/Wave_Train_Overlapping_Plot`](https://github.com/smarras79/Jexpresso/tree/yt/Laguerre_test_suite/problems/equations/AdvDiff/Wave_Train_Overlapping_Plot) and by default output will be written to `output/AdvDiff/Wave_Train_Overlapping_Plot`. To run this version of the problem execute the following from the Julia command line:

```bash
push!(empty!(ARGS), "AdvDiff", "Wave_Train_Overlapping_Plot");
include("./src/Jexpresso.jl")
```

<img src="assets/Wave_Train_overlap.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 3: 2D advection-diffusion equation

The problem is defined in [`problems/equations/AdvDiff/2D_laguerre`](https://github.com/smarras79/Jexpresso/tree/yt/Laguerre_test_suite/problems/equations/AdvDiff/2d_Laguerre) and by default output will be written to `output/AdvDiff/2D_laguerre`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "AdvDiff", "2D_laguerre");
include("./src/Jexpresso.jl")
```

<img src="assets/ad2d-4s-line.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 4: 2D Helmholtz equation

The problem is defined in [`problems/equations/Helmholtz/case1`](https://github.com/smarras79/Jexpresso/tree/yt/Laguerre_test_suite/problems/equations/Helmholtz/case1) and by default output will be written to `output/Helmholtz/case1`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "Helmholtz", "case1");
include("./src/Jexpresso.jl")
```

<img src="assets/Helmholtz_from_jexpresso-line.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 5: Rising thermal bubble

The problem is defined in [`problems/equations/CompEuler/theta_laguerre`](https://github.com/smarras79/Jexpresso/tree/yt/Laguerre_test_suite/problems/equations/CompEuler/theta_laguerre) and by default output will be written to `output/CompEuler/theta_laguerre`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "CompEuler", "theta_laguerre");
include("./src/Jexpresso.jl")
```

<img src="assets/48.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 6: Hydrostatic linear mountain waves

The problem is defined in [`problems/equations/CompEuler/HSmount_Lag_working`](https://github.com/smarras79/Jexpresso/tree/yt/Laguerre_test_suite/problems/equations/CompEuler/HSmount_Lag_working) and by default output will be written to `output/CompEuler/HSmount_Lag_working`. To solve this problem run the following commands from the Julia command line:

```bash      
push!(empty!(ARGS), "CompEuler", "HSmount_Lag_working");
include("./src/Jexpresso.jl")
```

<img src="assets/wvelo.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

## Plotting
Files can be written to VTK (recommended) or png. For the png plots, we use [Makie](https://github.com/MakieOrg/Makie.jl). If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

For non-periodic 2D tests, the output can also be written to VTK files by setting the value "vtk" for the usier_input key :outformat

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu)


## Manual

```@contents
Pages = ["Jexpresso.md"]
```
