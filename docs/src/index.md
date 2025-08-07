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

To install and run the code assume Julia 1.11.0

## Setup with CPUs on one core from the Julia REPL:

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

## Setup an MPI run:

1. Basic Execution
```bash
mpiexec -n <NPROCS> julia --project=. -e 'push!(empty!(ARGS), "<EQUATIONS>", "<CASE_NAME>"); include("./src/Jexpresso.jl")'
```

2. Implementation-Specific Examples
```bash
mpiexec -n 4 julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "3d"); include("./src/Jexpresso.jl")'
```

### Troubleshooting

- **Library conflicts:** Clear existing preferences:
  ```bash
  rm -f LocalPreferences.toml
  ```
- **Path issues:** Verify paths with:
  ```bash
  which mpiexec
  which mpirun
  ```
  You may have to use the full aboslute path to mpiexec or mpirun and to julia like this if necessary. For example, if your `mpirun` lives in `/opt/homebrew/Cellar/open-mpi/5.0.6/bin/mpirun`, then you may want to run the code with the full path like this:

  ```
  /opt/homebrew/Cellar/open-mpi/5.0.6/bin/mpirun -n 4 /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "theta"); include("./src/Jexpresso.jl")'
  ```

  Important notice: you need to use the Julia version that is correctly linked to the same MPI installation.
If you have multiple Julia installations on your machine, use the full path to the one linked to the mpirun at hand.

- **Version mismatches:** Ensure consistent versions:
  ```bash
  mpicc --version
  mpif90 --version
  ```


## Performance
Jexpresso leverages the properties of the Julia language to make it as fast as a compiled code. Some performance, measured against a legacy and massive Fortran 90/Modern Fortran code, are shown in page
``@contents
Pages = [
    "features/performance.md",
]

## Tutorials

The following tutorials will introduce you to the functionality of
Jexpresso.jl.

```@contents
Pages = [
    "tutorials/user_inputs.md",
    "tutorials/theta.md",
    "tutorials/define_output_variables.md",
    "tutorials/laguerre_paper.md",
    ]
Depth = 2
```

## Manual

```@contents
Pages = [
          Jexpresso.md",
     ]
```
