# <img src="./assets/logo-ext2.png" width="500" title="JEXPRESSO logo">

| **Documentation** |
|:------------ |
 [![](https://img.shields.io/badge/docs-stable-blue.svg)](https://smarras79.github.io/Jexpresso/dev/) [![](https://img.shields.io/badge/docs-dev-blue.svg)](https://smarras79.github.io/Jexpresso/dev/) |
|**Build Status** |
| [![CI](https://github.com/smarras79/Jexpresso/actions/workflows/CI.yml/badge.svg?branch=master)](https://github.com/smarras79/Jexpresso/actions?query=workflow%3ACI)
| **Contacts**  |
| [![Simone Marras](https://img.shields.io/badge/Simone%20Marras-smarras%40njit.edu-8e7cc3)](mailto:smarras@njit.edu) |
| [![Yassine Tissaoui](https://img.shields.io/badge/Yassine%20Tissaoui-yt277%40njit.edu-8e7cc3)](mailto:yt277@njit.edu) |
| [![Hang Wang](https://img.shields.io/badge/Hang%20Wang-hang.wang%40njit.edu-8e7cc3)](mailto:hang.wang@njit.edu) |
| **Citation** |
| [![DOI](https://img.shields.io/badge/article-arXiv:2401.05624-green)](https://doi.org/10.48550/arXiv.2401.05624) |

# JEXPRESSO
A CPU and GPU research software for the numerical solution of a system of arbitrary conservation laws using **continuous spectral elements** and finite differences in **1D, 2D, 3D**. DISCLAIMER: this will always be WIP! Contact us to join the team of developers!

Suggested Julia version: 1.11.2

If you use Jexpresso please drop us a line to let us know. We'd like to add a link to your paper or work on this page.

Please cite Jexpresso using:

```
@article{tissaoui2024,
  author = {Y. Tissaoui and J. F. Kelly and S. Marras}
  title = {Efficient Spectral Element Method for the Euler Equations on Unbounded Domains},
  volume ={487},
  pages={129080},
  year = {2024},
  journal = {App. Math. Comput.},
}
```

# Equations:
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
\end{bmatrix}\quad {\bf S}=\begin{bmatrix}
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

To install and run the code assume Julia 1.11.2

## Setup with CPUs

```bash
>> cd $JEXPRESSO_HOME
>> julia --project=. -e "using Pkg; Pkg.instantiate(); Pkg.API.precompile()"
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

Example of cloud simulations (please contact us to run this because its branch has not been merged into master yet)

<img src="assets/bomex.png"
     alt="Markdown icon"
     style="float: left; margin-right: 3.5px;" />


Examples available in this branch:

Example 1: to solve the 2D Euler equations with buoyancy and two passive tracers defined in `problems/equations/CompEuler/thetaTracers` you would do the following:
```bash
push!(empty!(ARGS), "CompEuler", "thetaTracers");
include("./src/Jexpresso.jl")
```

<img src="assets/thetaTracersMeshUnstr.png"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />


Example 2: to solve the 3D Euler equations with buoyancy defined in `problems/equations/CompEuler/3d` you would do the following:
```bash
push!(empty!(ARGS), "CompEuler", "3d");
include("./src/Jexpresso.jl")
```

<img src="assets/rtb3d.png"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />


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
  author = {Y. Tissaoui and J. F. Kelly and S. Marras}
  title = {Efficient Spectral Element Method for the Euler Equations on Unbounded Domains},
  volume ={487},
  pages={129080},
  year = {2024},
  journal = {App. Math. Comput.},
}
```

Test 1: 1D wave equation with Laguerre semi-infinite element absorbing layers

The problem is defined in [`problems/CompEuler/wave1d_lag`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/CompEuler/wave1d_lag) and by default output will be written to `output/CompEuler/wave1d_lag`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "CompEuler", "wave1d_lag");
include("./src/Jexpresso.jl")
```

<img src="assets/wave_v_4.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 2: 1D wave train for linearized shallow water equations

The problem is defined in [`problems/equations/AdvDiff/Wave_Train`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/AdvDiff/Wave_Train) and by default output will be written to `output/AdvDiff/Wave_Train`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "AdvDiff", "Wave_Train");
include("./src/Jexpresso.jl")
```

<img src="assets/Wave_Train_final.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />


Test 3: 2D advection-diffusion equation

The problem is defined in [`problems/equations/AdvDiff/2D_laguerre`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/AdvDiff/2d_Laguerre) and by default output will be written to `output/AdvDiff/2D_laguerre`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "AdvDiff", "2D_laguerre");
include("./src/Jexpresso.jl")
```

<img src="assets/ad2d-4s-line.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 4: 2D Helmholtz equation

The problem is defined in [`problems/equations/Helmholtz/case1`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/Helmholtz/case1) and by default output will be written to `output/Helmholtz/case1`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "Helmholtz", "case1");
include("./src/Jexpresso.jl")
```

<img src="assets/Helmholtz_from_jexpresso-line.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />

Test 5: Rising thermal bubble

The problem is defined in [`problems/equations/CompEuler/theta_laguerre`](https://github.com/smarras79/Jexpresso/tree/master/problems/equations/CompEuler/theta_laguerre) and by default output will be written to `output/CompEuler/theta_laguerre`. To solve this problem run the following commands from the Julia command line:

```bash
push!(empty!(ARGS), "CompEuler", "theta_laguerre");
include("./src/Jexpresso.jl")
```

<img src="assets/48.png"
     alt="Markdown icon"
     style="float: left; margin-right: 7px;" />


```bash
push!(empty!(ARGS), "CompEuler", "3d_bomex");
include("./src/Jexpresso.jl")
```
<img src="assets/bomex.png"
     alt="Markdown icon"
     style="float: left; margin-right: 3.5px;" />

## Setup and Run with MPI

JEXPRESSO supports parallel execution using either OpenMPI or MPICH. Follow these steps to configure and run with your preferred MPI implementation.

### 1. Install MPI Implementation

Choose either OpenMPI or MPICH:

#### OpenMPI Installation
```bash
# Ubuntu/Debian
sudo apt install libopenmpi-dev openmpi-bin

# macOS (Homebrew)
brew install open-mpi

# Verify installation
mpiexec --version
```

#### MPICH Installation
```bash
# Ubuntu/Debian
sudo apt install mpich libmpich-dev

# macOS (Homebrew) 
brew install mpich

# Verify installation
mpiexec --version
```

### 2. Configure MPI Preferences

#### Automatic Configuration (Default Path)
Use this command when MPI (OpenMPI/MPICH) is installed in standard system paths (`/usr/bin`, `/usr/local/bin`, etc.):
```bash
julia --project=. -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; MPIPreferences.use_system_binary()'
```

#### Manual Configuration (For Multiple MPI Installations or MPI not in Default Path)
For MPI installations in non-standard locations (e.g., /opt/openmpi, or custom paths):
```bash
julia --project=. -e 'using Pkg; Pkg.add("MPIPreferences"); using MPIPreferences; MPIPreferences.use_system_binary(;extra_paths = ["/where/your/mpi/lib"])'
```
If MPI is installed via homebrew on macOS, the MPI lib path is:
```bash
/opt/homebrew/lib
```

### 3. Running with MPI

#### Basic Execution
```bash
mpiexec -n <NPROCS> julia --project=. -e 'push!(empty!(ARGS), "<EQUATIONS>", "<CASE_NAME>"); include("./src/Jexpresso.jl")'
```

#### Implementation-Specific Examples
```bash
mpiexec -n 4 julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "3d"); include("./src/Jexpresso.jl")'
```
#### Script
You can simplify the run steps with a `runjexpresso` script like this:
```bash
#!/bin/bash

MPIRUN=/YOUR/PATH/TO/mpirun
JULIA=/YOUR/PATH/TO/julia

$MPIRUN -np $1 $JULIA --project=. -e 'push!(empty!(ARGS), "'"$2"'", "'"$3"'"); include("./src/Jexpresso.jl")' "$@"
```
and run it like this:
```bash
./runjexpresso 4 CompEuler theta
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
  You may have to use the full aboslute path to mpiexec or mpirun and to julia like this if necessary:
  ```
  /opt/homebrew/Cellar/open-mpi/5.0.6/bin/mpirun -n 4 /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia --project=. -e 'push!(empty!(ARGS), "CompEuler", "theta"); include("./src/Jexpresso.jl")'
  ```


- **Version mismatches:** Ensure consistent versions:
  ```bash
  mpicc --version
  mpif90 --version
  ```



<!-- <img src="assets/mpi_performance_comparison.png" width="700" alt="MPI Performance Comparison"> -->

## Plotting
Files can be written to VTK (recommended) or png (png is now only used for 1D results). For the png plots, we use [Makie](https://github.com/MakieOrg/Makie.jl). If you want to use a different package,
modify ./src/io/plotting/jplots.jl accordinly.

## Contacts
[Simone Marras](mailto:smarras@njit.edu), [Yassine Tissaoui](mailto:yt277@njit.edu), [Hang Wang](mailto:hang.wang@njit.edu)
