# Rising thermal bubble

```@contents
Pages = [
  "Jexpresso.md",
  ]
```

## Equations:

5. 2D Euler equations of compressible flows with gravity and N passive chemicals $c_i, \forall i=1,...,N$ 


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

Example 1: to solve the 2D Euler equations with buyoancy and two passive tracers defined in `problems/equations/CompEuler/thetaTracers` you would do the following:


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


```bash
push!(empty!(ARGS), "CompEuler", "thetaTracers");
include("./src/Jexpresso.jl")
```

<img src="assets/thetaTracersMesh.png"
     alt="Markdown icon"
     style="float: left; margin-right: 5px;" />


```@contents
Pages = ["Jexpresso.md"]
```
 