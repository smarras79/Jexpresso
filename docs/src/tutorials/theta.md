# Tutorial: Running the rising thermal bubble test case with Jexpresso.jl

This tutorial guides you through the process of creating a new physical problem from scratch. We will cover the necessary setup and command execution.

## Prerequisites

Before starting, ensure you have the following:

* **Julia Installation:** You need to have Julia installed on your system. You can download it from the official Julia website: [https://julialang.org/downloads/](https://julialang.org/downloads/)
* **Jexpresso Repository:** The `Jexpresso.jl` framework needs to be accessible. This tutorial assumes you have the repository cloned locally. If not, you can clone it using Git:
    ```bash
    git clone https://github.com/smarras79/Jexpresso.git
    git clone https://github.com/smarras79/JexpressoMeshes.git
    cd Jexpresso
    ln -s ../JexpressoMeshes/meshes .
    ```

# Creating a New Problem called Hello in, e.g., ```Jexpresso/problems/equations/CompEuler/```

## Step 1: Navigate to the Problems Directory

First, open your terminal and navigate to the `problems/equations/CompEuler` directory within your CompEuler repository.

```bash
cd problems/equations/CompEuler
```

## Step 2: Create a New Directory for Your Case
Next, create a new directory to house the files for your specific problem. Choose a descriptive name for Hello.

```
mkdir Hello
```

## Step 3: Enter the New Case Directory
Navigate into the newly created directory.

```
cd Hello
```

## Step 4: Copy Essential Files
Copy the generic solver files from the theta directory into your new case directory. These files provide the basic structure for your simulation.

```
cp ../theta/*.jl .
```

## Step 5: Configure Initial Conditions in initialize.jl
Open the initialize.jl file in a text editor. In this file, you will need to define and initialize the solution array q. The structure of q depends on the dimensionality of your problem and the number of conserved variables (e.g., density, momentum components, energy).

### Example of initialization for the 2D Euler equations density (rho), momentum (rho u), and potential temperature

```
function initialize(SD::NSD_2D, 
                    PT, 
                    mesh::St_mesh, 
                    inputs::Dict, 
                    OUTPUT_DIR::String, 
                    TFloat)
```

Define the solution variables and, optional, the array of output variables. 

The length of `qvars` will define the size of the problem. 
However, the optional array `qoutvars` can be longer or shorter and is only used by the output writing function.
If `qoutvars` is not defined, the the default output variables coincide with `qvars`.

```
    qvars    = ["ρ", "ρu", "ρv", "ρθ"]
    qoutvars = ["ρ", "u", "w", "θ", "p"]
```
Allocate space for the solution array:
```
    q = define_q(SD, 
                 mesh.nelem, 
                 mesh.npoin, 
                 mesh.ngl, 
                 qvars, 
                 TFloat, 
                 inputs[:backend]; 
                 neqs=length(qvars), 
                 qoutvars=qoutvars)
```

## Now initialize:

For example, a minimal version of [Jepresso/problems/equations/CompEuler/theta/initialize.jl](https://github.com/smarras79/Jexpresso/blob/master/problems/equations/CompEuler/theta/initialize.jl): may looks like this:

```
    if (inputs[:backend] == CPU())    
        PhysConst = PhysicalConst{Float64}()
        
            xc = 0.0; yc = 2500.0
        
            for ip = 1:mesh.npoin

                x, y = mesh.x[ip], mesh.y[ip]
                r = sqrt( (x - xc)^2 + (y - yc)^2 )
            
                θ = 300
                p = 100000
                ρ = 1.25

                u = 0.0
                v = 0.0

                q.qn[ip,1] = ρ
                q.qn[ip,2] = ρ*u
                q.qn[ip,3] = ρ*v
                q.qn[ip,4] = ρ*θ
                q.qn[ip,end] = p

            end
        end
    return q
end
```

WARNING: refer to a proper working [code](https://github.com/smarras79/Jexpresso/blob/master/problems/equations/CompEuler/theta/initialize.jl) rather than the simplified version above. The one above was given as an example of what an initialization file may look like.

## Add the fluxes and sources depending on the equations that you are solving

If we were to solve the 2D Euler equations of compressible flows with gravity, where `q` and the fluxes are defined as

$${\bf q}=\begin{bmatrix}
\rho \\
\rho u\\
\rho v\\
\rho \theta
\end{bmatrix}\quad {\bf F}1=\begin{bmatrix}
\rho u\\
\rho u^2 + p\\
\rho u v\\
\rho u \theta
\end{bmatrix}\quad {\bf F}2=\begin{bmatrix}
\rho v\\
\rho v u\\
\rho v^2 + p\\
\rho v \theta
\end{bmatrix}\quad {\bf S}=\begin{bmatrix}
0\\
0\\
-\rho g\\
0
\end{bmatrix}\quad \mu\nabla^2{\bf q}=\mu\begin{bmatrix}
0\\
u_{xx} + u_{zz}\\
v_{xx} + v_{zz}\\
\theta_{xx} + \theta_{zz}
\end{bmatrix},$$

then the function [user_flux.jl](https://github.com/smarras79/Jexpresso/blob/master/problems/equations/CompEuler/theta/user_flux.jl) is imply:

```
function user_flux!(F, 
                    G, 
                    SD::NSD_2D, 
                    q, qe,
                    mesh::St_mesh, 
                    ::CL, 
                    ::TOTAL; 
                    neqs=4, ip=1)

    PhysConst = PhysicalConst{Float64}()
    
    ρ  = q[1]
    ρu = q[2]
    ρv = q[3]
    ρθ = q[4]
    θ  = ρθ/ρ
    u  = ρu/ρ
    v  = ρv/ρ
    Pressure = perfectGasLaw_ρθtoP(PhysConst, ρ=ρ, θ=θ)
    
    F[1] = ρu
    F[2] = ρu*u .+ Pressure
    F[3] = ρv*u
    F[4] = ρθ*u

    G[1] = ρv
    G[2] = ρu*v
    G[3] = ρv*v .+ Pressure
    G[4] = ρθ*v
end
```
Notice how there are no loops and the `F` and `G` are exactly defined as you'd write them on paper.



