---
title: 'Jexpresso: a fast, modular and general solver of partial differential equations on CPUs and GPUs using Julia'
tags:
  - julia
  - pdes
  - partial differential equations
  - spectral elements
  - finite differences
  - finite volumes
authors:
  - name: Simone Marras
    orcid: 0000-0002-7498-049X
    affiliation: "1, 2"
  - name: Yassine Tissaoui
    orcid: 0000-0000-0000-0000
    affiliation: 3
  - name: Hang Wang
    orcid: 0000-0000-0000-0000
    affiliation: 1
affiliations:
 - name: Department of Mechanical and Industrial Engineering, New Jersey Institute of Technology, Newark (NJ), USA
   index: 1
 - name: Center for Applied Mathematics and Statistics, New Jersey Institute of Technology, Newark (NJ), USA
   index: 2
 - name: Department of Mathematics, University of Wisconsin, Madison (WI), USA
   index: 3
date: 19 Mau 2024
bibliography: paper.bib
---

# Summary

Jexpresso (for Julia EXascale-REady SOlver) is general purpose solver of arbitrary systems of partial differential equations expressed in the form of conservation laws. 

While requiring minimal user intervention, the equations are provided by the user via the 


(see [@tissaoui2024])

![](code.pdf)

Jexpresso is an open-source software hosted at Github and distributed under an MIT license.

# Introdution
Jexpresso was used for a scientific study for the first time 

# Statement of need


# State of the field

# Code structure and user-interaction


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

### Add the sources:
Similarly, you handle the source through [user_source.jl](https://github.com/smarras79/Jexpresso/blob/master/problems/equations/CompEuler/theta/user_source.jl).

### Outouts:
By default, the output is written for the solution variables. 
Follow the output tutorial [here](./define_output_variables.md)


# References