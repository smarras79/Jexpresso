# Performance of Jexpresso on CPU

Jexpresso was coded to minimize memory access and speed. With the goal of making Jexpresso a community solver of PDEs, it was mandatory to make it at least as fast as a compiled-language code.
Jexpresso was benchmarked against a legacy code for atmospheric modeling written in Fortran 90/95/modern Fortran. The two software packages were compared on the same CPU core.

The code speed was measured for the solution of the compressible Naver-Stokes equations with gravity

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
\end{bmatrix}.$$

## Speed

Table: Wall clock time of Jexpresso vs a legacy F90/Modern Fortran code for numerical weather prediction. Simulated 100 seconds of a rising-thermal-bubble test. The name of the time integrators may be different for the two codes so that the notation jexpresso/numa is used to indicate both. The wall clock times are to be taken with a $\pm 0.2$ due to a small variability from one simulation to the next one. 

Timing was measured using Julia 1.9.3 on a Macbook Air M1 2020, with macOS Big Sur Version 11.6.

| Time integrator           | max $\Delta t$ (s)        | Effective resolution (m)  | Order	         | ${\color{red}{Jexpresso}}$ (s)   | ${\color{blue}{F90}}$ (s)|
| :-------------------------| :-------------------------| :-------------------------| :-------------------------| :-------------------------| :------------------------|
| SSPRK53        	    | 0.3                       | "                         | "    			| 9.00  		    | 10.53  		       |
| SSPRK33		    | 0.2      		        | $$125\times 125$$         | 4   		        | 9.75			    | 9.2028		       |
| SSPRK54     	     	    | 0.4                       | "                         | "    			| 10.47 		    |       NA 		       |
| DP5 (Dormand-Prince RK54) | 0.6                       | "                         | "   			| 19.80 		    | 	    NA 		       |
| SSPRK73                   | 0.4                       | "                         | "    		 	| 12.95 		    | 	    NA  	       |
| SSPRK104    	            | 0.6                       | "                         | "                         | 12.50 		    | 	    NA		       |
| CarpenterKennedy2N54      | 0.4                       | "                         | "                         | 10.57 		    | 	    NA		       |
| Tsit5                     | 2.0 (adaptive)            | "                         | "   	                | 19.08 		    | 	    NA		       |


## Mass conservation

Table: Mass conservation of the **advective vs flux forms** of the equations and sensitivity to the time integrators of  DifferentialEquations.jl for the RTB at t = 1000 s viscous. Results are for inexact integration.

| Time integrator       | Advection form        | Flux from              |
| :---------------------| :---------------------| :----------------------|
| MSRK5                 | 7.818062181220379e-16 | 3.9090310906101895e-16 |
| SSPRK53               | 7.622610626689869e-15 | 1.9545155453050947e-16 | 
| SSPRK33               | 5.081740417793246e-15 | 1.1727093271830568e-15 | 

