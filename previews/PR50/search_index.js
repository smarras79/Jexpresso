var documenterSearchIndex = {"docs":
[{"location":"features/performance/#Performance-of-Jexpresso-on-CPU","page":"Performance","title":"Performance of Jexpresso on CPU","text":"","category":"section"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"Jexpresso was coded to minimize memory access and speed. With the goal of making Jexpresso a community solver of PDEs, it was mandatory to make it at least as fast as a compiled-language code. Jexpresso was benchmarked against a legacy code for atmospheric modeling written in Fortran 90/95/modern Fortran. The two software packages were compared on the same CPU core.","category":"page"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"The code speed was measured for the solution of the compressible Naver-Stokes equations with gravity","category":"page"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"bf q=beginbmatrix\nrho \nrho u\nrho v\nrho theta\nendbmatrixquad bf F1=beginbmatrix\nrho u\nrho u^2 + p\nrho u v\nrho u theta\nendbmatrixquad bf F2=beginbmatrix\nrho v\nrho v u\nrho v^2 + p\nrho v theta\nendbmatrixquad bf S=beginbmatrix\n0\n0\n-rho g\n0\nendbmatrixquad munabla^2bf q=mubeginbmatrix\n0\nu_xx + u_zz\nv_xx + v_zz\ntheta_xx + theta_zz\nendbmatrix","category":"page"},{"location":"features/performance/#Speed","page":"Performance","title":"Speed","text":"","category":"section"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"Table: Wall clock time of Jexpresso vs a legacy F90/Modern Fortran code for numerical weather prediction. Simulated 100 seconds of a rising-thermal-bubble test. The name of the time integrators may be different for the two codes so that the notation jexpresso/numa is used to indicate both. The wall clock times are to be taken with a pm 02 due to a small variability from one simulation to the next one. ","category":"page"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"Timing was measured using Julia 1.9.3 on a Macbook Air M1 2020, with macOS Big Sur Version 11.6.","category":"page"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"Time integrator max Delta t (s) Effective resolution (m) Order colorredJexpresso (s) colorblueF90 (s)\nSSPRK33 0.2 125times 125 4 9.75 9.2028\nSSPRK53 0.3 \" \" 9.00 10.53\nSSPRK54 0.4 \" \" 10.47 NA\nDP5 (Dormand-Prince RK54) 0.6 \" \" 19.80 NA\nSSPRK73 0.4 \" \" 12.95 NA\nSSPRK104 0.6 \" \" 12.50 NA\nCarpenterKennedy2N54 0.4 \" \" 10.57 NA\nTsit5 2.0 (adaptive) \" \" 19.08 NA","category":"page"},{"location":"features/performance/#Mass-conservation","page":"Performance","title":"Mass conservation","text":"","category":"section"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"Table: Mass conservation of the advective vs flux forms of the equations and sensitivity to the time integrators of  DifferentialEquations.jl for the RTB at t = 1000 s viscous. Results are for inexact integration.","category":"page"},{"location":"features/performance/","page":"Performance","title":"Performance","text":"Time integrator Advection form Flux from\nMSRK5 7.818062181220379e-16 3.9090310906101895e-16\nSSPRK53 7.622610626689869e-15 1.9545155453050947e-16\nSSPRK33 5.081740417793246e-15 1.1727093271830568e-15","category":"page"},{"location":"#Jexpresso.jl","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Documentation of Jexpresso.jl.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"note: Note\nThis documentation is and will always be WIP!","category":"page"},{"location":"#Introduction","page":"Jexpresso.jl","title":"Introduction","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Jexpresso is a CPU/GPU research software for the numerical solution of a system of arbitrary conservation laws in 1D, 2D, 3D using continuous spectral elements. Nevertheless, the code is built so that any other numerical method can be added. For example, the Jexpresso already contains a 1D finite difference implementation as well.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Jexpresso is written in the Julia programming language and was thought to be modular and allow any user to add any equations in any dimensions without knowing anything about numerical methods. ","category":"page"},{"location":"#Do-I-need-to-know-Julia-to-use-Jexpresso?","page":"Jexpresso.jl","title":"Do I need to know Julia to use Jexpresso?","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Yes and no. It depends how much you are interested in adding your own equation set in the code rather than using it as a black box. ","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The following are useful resources about Julia:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Julia webpage docs.julialang.org\nOfficial list of learning resources julialang.org/learning","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Pages = [\n  \"Jexpresso.md\",\n  ]","category":"page"},{"location":"#Equations:","page":"Jexpresso.jl","title":"Equations:","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Jexpresso uses arbitrarily high-order (3rd and above) continuous spectral elements to solve","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"fracpartial bf qpartial t + sum_i=1^ndnablacdotbf F_i(bf q) = munabla^2bf q + bf S(bf q) + rm bc","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"where the vectors bf q, bf F, and bf S are problem-dependent as shown below, and are taken to be zero vectors of the appropriate size when not explicitly stated otherwise.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The Julia package DifferentialEquations.jl is used for time discretization and stepping.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"In order, we provide tests and results for the following equations:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"1D wave equation:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"bf q=beginbmatrix\nu \nv\nendbmatrixquad bf F=beginbmatrix\nv\nu\nendbmatrix","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"2: 1D shallow water:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"bf q=beginbmatrix\nh \nu\nendbmatrixquad bf F=beginbmatrix\nUh + Hu\ngh + Uu\nendbmatrix","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"where H and U are a reference height and velocity, respectively.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"2D Helmholtz:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"bf S=beginbmatrix\nalpha^2 u + f(xz)\nendbmatrixquad munabla^2bf q=mubeginbmatrix\nu_xx + u_zz\nendbmatrix","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"for a constant value of alpha and mu, which are case-dependent.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"2D scalar advection-diffusion:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"bf q=beginbmatrix\nq\nendbmatrixquad bf F=beginbmatrix\nqu\nendbmatrixquad bf F=beginbmatrix\nqv\nendbmatrixquad munabla^2bf q=mubeginbmatrix\nq_xx + q_zz\nendbmatrix","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"2D Euler equations of compressible flows with gravity and N passive chemicals c_i forall i=1N ","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"bf q=beginbmatrix\nrho \nrho u\nrho v\nrho theta\nrho c1\n\nrho cN\nendbmatrixquad bf F1=beginbmatrix\nrho u\nrho u^2 + p\nrho u v\nrho u theta\nrho u c1\n\nrho u cN\nendbmatrixquad bf F2=beginbmatrix\nrho v\nrho v u\nrho v^2 + p\nrho v theta\nrho v c1\n\nrho v cN\nendbmatrixquad bf S=beginbmatrix\n0\n0\n-rho g\n0\n0\n\n0\nendbmatrixquad munabla^2bf q=mubeginbmatrix\n0\nu_xx + u_zz\nv_xx + v_zz\ntheta_xx + theta_zz\nc1_xx + c1_zz\n\ncN_xx + cN_zz\nendbmatrix","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The equation of state for a perfect gas is used to close the system.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"3D Euler equations of compressible flows with gravity","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"bf q=beginbmatrix\nrho \nrho u\nrho v\nrho w\nrho theta\nendbmatrixquad bf F1=beginbmatrix\nrho u\nrho u^2 + p\nrho u v\nrho u w\nrho u theta\nendbmatrixquad bf F2=beginbmatrix\nrho v\nrho v u\nrho v^2 + p\nrho v w\nrho v theta\nendbmatrixquad bf F3=beginbmatrix\nrho w\nrho w u\nrho w v\nrho w^2 + p\nrho w theta\nendbmatrixquad bf S=beginbmatrix\n0\n0\n0\n-rho g\n0\nendbmatrixquad munabla^2bf q=mubeginbmatrix\n0\nu_xx + u_yy + u_zz\nv_xx + v_yy + v_zz\nw_xx + w_yy + w_zz\ntheta_xx + theta_yy + theta_zz\nendbmatrix","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"If you are interested in contributing, please get in touch: Simone Marras, Yassine Tissaoui","category":"page"},{"location":"#Some-notes-on-using-JEXPRESSO","page":"Jexpresso.jl","title":"Some notes on using JEXPRESSO","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"To install and run the code assume Julia 1.10.0","category":"page"},{"location":"#Setup-with-CPUs","page":"Jexpresso.jl","title":"Setup with CPUs","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":">> cd $JEXPRESSO_HOME\n>> julia --project=. -e \"using Pkg; Pkg.instantiate(); Pkg.API.precompile()\"","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"followed by the following:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Push problem name to ARGS You need to do this only when you run a new problem","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), EQUATIONS::String, EQUATIONS_CASE_NAME::String);\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"PROBLEMNAME is the name of your problem directory as JEXPRESSO/problems/equations/problemname\nPROBLEMCASENAME is the name of the subdirectory containing the specific setup that you want to run: ","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The path would look like  $JEXPRESSO/problems/equations/PROBLEM_NAME/PROBLEM_CASE_NAME","category":"page"},{"location":"#Tutorials","page":"Jexpresso.jl","title":"Tutorials","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The following tutorials will introduce you to the functionality of Jexpresso.jl.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Pages = [\n    \"features/performance.md\",\n    \"tutorials/theta.md\",\n    \"tutorials/theta.md\",\n    ]\nDepth = 2","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Example 1: to solve the 2D Euler equations with buyoancy and two passive tracers defined in problems/equations/CompEuler/thetaTracers you would do the following:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"CompEuler\", \"thetaTracers\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/thetaTracersMesh.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 5px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Example 2: to solve the 2D Euler equations leading to a density current defined in problems/equations/CompEuler/dc you would do the following:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"CompEuler\", \"dc\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/dc.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Example 3: to solve the 1D wave equation  defined in problems/equations/CompEuler/wave1d you would do the following:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"CompEuler\", \"wave1d\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/wave1d-v.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"For ready to run tests, there are the currently available equations names:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"CompEuler (option with total energy and theta formulation)","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The code is designed to create any system of conservsation laws. See CompEuler/case1 to see an example of each file. Details will be given in the documentation (still WIP). Write us if you need help.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"More are already implemented but currently only in individual branches. They will be added to master after proper testing.","category":"page"},{"location":"#Laguerre-semi-infinite-element-test-suite","page":"Jexpresso.jl","title":"Laguerre semi-infinite element test suite","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"This section contains instructions to run all of the test cases presented in","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"@article{tissaoui2024,\n  doi = {},\n  url = {},\n  year = {2020},\n  volume = {},\n  number = {},\n  pages = {},\n  author = {Yassine Tissaoui and James F. Kelly and Simone Marras}\n  title = {Efficient Spectral Element Method for the Euler Equations on Unbounded Domains in Multiple Dimensions},\n  journal = {arXiv},\n}","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Test 1: 1D wave equation with Laguerre semi-infinite element absorbing layers","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The problem is defined in problems/CompEuler/wave1d_lag and by default output will be written to output/CompEuler/wave1d_lag. To solve this problem run the following commands from the Julia command line:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"CompEuler\", \"wave1d_lag\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/wavev4.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Test 2: 1D wave train for linearized shallow water equations","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The problem is defined in problems/equations/AdvDiff/Wave_Train and by default output will be written to output/AdvDiff/Wave_Train. To solve this problem run the following commands from the Julia command line:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"AdvDiff\", \"Wave_Train\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/WaveTrainfinal.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"A second version of this tests generate images with the solutions at different times overlapped.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"This version is defined in problems/equations/AdvDiff/Wave_Train_Overlapping_Plot and by default output will be written to output/AdvDiff/Wave_Train_Overlapping_Plot. To run this version of the problem execute the following from the Julia command line:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"AdvDiff\", \"Wave_Train_Overlapping_Plot\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/WaveTrainoverlap.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Test 3: 2D advection-diffusion equation","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The problem is defined in problems/equations/AdvDiff/2D_laguerre and by default output will be written to output/AdvDiff/2D_laguerre. To solve this problem run the following commands from the Julia command line:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"AdvDiff\", \"2D_laguerre\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/ad2d-4s-line.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Test 4: 2D Helmholtz equation","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The problem is defined in problems/equations/Helmholtz/case1 and by default output will be written to output/Helmholtz/case1. To solve this problem run the following commands from the Julia command line:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"Helmholtz\", \"case1\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/Helmholtzfromjexpresso-line.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Test 5: Rising thermal bubble","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The problem is defined in problems/equations/CompEuler/theta_laguerre and by default output will be written to output/CompEuler/theta_laguerre. To solve this problem run the following commands from the Julia command line:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"CompEuler\", \"theta_laguerre\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/48.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Test 6: Hydrostatic linear mountain waves","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"The problem is defined in problems/equations/CompEuler/HSmount_Lag_working and by default output will be written to output/CompEuler/HSmount_Lag_working. To solve this problem run the following commands from the Julia command line:","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"julia> push!(empty!(ARGS), \"CompEuler\", \"HSmount_Lag_working\");\njulia> include(\"./src/Jexpresso.jl\")","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"<img src=\"assets/wvelo.png\"      alt=\"Markdown icon\"      style=\"float: left; margin-right: 7px;\" />","category":"page"},{"location":"#Plotting","page":"Jexpresso.jl","title":"Plotting","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Files can be written to VTK (recommended) or png. For the png plots, we use Makie. If you want to use a different package, modify ./src/io/plotting/jplots.jl accordinly.","category":"page"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"For non-periodic 2D tests, the output can also be written to VTK files by setting the value \"vtk\" for the usier_input key :outformat","category":"page"},{"location":"#Contacts","page":"Jexpresso.jl","title":"Contacts","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Simone Marras, Yassine Tissaoui","category":"page"},{"location":"#Manual","page":"Jexpresso.jl","title":"Manual","text":"","category":"section"},{"location":"","page":"Jexpresso.jl","title":"Jexpresso.jl","text":"Pages = [\"Jexpresso.md\"]","category":"page"},{"location":"Jexpresso/#Jexpresso","page":"Jexpresso","title":"Jexpresso","text":"","category":"section"},{"location":"Jexpresso/","page":"Jexpresso","title":"Jexpresso","text":"Jexpresso","category":"page"},{"location":"Jexpresso/#Jexpresso","page":"Jexpresso","title":"Jexpresso","text":"A research software for the numerical solution of a system of an arbitrary number of conservation laws using continuous spectral elements. DISCLAIMER: this is WIP and only 2D is being maintained until parallelization is complete.\n\nIf you are interested in contributing, please get in touch. Simone Marras, Yassine Tissaoui\n\n\n\n\n\n","category":"module"}]
}
