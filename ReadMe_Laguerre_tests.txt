This file provides instructions for running our Laguerre semi-infinite element test suite.
The instructions are to be executed through the julia command line.

We first provide the instruction to load each specific test case:

-1D wave equation: push!(empty!(ARGS), "CompEuler", "wave1d_lag");
-1D wave train for linearized shallow water equations: push!(empty!(ARGS), "AdvDiff", "Wave_Train");
To generate a plot of solutions overlapped over time use push!(empty!(ARGS), "AdvDiff", "Wave_Train_Overlapping_Plot"); instead
-2D advection-diffusion equation: push!(empty!(ARGS), "AdvDiff", "2D_laguerre");
-2D Helmholtz equation: push!(empty!(ARGS), "Helmholtz", "case1");
-Rising thermal bubble: push!(empty!(ARGS), "CompEuler", "theta_laguerre");
-Hydrostatic Linear mountain waves: push!(empty!(ARGS), "CompEuler", "HSmount_Lag_working");

To run a case execute the following command in the julia command line: include("./src/Jexpresso.jl")

The output files will be saved to your Jexpresso directory at "Jexpresso/output/ProblemType/Case".
ProblemType is the first argument of the loading command and Case is the second argument.
For example the outputs for the rising thermal bubble will be found at "Jexpresso/output/CompEuler/theta_laguerre"

Note that output for all except the last two tests is generated as image files.
The last two tests provide their output as ".vtu" files which can be read with software such as "Paraview"
