#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using Gridap
using GridapGmsh
using MPI
using Revise
using WriteVTK

#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../src/io/mod_inputs.jl")
include("../src/problems/AdvDiff/drivers.jl") #automate this based on input
#--------------------------------------------------------

#MPI.Init()
#comm = MPI.COMM_WORLD

#if MPI.Comm_rank(comm) == 0
mod_inputs_print_welcome()

#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
inputs        = Dict{}()
inputs, nvars = mod_inputs_user_inputs()

#--------------------------------------------------------
# Problem setup
# !!!!!!
# !!!!!! WARNING: MOVE all the setup parameters to user_input.jl
# !!!!!!
#--------------------------------------------------------
PROBLEM_EQUATIONS = Adv2D()
#PROBLEM_EQUATIONS = Wave1D()
driver(CG(),   # Space discretization type    
       PROBLEM_EQUATIONS, # Equation subtype
       inputs, # input parameters from src/user_input.jl
       TFloat)
