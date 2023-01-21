const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../src/io/mod_inputs.jl")
include("../src/problems/AdvDiff/drivers.jl") #automate this based on input

#--------------------------------------------------------
return 
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
PROBLEM_EQUATIONS = inputs[:problem]

driver(CG(),   # Space discretization type    
       PROBLEM_EQUATIONS, # Equation subtype
       inputs, # input parameters from src/user_input.jl
       TFloat)
