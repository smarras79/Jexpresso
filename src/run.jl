const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("../src/io/mod_inputs.jl")

#--------------------------------------------------------
# The problem name is a command line argument:
#
# >> julia --project=. ./src/run.jl problem_name
#
# The problem name must be the name of your problem directory
# as $JEXPRESSO/src/problems/problem_name
#
# Ex.
# If you build the Advection Diffusion problem in $JEXPRESSO/src/problems/AdvDiff
# then `AdvDiff` should be your problem name that is passed to jexpresso from
# the command line.
#--------------------------------------------------------
parsed_args  = parse_commandline()
problem_name = string(parsed_args["arg1"])
problem_dir  = string("../src/problems/")
driver_dir   = string(problem_dir, problem_name, "/drivers.jl")
include(driver_dir)
#--------------------------------------------------------
mod_inputs_print_welcome()

#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
inputs        = Dict{}()
inputs, nvars = mod_inputs_user_inputs!(problem_name, problem_dir)

#--------------------------------------------------------
# Problem setup
# !!!!!!
# !!!!!! WARNING: MOVE all the setup parameters to user_input.jl
# !!!!!!
#--------------------------------------------------------
driver(CG(),   # Space discretization type    
       inputs, # input parameters from src/user_input.jl
       TFloat)
