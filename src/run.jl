using Dates
using Revise

const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# USER INPUT ARGUMENT:
#--------------------------------------------------------
problem_name = "advdiff"
#--------------------------------------------------------
# END USER INPUT ARGUMENT:
#--------------------------------------------------------


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
if !isempty(ARGS)
    empty!(ARGS)
    push!(ARGS, problem_name)
else
    push!(ARGS, problem_name)
end

include("../src/io/mod_inputs.jl")
parsed_args  = parse_commandline()
problem_name = string(parsed_args["arg1"])
problem_dir  = string("problems")
driver_dir   = string("./", problem_dir, "/", problem_name, "/drivers.jl")
include(driver_dir)

#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
mod_inputs_print_welcome()
inputs        = Dict{}()
inputs, nvars = mod_inputs_user_inputs!(problem_name, problem_dir)

#--------------------------------------------------------
#Create output directory if it doesn't exist:
#--------------------------------------------------------
OUTPUT_DIR = string(dirname(@__DIR__()), "/src/", problem_dir, "/", problem_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
if !isdir(OUTPUT_DIR)
    mkdir(OUTPUT_DIR)
end

#--------------------------------------------------------
# Problem setup
# !!!!!!
# !!!!!! WARNING: MOVE all the setup parameters to user_input.jl
# !!!!!!
#--------------------------------------------------------
driver(ContGal(),   # Space discretization type    
       inputs, # input parameters from src/user_input.jl
       OUTPUT_DIR,
       TFloat)
