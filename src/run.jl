using Dates
using Revise

const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# The problem name is a command line argument:
#
# 1. Launch Julia:
# >> julia --project=.
#
# 2. Push problem name to ARGS
#    You need this only when you run a new problem
#
#    julia > push!(empty!(ARGS), PROBLEM_NAME::String);
#    julia > include(./src/run.jl)
#
# PROBLEM_NAME is the name of your problem directory
# as $JEXPRESSO/src/problems/problem_name
#
# Ex. If you run the Advection Diffusion problem in $JEXPRESSO/src/problems/AdvDiff
# 
#  julia > push!(empty!(ARGS), "AdvDiff");
#  julia > include(./src/run.jl)
#
#--------------------------------------------------------
if isempty(ARGS)
    s = """
            
            Please, run the following every time that PROBLEM_NAME changes:
                julia> push!(empty!(ARGS), PROBLEM_NAME::String, , PROBLEM_CASE_NAME::String);

            and only then run jexpresso with:
                julia> include(./src/run.jl)

            Currently avaiable PROBLEM_NAME options:
            - Elliptic
            - AdvDiff

            PROBLEM_CASE_NAME is user defined and must be the name of the case directory inside PROBLEM_NAME:
                For example, if "AdvDiff" contains a directory called "Case1", you would do the following:
                    julia> push!(empty!(ARGS), "AdvDiff", "Case1");
                but if you don't have any CASE DIRECTORY inside "AdvDiff", then you simply 
                    julia> push!(empty!(ARGS), "AdvDiff");
        """
    error(s)
end

include("../src/io/mod_inputs.jl")
parsed_args  = parse_commandline()
problem_name = string(parsed_args["arg1"])
if (parsed_args["arg2"] === nothing)
    problem_case_name = ""
else
    problem_case_name = string(parsed_args["arg2"])
end
problem_dir  = string("problems")
driver_dir   = string("./", problem_dir, "/", problem_name, "/", problem_case_name, "/drivers.jl")
include(driver_dir)

#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
mod_inputs_print_welcome()
inputs        = Dict{}()
inputs        = mod_inputs_user_inputs!(problem_name, problem_case_name, problem_dir)

#--------------------------------------------------------
#Create output directory if it doesn't exist:
#--------------------------------------------------------
user_defined_output_dir = inputs[:output_dir]
if isempty(user_defined_output_dir)
    OUTPUT_DIR = string(dirname(@__DIR__()), "/src/", problem_dir, "/", problem_name, "/", problem_case_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
else
    @info user_defined_output_dir
    OUTPUT_DIR = string(dirname(user_defined_output_dir), "/", problem_name, "/", problem_case_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
end
if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
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
