#--------------------------------------------------------
# The problem name is a command line argument:
#
# 1. Launch Julia:
# >> julia --project=.
#
# 2. Push equations name to ARGS
#    You need this only when you run a new equations
#
#    julia > push!(empty!(ARGS), EQUATIONS::String);
#    julia > include(./src/run.jl)
#
# EQUATIONS is the name of your equations directory
# as $JEXPRESSO/src/equations/equations
#
# Ex. If you run the Advection Diffusion equations in $JEXPRESSO/src/equations/AdvDiff
# 
#  julia > push!(empty!(ARGS), "AdvDiff");
#  julia > include(./src/run.jl)
#
#--------------------------------------------------------
if isempty(ARGS)
    s = """
            
            Please, run the following every time that EQUATIONS changes:
                julia> push!(empty!(ARGS), EQUATIONS::String, , EQUATIONS_CASE_NAME::String);

            and only then run jexpresso with:
                julia> include(./src/run.jl)

            Currently avaiable EQUATIONS options:
            - Elliptic
            - AdvDiff

            EQUATIONS_CASE_NAME is user defined and must be the name of the case directory inside EQUATIONS:
                For example, if "AdvDiff" contains a directory called "Case1", you would do the following:
                    julia> push!(empty!(ARGS), "AdvDiff", "Case1");
                but if you don't have any CASE DIRECTORY inside "AdvDiff", then you simply 
                    julia> push!(empty!(ARGS), "AdvDiff");
        """
    error(s)
end

include("./io/mod_inputs.jl")
parsed_args  = parse_commandline()
equations = string(parsed_args["arg1"])
if (parsed_args["arg2"] === nothing)
    equations_case_name = ""
else
    equations_case_name = string(parsed_args["arg2"])
end
equations_dir  = string("equations")
driver_dir   = string("./", equations_dir, "/", equations, "/", equations_case_name, "/drivers.jl")
include(driver_dir)

#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
mod_inputs_print_welcome()
inputs        = Dict{}()
inputs        = mod_inputs_user_inputs!(equations, equations_case_name, equations_dir)

#--------------------------------------------------------
#Create output directory if it doesn't exist:
#--------------------------------------------------------
user_defined_output_dir = inputs[:output_dir]
if isempty(user_defined_output_dir)
    OUTPUT_DIR = string(dirname(@__DIR__()), "/src/", equations_dir, "/", equations, "/", equations_case_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
else
    @info user_defined_output_dir
    OUTPUT_DIR = string(dirname(user_defined_output_dir), "/", equations, "/", equations_case_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
end
if !isdir(OUTPUT_DIR)
    mkpath(OUTPUT_DIR)
end

#--------------------------------------------------------
# Equations setup
# !!!!!!
# !!!!!! WARNING: MOVE all the setup parameters to user_input.jl
# !!!!!!
#--------------------------------------------------------
driver(ContGal(),   # Space discretization type    
       inputs, # input parameters from src/user_input.jl
       OUTPUT_DIR,
       TFloat)
