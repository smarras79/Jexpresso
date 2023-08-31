#--------------------------------------------------------
# The problem name is a command line argument:
#
# 1. Launch Julia:
# >> julia --project=.
#
# 2. Push equations name to ARGS
#    You need this only when you run a new equations
#
#    julia > push!(empty!(ARGS), EQUATIONS::String, EQUATIONS_CASE_NAME::String);
#    julia > include(./src/Jexpresso.jl)
#
#    EQUATIONS is the name of your equations directory as $JEXPRESSO/src/equations/EQUATIONS
#    EQUATIONS_CASE_NAME is the name of the subdirectory $JEXPRESSO/src/equations/EQUATIONS_CASE_NAME
#
# Ex. To run the Compressible Euler equations in $JEXPRESSO/src/equations/CompEuler/theta
# 
#  julia > push!(empty!(ARGS), "CompEuler", "theta");
#  julia > include(./src/Jexpresso.jl)
#
#--------------------------------------------------------
function myparsing()
    if size(ARGS, 1) < 2
        println("")
        println(" ERROR: Missing command line arguments")
        println(" Usage example:")
        println("     julia> push!(empty!(ARGS), \"CompEuler\", \"theta\")")
        println("     julia> include(./src/Jexpresso.jl)")
        println(" -- See detailed instructions in the README file on github.")
        println("")
    else
        include("./io/mod_inputs.jl")

        #parsed_args  = parse_commandline()
        s = ArgParseSettings()
        parsed_args = parse_args(s)
        equations = string(parsed_args["arg1"])
        if (parsed_args["arg2"] === nothing)
            equations_case_name = ""
        else
            equations_case_name = string(parsed_args["arg2"])
        end
        equations_dir  = string("equations")
        driver_dir   = string("./", equations_dir, "/", equations, "/", equations_case_name, "/drivers.jl")
        include(driver_dir)
    end
end
#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
myparsing()
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
