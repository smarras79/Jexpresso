using ArgParse

#using Profile
#using PProf

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
function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "eqs"
        help = "equations"
        default = "CompEuler"
        required = false
        
        "eqs_case"
        help = "case name in equations directory"
        default = "theta"
        required = false
        
    end

    return parse_args(s)
end

#--------------------------------------------------------
#Parse command line args:
#--------------------------------------------------------
parsed_args         = parse_commandline()
equations           = string(parsed_args["eqs"])
equations_case_name = string(parsed_args["eqs_case"])
equations_dir       = string("equations")

driver_dir          = string(dirname(@__DIR__()), "/src/", equations_dir, "/", equations, "/", equations_case_name, "/drivers.jl")
user_flux_dir       = string(dirname(@__DIR__()), "/src/", equations_dir, "/", equations, "/", equations_case_name, "/user_flux.jl")
user_source_dir     = string(dirname(@__DIR__()), "/src/", equations_dir, "/", equations, "/", equations_case_name, "/user_source.jl")
user_bc_dir         = string(dirname(@__DIR__()), "/src/", equations_dir, "/", equations, "/", equations_case_name, "/user_bc.jl")

include(driver_dir)
include(user_flux_dir)
include(user_source_dir)
include(user_bc_dir)

#--------------------------------------------------------
#Read User Inputs:
#--------------------------------------------------------
mod_inputs_print_welcome()
inputs = Dict{}()
inputs = mod_inputs_user_inputs!(equations, equations_case_name, equations_dir)

#--------------------------------------------------------
#Create output directory if it doesn't exist:
#--------------------------------------------------------
user_defined_output_dir = inputs[:output_dir]
@info user_defined_output_dir
if user_defined_output_dir == "none"
    OUTPUT_DIR = string(dirname(@__DIR__()), "/src/", equations_dir, "/", equations, "/", equations_case_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
else
    OUTPUT_DIR = string(user_defined_output_dir, "/", equations, "/", equations_case_name, "/output-",  Dates.format(now(), "dduyyyy-HHMMSS/"))
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
#Profile.clear()

driver(ContGal(),   # Space discretization type    
       inputs, # input parameters from src/user_input.jl
       OUTPUT_DIR,
       TFloat)

# Export pprof profile and open interactive profiling web interface.
#pprof()
