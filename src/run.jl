#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using DifferentialEquations
using Gridap
using GridapGmsh
using MPI
using Revise
using WriteVTK

#Plots
using Plots; gr()
plotlyjs()

#Constants
const TInt   = Int64
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("./test.jl")
include("./IO/mod_inputs.jl")
include("./Mesh/mod_mesh.jl")
include("./drivers/drivers.jl")
include("../tests/plot_lagrange_polynomial.jl")
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
driver(CG(),   # Space discretization type    
       AD1D(), # Equation subtype
       inputs, # input parameters from src/user_input.jl
       TFloat)

#test_driver(NSD_1D(),        # Number of Space Dimensions
#            INTERPOLATION(), # Problem Type
#            inputs,          # input parameters from src/user_input.jl
#            TFloat)
