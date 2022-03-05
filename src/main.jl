#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using DifferentialEquations
using Gridap
using GridapGmsh
using MPI
using Revise

#Plots
using Plots; gr()
plotlyjs()

#Constants
const TInt   = Int8
const TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("./IO/mod_inputs.jl")
include("./Mesh/mod_mesh.jl")
include("./solver/mod_solution.jl")
#--------------------------------------------------------

mod_inputs_print_welcome()

#--------------------------------------------------------
#User inputs:
#--------------------------------------------------------
inputs, nvars = mod_inputs_user_inputs()

#--------------------------------------------------------
#Print inputs to screen
#--------------------------------------------------------
mod_inputs_print(inputs; nvars)

#--------------------------------------------------------
# Build mesh    
#--------------------------------------------------------
#
# Initialize mesh struct
mesh = St_mesh{TInt,TFloat}(zeros(inputs[:npx]), zeros(inputs[:npy]), zeros(inputs[:npz]),
                            inputs[:xmin], inputs[:xmax],
                            inputs[:ymin], inputs[:ymax],
                            inputs[:zmin], inputs[:zmax],
                            inputs[:npx], inputs[:npy], inputs[:npz], 1, 1)

# Create mesh
mod_mesh_build_mesh!(mesh)
#--------------------------------------------------------
# END Build mesh    
#--------------------------------------------------------

#--------------------------------------------------------
# Initialize solution struct
#--------------------------------------------------------
#
# Initlaize solution struct
qsol = St_solution{TInt,TFloat}(zeros(TFloat, nvars, inputs[:npx]*inputs[:npy]*inputs[:npz]),
                                zeros(TFloat, nvars, inputs[:npx]*inputs[:npy]*inputs[:npz]),
                                zeros(TFloat, nvars, inputs[:npx]*inputs[:npy]*inputs[:npz]),
                                zeros(TFloat, nvars, inputs[:npx]*inputs[:npy]*inputs[:npz]),
                                zeros(TFloat, nvars, inputs[:npx]*inputs[:npy]*inputs[:npz]),
                                zeros(TFloat, inputs[:nsd]*inputs[:nsd]))

# Build initial conditions
mod_solution_initial_conditions!(mesh,
                                 qsol,
                                 inputs[:problem])



plt = plot(mesh.x, qsol.q[1,:], w = 3)
plot(scatter!(mesh.x, zeros(length(mesh.x)), x=:sepal_width, y=:sepal_length, mode="markers"))
savefig("~/Work/Codes/jexpresso/figs/initial_conditions.png")

#=mod_solution_solveODE!(mesh,
                       qsol,
                       inputs[:nsd],
                       inputs[:npx], inputs[:npy], inputs[:npz],
                       inputs[:problem])
=#

