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
include("./IO/mod_inputs.jl")
include("./Mesh/mod_mesh.jl")
#include("./basis/basis_structs.jl")
include("./drivers/cg_driver.jl")
#include("./Infrastructure/Kopriva_functions.jl")
#include("./Infrastructure/2D_3D_structures.jl")
#include("./element_matrices.jl")
include("../tests/plot_lagrange_polynomial.jl")
#--------------------------------------------------------

#MPI.Init()
#comm = MPI.COMM_WORLD

#if MPI.Comm_rank(comm) == 0
mod_inputs_print_welcome()

#--------------------------------------------------------
#User inputs:
#--------------------------------------------------------
inputs        = Dict{}()
inputs, nvars = mod_inputs_user_inputs()

#--------------------------------------------------------
# Create/read mesh
# return mesh::St_mesh
#--------------------------------------------------------
mesh = mod_mesh_mesh_driver(inputs)

#--------------------------------------------------------
# Problem setup
# !!!!!!
# !!!!!! WARNING: MOVE all the setup parameters to user_input.jl
# !!!!!!
#--------------------------------------------------------
cg_driver(AD1D(),
          inputs[:lexact_integration],
          inputs[:nop],
          TFloat)
