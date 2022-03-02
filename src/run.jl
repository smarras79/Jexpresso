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
include("./basis/basis_structs.jl")
include("./solver/mod_solution.jl")
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
# Build mesh    
#--------------------------------------------------------

if (haskey(inputs, :lread_gmsh) && inputs[:lread_gmsh]==true)
    
    println(" # Read gmsh grid")
    # Initialize mesh struct: the arrays length will be increased in mod_mesh_read_gmsh
    mesh = St_mesh{TInt,TFloat}(zeros(2), zeros(2), zeros(2),
                                0.0, 0.0,
                                0.0, 0.0,
                                0.0, 0.0,
                                1, 1, 1, 1, 1,
                                inputs[:nop])
    
    # Read gmsh grid using the GridapGmsh reader
    mod_mesh_read_gmsh!(mesh, inputs[:gmsh_filename])
    
    println(" # Read gmsh grid ........................ DONE")
else
    println(" # Build grid")
    # Initialize mesh struct
    mesh = St_mesh{TInt,TFloat}(zeros(inputs[:npx]), zeros(inputs[:npy]), zeros(inputs[:npz]),
                                inputs[:xmin], inputs[:xmax],
                                inputs[:ymin], inputs[:ymax],
                                inputs[:zmin], inputs[:zmax],
                                inputs[:npx], inputs[:npy], inputs[:npz], 1, 1,
                                inputs[:nop])
    
    mod_mesh_build_mesh!(mesh)
    
    println(" # Build grid ........................ DONE")
end

#--------------------------------------------------------
# END Build mesh    
#--------------------------------------------------------

#--------------------------------------------------------
# Initialize solution struct
#--------------------------------------------------------
#
#= Initlaize solution struct
qsol = St_solution{TInt,TFloat}(zeros(TFloat, nvars, mesh.npx*mesh.npy*mesh.npz),
                                zeros(TFloat, nvars, mesh.npx*mesh.npy*mesh.npz),
                                zeros(TFloat, nvars, mesh.npx*mesh.npy*mesh.npz),
                                zeros(TFloat, nvars, mesh.npx*mesh.npy*mesh.npz),
                                zeros(TFloat, nvars, mesh.npx*mesh.npy*mesh.npz),
                                zeros(TFloat, mesh.nsd*mesh.nsd))

# Build initial conditions
mod_solution_initial_conditions!(mesh,
                                 qsol,
                                 inputs[:problem])

=#

Legendre = St_legendre{TFloat}(0.0, 0.0, 0.0, 0.0)
lgl      = St_lgl{TFloat}(zeros(TFloat, mesh.nop+1),
                     zeros(TFloat, mesh.nop+1))
build_lgl!(Legendre, lgl, mesh.nop)

#
# Plot to file:
#
#plt = plot(mesh.x, qsol.q[1,:], w = 3)
#plot(scatter!(mesh.x, zeros(length(mesh.x)), x=:sepal_width, y=:sepal_length, mode="markers"))
#savefig("~/Work/Codes/jexpresso/figs/initial_conditions.png")
