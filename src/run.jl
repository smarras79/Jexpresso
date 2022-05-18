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
include("./basis/basis_structs.jl")
include("./solver/mod_solution.jl")
#--------------------------------------------------------

struct EDGES <:At_geo_entity end
struct FACES <:At_geo_entity end
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
    
    println(" # Read gmsh grid and populate with high-order points ")
    
    # Initialize mesh struct: the arrays length will be increased in mod_mesh_read_gmsh
    mesh = St_mesh{TInt,TFloat}(nsd=Int8(inputs[:nsd]),
                                nop=Int8(inputs[:nop]))
    
    # Read gmsh grid using the GridapGmsh reader
    
    mod_mesh_read_gmsh!(mesh, inputs[:gmsh_filename])
    
    println(" # Read gmsh grid and populate with high-order points ........................ DONE")
else
    println(" # Build grid")
    
    # Initialize mesh struct for native structured grid:
    mesh = St_mesh{TInt,TFloat}(x = zeros(Int8(inputs[:npx])),
                                y = zeros(Int8(inputs[:npy])),
                                z = zeros(Int8(inputs[:npz])),
                                npx  = Int8(inputs[:npx]),
                                npy  = Int8(inputs[:npy]),
                                npz  = Int8(inputs[:npz]), 
                                xmin = Float64(inputs[:xmin]), xmax = Float64(inputs[:xmax]),
                                ymin = Float64(inputs[:ymin]), ymax = Float64(inputs[:ymax]),
                                zmin = Float64(inputs[:zmin]), zmax = Float64(inputs[:zmax]),
                                nsd=Int8(inputs[:nsd]),
                                nop=Int8(inputs[:nop]))

    
    #@info mesh
    

    
    #Write structured grid to VTK
    vtkfile = vtk_grid("mySTRUCTURED_GRID", mesh.x, mesh.y, mesh.z) # 3-D
    outfiles = vtk_save(vtkfile)
    
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


#
# Plot to file:
#
#plt = plot(mesh.x, qsol.q[1,:], w = 3)
#plot(scatter!(mesh.x, zeros(length(mesh.x)), x=:sepal_width, y=:sepal_length, mode="markers"))
#savefig("~/Work/Codes/jexpresso/figs/initial_conditions.png")
