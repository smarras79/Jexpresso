#--------------------------------------------------------
# external packages
#--------------------------------------------------------
using Crayons.Box
using DifferentialEquations
using Gridap
using GridapGmsh
using MPI
using PlotlyJS

#Constants
TInt   = Int64
TFloat = Float64

#--------------------------------------------------------
# jexpresso modules
#--------------------------------------------------------
include("./user_inputs.jl")
include("./mesh/mod_mesh.jl")
include("./solver/mod_solution.jl")
#--------------------------------------------------------

MPI.Init()
comm = MPI.COMM_WORLD

if MPI.Comm_rank(comm) == 0
    print(BLUE_FG("\n #------------------------------------------------------------\n"))
    print(BLUE_FG(" # Welcome to ", RED_FG("jExpresso !!!\n")))
    print(BLUE_FG(" # A Julia code to solve turbulence problems in the atmosphere\n"))
    print(BLUE_FG(" #------------------------------------------------------------\n\n"))

    
    #npx, npy, npz = 100, 100, 100    
    inputs = user_inputs()
#=
    equation_set, problem,
    nsd,
    npx,
    xmin, xmax,
    ymin, ymax,
    zmin, zmax = inputs[1],inputs[2], inputs[3], inputs[4], inputs[5], inputs[6], inputs[7], inputs[8], inputs[9], inputs[10]
=#
    println( " #--------------------------------------------------------------------------------\n")
    println( " # User inputs:\n # Equation set: ", equation_set)
    println( "\n # Problem:       ", problem)
    println( "\n # N. space dims: ", nsd)
    println( "\n N. x-points:     ", npx)
    println( "\n [xmin, xmax]:    ", xmin, " ", xmax)
    if (nsd > 1)
        println( "\n N. x-points:     ", npy)
        println( "\n [ymin, ymax]:    ", ymin, " ", ymax)
    end
    if (nsd == 3)
        println( "\n N. x-points:     ", npz)
        println( "\n [zmin, zmax]:    ", zmin, " ", zmax)
    end
    println( "\n # End user inputs.\n")
    println( " #--------------------------------------------------------------------------------\n")
    
    #=    mesh = St_mesh{TFloat, TInt}(zeros(npx), zeros(npy), zeros(npz),
    xmin, xmax, ymin, ymax, zmin, zmax,
                                 npx, npy, npz)
    build_mesh2d!(mesh)

    qs = St_solution{TInt, TFloat}(zeros(npx,nvars))
    #initial_conditions!(mesh,
    #                    qinit,
    #                    npx, npy, npz;
    #                    problem)
    =#
end

MPI.Barrier(comm)
MPI.Finalize()
