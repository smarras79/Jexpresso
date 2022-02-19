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
    
    println( " #--------------------------------------------------------------------------------\n")
    println( " # User inputs:\n # Equation set: ", inputs.equation_set)
    println( " # Problem:       ", inputs.problem)
    println( " # N. space dims: ", inputs.nsd)
    println( " # N. x-points:     ", inputs.npx)
    println( " # [xmin, xmax]:    ", inputs.xmin, " ", inputs.xmax)
    if (inputs.nsd > 1)
        println( " # N. x-points:     ", inputs.npy)
        println( " # [ymin, ymax]:    ", inputs.ymin, " ", inputs.ymax)
    end
    if (inputs.nsd == 3)
        println( " # N. x-points:     ", inputs.npz)
        println( " # [zmin, zmax]:    ", inputs.zmin, " ", inputs.zmax)
    end
    println( " # End user inputs.\n")
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
