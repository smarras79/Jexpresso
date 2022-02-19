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
    print(BLUE_FG(" #--------------------------------------------------------------------------------\n"))
    print(BLUE_FG(" # Welcome to ", RED_FG("jexpresso\n")))
    print(BLUE_FG(" # A Julia code to solve turbulence problems in the atmosphere\n"))

    #Load user inputs (NamedTuple defined in user_inputs.jl
    inputs = user_inputs()

    #=
    NOTE: SM
    ADD HERE A FUNCTION TO CHECK IF SOME NECESSARY INPUTS WERE NOT DEFINED
    THINK OF A WAY TO CREATE A DYNAMIC INPUT SETUP to make the inputs list flexible.
    =#
    println( " #--------------------------------------------------------------------------------")
    print(GREEN_FG(" # User inputs:\n"))
    println( " # Equation set: ", inputs.equation_set)
    println( " # Problem:       ", inputs.problem)
    println( " # N. space dims: ", inputs.nsd)
    println( " # N. x-points:   ", inputs.npx)
    println( " # [xmin, xmax]:  ", inputs.xmin, " ", inputs.xmax)
    if (inputs.nsd > 1)
        println( " # N. y-points:   ", inputs.npy)
        println( " # [ymin, ymax]:  ", inputs.ymin, " ", inputs.ymax)
    end
    if (inputs.nsd == 3)
        println( " # N. z-points:   ", inputs.npz)
        println( " # [zmin, zmax]:  ", inputs.zmin, " ", inputs.zmax)
    end

    #
    # Build mesh
    #
    mesh = St_mesh{TFloat, TInt}(zeros(npx), zeros(npy), zeros(npz),
                                 inputs.xmin, inputs.xmax,
                                 inputs.ymin, inputs.ymax,
                                 inputs.zmin, inputs.zmax,
                                 inputs.npx, inputs.npy, inputs.npz)
    
    build_mesh2d!(mesh)
    #
    # END Build mesh
    #

    qs = St_solution{TInt, TFloat}(zeros(npx,nvars))
    #initial_conditions!(mesh,
    #                    qinit,
    #                    npx, npy, npz;
    #                    problem)
    
end

MPI.Barrier(comm)
MPI.Finalize()
