using Crayons.Box
using Gridap
using GridapGmsh
using MPI
using Mod_grid

MPI.Init()
comm = MPI.COMM_WORLD

if MPI.Comm_rank(comm) == 0
    print(WHITE_FG("\n #------------------------------------------------------------\n"))
    print(WHITE_FG(" # Welcome to ", YELLOW_FG("jExpresso !!!\n")))
    print(WHITE_FG(" # A Julia code to solve turbulence problems in the atmosphere\n"))
    print(WHITE_FG(" #------------------------------------------------------------\n\n"))

    
    
end

MPI.Barrier(comm)
MPI.Finalize()
