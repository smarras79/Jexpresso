using MPI
using Crayons.Box

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
