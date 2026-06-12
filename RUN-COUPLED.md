
#
# Run coupled Jexpresso + Alya directly:
#
#2D
mpirun -np 2 ./AlyaProxy/Alya.x : -x JEXPRESSO_COUPLED=1 -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya

#3D
mpirun -np 2 ./AlyaProxy/Alya.x \
     : -np 2 -x JEXPRESSO_COUPLED=1 "$JULIA_BIN" --project=. ./src/Jexpresso.jl CompEuler 3dAlya

#If you get a `prterun` error, then try this:

JULIA_BIN=$(julia -e 'print(joinpath(Sys.BINDIR, "julia"))')
echo "$JULIA_BIN"   # sanity check, e.g. ~/.julia/juliaup/julia-1.11.x+.../bin/julia

mpirun -np 2 ./AlyaProxy/Alya.x \
     : -np 2 -x JEXPRESSO_COUPLED=1 "$JULIA_BIN" --project=. ./src/Jexpresso.jl CompEuler 3dAlya


#
# WITH MPICH on OSX
#
julia --project=. -e '
  using MPI
  run(`$(mpiexec()) -n 2 ./AlyaProxy/Alya.x : -n 2 -env JEXPRESSO_COUPLED 1 $(Base.julia_cmd()) --project=. src/Jexpresso.jl CompEuler 3dAlya`)'