
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

Two things to watch with the MPICH call above:
The bundled launcher is MPICH's hydra, so the Open MPI flag -x VAR=VAL doesn't exist — use -env VAR VAL (or simply export JEXPRESSO_COUPLED=1 before launching; Alya ignores it, and Jexpresso reads it in src/run.jl).
Alya.x must be recompiled against that same MPICH. An Alya.x built with system Open MPI's mpif90 will fail to join MPI_COMM_WORLD under the bundled MPICH mpiexec (the PMI handshake between the two implementations is incompatible). If recompiling Alya against MPICH is inconvenient, use Option A.
One caveat on the in-repo helpers: RUN-COUPLED.md still documents the old system-mpirun recipe (your original command), and run_coupled.sh is stale — it references run_jexpresso.jl and scripts/make_sysimage.jl, neither of which exists in the current tree — so don't rely on that script as-is. jexp_mpich.sh only launches Jexpresso standalone, but Option B above is exactly its launch pattern extended with the Alya half of the MPMD command line.