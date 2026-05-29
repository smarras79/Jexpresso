
#
# Run coupled Jexpresso + Alya directly:
#

#2D
mpirun -np 2 ./AlyaProxy/Alya.x : -x JEXPRESSO_COUPLED=1 -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya

#3D
mpirun -np 2 ./AlyaProxy/Alya.x : -x JEXPRESSO_COUPLED=1 -np 2 julia --project=. ./src/Jexpresso.jl CompEuler 3dAlya
