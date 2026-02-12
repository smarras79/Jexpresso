
#
# Run coupled Jexpresso + Alya as:
#
mpirun --tag-output -np 2 ./alya/Alya_enhanced.x : -np 2 julia --project=. ./src/Jexpresso-mini-coupled.jl

#
# Or equivalently, running Jexpresso.jl directly with problem args:
#
mpirun --tag-output -np 2 ./alya/Alya_enhanced.x : -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
