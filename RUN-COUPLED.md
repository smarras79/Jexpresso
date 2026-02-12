
#
# Run coupled Jexpresso + Alya directly:
#
mpirun --tag-output -np 2 ./alya/Alya_enhanced.x : -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya

#
# Alternatively, using the mini-coupled wrapper (for backward compatibility):
#
mpirun --tag-output -np 2 ./alya/Alya_enhanced.x : -np 2 julia --project=. ./src/Jexpresso-mini-coupled.jl
