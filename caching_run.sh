#!/bin/bash

# Check if the first argument is "cache"
if [ "$1" == "cache" ]; then
    echo " # Clearing JLD2 cache before new caching run ..........."
    rm -f ./meshes/gmsh_grids/*.jld2
    echo " # Clearing JLD2 cache before new caching run ........... DONE"
    mpirun -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
else
    mpirun -np 2 ./AlyaProxy/Alya.x : -np 2 julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
fi



