#!/bin/bash                                                                     
MPIRUN=/app/mpich-3.1/bin/mpirun
JULIA=../julia-1.12.5/bin/julia
 
$MPIRUN -np $1 $JULIA --project=. -e 'push!(empty!(ARGS), "'"$2"'", "'"$3"'"); include("./src/Jexpresso.jl")' "$@"
