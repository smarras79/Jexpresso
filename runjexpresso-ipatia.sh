#!/bin/bash                                                                     
MPIRUN=/opt/local/bin/mpirun
JULIA= /Applications/Julia-1.12.app/Contents/Resources/julia/bin/julia
 
$MPIRUN -np $1 $JULIA --project=. -e 'push!(empty!(ARGS), "'"$2"'", "'"$3"'"); include("./src/Jexpresso.jl")' "$@"
