#!/bin/bash

JULIA=/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia
jexp_mpich() {
    $JULIA --project=. -e "
      using MPI
      run(\`\$(mpiexec()) -n $1 \$(Base.julia_cmd()) --project=. src/Jexpresso.jl $2 $3\`)"
}
jexp_mpich "$@"
