#!/usr/bin/env julia
#
# Launcher script for Jexpresso with proper argument handling
#
# Usage:
#   julia --project=. run_jexpresso.jl [OPTIONS] EQS EQS_CASE
#
# Example:
#   julia --project=. run_jexpresso.jl CompEuler theta
#
# With MPI:
#   mpirun -np 4 julia --project=. run_jexpresso.jl CompEuler theta
#
# With coupling:
#   mpirun -np 4 julia --project=. run_jexpresso.jl --coupling --code-id 1 --n-codes 2 CompEuler theta
#

# Arguments are already in ARGS from command line
# Just load the Jexpresso module, which will execute run.jl
using Jexpresso
