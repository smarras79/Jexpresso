#!/bin/bash

# USER: change to your julia path: ---------------------------------------
JULIA=/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia
# END USER ---------------------------------------------------------------

jexp_mpich() {
    local msg="I am starting Julia; please be patient. Jexpresso hasn't started yet!"
    local border
    border=$(printf '═%.0s' $(seq 1 $(( ${#msg} + 2 ))))
    printf '\033[1;31m╔%s╗\n║ %s ║\n╚%s╝\033[0m\n' "$border" "$msg" "$border"

    # Defaults — silence the heartbeat and the alloc summary unless the
    # caller already set them. The `: "${VAR:=default}"` pattern only
    # assigns when VAR is unset or empty, so a command-line override like
    #   JEXPRESSO_STEP_HEARTBEAT=1 ./jexp_mpich.sh 4 CompEuler city2d
    # still wins.
    : "${JEXPRESSO_STEP_HEARTBEAT:=0}"
    : "${JEXPRESSO_ALLOC_SUMMARY:=0}"
    : "${JEXPRESSO_PRECOMPILE_WARMUP:=1}"
    export JEXPRESSO_STEP_HEARTBEAT JEXPRESSO_ALLOC_SUMMARY JEXPRESSO_PRECOMPILE_WARMUP

    $JULIA --project=. -e "
      using MPI
      run(\`\$(mpiexec()) -n $1 \$(Base.julia_cmd()) --project=. src/Jexpresso.jl $2 $3\`)"
}

jexp_mpich "$@"
