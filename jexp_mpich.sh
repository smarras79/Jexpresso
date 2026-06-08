#!/bin/bash

# USER: change to your julia path: ---------------------------------------
JULIA=/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia
# END USER ---------------------------------------------------------------

#!/bin/bash
jexp_mpich() {
    local msg="I am starting Julia; please be patient. Jexpresso hasn't started yet!"
    local border
    border=$(printf '═%.0s' $(seq 1 $(( ${#msg} + 2 ))))
    printf '\033[1;31m╔%s╗\n║ %s ║\n╚%s╝\033[0m\n' "$border" "$msg" "$border"

    $JULIA --project=. -e "
      using MPI
      run(\`\$(mpiexec()) -n $1 \$(Base.julia_cmd()) --project=. src/Jexpresso.jl $2 $3\`)"
}
jexp_mpich "$@"

