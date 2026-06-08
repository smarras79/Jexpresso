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

    # Launch Julia in the background so we can animate while it starts
    $JULIA --project=. -e "
      using MPI
      run(\`\$(mpiexec()) -n $1 \$(Base.julia_cmd()) --project=. src/Jexpresso.jl $2 $3\`)" &
    local pid=$!

    # Ctrl-C should kill Julia and restore the cursor
    trap 'kill "$pid" 2>/dev/null; tput cnorm 2>/dev/null; printf "\r\033[K"; trap - INT; return 130' INT

    local frames='⠋⠙⠹⠸⠼⠴⠦⠧⠇⠏'
    local n=${#frames}
    local i=0
    tput civis 2>/dev/null                  # hide cursor
    while kill -0 "$pid" 2>/dev/null; do
        printf '\r\033[1;31m%s waiting…\033[0m\033[K' "${frames:i++%n:1}"
        sleep 0.1
    done
    tput cnorm 2>/dev/null                   # restore cursor
    printf '\r\033[K'                        # erase the spinner line
    trap - INT
    wait "$pid"                              # propagate Julia's exit status
}
jexp_mpich "$@"
}
