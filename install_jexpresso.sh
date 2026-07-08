#!/bin/bash

JULIA=$(command -v julia || "${SHELL:-/bin/zsh}" -ic "which julia" 2>/dev/null | tail -n1 | awk '{print $NF}')

if [[ ! -x "$JULIA" ]]; then
    echo "Error: julia not found or not executable: '$JULIA'" >&2
    exit 1
fi

if [[ "$(uname)" == "Darwin" ]]; then
    #echo "$(uname)"
    if [ ! -d "$HOME/GridapP4est.jl" ]; then
        git clone -b arm64-cfunction-fix git@github.com:Hwang1229/GridapP4est.jl.git ~/GridapP4est.jl
    fi
    "$JULIA" --project=. -e 'using Pkg; Pkg.develop(path=expanduser("~/GridapP4est.jl"))'
fi

"$JULIA" --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'

"$JULIA" --project=. -e 'using Pkg; Pkg.precompile()'
