#!/bin/bash

JULIA=$(command -v julia || "${SHELL:-/bin/zsh}" -ic "which julia" 2>/dev/null | tail -n1 | awk '{print $NF}')

if [[ ! -x "$JULIA" ]]; then
    echo "Error: julia not found or not executable: '$JULIA'" >&2
    exit 1
fi

"$JULIA" --project=. -e 'ENV["JULIA_PKG_PRECOMPILE_AUTO"]=0; using Pkg; Pkg.instantiate()'

"$JULIA" --project=. -e 'using Pkg; Pkg.precompile()'

# Optional AMR support (GridapP4est / P4est_wrapper, installed into envs/amr/;
# the patched GridapP4est fork is pinned there). Skipped on Windows, where
# P4est_wrapper does not build, and when JEXPRESSO_AMR=0.
case "$(uname -s)" in
    MINGW*|MSYS*|CYGWIN*) JEXPRESSO_AMR=0 ;;
esac
if [[ "${JEXPRESSO_AMR:-1}" != "0" ]]; then
    "$JULIA" --project=. tools/setup_amr.jl
else
    echo "Skipping AMR setup (run 'julia --project=. tools/setup_amr.jl' to add it later)."
fi
