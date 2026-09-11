#!/usr/bin/env bash
#
# Convergence history of the smooth MHD vortex (Dao & Nazarov 2022, §5.1.1
# and Fig. 1): every polynomial order on every mesh, which
# problems/MHD/smoothVortex/user_plot.jl turns into
# convergence_dsgs-it<n>.png and convergence_galerkin-it<n>.png.
#
#   tools/smooth_vortex_scan.sh                    # orders 2-4, meshes 4-32
#   SV_NOPS="3 4" SV_NELX="8 16 32" tools/smooth_vortex_scan.sh
#   SV_VISC="dsgs none" tools/smooth_vortex_scan.sh    # both panels of Fig. 1
#
# The errors accumulate, so a partial sweep can be completed later; delete
# problems/MHD/smoothVortex/errors to start a fresh one.
set -u
cd "$(dirname "$0")/.."

NOPS=${SV_NOPS:-"2 3 4"}
NELX=${SV_NELX:-"4 8 16 32"}
VISCS=${SV_VISC:-"dsgs"}
JULIA=${JULIA:-julia}

[ "${SV_FRESH:-0}" = "1" ] && rm -rf problems/MHD/smoothVortex/errors

SV_NELX="$NELX" tools/smooth_vortex_mesh.sh

for V in $VISCS; do
    for N in $NOPS; do
        for M in $NELX; do
            echo "=== smooth vortex: nop $N, ${M}x${M} elements, visc $V"
            env JEXPRESSO_SV_NOP="$N" JEXPRESSO_SV_NELX="$M" JEXPRESSO_SV_VISC="$V" \
                "$JULIA" --project=. src/Jexpresso.jl MHD smoothVortex || exit 1
        done
    done
done
echo "=== done; the figures of the last run hold the whole sweep"
