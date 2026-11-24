#!/usr/bin/env bash
#
# Convergence history of the smooth MHD vortex (Dao & Nazarov 2022, §5.1.1
# and Fig. 1): every polynomial order on every mesh, which
# problems/MHD/smoothVortex/user_plot.jl turns into
# convergence_dsgs-it<n>.png and convergence_galerkin-it<n>.png.
#
#   tools/smooth_vortex_scan.sh                    # orders 4-7, meshes 4-32
#   SV_NOPS="3 4" SV_NELX="8 16 32" tools/smooth_vortex_scan.sh
#   SV_VISC=dsgs tools/smooth_vortex_scan.sh       # only the DynSGS panel
#
# BOTH panels of the paper's Fig. 1 are run by default: the DynSGS one and the
# plain Galerkin one (SV_VISC="dsgs none"), since the point of the test is the
# comparison — the residual viscosity must not cost the high-order solution its
# accuracy where the solution is smooth.
#
# THE ERROR STORE IS CLEARED FIRST. Every run redraws the convergence figure
# from every error in problems/MHD/smoothVortex/errors, so a leftover sweep
# (a different set of orders, other coefficients) would appear on the figure
# of this one from its very first run. Set SV_KEEP=1 to accumulate instead,
# to finish a partial sweep.
set -u
cd "$(dirname "$0")/.."

NOPS=${SV_NOPS:-"4 5 6 7"}
NELX=${SV_NELX:-"4 8 16 32"}
VISCS=${SV_VISC:-"dsgs none"}
JULIA=/Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia

[ "${SV_KEEP:-0}" = "1" ] || rm -rf problems/MHD/smoothVortex/errors

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
