#!/usr/bin/env bash
#
# Convergence history of the Brio-Wu shock tube: every polynomial order at
# every resolution, which is what problems/MHD/brioWu1d/user_plot.jl turns
# into density-it<n>.png, convergence-it<n>.png and
# convergence_smooth-it<n>.png (Dao & Nazarov 2022, Fig. 1 and Fig. 2).
#
#   tools/brio_wu_order_scan.sh                 # orders 4-7, 150-1200 DOFs
#   BW_DOFS="300 600" BW_NOPS="4 7" tools/brio_wu_order_scan.sh
#   BW_CMIN=0.0 tools/brio_wu_order_scan.sh     # with the DynSGS floor off
#
# The store accumulates, so a partial sweep can be completed later; delete
# problems/MHD/brioWu1d/curves to start a fresh one.
set -u
cd "$(dirname "$0")/.."

DOFS=${BW_DOFS:-"150 300 600 1200"}
NOPS=${BW_NOPS:-"4 5 6 7"}
CMIN=${BW_CMIN:-""}
JULIA=${JULIA:-julia}

[ "${BW_FRESH:-0}" = "1" ] && rm -rf problems/MHD/brioWu1d/curves

for D in $DOFS; do
    for N in $NOPS; do
        echo "=== Brio-Wu: nop $N, ~$D DOFs ${CMIN:+(Cmin $CMIN)}"
        env JEXPRESSO_BW_NOP="$N" JEXPRESSO_BW_DOFS="$D" \
            ${CMIN:+JEXPRESSO_BW_CMIN="$CMIN"} \
            "$JULIA" --project=. src/Jexpresso.jl MHD brioWu1d || exit 1
    done
done
echo "=== done; the figures of the last run hold the whole sweep"
