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
#   BW_HSCALE=ngl tools/brio_wu_order_scan.sh   # the kernels' Δ_K/(k+1) instead
#   BW_NODAL=1 tools/brio_wu_order_scan.sh      # the nodal coefficient
#
# THE STORE IS CLEARED FIRST. Every run redraws the figure from every curve
# in problems/MHD/brioWu1d/curves, so a leftover sweep would appear on the
# comparison of this one from its very first run. Set BW_KEEP=1 to accumulate
# instead, to finish a partial sweep.
set -u
cd "$(dirname "$0")/.."

DOFS=${BW_DOFS:-"150 300 600 1200"}
NOPS=${BW_NOPS:-"4 5 6 7"}
CMIN=${BW_CMIN:-""}
HSCALE=${BW_HSCALE:-nop}
NODAL=${BW_NODAL:-""}
JULIA=${JULIA:-julia}

[ "${BW_KEEP:-0}" = "1" ] || rm -rf problems/MHD/brioWu1d/curves

echo "=== Brio-Wu order scan: nops [$NOPS] x DOFs [$DOFS], length scale \
$([ "$HSCALE" = nop ] && echo 'Delta_K/k' || echo 'Delta_K/(k+1)')\
${NODAL:+, nodal coefficient}${CMIN:+, Cmin $CMIN}"

for D in $DOFS; do
    for N in $NOPS; do
        echo "=== Brio-Wu: nop $N, ~$D DOFs ${CMIN:+(Cmin $CMIN)}"
        env JEXPRESSO_BW_NOP="$N" JEXPRESSO_BW_DOFS="$D" \
            JEXPRESSO_BW_HSCALE="$HSCALE" \
            ${NODAL:+JEXPRESSO_BW_NODAL="$NODAL"} \
            ${CMIN:+JEXPRESSO_BW_CMIN="$CMIN"} \
            "$JULIA" --project=. src/Jexpresso.jl MHD brioWu1d || exit 1
    done
done
echo "=== done; the figures of the last run hold the whole sweep"
