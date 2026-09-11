#!/usr/bin/env bash
#
# Convergence history of the smooth MHD vortex on MULTIPLE CORES: every
# polynomial order on every mesh, each case run under MPI, which
# problems/MHD/smoothVortex/user_plot.jl turns into the panels of
# Dao & Nazarov (2022), Fig. 1.
#
# Two kinds of parallelism, and this script does both:
#
#   SV_NP   ranks per case (MPI, one simulation split over SV_NP cores) —
#           what a big grid needs, since it is one run that must fit and
#           finish;
#   SV_JOBS cases run at the same time (they are independent), so a sweep of
#           small cases fills the machine as well.
#
# The product SV_NP × SV_JOBS is what you are actually asking the machine for;
# keep it at or below the cores you have.
#
#   tools/smooth_vortex_mpi_scan.sh                       # 4 ranks, orders 1-4
#   SV_NP=8 SV_NELX="8 16 32 64" tools/smooth_vortex_mpi_scan.sh
#   SV_NP=16 SV_JOBS=4 SV_NOPS="1 3" tools/smooth_vortex_mpi_scan.sh
#   SV_NP=1 SV_JOBS=8 tools/smooth_vortex_mpi_scan.sh     # serial cases, 8 at once
#
# THE TIME ERROR. CarpenterKennedy2N54 is fourth order, so with Δt ∝ h the
# measured rate saturates at 4 whatever the polynomial order (see the note by
# _sv_solver in the deck). This script therefore
#   * runs every case at ONE Δt, computed from the finest (nelx, nop) of the
#     sweep, so the time error is the same constant everywhere, and
#   * uses Vern9 (ninth order) by default, so that constant is far below the
#     spatial error.
# SV_DT and SV_SOLVER override both; SV_SOLVER=ck54 restores the low-storage
# integrator of the production decks.
#
# THE ERROR STORE IS CLEARED FIRST (SV_KEEP=1 to add to it instead): every run
# redraws the figures from the whole store, so a leftover sweep would appear
# on the comparison of this one.
set -u
cd "$(dirname "$0")/.."

NOPS=${SV_NOPS:-"1 2 3 4"}
NELX=${SV_NELX:-"4 8 16 32"}
VISCS=${SV_VISC:-"dsgs none"}
NP=${SV_NP:-4}
JOBS=${SV_JOBS:-1}
SOLVER=${SV_SOLVER:-vern9}
TEND=${SV_TEND:-1.0}
JULIA=${JULIA:-julia}
MPIEXEC=${MPIEXEC:-mpiexec}
PLOT_NOPS=${SV_PLOT_NOPS:-}

# One Δt for the sweep: the deck's rule at the finest (nelx, nop) of it.
# 2.0e-3 * 64 / (nelx*nop) is that rule; keep the two in step.
if [ -n "${SV_DT:-}" ]; then
    DT=$SV_DT
else
    MAXN=0; for N in $NOPS;  do [ "$N" -gt "$MAXN" ] && MAXN=$N; done
    MAXM=0; for M in $NELX; do [ "$M" -gt "$MAXM" ] && MAXM=$M; done
    DT=$(awk -v n="$MAXN" -v m="$MAXM" 'BEGIN{printf "%.6g", 2.0e-3*64.0/(n*m)}')
fi

[ "${SV_KEEP:-0}" = "1" ] || rm -rf problems/MHD/smoothVortex/errors

SV_NELX="$NELX" tools/smooth_vortex_mesh.sh

echo "=== smooth vortex: nops [$NOPS] x nelx [$NELX] x [$VISCS]"
echo "=== $NP rank(s) per case, $JOBS case(s) at a time, dt = $DT, solver $SOLVER, tend $TEND"

run_one() {   # $1 visc  $2 nop  $3 nelx
    echo "--- visc $1, nop $2, ${3}x${3} elements   $(date +%T)"
    env JEXPRESSO_SV_NOP="$2" JEXPRESSO_SV_NELX="$3" JEXPRESSO_SV_VISC="$1" \
        JEXPRESSO_SV_DT="$DT" JEXPRESSO_SV_SOLVER="$SOLVER" JEXPRESSO_SV_TEND="$TEND" \
        ${PLOT_NOPS:+JEXPRESSO_SV_PLOT_NOPS="$PLOT_NOPS"} \
        $([ "$NP" -gt 1 ] && echo "$MPIEXEC -n $NP") \
        "$JULIA" --project=. src/Jexpresso.jl MHD smoothVortex \
        > "logs/sv_${1}_nop$2_nelx$3.log" 2>&1 \
        || echo "    FAILED: visc $1, nop $2, nelx $3 — see logs/sv_${1}_nop$2_nelx$3.log"
}

mkdir -p logs
# The runs of one (visc, mesh) group are independent; JOBS of them at a time.
for V in $VISCS; do
    for M in $NELX; do
        n=0
        for N in $NOPS; do
            run_one "$V" "$N" "$M" &
            n=$((n + 1))
            if [ "$n" -ge "$JOBS" ]; then wait; n=0; fi
        done
        wait
    done
done

echo "=== done. The figures of the last run hold the whole sweep:"
echo "    output/MHD/smoothVortex/output/convergence_{dsgs,galerkin}[-<subset>]-it<n>.png"
echo "    (replot any subset without rerunning: tools/smooth_vortex_plot.jl 1,3)"
