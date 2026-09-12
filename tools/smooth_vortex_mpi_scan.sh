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
#     sweep, so the time error is the same constant everywhere rather than
#     something that shrinks with h and contaminates the slope.
#
# The integrator is the deck's CarpenterKennedy2N54. Vern9 was the default
# here for a while, on the theory that a 4th-order integrator caps the
# measurable rate at 4 — MEASURED, at 32x32 elements and nop 4, it does not:
#
#   ck54,  dt = 1.0e-3     2.318e-5
#   vern9, dt = 1.0e-3     2.320e-5
#   vern9, dt = 3.3e-4     2.314e-5
#
# The error there is purely spatial: neither the integrator nor the step
# moves it, and Vern9's 16 stages per step cost 3x for nothing.
# SV_SOLVER=vern9 brings it back; the honest check on any sweep is to re-run
# its FINEST case at half SV_DT and see that the error does not move.
#
# THE ERROR STORE IS CLEARED FIRST (SV_KEEP=1 to add to it instead): every run
# redraws the figures from the whole store, so a leftover sweep would appear
# on the comparison of this one.
set -u
cd "$(dirname "$0")/.."

# WHICH CASE. The MHD vortex by default; SV_CASE=CompEuler/smoothVortex runs
# the hydrodynamic control (the classical isentropic vortex: same box, same
# meshes, same figures, no magnetic field), whose deck reads JEXPRESSO_EV_*
# instead of JEXPRESSO_SV_*.
CASE=${SV_CASE:-MHD/smoothVortex}
EQNS=${CASE%%/*}
CNAME=${CASE##*/}
if [ "${SV_PREFIX:-}" != "" ]; then PFX=$SV_PREFIX
elif [ "$EQNS" = "CompEuler" ];  then PFX=EV
else                                  PFX=SV
fi

NOPS=${SV_NOPS:-"1 2 3 4"}
NELX=${SV_NELX:-"4 8 16 32"}
VISCS=${SV_VISC:-"dsgs none"}
NP=${SV_NP:-4}
JOBS=${SV_JOBS:-1}
SOLVER=${SV_SOLVER:-ck54}
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

[ "${SV_KEEP:-0}" = "1" ] || rm -rf "problems/$CASE/errors"

SV_NELX="$NELX" tools/smooth_vortex_mesh.sh || true
for M in $NELX; do
    [ -f "problems/MHD/smoothVortex/vortex_${M}x${M}.msh" ] || {
        echo "MISSING problems/MHD/smoothVortex/vortex_${M}x${M}.msh — generate it with"
        echo "    SV_NELX=\"$NELX\" tools/smooth_vortex_mesh.sh   (needs gmsh)"
        exit 1
    }
done

# Does the launcher work at all? One 2-rank hello before committing hours to
# it: a wrong --mpi= flag or a launcher the cluster will not run inside a job
# step fails here, in one line, instead of eight times in eight log files.
if [ "$NP" -gt 1 ] && [ "${SV_SKIP_LAUNCHER_CHECK:-0}" != "1" ]; then
    echo "=== launcher check: $MPIEXEC -n 2 $JULIA -e 'using MPI; ...'"
    if ! $MPIEXEC -n 2 "$JULIA" --project=. -e \
            'using MPI; MPI.Init(); r=MPI.Comm_rank(MPI.COMM_WORLD); n=MPI.Comm_size(MPI.COMM_WORLD); println("    rank $r of $n on ", gethostname()); MPI.Finalize()' \
            2>&1 | sed 's/^/    /'; then
        echo "=== the launcher FAILED. Nothing else will run. Check:"
        echo "      srun --mpi=list                 (pmi2 for MPICH, pmix for OpenMPI)"
        echo "      MPIEXEC=\"mpirun\" SV_NP=$((NP*JOBS)) SV_JOBS=1   as a fallback"
        exit 1
    fi
fi

echo "=== $CASE: nops [$NOPS] x nelx [$NELX] x [$VISCS]"
echo "=== $NP rank(s) per case, $JOBS case(s) at a time, dt = $DT, solver $SOLVER, tend $TEND"

run_one() {   # $1 visc  $2 nop  $3 nelx
    log="logs/${CNAME}_${1}_nop$2_nelx$3.log"
    launcher=$([ "$NP" -gt 1 ] && echo "$MPIEXEC -n $NP")
    echo "--- START visc $1, nop $2, ${3}x${3} elements   $(date +%T)"
    echo "    $launcher $JULIA --project=. src/Jexpresso.jl $EQNS $CNAME"
    echo "    everything this case prints goes to $log"
    if [ "${DRYRUN:-0}" = "1" ]; then echo "    (dry run)"; return 0; fi
    env JEXPRESSO_${PFX}_NOP="$2" JEXPRESSO_${PFX}_NELX="$3" JEXPRESSO_${PFX}_VISC="$1" \
        JEXPRESSO_${PFX}_DT="$DT" JEXPRESSO_${PFX}_SOLVER="$SOLVER" JEXPRESSO_${PFX}_TEND="$TEND" \
        ${PLOT_NOPS:+JEXPRESSO_${PFX}_PLOT_NOPS="$PLOT_NOPS"} \
        $launcher "$JULIA" --project=. src/Jexpresso.jl "$EQNS" "$CNAME" \
        > "$log" 2>&1
    rc=$?
    err="problems/$CASE/errors/nop$2_nelx$3_$([ "$1" = none ] && echo galerkin || echo dsgs).dat"
    if [ "$rc" -ne 0 ]; then
        echo "--- FAILED (exit $rc) visc $1, nop $2, nelx $3   $(date +%T)"
        echo "    last lines of $log:"
        tail -n 15 "$log" | sed 's/^/      /'
    elif [ -f "$err" ]; then
        echo "--- done visc $1, nop $2, ${3}x${3}   $(date +%T)   error: $(tail -n 1 "$err")"
    else
        echo "--- done visc $1, nop $2, ${3}x${3}   $(date +%T)   BUT NO ERROR FILE at $err"
        echo "    last lines of $log:"
        tail -n 8 "$log" | sed 's/^/      /'
    fi
}

mkdir -p logs
# EVERY case is independent of every other, so the slots are filled from one
# flat list rather than from each (visc, mesh) group in turn: a sweep over two
# orders would otherwise never run more than two at a time, whatever SV_JOBS
# said, and half of a 64-core allocation would sit idle.
CASES=""
for V in $VISCS; do
    for M in $NELX; do
        for N in $NOPS; do
            CASES="$CASES $V:$N:$M"
        done
    done
done

for c in $CASES; do
    V=${c%%:*}; rest=${c#*:}; N=${rest%%:*}; M=${rest##*:}
    # keep SV_JOBS of them in flight
    while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do sleep 5; done
    run_one "$V" "$N" "$M" &
done
wait

echo "=== sweep finished $(date +%T). What is in the store:"
n_have=0; n_want=0
for V in $VISCS; do
    tag=$([ "$V" = none ] && echo galerkin || echo dsgs)
    for M in $NELX; do
        for N in $NOPS; do
            n_want=$((n_want + 1))
            f="problems/$CASE/errors/nop${N}_nelx${M}_${tag}.dat"
            if [ -f "$f" ]; then n_have=$((n_have + 1)); else echo "    MISSING $f"; fi
        done
    done
done
echo "=== $n_have of $n_want cases stored an error"
echo "=== done. The figures of the last run hold the whole sweep:"
echo "    \$PWD/output/$CASE/output/convergence_{dsgs,galerkin}[_<subset>]-it<n>.png"
echo "    (a case writes them only when it REACHES ITS FINAL TIME, so nothing"
echo "     appears until the first case completes; watch logs/sv_*.log meanwhile)"
echo "    (replot any subset without rerunning:"
echo "         julia --project=. tools/smooth_vortex_plot.jl 1,3 --case=$CASE)"
