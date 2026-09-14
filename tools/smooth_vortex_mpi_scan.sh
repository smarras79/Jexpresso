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
# The box WIDTH. The default 20 is [-10,10]^2, the box Dao & Nazarov run this
# accuracy test on. SV_L=10 is [-5,5]^2, where the vortex tail at the periodic
# seam floors every order at ~1e-5 (see the header of
# tools/smooth_vortex_mesh.sh).
LBOX=${SV_L:-20}
LTAG="L$(printf '%s' "$LBOX" | sed 's/\.0*$//')_"

# One Δt for the sweep: the deck's rule at the finest (nelx, nop) of it.
# 2.0e-3 * 64 / (nelx*nop) is that rule; keep the two in step.
if [ "${SV_DT:-}" = "auto" ] || [ "${SV_DT:-}" = "percase" ]; then
    # Δt PER CASE, from the deck's own CFL rule (DT_REF*64/(nelx*nop)): the
    # step is refined WITH the mesh instead of being held fixed for the whole
    # sweep. That is the study to run when the residual viscosity is on: nu
    # is C_R h^2 R, and at a fixed step R stops falling once it reaches the
    # time-discretization error, so nu stalls and every RV curve eventually
    # bends to slope h^2. A step that shrinks with h keeps that floor falling.
    DT=percase
elif [ -n "${SV_DT:-}" ]; then
    DT=$SV_DT
else
    MAXN=0; for N in $NOPS;  do [ "$N" -gt "$MAXN" ] && MAXN=$N; done
    MAXM=0; for M in $NELX; do [ "$M" -gt "$MAXM" ] && MAXM=$M; done
    DT=$(awk -v n="$MAXN" -v m="$MAXM" 'BEGIN{printf "%.6g", 2.0e-3*64.0/(n*m)}')
fi

# A sweep is self-contained: it starts from an empty store unless SV_KEEP=1
# says to add to what is there. The old store is MOVED ASIDE, not deleted —
# a sweep that took hours should not be lost to a mistyped command — and a
# dry run touches nothing at all.
if [ "${SV_KEEP:-0}" != "1" ] && [ "${DRYRUN:-0}" != "1" ] && [ -d "problems/$CASE/errors" ]; then
    rm -rf "problems/$CASE/errors.prev"
    mv "problems/$CASE/errors" "problems/$CASE/errors.prev"
    echo "=== the previous store is in problems/$CASE/errors.prev"
fi

SV_NELX="$NELX" SV_L="$LBOX" tools/smooth_vortex_mesh.sh || true
for M in $NELX; do
    [ -f "problems/MHD/smoothVortex/vortex_${LTAG}${M}x${M}.msh" ] || {
        echo "MISSING problems/MHD/smoothVortex/vortex_${LTAG}${M}x${M}.msh — generate it with"
        echo "    SV_NELX=\"$NELX\" SV_L=$LBOX tools/smooth_vortex_mesh.sh   (needs gmsh)"
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

echo "=== $CASE: nops [$NOPS] x nelx [$NELX] x [$VISCS] on the [-$(awk -v l="$LBOX" 'BEGIN{printf "%g", l/2}'), $(awk -v l="$LBOX" 'BEGIN{printf "%g", l/2}')]^2 box"
echo "=== $NP rank(s) per case, $JOBS case(s) at a time, dt = $([ "$DT" = percase ] && echo "per case (the deck's CFL rule)" || echo "$DT"), solver $SOLVER, tend $TEND"

run_one() {   # $1 visc  $2 nop  $3 nelx
    log="logs/${CNAME}_${1}_nop$2_nelx$3.log"
    launcher=$([ "$NP" -gt 1 ] && echo "$MPIEXEC -n $NP")
    # SV_KEEP=1 adds to the store instead of starting a fresh one, so a case
    # whose record is already there has nothing to add: skip it. That makes a
    # sweep RESUMABLE — a machine that goes down, or a job that hits its wall
    # clock, costs only the cases that were in flight, not the whole run.
    if [ "${SV_KEEP:-0}" = "1" ] && [ "${SV_REDO:-0}" != "1" ]; then
        _e=$(err_file "$2" "$3" "$(_tag_of "$1")")
        if [ -f "$_e" ]; then
            echo "--- have visc $1, nop $2, ${3}x${3} already: $(tail -n 1 "$_e")"
            return 0
        fi
    fi
    echo "--- START visc $1, nop $2, ${3}x${3} elements   $(date +%T)"
    echo "    $launcher $JULIA --project=. src/Jexpresso.jl $EQNS $CNAME"
    echo "    everything this case prints goes to $log"
    if [ "${DRYRUN:-0}" = "1" ]; then echo "    (dry run)"; return 0; fi
    env JEXPRESSO_${PFX}_NOP="$2" JEXPRESSO_${PFX}_NELX="$3" JEXPRESSO_${PFX}_VISC="$1" \
        $([ "$DT" = percase ] || echo JEXPRESSO_${PFX}_DT="$DT") \
        JEXPRESSO_${PFX}_SOLVER="$SOLVER" JEXPRESSO_${PFX}_TEND="$TEND" \
        JEXPRESSO_${PFX}_L="$LBOX" \
        ${SV_CUTOFF:+JEXPRESSO_${PFX}_CUTOFF="$SV_CUTOFF"} \
        ${PLOT_NOPS:+JEXPRESSO_${PFX}_PLOT_NOPS="$PLOT_NOPS"} \
        $launcher "$JULIA" --project=. src/Jexpresso.jl "$EQNS" "$CNAME" \
        > "$log" 2>&1
    rc=$?
    err=$(err_file "$2" "$3" "$(_tag_of "$1")")
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

# Where a case stores its error. The Euler deck tags the record with the vortex
# strength (nop4_nelx8_b1_dsgs.dat), the MHD deck does not (nop4_nelx8_dsgs.dat),
# so match either and report the plain name when nothing is there yet.
# What the deck will call this run's record: a run with the smoothness cutoff
# on is its own experiment and stores its own curve (SV_CUTOFF > 0).
_tag_of() {
    [ "$1" = none ] && { echo galerkin; return 0; }
    case "${SV_CUTOFF:-0}" in 0|0.0|"") echo dsgs ;; *) echo dsgs_cut ;; esac
}

err_file() {
    _plain="problems/$CASE/errors/nop$1_nelx$2_$3.dat"
    for _f in "$_plain" problems/"$CASE"/errors/nop"$1"_nelx"$2"_b*_"$3".dat; do
        [ -f "$_f" ] && { echo "$_f"; return 0; }
    done
    echo "$_plain"
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
    tag=$(_tag_of "$V")
    for M in $NELX; do
        for N in $NOPS; do
            n_want=$((n_want + 1))
            f=$(err_file "$N" "$M" "$tag")
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
