#!/usr/bin/env bash
#
# THE BDF-k STUDY of the DynSGS residual sensor on the two smooth vortices
# (problems/MHD/smoothVortex and its hydrodynamic control
# problems/CompEuler/smoothVortex): is the excess of DynSGS over Galerkin from
# N = 5 on the floor of the sensor's BDF2 time difference, -(dt^2/3) d3q/dt3,
# as problems/MHD/setup_notes/dsgsCutoff.pdf found on the Alfven wave? With
# :dsgs_bdf_order = k that floor is -dt^k/(k+1) d^(k+1)q/dt^(k+1) instead
# (README of the MHD case, "The order of the sensor's time difference").
#
# For each case, through tools/smooth_vortex_mpi_scan.sh (so the records, tags
# and figures are the ones the sweep driver always writes):
#
#   dtcheck  the decisive one: the finest mesh at dt and at dt/2, Galerkin and
#            DynSGS without the cutoff (BDF2). A BDF2 time-difference floor
#            falls 4x when dt halves; a spatial residual does not move.
#            The dt/2 records are kept apart, in errors.dt2/ beside the store.
#   sweep    Galerkin and DynSGS without the cutoff with BDF2 (the sensor as it
#            has always been), then with BDF5 and BDF4: every order of SV_NOPS
#            on every mesh of SV_NELX, all at ONE dt (the driver's rule at the
#            finest point), records dsgs_nocut[_bdf<k>] (MHD) and
#            dsgs[_bdf<k>] (Euler, whose deck has no cutoff).
#
# and then prints, per case, every error against Galerkin's and draws the
# convergence figures of Galerkin, BDF2 and the highest BDF order alone.
#
#   tools/smooth_vortex_bdf_study.sh                                  # both cases, everything
#   SV_NP=1 SV_JOBS=8 tools/smooth_vortex_bdf_study.sh                # 8 serial cases at a time
#   STUDY_STEPS=dtcheck tools/smooth_vortex_bdf_study.sh              # the dt check alone (12 runs)
#   STUDY_CASES=MHD/smoothVortex STUDY_BDFS=5 tools/smooth_vortex_bdf_study.sh
#   SV_NELX="8 16 32 64" SV_NP=4 SV_JOBS=4 tools/smooth_vortex_bdf_study.sh   # a 64^2 point too
#   DRYRUN=1 tools/smooth_vortex_bdf_study.sh                         # print the plan only
#
# | variable          | default                                  |
# |-------------------|------------------------------------------|
# | STUDY_CASES       | MHD/smoothVortex CompEuler/smoothVortex  |
# | STUDY_STEPS       | dtcheck sweep                            |
# | STUDY_BDFS        | 5 4      (BDF2 always runs)              |
# | STUDY_CHECK_NOPS  | 5 6 7    (orders of the dt check)        |
# | SV_NOPS, SV_NELX  | 4 5 6 7, 4 8 16 32                       |
# | SV_NP, SV_JOBS    | 1, 4     (ranks per case, cases at once) |
# | SV_TEND           | 1.0                                      |
# | JULIA, MPIEXEC, JEXPRESSO_EV_BETA (5)  passed through            |
#
# ONE (mesh, order) AT A TIME. The first case of a (mesh, order) writes its SEM
# cache (.jexpresso_cache/*.jld2) and the others read it, memory-mapped; one
# that reads while another writes dies with a bus error (seen here, two cases
# started together on a mesh with no cache yet). The driver lists its cases
# variant by variant, so each call below runs ONE variant: every (mesh, order)
# is in it once, and the first call builds the caches the later ones read —
# whatever SV_JOBS is.
#
# A store already there is MOVED ASIDE first, never deleted
# (problems/<case>/errors.before-bdf-study-<date>), and so is an earlier
# errors.dt2: a study starts from an empty store, which every figure of it is
# drawn from.
set -u
cd "$(dirname "$0")/.."

CASES=${STUDY_CASES:-"MHD/smoothVortex CompEuler/smoothVortex"}
STEPS=${STUDY_STEPS:-"dtcheck sweep"}
BDFS=${STUDY_BDFS:-"5 4"}
CHECK_NOPS=${STUDY_CHECK_NOPS:-"5 6 7"}
NOPS=${SV_NOPS:-"4 5 6 7"}
NELX=${SV_NELX:-"4 8 16 32"}
export SV_NP=${SV_NP:-1} SV_JOBS=${SV_JOBS:-4} SV_TEND=${SV_TEND:-1.0}
DRY=${DRYRUN:-0}

# One dt for everything: the driver's rule, 2e-3*64/(nelx*nop), at the finest
# point of the sweep; the dt check halves it.
MAXN=0; for N in $NOPS; do [ "$N" -gt "$MAXN" ] && MAXN=$N; done
MAXM=0; for M in $NELX; do [ "$M" -gt "$MAXM" ] && MAXM=$M; done
DT=$(awk -v n="$MAXN" -v m="$MAXM" 'BEGIN{printf "%.17g", 2.0e-3*64.0/(n*m)}')
DT2=$(awk -v d="$DT" 'BEGIN{printf "%.17g", d/2}')
STAMP=$(date +%Y%m%d-%H%M%S)

# The julia the figures are drawn with: as the sweep driver finds it (JULIA,
# else the macOS app's, else the PATH's).
JULIA_BIN=${JULIA:-}
if [ -z "$JULIA_BIN" ]; then
    for j in /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia "$(command -v julia 2>/dev/null)"; do
        [ -n "$j" ] && "$j" --version >/dev/null 2>&1 && { JULIA_BIN=$j; break; }
    done
fi
JULIA_BIN=${JULIA_BIN:-julia}
mkdir -p logs

echo "=== BDF-k study: cases [$CASES], steps [$STEPS], BDF [2 $BDFS]"
echo "=== nops [$NOPS] x nelx [$NELX], dt = $DT (dt check at $DT2 on ${MAXM}^2, nops [$CHECK_NOPS])"
echo "=== $SV_NP rank(s) per case, $SV_JOBS case(s) at a time, tend $SV_TEND"

scan() {   # the sweep driver, adding to the store; arguments are VAR=value pairs
    if [ "$DRY" = 1 ]; then echo "    would run: $* tools/smooth_vortex_mpi_scan.sh"; return 0; fi
    env SV_KEEP=1 SV_CUTOFF=0 "$@" tools/smooth_vortex_mpi_scan.sh
}

# what the decks call the DynSGS-without-cutoff record (the Euler deck has no
# cutoff, so its run is plain "dsgs")
dsgs_tag() { case "$1" in MHD/*) echo dsgs_nocut ;; *) echo dsgs ;; esac; }

# The last line of a record: "L1 L2 Linf". The Euler records carry the vortex
# strength in their name (nop5_nelx32_b5_dsgs.dat).
err_of() {   # $1 store, $2 nop, $3 nelx, $4 tag  ->  L1 (or nothing)
    for f in "$1/nop$2_nelx$3_$4.dat" "$1"/nop"$2"_nelx"$3"_b*_"$4".dat; do
        [ -f "$f" ] && { tail -n 1 "$f" | awk '{print $1}'; return 0; }
    done
}

num() { [ -n "$1" ] && printf '%.4e' "$1" || printf -- '-'; }

table() {   # $1 case, $2 store: every L1 error against Galerkin's
    local case=$1 st=$2 d=$(dsgs_tag "$1")
    echo "=== $case: velocity L1 error at t = $SV_TEND and its ratio to Galerkin's (store $st)"
    printf '    %3s %5s %12s' N nelx Galerkin
    for k in 2 $BDFS; do printf ' %12s %7s' "BDF$k" ratio; done; echo
    for N in $NOPS; do for M in $NELX; do
        g=$(err_of "$st" "$N" "$M" galerkin)
        printf '    %3s %5s %12s' "$N" "$M" "$(num "$g")"
        for k in 2 $BDFS; do
            t=$d; [ "$k" = 2 ] || t="${d}_bdf$k"
            e=$(err_of "$st" "$N" "$M" "$t")
            r=$( [ -n "$e" ] && [ -n "$g" ] && awk -v e="$e" -v g="$g" 'BEGIN{printf "%.3f", e/g}')
            printf ' %12s %7s' "$(num "$e")" "${r:--}"
        done; echo
    done; done
}

for CASE in $CASES; do
    ST="problems/$CASE/errors"
    echo; echo "################ $CASE"
    if [ "$DRY" != 1 ]; then
        [ -d "$ST" ]     && mv "$ST"     "$ST.before-bdf-study-$STAMP"     && echo "=== the store was moved to $ST.before-bdf-study-$STAMP"
        [ -d "$ST.dt2" ] && mv "$ST.dt2" "$ST.dt2.before-bdf-study-$STAMP" && echo "=== errors.dt2 was moved to $ST.dt2.before-bdf-study-$STAMP"
    fi
    D=$(dsgs_tag "$CASE")

    if [[ " $STEPS " == *" dtcheck "* ]]; then
        echo "=== dt check: ${MAXM}^2, nops [$CHECK_NOPS], Galerkin and DynSGS-BDF2 at dt/2 = $DT2"
        # one variant per driver call: see "ONE (mesh, order) AT A TIME" above
        scan SV_CASE="$CASE" SV_NELX="$MAXM" SV_NOPS="$CHECK_NOPS" SV_VISC=none SV_DT="$DT2"
        scan SV_CASE="$CASE" SV_NELX="$MAXM" SV_NOPS="$CHECK_NOPS" SV_VISC=dsgs SV_DT="$DT2"
        if [ "$DRY" != 1 ]; then
            mkdir -p "$ST.dt2" && mv "$ST"/*.dat "$ST.dt2"/ 2>/dev/null
            echo "=== the dt/2 records are in $ST.dt2"
        fi
        # ... and the same at dt, the reference of the check (the sweep would
        # run these too; they stay in the store for it)
        scan SV_CASE="$CASE" SV_NELX="$MAXM" SV_NOPS="$CHECK_NOPS" SV_VISC="dsgs none" SV_DT="$DT"
    fi

    if [[ " $STEPS " == *" sweep "* ]]; then
        echo "=== sweep: Galerkin"
        scan SV_CASE="$CASE" SV_NELX="$NELX" SV_NOPS="$NOPS" SV_VISC=none SV_DT="$DT"
        echo "=== sweep: DynSGS without the cutoff, BDF2"
        scan SV_CASE="$CASE" SV_NELX="$NELX" SV_NOPS="$NOPS" SV_VISC=dsgs SV_DT="$DT"
        for k in $BDFS; do
            echo "=== sweep: DynSGS without the cutoff, BDF$k"
            scan SV_CASE="$CASE" SV_NELX="$NELX" SV_NOPS="$NOPS" SV_VISC=dsgs SV_BDF="$k" SV_DT="$DT"
        done
    fi
    [ "$DRY" = 1 ] && continue

    echo
    table "$CASE" "$ST"
    if [ -d "$ST.dt2" ]; then
        echo "=== dt check on ${MAXM}^2: L1 at dt and dt/2; the DynSGS excess over Galerkin falls 4x if it is the BDF2 floor"
        printf '    %3s %12s %12s %12s %12s %10s %10s\n' N "Gal dt" "Gal dt/2" "BDF2 dt" "BDF2 dt/2" "excess dt" "dt/2"
        for N in $CHECK_NOPS; do
            g1=$(err_of "$ST" "$N" "$MAXM" galerkin);  g2=$(err_of "$ST.dt2" "$N" "$MAXM" galerkin)
            d1=$(err_of "$ST" "$N" "$MAXM" "$D");      d2=$(err_of "$ST.dt2" "$N" "$MAXM" "$D")
            awk -v n="$N" -v g1="${g1:-nan}" -v g2="${g2:-nan}" -v d1="${d1:-nan}" -v d2="${d2:-nan}" \
                'BEGIN{printf "    %3s %12.4e %12.4e %12.4e %12.4e %10.3e %10.3e\n", n, g1, g2, d1, d2, d1-g1, d2-g2}'
        done
    fi
    # the comparison figures: Galerkin, BDF2 and the highest BDF order run
    top=0; for k in $BDFS; do [ "$k" -gt "$top" ] && top=$k; done
    out="output/$CASE/output/bdf_study"
    vis="galerkin,$D"; [ "$top" -gt 0 ] && vis="$vis,${D}_bdf$top"
    if "$JULIA_BIN" --project=. tools/smooth_vortex_plot.jl "$(echo $NOPS | tr ' ' ',')" --case="$CASE" \
           --visc="$vis" --out="$out" > "logs/bdf_study_plot_${CASE%%/*}.log" 2>&1; then
        echo "=== figures ($vis): $PWD/$out/convergence[_L1|_L2|_Linf]-it0.png"
    else
        echo "=== the figures failed: see logs/bdf_study_plot_${CASE%%/*}.log"
    fi
done
