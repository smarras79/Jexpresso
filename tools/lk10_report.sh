#!/bin/bash
# tools/lk10_report.sh — everything a diagnosis of a run of problems/MHD/fluxEmergenceLeake2010 needs, in ONE
# text file to send. Run it from the top of the repository, on the machine that holds the output.
#
#   bash tools/lk10_report.sh RUN [--base BASE] [--log LOGFILE]... [--py PYTHON] [--times T1,T2,...]
#
#   RUN      the output directory of the run (output/MHD/fluxEmergenceLeake2010/output-<date>);
#            DIR1,DIR2 for a run made in several jobs
#   --base   a run whose FIRST snapshot is the unit of mass of the tables and the pressure they are compared
#            with, so that the tables of several restarts can be compared: the run from t = 0, or the earliest
#            run of the chain that is still on disk. Without it, the first snapshot of RUN itself.
#   --log    the job's standard output (the file the queue wrote); may be given several times
#   --py     the Python to use (default python3; it needs numpy and matplotlib, as the figures do)
#   --times  the times [s] of the vertical profiles (default 1380,1610,1955,2185: those the run has)
#
# It writes lk10_report_<last directory of RUN>.txt in the current directory and reads the snapshots
# only: it changes nothing in the run. A few minutes on a full REF 4 run.
set -u
RUN=""; BASE=""; PY=python3; TIMES="1380,1610,1955,2185"; LOGS=()
while [ $# -gt 0 ]; do
  case "$1" in
    --base)  BASE="$2"; shift 2;;
    --log)   LOGS+=("$2"); shift 2;;
    --py)    PY="$2"; shift 2;;
    --times) TIMES="$2"; shift 2;;
    -h|--help) sed -n 2,18p "$0"; exit 0;;
    *) if [ -z "$RUN" ]; then RUN="$1"; shift; else echo "unknown argument: $1" >&2; exit 2; fi;;
  esac
done
[ -n "$RUN" ] || { sed -n 2,18p "$0"; exit 2; }
[ -f tools/plot_leake2010.py ] || { echo "run this from the top of the repository (tools/plot_leake2010.py not found here)" >&2; exit 2; }

LAST=${RUN##*,}; LAST=${LAST%%@*}; LAST=${LAST%/}
OUT="lk10_report_$(basename "$LAST").txt"
FIGS="lk10_report_files"; mkdir -p "$FIGS"
REF=(); [ -n "$BASE" ] && REF=(--ref "$BASE")
P="$PY tools/plot_leake2010.py"

sec() { printf '\n===== %s =====\n' "$1"; }
run() { echo "\$ $*"; "$@" 2>&1; echo "[exit $?]"; }

{
  echo "lk10_report: $(date -u '+%Y-%m-%d %H:%M UTC')"
  echo "source: $(git log -1 --format='%h %ad %s' --date=short 2>/dev/null || echo 'not a git checkout')"
  echo "run: $RUN"
  echo "base: ${BASE:-none given}"
  echo "python: $($PY --version 2>&1)"

  sec "configuration of each directory (lk10_case.txt) and its snapshots"
  for spec in ${RUN//,/ } ${BASE:+$BASE}; do
    d=${spec%%@*}
    echo "--- $d: $(ls "$d" 2>/dev/null | grep -c '\.pvtu$') snapshots"
    if [ -f "$d/lk10_case.txt" ]; then cat "$d/lk10_case.txt"; else echo "(no lk10_case.txt)"; fi
  done

  sec "drain: why the lifted plasma drains or does not"
  run $P drain "$RUN" ${REF[@]+"${REF[@]}"} --outdir "$FIGS"

  sec "mass: the plasma above 3 Mm"
  run $P mass "$RUN" ${REF[@]+"${REF[@]}"} --outdir "$FIGS"

  sec "null: the top of the emerged flux on x = 0"
  run $P null "$RUN" --outdir "$FIGS"

  sec "dome: top of log rho > -13 on x = 0"
  run $P dome "$RUN" --outdir "$FIGS"
  for f in "$FIGS"/dome_*.csv; do [ -f "$f" ] && { echo "--- $f"; cat "$f"; }; done

  for x in 0 5 10; do
    sec "column: the state on the vertical line x = $x Mm"
    run $P column "$RUN" --times "$TIMES" --x $x --outdir "$FIGS"
  done

  for L in ${LOGS[@]+"${LOGS[@]}"}; do
    sec "log $L"
    if [ ! -f "$L" ]; then echo "(not found)"; continue; fi
    echo "--- how the run was set up"
    grep -E 'DynSGS: C_R|Restart:|Subnormal|Dissipation:|Fastest wave|GLM: c_h|Δt *=>|dsgs_|lk10_|mhd_anomalous|Time loop running' "$L" | grep -v '^ *#' | head -80
    echo "--- positivity floors (first 15 and last 15 reports)"
    grep 'positivity floors' "$L" | head -15
    echo "   ..."
    grep 'positivity floors' "$L" | tail -15
    echo "--- stability cap, anomalous resistivity, failures"
    grep -E 'held at the stability cap|anomalous resistivity|non-finite|ERROR|Killed|BAD TERMINATION|:dsgs_stages' "$L" | head -40
    echo "--- CFL numbers at each snapshot (time, fast, viscous)"
    grep -E 'writing .*pvtu|Fast +CFL|Viscous +CFL|GLM +CFL' "$L" | sed -E 's/.*(iter_[0-9]+\.pvtu at t=[0-9.e+-]+).*/\1/' | paste -sd' ' | sed -E 's/ (iter_)/\n\1/g' | cut -c1-240 | tail -60
    echo "--- progress and the end of the run"
    grep -E 'progress:|DSGS over the run|Solving ODE|seconds \(' "$L" | tail -8
  done
} > "$OUT" 2>&1

echo "wrote $OUT ($(wc -l < "$OUT") lines): send this file. Figures and CSV tables are in $FIGS/."
