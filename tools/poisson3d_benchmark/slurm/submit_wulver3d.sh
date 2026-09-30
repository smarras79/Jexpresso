#!/bin/bash
# =============================================================================
#  Submit the 3D periodic Poisson solver comparison to SLURM on NJIT Wulver:
#  one job per configuration (solver, ne, N), each on its own cores, with
#  cores / memory / time sized from the problem; then one merge job.
#
#      bash tools/poisson3d_benchmark/slurm/submit_wulver3d.sh [env-file]
#
#  env-file defaults to wulver3d.env next to this script. Environment options:
#      DRY_RUN=1          print the job table and the sbatch commands only
#      SKIP_PRECOMPILE=1  do not submit the precompile/verification job
#
#  Job chain:  precompile (+ verify_2d, verify_3d)
#                -> one job per configuration  (afterok)
#                -> merge + figures            (afterany: when all have ended)
#  A configuration whose $OUTDIR/parts/<solver>_ne<ne>_N<N>/results.csv says
#  "ok" is not submitted again, so re-running this script resubmits only what
#  is missing (failed, timed out, preempted, or not yet run).
# =============================================================================
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVFILE="$(realpath "${1:-$HERE/wulver3d.env}")"
source "$ENVFILE"
DRY_RUN="${DRY_RUN:-0}"; SKIP_PRECOMPILE="${SKIP_PRECOMPILE:-0}"
mkdir -p "$OUTDIR/logs" "$OUTDIR/parts"
cp "$ENVFILE" "$OUTDIR/logs/settings.env.used"

common=(--account="$ACCOUNT" --qos="$QOS" --output="$OUTDIR/logs/%x.%j.out")
[ "$QOS" = "low" ] && common+=(--requeue)
[ "$EXCLUSIVE" = "1" ] && common+=(--exclusive)
ENVX="ALL,P3D_ENV=$ENVFILE,P3D_DIR=$HERE"   # batch scripts run from a spool copy: pass our directory

submit() {
    if [ "$DRY_RUN" = "1" ]; then echo "sbatch $*" >&2; echo "DRY$RANDOM"; else sbatch --parsable "$@"; fi
}
hms() { printf "%02d:%02d:00" $(($1/60)) $(($1%60)); }

# ---- resource model ------------------------------------------------------------
# Calibrated in 3D (N = 2..6, n up to 1.1e5; see README.md):
#   nnz(K) = (3N+4) n;  METIS Cholesky factor of K: nnz ~ 25 n^(4/3);
#   skeleton: n_s = n (1 - ((N-1)/N)^3), nnz(B) ~ 0.7 (N+1)^3 n_s,
#   its Cholesky factor ~ 5 (N+1)^1.5 n_s^(4/3)   (all ×(size/1e5)^0.05 drift).
# Times (1 core, conservative): direct factorisation ~ 10 s (n/1e5)^2 (3D
# nested dissection is O(n^2)); AMG-CG ~ 3e-5 s n; Jacobi-CG ~ 1e-5 s n (n/1.4e4)^(1/3)
# (iterations grow like h^-1); spectral solvers are seconds. Direct and
# condensation times are divided by threads^0.6.
resources() {  # solver ne N -> "n thr cpus memGB minutes partition status"
    awk -v s="$1" -v ne="$2" -v N="$3" -v d="$D" \
        -v b1="$THREADS_BREAK1" -v b2="$THREADS_BREAK2" -v ts="$THREADS_SMALL" -v tm="$THREADS_MEDIUM" -v tl="$THREADS_LARGE" \
        -v gp="$PARTITION" -v bp="$BIGMEM_PARTITION" -v gmax="$GENERAL_MAX_GB" -v bmax="$BIGMEM_MAX_GB" \
        -v hmax="$MAX_HOURS" -v safety="$MEM_SAFETY" '
    function ceil(x) { return (x == int(x)) ? x : int(x) + 1 }
    BEGIN {
        ng = ne * N; n = ng^d; nnzK = (3 * N + 4) * n
        drift = (n / 1e5)^0.05; if (drift < 1) drift = 1
        ns = n * (1 - ((N - 1) / N)^3); nnzB = 0.7 * (N + 1)^3 * ns
        thr = (n < b1) ? ts : ((n < b2) ? tm : tl); sp = thr^0.6
        base = 3.0
        if (s == "sem")        { mem = base + (40*nnzK + 12*25*n^(4/3)*drift + 200*n) / 1e9;          t = 10 * (n/1e5)^2 / sp + 1e-5*n }
        else if (s == "sem_amg")    { mem = base + (100*nnzK + 400*n) / 1e9;                           t = 3e-5 * n }
        else if (s == "sem_jacobi") { mem = base + (40*nnzK + 100*n) / 1e9;                            t = 1e-5 * n * (n/1.4e4)^(1/3) }
        else if (s == "sc_direct")  { mem = base + (40*nnzK + 40*nnzB + 12*5*(N+1)^1.5*ns^(4/3)*drift + 0.7e9) / 1e9
                                      t = (2e-4*(N/4)^6*ne^d + 30*(ns/1e5)^2) / sp + 1e-5*n }
        else if (s == "sc_amg")     { mem = base + (40*nnzK + 100*nnzB + 0.7e9) / 1e9;                 t = 2e-4*(N/4)^6*ne^d/sp + 5e-5*n }
        else if (s == "ps")         { mem = base + 100*n / 1e9;                                        t = 1e-9 * 2 * d * ng^(d+1) + 1e-6*n }
        else                        { mem = base + 80*n / 1e9;                                         t = 1e-6 * n }
        memgb = ceil(mem * safety) + 1
        minutes = ceil(20 + 3 * t / 60)
        cpus = thr; if (ceil(memgb / 4) > cpus) cpus = ceil(memgb / 4); if (cpus > 128) cpus = 128
        part = gp; st = "ok"
        if (memgb > gmax) part = bp
        if (memgb > bmax) st = sprintf("skip:memory(%dGB)", memgb)
        if (minutes > hmax * 60) st = sprintf("skip:time(%dh)", minutes / 60)
        if (minutes > hmax * 60) minutes = hmax * 60
        printf "%d %d %d %d %d %s %s\n", n, thr, cpus, memgb, minutes, part, st
    }'
}

# ---- precompile + verification ---------------------------------------------------
dep=()
if [ "$SKIP_PRECOMPILE" != "1" ]; then
    pre=$(submit "${common[@]}" --partition="$PARTITION" --export="$ENVX" "$HERE/precompile.sbatch")
    echo "precompile + verification job: $pre"
    dep=(--dependency=afterok:$pre)
fi

# ---- one job per configuration ------------------------------------------------------
ids=()
printf "%-10s %3s %2s %12s %4s %5s %7s %9s %-8s %s\n" solver ne N unknowns thr cpus mem time partition job
for sweep in $SWEEPS; do
    N="${sweep%%:*}"
    for ne in $(echo "${sweep#*:}" | tr ',' ' '); do
        for S in $SOLVERS; do
            read -r n thr cpus mem mins part st <<< "$(resources "$S" "$ne" "$N")"
            done_csv="$OUTDIR/parts/${S}_ne${ne}_N${N}/results.csv"
            if [ -s "$done_csv" ] && grep -q ",ok$" "$done_csv"; then
                printf "%-10s %3s %2s %12s  done\n" "$S" "$ne" "$N" "$n"; continue
            fi
            if [ "$st" != "ok" ]; then
                printf "%-10s %3s %2s %12s  %s\n" "$S" "$ne" "$N" "$n" "$st"; continue
            fi
            id=$(submit "${common[@]}" "${dep[@]}" --partition="$part" \
                        --job-name="p3d-${S}-ne${ne}-N${N}" --cpus-per-task="$cpus" --mem="${mem}G" \
                        --time="$(hms "$mins")" \
                        --export="$ENVX,P3D_SOLVER=$S,P3D_NE=$ne,P3D_NOP=$N,P3D_THREADS=$thr" \
                        "$HERE/job.sbatch")
            ids+=("$id")
            printf "%-10s %3s %2s %12s %4s %5s %6sG %9s %-8s %s\n" "$S" "$ne" "$N" "$n" "$thr" "$cpus" "$mem" "$(hms "$mins")" "$part" "$id"
        done
    done
done

# ---- merge + figures when every job has ended -------------------------------------------
if [ ${#ids[@]} -gt 0 ]; then
    mdep=$(IFS=:; echo "${ids[*]}")
    mid=$(submit "${common[@]}" --partition="$PARTITION" --dependency=afterany:$mdep --export="$ENVX" "$HERE/merge.sbatch")
    echo "merge job: $mid  (after ${#ids[@]} benchmark jobs)"
else
    echo "nothing to run; merge with:  sbatch ${common[*]} --partition=$PARTITION --export=$ENVX $HERE/merge.sbatch"
fi
echo "results: $OUTDIR  (logs: $OUTDIR/logs)"
