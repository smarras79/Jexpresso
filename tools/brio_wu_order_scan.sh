#!/usr/bin/env bash
#
# The Brio-Wu shock tube at every polynomial order and resolution asked for,
# and the figures of the whole sweep. At the end tools/brio_wu_plot.jl draws,
# from the curves the runs stored, in output/MHD/brioWu1d/output/:
#
#   density_nu_dof<M>-it5.svg   the density of every order at ~M DOFs with the
#                               zoom boxes, and their DynSGS coefficient ν: the
#                               layout of Dao & Nazarov (2024), Fig. 2
#   density_nu-it5.svg          the same at the finest resolution
#   density_dof<M>-it5.svg      the density alone (Dao & Nazarov 2022, Fig. 2)
#   convergence-it5.svg, convergence_smooth-it5.svg   error against 1/#DOFs
#
#   tools/brio_wu_order_scan.sh                       # N = 4 ... 7 at 600 and 1200 DOFs
#   BW_JOBS=4 tools/brio_wu_order_scan.sh             # the four orders at once
#   BW_DOFS=600 tools/brio_wu_order_scan.sh           # the paper's resolution alone
#   BW_DOFS="150 300 600 1200" tools/brio_wu_order_scan.sh   # for the convergence figures
#   BW_NOPS="4 7" tools/brio_wu_order_scan.sh
#   BW_CMIN=0.0 tools/brio_wu_order_scan.sh           # with the DynSGS floor off
#   BW_HSCALE=nop tools/brio_wu_order_scan.sh         # Dao & Nazarov's Δ_K/k instead of h_e/(N+1)
#   BW_NODAL=1 tools/brio_wu_order_scan.sh            # the nodal coefficient
#   JEXPRESSO_BW_FREEZE=0 tools/brio_wu_order_scan.sh # DynSGS per stage, as before b956ccd
#   JEXPRESSO_BW_FIGFMT=pdf tools/brio_wu_order_scan.sh   # PDF figures (default SVG)
#   JULIA=/path/to/julia tools/brio_wu_order_scan.sh  # default: the macOS app's julia, else the PATH's
#
# At ~M DOFs order N runs on M/N elements (150, 120, 100, 86 at 600 DOFs) with
# Δt = 5e-5·600/M, so every order costs about the same.
#
# PARALLEL. BW_JOBS=n runs up to n orders at once, one Julia session per order
# (each compiles the case once and then runs its resolutions one after the
# other); the log of order N goes to output/MHD/brioWu1d/scan-nop<N>.log.
# BW_JOBS=1, the default, is a single session for the whole sweep.
#
# THE STORE IS CLEARED FIRST, with the comparison figures an earlier sweep drew
# from it. Every run redraws the figures from every curve in
# problems/MHD/brioWu1d/curves, so a leftover sweep would appear on the
# comparison of this one. BW_KEEP=1 keeps both instead, to finish a partial
# sweep. A run that fails is reported and the sweep goes on.
set -u
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$0")/.."

CASE=problems/MHD/brioWu1d
OUT=output/MHD/brioWu1d
DOFS=${BW_DOFS:-"600 1200"}
NOPS=${BW_NOPS:-"4 5 6 7"}
JOBS=${BW_JOBS:-1}
export BW_DOFS="$DOFS"     # for the per-order children of a parallel sweep
# the length scale: the deck's default (the paper's h_e/(N+1)) unless BW_HSCALE says otherwise
[ -n "${BW_HSCALE:-}" ] && export JEXPRESSO_BW_HSCALE="$BW_HSCALE"
[ -n "${BW_CMIN:-}" ]  && export JEXPRESSO_BW_CMIN="$BW_CMIN"
[ -n "${BW_NODAL:-}" ] && export JEXPRESSO_BW_NODAL="$BW_NODAL"
# JULIA: as given, or else the first that starts of the macOS app's julia and
# the julia on the PATH. On a Mac the one on the PATH can be a juliaup launcher
# with no default channel, which fails every call ("failed to figure out which
# juliaup channel to use"). Whichever it is, the runs then call the binary
# itself, not a launcher in front of it: several launchers started at once can
# also fail on juliaup's lock.
if [ -z "${JULIA:-}" ]; then
    for j in /Applications/Julia-1.11.app/Contents/Resources/julia/bin/julia "$(command -v julia 2>/dev/null)"; do
        [ -n "$j" ] && "$j" --version >/dev/null 2>&1 && { JULIA=$j; break; }
    done
fi
JULIA_IN=${JULIA:-julia}
JULIA=$("$JULIA_IN" --startup-file=no -e 'print(joinpath(Sys.BINDIR, Base.julia_exename()))' 2>/dev/null) || {
    echo "=== $JULIA_IN does not start: set JULIA=/path/to/julia (or run 'juliaup default 1.11')"; exit 1; }
export JULIA
export GKSwstype=${GKSwstype:-100}

# Every run of the orders given, in ONE Julia session: the case is compiled
# once and each run only re-reads its JEXPRESSO_BW_* settings.
run_orders() {
    runs=""
    for n in "$@"; do
        for d in $DOFS; do runs="$runs $n:$d"; done
    done
    "$JULIA" --project=. -e '
using Jexpresso
runs   = copy(ARGS)      # run_case empties ARGS and puts the case name in it
failed = String[]
for r in runs
    n, d = split(r, ":")
    ENV["JEXPRESSO_BW_NOP"] = n; ENV["JEXPRESSO_BW_DOFS"] = d
    println("=== Brio-Wu: N = ", n, ", ~", d, " DOFs")
    t0 = time()
    try
        Jexpresso.run_case("MHD", "brioWu1d")
        println("=== N = ", n, ", ", d, " DOFs: ", round(time() - t0, digits = 1), " s")
    catch err
        push!(failed, string("N = ", n, ", ", d, " DOFs"))
        println("=== N = ", n, ", ", d, " DOFs: FAILED, ", first(sprint(showerror, err), 300))
    end
end
isempty(failed) || println("=== FAILED: ", join(failed, "; "))
' $runs
}

# one order of a parallel sweep (BW_JOBS > 1); xargs calls the script back
if [ "${1:-}" = "--order" ]; then
    log="$OUT/scan-nop$2.log"
    echo "=== N = $2: started, log $log"
    run_orders "$2" > "$log" 2>&1
    echo "=== N = $2: $(grep -c '^=== N = .* s$' "$log") of $(echo $DOFS | wc -w) runs done$(grep '^=== FAILED' "$log" | sed 's/^===/,/')"
    exit 0
fi

if [ "${BW_KEEP:-0}" != "1" ]; then
    # the store, and the comparison figures drawn from it by an earlier sweep
    rm -rf "$CASE/curves"
    rm -f "$OUT"/output/density_dof*-it* "$OUT"/output/density_nu*-it* "$OUT"/output/convergence*-it*
fi
mkdir -p "$OUT"
echo "=== Brio-Wu order scan: N [$NOPS] x DOFs [$DOFS], length scale \
$([ "${JEXPRESSO_BW_HSCALE:-ngl}" = nop ] && echo 'Delta_K/k' || echo 'Delta_K/(k+1)')\
${BW_NODAL:+, nodal coefficient}${BW_CMIN:+, Cmin $BW_CMIN}, $JOBS session(s)"

if [ "$JOBS" -gt 1 ]; then
    printf '%s\n' $NOPS | xargs -P "$JOBS" -I{} "$SELF" --order {}
else
    run_orders $NOPS
fi

echo "=== the whole sweep:"
"$JULIA" --project=. tools/brio_wu_plot.jl
