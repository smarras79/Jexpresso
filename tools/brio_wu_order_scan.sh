#!/usr/bin/env bash
#
# The Brio-Wu shock tube at every polynomial order and resolution asked for,
# with the DynSGS settings of the deck, and the figures of the sweep. The runs
# store their curves in
# problems/MHD/brioWu1d/curves/<store>/, and when the last run is done
# tools/brio_wu_plot.jl draws from them, in output/MHD/brioWu1d/output/<store>/:
#
#   rho_By_dof<M>-it5.svg       every order at ~M DOFs, one figure per resolution:
#                               ρ (left) and B_y (right) against the reference,
#                               with zoom boxes; rho_By-it5.svg the finest
#   nu_all-it5.svg, nu_dof<M>-it5.svg, nu-it5.svg   the DynSGS coefficient ν alone
#   density_dof<M>-it5.svg, density-it5.svg    the density alone
#   convergence-it5.svg, convergence_smooth-it5.svg   error against 1/#DOFs
#
#   tools/brio_wu_order_scan.sh                       # N = 4 ... 7 at 600 and 1200 DOFs
#   BW_JOBS=4 tools/brio_wu_order_scan.sh             # the four orders at once
#   BW_DOFS=600 tools/brio_wu_order_scan.sh           # the paper's resolution alone
#   BW_DOFS="150 300 600 1200" tools/brio_wu_order_scan.sh   # for the convergence figures
#   BW_NOPS="4 7" tools/brio_wu_order_scan.sh
#   JEXPRESSO_BW_FIGFMT=pdf tools/brio_wu_order_scan.sh   # PDF figures (default SVG)
#   JULIA=/path/to/julia tools/brio_wu_order_scan.sh  # default: the macOS app's julia, else the PATH's
#
# THE DYNSGS SETTINGS ARE THE DECK'S. C_min, the cutoff, C_1, the length-scale
# convention and the element/nodal form are literal entries of
# problems/MHD/brioWu1d/user_inputs.jl (:dsgs_Cmin, :dsgs_cutoff, :dsgs_C1,
# :bw_hscale, :ldsgs_nodal); to compare values, edit the deck and run the scan
# again, once per value. The script reads them from the deck only to name the
# store, <store> = Cmin<C_min>, plus _nop (:bw_hscale => "nop"), _nodal,
# _C1_<C_1> (C_1 ≠ 1) and _cut<c> (c > 0): Cmin0.06, Cmin0.03_cut1e-3, ...
# so that different methods never share a store or a figure.
#
# At ~M DOFs order N runs on M/N elements (150, 120, 100, 86 at 600 DOFs) with
# Δt = 5e-5·600/M, so every order costs about the same.
#
# PARALLEL. BW_JOBS=n runs up to n orders at once, one Julia session per order
# (each compiles the case once and then runs its resolutions one after the
# other); the log of order N goes to output/MHD/brioWu1d/scan-<store>-nop<N>.log.
# BW_JOBS=1, the default, is a single session.
#
# EACH SWEEP STARTS FROM AN EMPTY STORE: the store and the figures of the same
# settings from an earlier sweep are deleted first, so that a leftover curve can
# never appear on the comparison. BW_KEEP=1 keeps them, to finish a partial
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
# The DynSGS knobs of this script before Nov 2026 are refused: those settings
# are the deck's now (see the top of this file).
for v in BW_CMIN BW_CUTOFF BW_C1 BW_CR BW_HSCALE BW_NODAL; do
    [ -n "${!v:-}" ] && { echo "=== $v is no longer read: set the DynSGS settings in $CASE/user_inputs.jl"; exit 1; }
done
# The literal value of a key of the deck's Dict (an uncommented `:key => value,` line).
deck_val() {
    sed -n -E "s/^[[:space:]]*:$1[[:space:]]*=>[[:space:]]*\"?([^\",[:space:]]+)\"?[[:space:]]*,.*/\1/p" \
        "$CASE/user_inputs.jl" | head -n 1
}
CMIN=$(deck_val dsgs_Cmin); CUT=$(deck_val dsgs_cutoff); C1=$(deck_val dsgs_C1)
HSCALE=$(deck_val bw_hscale); NODAL=$(deck_val ldsgs_nodal)
for v in CMIN CUT C1 HSCALE NODAL; do
    [ -n "${!v}" ] || { echo "=== no literal value of $v in $CASE/user_inputs.jl"; exit 1; }
done
# the runs only store their curves; each sweep draws its figures once, at the end
export JEXPRESSO_BW_DRAWSTORE=0
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
    println("=== Brio-Wu: N = ", n, ", ~", d, " DOFs, store ", ENV["JEXPRESSO_BW_STORE"])
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
    log="$OUT/scan-$JEXPRESSO_BW_STORE-nop$2.log"
    echo "=== N = $2: started, log $log"
    run_orders "$2" > "$log" 2>&1
    echo "=== N = $2: $(grep -c '^=== N = .* s$' "$log") of $(echo $DOFS | wc -w) runs done$(grep '^=== FAILED' "$log" | sed 's/^===/,/')"
    exit 0
fi

mkdir -p "$OUT"
# comparison figures the runs drew in the case's output directory before the
# sweeps drew theirs in output/Cmin<C_min>: stale, whatever they show
rm -f "$OUT"/output/density_dof*-it* "$OUT"/output/density_nu*-it* "$OUT"/output/convergence*-it* \
      "$OUT"/output/rho_By*-it* "$OUT"/output/nu-it* "$OUT"/output/nu_*-it*
case "$HSCALE" in nop) suffix="_nop" ;; *) suffix="" ;; esac
case "$NODAL" in true) suffix="${suffix}_nodal" ;; esac
case "$CUT" in 0|0.0|0e0|0e-0) cut="" ;; *) cut="_cut$CUT" ;; esac
case "$C1" in 1|1.0) crs="" ;; *) crs="_C1_$C1" ;; esac
export JEXPRESSO_BW_STORE="Cmin$CMIN$suffix$crs$cut"
figs="$OUT/output/$JEXPRESSO_BW_STORE"
[ "${BW_KEEP:-0}" = "1" ] || rm -rf "$CASE/curves/$JEXPRESSO_BW_STORE" "$figs"
echo "=== Brio-Wu order scan, the deck's DynSGS (C_min = $CMIN, cutoff $CUT, C_1 $C1, \
length scale $([ "$HSCALE" = nop ] && echo 'Delta_K/k' || echo 'Delta_K/(k+1)')\
$([ "$NODAL" = true ] && echo ', nodal coefficient')): N [$NOPS] x DOFs [$DOFS], $JOBS session(s)"
if [ "$JOBS" -gt 1 ]; then
    printf '%s\n' $NOPS | xargs -P "$JOBS" -I{} "$SELF" --order {}
else
    run_orders $NOPS
fi
echo "=== the whole sweep, in $figs:"
"$JULIA" --project=. tools/brio_wu_plot.jl --out="$figs"
