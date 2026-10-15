#!/usr/bin/env bash
#
# Convergence history of the smooth Alfvén wave (Babbar & Ranocha 2026, §5.5.1
# and Fig. 13) for the polynomial orders N = 3 ... 7: every order on its own
# range of meshes, with DynSGS and with plain Galerkin. At the end
# tools/alfven_wave_plot.jl draws the figure of the whole sweep,
# output/MHD/alfvenWaveBabbarRanocha2026/output/convergence-it0.png, and prints
# the table of errors and observed orders.
#
#   tools/alfven_wave_scan.sh                       # N = 3 ... 7, DSGS and Galerkin
#   AW_JOBS=5 tools/alfven_wave_scan.sh             # the five orders at once
#   AW_NOPS="3 4" AW_NELX="16 32 64 128" tools/alfven_wave_scan.sh   # the paper's meshes
#   AW_NELX_7="4 8 16" tools/alfven_wave_scan.sh    # other meshes for one order
#   AW_VISC=dsgs tools/alfven_wave_scan.sh          # only the DynSGS runs
#   JULIA=/path/to/julia tools/alfven_wave_scan.sh  # default: the macOS app's julia, else the PATH's
#
# MESHES. Each order gets the meshes that take its Galerkin error from the
# coarse-mesh range down to about 1e-12, short of the round-off level, so that
# the five curves cover about the same span of #DOFs:
#   N = 3, 4    8², 16², 32², 64² elements   (64²: 4e-7 and 2e-9)
#   N = 5, 6    4², 8², 16², 32²             (32²: 5e-10 and 4e-12)
#   N = 7       4², 8², 16²                  (16²: 6e-12; 32² would be round-off)
# AW_NELX sets the meshes of every order, AW_NELX_<N> those of order N alone.
#
# TIME INTEGRATION follows the deck (user_inputs.jl): CarpenterKennedy2N54 up
# to N = 5 and Vern7 above, at a Courant number of 0.1 for every run.
#
# COST: about 3.5 h one order after the other (AW_JOBS=1, the default: one
# Julia session, the case compiled once); 83 min with AW_JOBS=4 on four cores,
# the time of the N = 4 session. The longest run is N = 4 on 64² with DynSGS,
# about 40 min. AW_JOBS=n runs up to n orders at once, one Julia session per
# order (each compiles the case, 2-3 min, and holds about 2 GB). The log of
# order N then goes to
# output/MHD/alfvenWaveBabbarRanocha2026/scan-nop<N>.log. The order is the unit
# because two sessions on the same mesh AND order would write the same
# .jexpresso_cache file at the same time, and JLD2 then dies with a bus error.
#
# THE ERROR STORE IS CLEARED FIRST, since every run redraws the figure from all
# the errors stored for the same final time; AW_KEEP=1 accumulates instead, to
# finish a partial sweep. A run that fails is reported and the sweep goes on.
set -u
SELF="$(cd "$(dirname "$0")" && pwd)/$(basename "$0")"
cd "$(dirname "$0")/.."

CASE=problems/MHD/alfvenWaveBabbarRanocha2026
OUT=output/MHD/alfvenWaveBabbarRanocha2026
NOPS=${AW_NOPS:-"3 4 5 6 7"}
VISCS=${AW_VISC:-"dsgs none"}
JOBS=${AW_JOBS:-1}
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

# the meshes of order $1
nelx_of() {
    eval "v=\${AW_NELX_$1:-}"
    if   [ -n "$v" ];             then echo "$v"
    elif [ -n "${AW_NELX:-}" ];   then echo "$AW_NELX"
    elif [ "$1" -le 4 ];          then echo "8 16 32 64"
    elif [ "$1" -le 6 ];          then echo "4 8 16 32"
    else                               echo "4 8 16"
    fi
}

# Every run of the orders given, in ONE Julia session: the case is compiled
# once and each run only re-reads its JEXPRESSO_AW_* settings.
run_orders() {
    runs=""
    for n in "$@"; do
        for v in $VISCS; do
            for m in $(nelx_of "$n"); do runs="$runs $v:$n:$m"; done
        done
    done
    "$JULIA" --project=. -e '
using Jexpresso
runs   = copy(ARGS)      # run_case empties ARGS and puts the case name in it
failed = String[]
for r in runs
    v, n, m = split(r, ":")
    ENV["JEXPRESSO_AW_VISC"] = v; ENV["JEXPRESSO_AW_NOP"] = n; ENV["JEXPRESSO_AW_NELX"] = m
    println("=== smooth Alfvén wave: N = ", n, ", ", m, "x", m, " elements, visc ", v)
    t0 = time()
    try
        Jexpresso.run_case("MHD", "alfvenWaveBabbarRanocha2026")
        println("=== N = ", n, ", ", m, "x", m, ", ", v, ": ", round(time() - t0, digits = 1), " s")
    catch err
        push!(failed, string("N = ", n, ", ", m, "x", m, ", ", v))
        println("=== N = ", n, ", ", m, "x", m, ", ", v, ": FAILED, ", first(sprint(showerror, err), 300))
    end
end
isempty(failed) || println("=== FAILED: ", join(failed, "; "))
' $runs
}

# one order of a parallel sweep (AW_JOBS > 1); xargs calls the script back
if [ "${1:-}" = "--order" ]; then
    log="$OUT/scan-nop$2.log"
    echo "=== N = $2: started, log $log"
    run_orders "$2" > "$log" 2>&1
    nrun=$(( $(echo $VISCS | wc -w) * $(nelx_of "$2" | wc -w) ))
    echo "=== N = $2: $(grep -c '^=== N = .* s$' "$log") of $nrun runs done$(grep '^=== FAILED' "$log" | sed 's/^===/,/')"
    exit 0
fi

[ "${AW_KEEP:-0}" = "1" ] || rm -rf "$CASE/errors"
meshes=$(for n in $NOPS; do nelx_of "$n"; done | tr ' ' '\n' | sort -n | uniq | tr '\n' ' ')
AW_NELX="$meshes" tools/alfven_wave_mesh.sh || exit 1
mkdir -p "$OUT"

if [ "$JOBS" -gt 1 ]; then
    printf '%s\n' $NOPS | xargs -P "$JOBS" -I{} "$SELF" --order {}
else
    run_orders $NOPS
fi

echo "=== the whole sweep:"
"$JULIA" --project=. tools/alfven_wave_plot.jl
