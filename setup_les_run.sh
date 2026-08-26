#!/bin/bash
# ===========================================================================
# setup_les_run.sh -- one-shot setup for a Jexpresso LES run.
#
# Generates the mesh, orders it so the decomposition is balanced, patches the
# case's user_inputs.jl, writes a matching SLURM script, and optionally
# submits. Everything that has cost a queue cycle in practice is derived here
# rather than left to be edited by hand and kept consistent by memory:
#
#   * mesh gmsh generates       = TARGET_GRID / 2^REFINE_LVL
#   * rank count                = largest valid one that keeps ranks fed
#   * --nodes / --ntasks-per-node
#   * :lxy_partition / :linitial_refine / :init_refine_lvl / :gmsh_filename
#
# EDIT ONLY THE BLOCK BELOW.
# ===========================================================================
MODULES=(Julia/1.11.9 GCC MPICH)

CASE="LESICP2-64x64x60-imex"        # problems/CompEuler/<CASE>

for m in "${MODULES[@]}"; do module load "$m" || exit 1; done
# ---------------------------------------------------------------------------
# MESH: supply your own (the default), or have gmsh build one.
#
# Wulver has no gmsh, so the normal workflow is to run gmsh on your laptop and
# copy the .msh over. Set COARSE_MESH to that file and this script takes it
# from there -- it reads the element counts and extents OUT OF THE FILE, so
# there is nothing to keep in sync by hand.
#
# RUN_GMSH is OFF by default and must be turned on deliberately. Meshing is the
# one step here that overwrites a file you may have spent effort producing, and
# an empty COARSE_MESH is a likely typo -- so an unset mesh is now an error
# that says so, rather than a silent decision to build a different one.
# ---------------------------------------------------------------------------
RUN_GMSH="no"               # "yes" to write a .geo and call gmsh (needs gmsh on PATH)
COARSE_MESH=""              # e.g. "meshes/LESICP_16x16x15_10240mX10240mX5000m.msh"

# TARGET_GRID and DOMAIN are used ONLY when COARSE_MESH is empty, i.e. when
# this script has to generate the mesh. When you supply a mesh, both are read
# from it and these are ignored.
#
# TARGET_GRID is the FINAL resolution -- what the simulation solves on, NOT
# what gmsh produces:
#
#     mesh gmsh generates  =  TARGET_GRID / 2^REFINE_LVL
#
#     TARGET_GRID   REFINE_LVL   gmsh generates   simulation solves on
#      32x32x30          0          32x32x30          32x32x30
#      32x32x30          1          16x16x15          32x32x30
#     128x128x120        3          16x16x15         128x128x120
#
TARGET_GRID="64x64x60"      # final resolution; also used by REFINE_LVL="auto"
DOMAIN="10240x10240x5000"   # metres             (ignored if COARSE_MESH set)

# REFINE_LVL applies either way. With a supplied mesh, the resolution you end
# up solving on is  <mesh dimensions> x 2^REFINE_LVL.
# REFINE_LVL may be an integer, or "auto".
#   integer : refine the mesh this many times; the resolution you solve on is
#             <mesh dimensions> x 2^REFINE_LVL.
#   "auto"  : derive it from TARGET_GRID -- the script works out how many
#             refinements take your mesh to that resolution, and refuses if no
#             whole number of levels does.
REFINE_LVL=0                # integer, or "auto" (needs TARGET_GRID)
MAX_CORES=256               # upper bound on ranks
MEM_PER_CPU="7000M"         # per rank; setup peaks well above the time loop
SUBMIT="no"                 # "yes" to sbatch at the end

# ===========================================================================
# Nothing below here needs editing.
# ===========================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$ROOT"

EQS="CompEuler"
NOP=$(awk -F'=>' '/^[[:space:]]*:nop[[:space:]]*=>/ {if (match($2,/[0-9]+/)) {print substr($2,RSTART,RLENGTH); exit}}' \
      "problems/$EQS/$CASE/user_inputs.jl" 2>/dev/null || echo 4)
: "${NOP:=4}"

CASEDIR="problems/$EQS/$CASE"
[ -d "$CASEDIR" ] || { echo "ERROR: no such case directory: $CASEDIR" >&2; exit 1; }

case "$RUN_GMSH" in
    yes|YES|true|1)  RUN_GMSH=yes ;;
    no|NO|false|0|"") RUN_GMSH=no ;;
    *) echo "ERROR: RUN_GMSH must be yes or no; got '$RUN_GMSH'." >&2; exit 1 ;;
esac

# The two mesh sources are mutually exclusive, and BOTH ways of getting it
# wrong are silent if not caught here: with neither set the script used to
# generate a mesh nobody asked for, and with both set it would generate one and
# then ignore it in favour of COARSE_MESH.
if [ "$RUN_GMSH" = yes ] && [ -n "$COARSE_MESH" ]; then
    echo "ERROR: RUN_GMSH=yes and COARSE_MESH are both set." >&2
    echo "       Pick one: unset COARSE_MESH to build a mesh from TARGET_GRID" >&2
    echo "       and DOMAIN, or set RUN_GMSH=no to use the mesh you supplied." >&2
    exit 1
fi
if [ "$RUN_GMSH" = no ] && [ -z "$COARSE_MESH" ]; then
    echo "ERROR: no mesh. COARSE_MESH is empty and RUN_GMSH is no, so there is" >&2
    echo "       nothing to run on and nothing is being generated." >&2
    echo "       Either set COARSE_MESH to a .msh you already have, or set" >&2
    echo "       RUN_GMSH=\"yes\" to build one from TARGET_GRID/DOMAIN (needs" >&2
    echo "       gmsh on PATH -- Wulver has none)." >&2
    exit 1
fi

POW=1
[ "$REFINE_LVL" = "auto" ] || POW=$(( 2 ** REFINE_LVL ))
if [ "$RUN_GMSH" = no ]; then
    # Supplied mesh is the source of truth. Read its shape rather than trusting
    # a restatement of it -- a mesh whose real dimensions differ from what the
    # config claims is exactly the failure that has been hardest to see.
    [ -f "$COARSE_MESH" ] || { echo "ERROR: COARSE_MESH not found: $COARSE_MESH" >&2; exit 1; }
    # Check the inspection SUCCEEDED before reading what it was meant to set.
    # `eval` of a failed command sets nothing, and the next line then trips
    # `set -u` with "MESH_STRUCTURED: unbound variable" -- which says nothing
    # about julia being missing, the mesh being unreadable, or the parse having
    # failed, all of which land here identically.
    MESH_FACTS="$(julia --project=. --startup-file=no tools/inspect_msh.jl "$COARSE_MESH" --machine)" || {
        echo "ERROR: could not inspect $COARSE_MESH." >&2
        echo "       tools/inspect_msh.jl failed -- is julia on PATH, and is the" >&2
        echo "       file a readable gmsh mesh? Its own output is above." >&2
        exit 1; }
    eval "$MESH_FACTS"
    if [ "${MESH_STRUCTURED:-}" != "1" ]; then
        echo "ERROR: $COARSE_MESH is not a regular box lattice." >&2
        echo "       The rank planning here assumes one. Run tools/inspect_msh.jl for detail." >&2
        exit 1
    fi
    if [ "$REFINE_LVL" = "auto" ]; then
        # Derive the level from the target instead of making the user divide by
        # powers of two and get it right. Every axis must agree on the SAME
        # power, otherwise no uniform refinement produces that grid -- p4est
        # refines all three directions together.
        IFS=x read -r TNX TNY TNZ <<< "$TARGET_GRID"
        REFINE_LVL=""
        for l in 0 1 2 3 4 5 6; do
            q=$(( 2 ** l ))
            if [ $(( MESH_NX * q )) -eq "$TNX" ] && \
               [ $(( MESH_NY * q )) -eq "$TNY" ] && \
               [ $(( MESH_NZ * q )) -eq "$TNZ" ]; then
                REFINE_LVL=$l; break
            fi
        done
        if [ -z "$REFINE_LVL" ]; then
            echo "ERROR: no whole number of refinements takes" >&2
            echo "         ${MESH_NX} x ${MESH_NY} x ${MESH_NZ}   (your mesh)" >&2
            echo "       to  ${TNX} x ${TNY} x ${TNZ}   (TARGET_GRID)" >&2
            echo "       Uniform refinement doubles all three axes at once, so the" >&2
            echo "       ratio must be the same power of 2 in x, y and z. Reachable" >&2
            echo "       targets from this mesh:" >&2
            for l in 0 1 2 3 4; do
                q=$(( 2 ** l ))
                echo "         lvl $l -> $(( MESH_NX * q )) x $(( MESH_NY * q )) x $(( MESH_NZ * q ))" >&2
            done
            exit 1
        fi
        POW=$(( 2 ** REFINE_LVL ))
        echo "    REFINE_LVL=auto  ->  derived $REFINE_LVL"
    fi
    NX=$(( MESH_NX * POW )); NY=$(( MESH_NY * POW )); NZ=$(( MESH_NZ * POW ))
    LX=$MESH_LX; LY=$MESH_LY; LZ=$MESH_LZ
else
    IFS=x read -r NX NY NZ <<< "$TARGET_GRID"
    IFS=x read -r LX LY LZ <<< "$DOMAIN"
fi

echo "=== 1/5  planning ==========================================="
if [ "$RUN_GMSH" = no ]; then
  echo "    using YOUR mesh          : $COARSE_MESH"
  echo "    which is                 : ${MESH_NX} x ${MESH_NY} x ${MESH_NZ} elements over ${LX} x ${LY} x ${LZ} m"
  echo "    REFINE_LVL=${REFINE_LVL}  =>  simulation will SOLVE on : ${NX} x ${NY} x ${NZ}"
else
  echo "    simulation will SOLVE on : ${NX} x ${NY} x ${NZ}   (TARGET_GRID)"
  echo "    gmsh will GENERATE       : $((NX/POW)) x $((NY/POW)) x $((NZ/POW))   (TARGET_GRID / 2^${REFINE_LVL})"
fi
echo
# The planner replicates the real partitioners, so a rank count it lists is one
# the code can actually build. See tools/plan_run.jl for why the two paths have
# different rules.
eval "$(julia --project=. --startup-file=no tools/plan_run.jl \
        "$NX" "$NY" "$NZ" "$NOP" "$REFINE_LVL" "$MAX_CORES" "$LX" "$LY" --machine)"
julia --project=. --startup-file=no tools/plan_run.jl \
        "$NX" "$NY" "$NZ" "$NOP" "$REFINE_LVL" "$MAX_CORES" "$LX" "$LY"

echo "=== 2/5  mesh =============================================="
if [ "$RUN_GMSH" = no ]; then
    MESHBASE="$(basename "${COARSE_MESH%.msh}")"
    RAW="$CASEDIR/${MESHBASE}.msh"
    if [ "$(cd "$(dirname "$COARSE_MESH")" && pwd)/$(basename "$COARSE_MESH")" != "$ROOT/$RAW" ]; then
        cp "$COARSE_MESH" "$RAW"
        echo "    copied your mesh into the case: $RAW"
    else
        echo "    using in place: $RAW"
    fi
else
    MESHBASE="LESICP_$((NX/POW))x$((NY/POW))x$((NZ/POW))_${LX}mX${LY}mX${LZ}m"
    GEO="$CASEDIR/${MESHBASE}.geo"
    RAW="$CASEDIR/${MESHBASE}.msh"
    cat > "$GEO" <<GEOEOF
// Generated by setup_les_run.sh -- do not hand-edit; re-run the script.
// $((NX/POW)) x $((NY/POW)) x $((NZ/POW)) elements over ${LX} x ${LY} x ${LZ} m.
// REFINE_LVL=${REFINE_LVL} takes this to the ${NX} x ${NY} x ${NZ} grid the run solves on.
nelemx = $((NX/POW));
nelemy = $((NY/POW));
nelemz = $((NZ/POW));

xmin = 0.0;  xmax = ${LX};
ymin = 0.0;  ymax = ${LY};
zmin = 0.0;  zmax = ${LZ};

gridsize = (xmax-xmin) / nelemx;

Point(1) = {xmin, ymin, zmin, gridsize};
Point(2) = {xmax, ymin, zmin, gridsize};
Point(3) = {xmax, ymin, zmax, gridsize};
Point(4) = {xmin, ymin, zmax, gridsize};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

npx = nelemx + 1;
npz = nelemz + 1;

Transfinite Line {1, 3} = npx;
Transfinite Line {4, -2} = npz Using Progression 1.0;

Line Loop(11) = {4, 1, 2, 3};
Plane Surface(12) = {11};
Transfinite Surface {12};
Recombine Surface {12};

surfaceVector = Extrude {0,(ymax-ymin),0} {
  Surface{12};
  Layers{nelemy};
  Recombine;
};

Physical Surface("periodicy") = {12,34};
Physical Volume("internal")   = {1};
Physical Surface("top_wall")  = {33};
Physical Surface("MOST")      = {25};
Physical Surface("periodicx") = {21,29};
GEOEOF
    echo "    wrote $GEO"
    if ! command -v gmsh >/dev/null 2>&1; then
        echo "ERROR: gmsh is not on PATH (Wulver has none)." >&2
        echo "       Run this on your laptop:  gmsh -3 $GEO -o $(basename "$RAW")" >&2
        echo "       then copy the .msh over and set COARSE_MESH to it." >&2
        exit 1
    fi
    gmsh -3 "$GEO" -o "$RAW" 2>&1 | tail -3
    echo "    wrote $RAW"
fi
FINAL_MESH="$RAW"
[ "$REFINE_LVL" -gt 0 ] && FINAL_MESH="$CASEDIR/${MESHBASE}_cols.msh"

echo
echo "=== 3/5  element ordering ==================================="
if [ "$REFINE_LVL" -gt 0 ]; then
    # p4est hands each rank a CONTIGUOUS range of the mesh's element order.
    # gmsh writes horizontal layers, so without this every rank would get a
    # horizontal slab and only the ones holding the bottom layer would own any
    # ground surface -- the wall model would run on a handful of ranks.
    julia --project=. --startup-file=no tools/reorder_msh_columns.jl \
          "$RAW" "$FINAL_MESH" "$RANKS"
else
    echo "    direct path: lxy_partition bins by coordinate, so element order"
    echo "    does not matter. Skipping reorder."
fi

echo
echo "=== 4/5  patching $CASEDIR/user_inputs.jl ==================="
UI="$CASEDIR/user_inputs.jl"
[ -f "$UI.orig" ] || cp "$UI" "$UI.orig"
# Replace the value of an ACTIVE (uncommented) key. Commented-out variants are
# left alone -- these decks carry several dead :gmsh_filename lines and a blind
# sed would rewrite the wrong one.
set_key() {
    local file="$1" key="$2" val="$3"
    awk -v k=":$key" -v v="$val" '
        {
          line=$0
          if (line ~ /^[[:space:]]*#/) { print line; next }
          if (index(line, k)==0)       { print line; next }
          if (line !~ ("^[[:space:]]*" k "[[:space:]]*=>")) { print line; next }
          match(line, /^[[:space:]]*/); ind=substr(line, 1, RLENGTH)
          printf "%s%s => %s,\n", ind, k, v
          found=1; next
        }
        END { if (!found) exit 3 }
    ' "$file" > "$file.tmp" || { echo "ERROR: key :$key not found in $file" >&2; rm -f "$file.tmp"; exit 1; }
    mv "$file.tmp" "$file"
    printf '    :%-18s => %s\n' "$key" "$val"
}
set_key "$UI" gmsh_filename   "\"./$FINAL_MESH\""
set_key "$UI" linitial_refine "$([ "$REFINE_LVL" -gt 0 ] && echo true || echo false)"
set_key "$UI" init_refine_lvl "$REFINE_LVL"
set_key "$UI" lxy_partition   "$([ "$REFINE_LVL" -gt 0 ] && echo false || echo true)"

julia --project=. --startup-file=no -e "
  f=\"$UI\"; Meta.parseall(read(f,String); filename=f); println(\"    user_inputs.jl parses OK\")"

echo
echo "=== 5/5  SLURM script ======================================="
SB="submit_${CASE}.sh"
sed -e "s|^#SBATCH --nodes=.*|#SBATCH --nodes=${NODES}|" \
    -e "s|^#SBATCH --ntasks-per-node=.*|#SBATCH --ntasks-per-node=${NTASKS_PER_NODE}|" \
    -e "s|^#SBATCH --job-name=.*|#SBATCH --job-name=${CASE}|" \
    -e "s|^#SBATCH --mem-per-cpu=.*|#SBATCH --mem-per-cpu=${MEM_PER_CPU}|" \
    -e "s|^CASE=\"\${CASE:-.*|CASE=\"\${CASE:-${CASE}}\"|" \
    submit_Jexpresso_precompile_pass.sh > "$SB"
chmod +x "$SB"
echo "    wrote $SB   (nodes=${NODES}, ntasks-per-node=${NTASKS_PER_NODE}, ranks=${RANKS}, mem-per-cpu=${MEM_PER_CPU})"

cat <<SUMEOF

-----------------------------------------------------------------
 case            : $EQS / $CASE
 solves on       : ${NX} x ${NY} x ${NZ}   (nop=${NOP})   <- $([ "$RUN_GMSH" = no ] && echo "derived from your mesh x2^${REFINE_LVL}" || echo "TARGET_GRID")
 mesh on disk    : $((NX/POW)) x $((NY/POW)) x $((NZ/POW))  ->  refined x2^${REFINE_LVL}
 path            : $([ "$REFINE_LVL" -gt 0 ] && echo "p4est, ${REFINE_LVL} refinement level(s)" || echo "direct read + lxy_partition")
 mesh            : $FINAL_MESH
 ranks           : ${RANKS}  (${NODES} x ${NTASKS_PER_NODE})
 gridpoints/rank : ${PTS_PER_RANK}
 ground/rank     : ${GROUND_PER_RANK}     <- must be > 0 on EVERY rank
 other valid rank counts: ${VALID_RANKS}
-----------------------------------------------------------------

Check these lines in the output when it runs:
   # Bottom faces (min/avg/max) : ${GROUND_PER_RANK} / ${GROUND_PER_RANK}.0 / ${GROUND_PER_RANK}
   # [mem] ... after ... bcast          <- should stay near baseline
SUMEOF

if [ "$SUBMIT" = "yes" ]; then
    echo; echo "submitting..."; sbatch "$SB"
else
    echo; echo "To submit:  sbatch $SB"
fi
