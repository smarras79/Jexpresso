#!/usr/bin/env bash
#
# Generate the periodic meshes the smooth-vortex convergence sweep needs,
# from the single parameterized problems/MHD/smoothVortex/vortex_periodic.geo.
#
#   tools/smooth_vortex_mesh.sh            # 4, 8, 16, 32 elements per side
#   SV_NELX="8 16 24" tools/smooth_vortex_mesh.sh
set -eu
cd "$(dirname "$0")/.."
DIR=problems/MHD/smoothVortex
NELX=${SV_NELX:-"4 8 16 32"}
GMSH=${GMSH:-gmsh}
for n in $NELX; do
    out="$DIR/vortex_${n}x${n}.msh"
    [ -f "$out" ] && { echo "have $out"; continue; }
    "$GMSH" -2 -setnumber nx "$n" "$DIR/vortex_periodic.geo" -o "$out" -format msh41 -v 1
    echo "wrote $out"
done
