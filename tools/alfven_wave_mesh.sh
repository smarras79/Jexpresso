#!/usr/bin/env bash
#
# Generate the periodic meshes of the smooth Alfvén wave convergence test
# (Babbar & Ranocha, J. Comput. Phys. 563 (2026) 115060, §5.5.1) from the single
# parameterized problems/MHD/alfvenWaveBabbarRanocha2026/alfven_periodic.geo:
# nx x nx quads on [0, 1/cos α] x [0, 1/sin α] = [0, √2]².
#
#   tools/alfven_wave_mesh.sh                  # 4, 8, 16, 32, 64, 128 per side
#   AW_NELX="256" tools/alfven_wave_mesh.sh
#
# 4² to 128² are in the repository; 16² to 128² are the meshes of the paper's
# Fig. 13. A mesh that exists is left alone.
set -eu
cd "$(dirname "$0")/.."
DIR=problems/MHD/alfvenWaveBabbarRanocha2026
NELX=${AW_NELX:-"4 8 16 32 64 128"}
GMSH=${GMSH:-gmsh}

for n in $NELX; do
    out="$DIR/alfven_${n}x${n}.msh"
    [ -f "$out" ] && { echo "have $out"; continue; }
    "$GMSH" -2 -setnumber nx "$n" "$DIR/alfven_periodic.geo" -o "$out" -format msh41 -v 1
    echo "wrote $out"
done
