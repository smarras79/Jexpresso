#!/usr/bin/env bash
#
# Generate the periodic meshes the smooth-vortex convergence sweep needs,
# from the single parameterized problems/MHD/smoothVortex/vortex_periodic.geo.
#
#   tools/smooth_vortex_mesh.sh            # 4, 8, 16, 32 elements per side
#   SV_NELX="8 16 24" tools/smooth_vortex_mesh.sh
#
# SV_L is the BOX WIDTH, and it decides how far an accuracy study can go.
# The vortex is a Gaussian, not compactly supported, so on [-L/2, L/2]² the
# exact solution is not periodic: its velocity perturbation at the middle of
# an edge is (L/2)exp((1 - (L/2)²)/2)/2π, with OPPOSITE SIGN on the two
# opposite edges. The initial condition therefore jumps across the periodic
# seam by twice that, and no scheme of any order converges below it:
#
#     L = 10 (the box of Balsara and of Dao & Nazarov)   ~ 4.9e-06   <-- default
#     L = 15                                             ~ 1.2e-12
#     L = 20                                             ~ 5.1e-22   (machine zero)
#
# So L = 10 is right for reproducing the published figure and floors out at
# ~1e-5; a study that must show the decay of a 6th-order element down to
# 1e-8 or below needs SV_L=15 or SV_L=20 (and, for the same resolution, the
# proportionally larger element count).
#
#   SV_L=20 SV_NELX="16 32 64" tools/smooth_vortex_mesh.sh
#
# Meshes of the default box keep their plain name, vortex_<n>x<n>.msh; any
# other box is tagged, vortex_L20_<n>x<n>.msh, so the two never collide.
set -eu
cd "$(dirname "$0")/.."
DIR=problems/MHD/smoothVortex
NELX=${SV_NELX:-"4 8 16 32"}
L=${SV_L:-10}
GMSH=${GMSH:-gmsh}

# vortex_16x16.msh for the default box, vortex_L20_16x16.msh otherwise.
LTAG=""
[ "$(printf '%s' "$L" | sed 's/\.0*$//')" = "10" ] || LTAG="L$(printf '%s' "$L" | sed 's/\.0*$//')_"

for n in $NELX; do
    out="$DIR/vortex_${LTAG}${n}x${n}.msh"
    [ -f "$out" ] && { echo "have $out"; continue; }
    "$GMSH" -2 -setnumber nx "$n" -setnumber L "$L" "$DIR/vortex_periodic.geo" \
            -o "$out" -format msh41 -v 1
    echo "wrote $out"
done
