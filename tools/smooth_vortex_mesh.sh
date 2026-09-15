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
#     L = 10  ([-5,5]^2, the box of Balsara)            ~ 4.9e-06
#     L = 15                                             ~ 1.2e-12
#     L = 20  ([-10,10]^2, the box of Dao & Nazarov)     ~ 5.1e-22   <-- default
#
# So the default L = 20 is the box the published accuracy test is run on, and
# the one an error study needs: on L = 10 every order flattens at ~1e-5, the
# high orders first, because the DATA is discontinuous at that size.
#
#   SV_L=20 SV_NELX="16 32 64" tools/smooth_vortex_mesh.sh
#
# Meshes of the default box keep their plain name, vortex_<n>x<n>.msh; any
# other box is tagged, vortex_L20_<n>x<n>.msh, so the two never collide.
set -eu
cd "$(dirname "$0")/.."
DIR=problems/MHD/smoothVortex
NELX=${SV_NELX:-"4 8 16 32"}
L=${SV_L:-20}
GMSH=${GMSH:-gmsh}

# The box is ALWAYS in the name: vortex_L20_16x16.msh. A mesh whose name does
# not say which box it is cannot be told from one that is a different box.
LTAG="L$(printf '%s' "$L" | sed 's/\.0*$//')_"

for n in $NELX; do
    out="$DIR/vortex_${LTAG}${n}x${n}.msh"
    [ -f "$out" ] && { echo "have $out"; continue; }
    "$GMSH" -2 -setnumber nx "$n" -setnumber L "$L" "$DIR/vortex_periodic.geo" \
            -o "$out" -format msh41 -v 1
    echo "wrote $out"
done
