#!/bin/bash
# Compile the Alya proxy (myAlya.f90 -> alya_all2all_time_loop.f90) against MPI.
#
# IMPORTANT: the proxy MUST be built with the SAME MPI that Julia's MPI.jl is
# bound to — same implementation, version, and ABI (see ../RUN-COUPLED.md).
#
# By default this uses the `mpif90` first on your PATH. To pin a specific MPI
# (recommended when several are installed), set MPIF90 to its wrapper, e.g.:
#
#   MPIF90=/opt/homebrew/Cellar/open-mpi/5.0.8/bin/mpif90 bash compilef90.sh
#
set -e

MPIF90="${MPIF90:-mpif90}"

echo "==> Fortran MPI wrapper: $MPIF90"
echo "==> It links against:"
"$MPIF90" -show        # shows the include/lib paths — confirm they match MPI.jl
echo "==> Compiling Alya.x ..."
"$MPIF90" -cpp -DUSEMPIF08 myAlya.f90 -o Alya.x
echo "==> Built ./Alya.x"
