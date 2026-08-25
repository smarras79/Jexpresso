#!/usr/bin/env bash
#
# run_rtb2d.sh -- the desktop A/B for the 2D (slab) rising thermal bubble.
#
#     ./problems/CompEuler/rtb2d_schur/run_rtb2d.sh full        five fields, 5*Np
#     ./problems/CompEuler/rtb2d_schur/run_rtb2d.sh schur       scalar Schur, fast H
#     ./problems/CompEuler/rtb2d_schur/run_rtb2d.sh schur-ref   scalar Schur, reference H
#     ./problems/CompEuler/rtb2d_schur/run_rtb2d.sh all         all three, in order
#
# No SLURM, no cluster. This is submit_Jexpresso_profile.sh's job, done with
# mpiexec on whatever machine you are sitting at.
#
# THE ARGUMENT IS REQUIRED AND A TYPO EXITS 2. Running the baseline by accident
# and calling it the Schur arm is the exact failure this guards: half an A/B
# that turns out to be two copies of the same thing looks like a result.
#
# Environment knobs, all optional:
#   NR=4            ranks. See below -- 4 is the default for a REASON.
#   DBG_MESH        20x1x80 (default) or 10x1x40
#   DBG_TEND        model seconds (default 60, enough for one profile block)
#   DBG_DT          time step (default 0.6)
#   DBG_RTOL        Krylov tolerance (default 1e-8)
#
# WHY NR DEFAULTS TO 4 AND NOT TO YOUR CORE COUNT
# -----------------------------------------------
# :lxy_partition never cuts z, and this slab is one element thick in y, so it
# has 20 element columns and that is the hard ceiling on ranks. Of the counts
# that leave no rank empty (1, 2, 4, 5, 10 -- tools/pick_nranks.jl), 4 gives
# 32,501 points and 5 columns per rank at halo/elem 2.4. At 10 ranks each rank
# owns 2 columns and 13,000 points at halo/elem 3.0: the run goes
# communication-bound and the profile stops being about the stage solve, which
# is the thing being measured. More cores than 4 are better spent on threads
# (JULIA_NUM_THREADS) or on the 3D case.
set -u

usage() {
    echo "usage: $0 {full|schur|schur-ref|all}" >&2
    exit 2
}
[ $# -ge 1 ] || usage

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../../.." && pwd)"
cd "$ROOT" || exit 1

NR="${NR:-4}"
export DBG_MESH="${DBG_MESH:-20x1x80}"
export DBG_TEND="${DBG_TEND:-60.0}"
export DBG_DT="${DBG_DT:-0.6}"
export DBG_RTOL="${DBG_RTOL:-1.0e-8}"
# Ask for the per-step breakdown. Without this the run prints a step time and
# nothing about where it went, and the SPLIT is the part that transfers between
# this case and the 3D one.
export JEXPRESSO_HEVI_PROFILE="${JEXPRESSO_HEVI_PROFILE:-1}"

# LAUNCHING. mpiexec is not always on PATH even where MPI.jl works, because
# MPI.jl may be using a launcher from its own artifact.
#
# DO NOT try to turn MPI.mpiexec() into a path. It returns a Cmd that carries an
# ENVIRONMENT (the artifact's LD_LIBRARY_PATH among it), so printing it yields
# the whole `setenv(\`...\`, [...])` form -- several kB of environment that the
# shell then tries to exec as a filename ("File name too long"). Run the Cmd
# from inside Julia instead, which is also the only way that environment
# survives. The exported DBG_* variables reach the ranks either way: the Cmd
# captures the current environment.
launch() {
    if command -v mpiexec >/dev/null 2>&1; then
        mpiexec -n "$NR" julia --project=. src/Jexpresso.jl CompEuler rtb2d_schur
    else
        julia --project=. -e '
            using MPI
            run(`$(MPI.mpiexec()) -n '"$NR"' $(Base.julia_cmd()[1]) --project=. src/Jexpresso.jl CompEuler rtb2d_schur`)'
    fi
}

one_arm() {
    local name="$1" schur="$2" kern="$3"
    echo
    echo "=================================================================="
    echo " $name   ($NR ranks, mesh $DBG_MESH, tend $DBG_TEND s)"
    echo "   DBG_SCHUR=$schur  DBG_SCHUR_KERN=$kern"
    echo "=================================================================="
    export DBG_SCHUR="$schur" DBG_SCHUR_KERN="$kern"
    launch
    local rc=$?
    # Report the launcher's status. Without this the script looks identical
    # whether the run finished or every rank was killed.
    [ $rc -eq 0 ] || echo "*** $name FAILED (exit $rc) ***" >&2
    return $rc
}

case "$1" in
    full)      one_arm "five-field baseline"      0 1 ;;
    schur)     one_arm "scalar Schur, kernel H"   1 1 ;;
    schur-ref) one_arm "scalar Schur, reference H" 1 0 ;;
    all)
        # Sequential, never concurrent: two runs sharing the cores would time
        # each other's contention rather than the solvers.
        one_arm "five-field baseline"       0 1
        one_arm "scalar Schur, reference H" 1 0
        one_arm "scalar Schur, kernel H"    1 1
        ;;
    *) usage ;;
esac
