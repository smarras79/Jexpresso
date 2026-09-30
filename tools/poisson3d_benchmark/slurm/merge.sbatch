#!/bin/bash -l
#SBATCH --job-name=p3d-merge
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
# Runs after every benchmark job has ended (afterany): merges the parts into
# $OUTDIR/results.csv and results.md and draws the figures into $OUTDIR/assets.
source "$P3D_DIR/common.sh"      # P3D_DIR: this directory, exported by submit_wulver3d.sh
$PYTHON tools/poisson3d_benchmark/plot3d.py "$OUTDIR"
