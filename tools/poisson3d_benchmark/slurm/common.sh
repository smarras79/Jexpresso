# Environment of the 3D benchmark jobs (sourced; needs P3D_ENV and P3D_DIR).
set -eo pipefail          # no -u: module scripts use unset variables
source "$P3D_ENV"
if [ -n "${JULIA_MODULE}" ]; then
    module load wulver 2>/dev/null || true
    module load ${JULIA_MODULE}
fi
if [ -n "${JULIA_DEPOT}" ]; then export JULIA_DEPOT_PATH="${JULIA_DEPOT}"; fi
export JULIA_PKG_OFFLINE=true             # compute nodes: never touch the registry
# SLURM >= 22.05: srun no longer inherits --cpus-per-task from sbatch
if [ -n "${SLURM_CPUS_PER_TASK:-}" ]; then export SRUN_CPUS_PER_TASK="$SLURM_CPUS_PER_TASK"; fi
cd "$REPO"
