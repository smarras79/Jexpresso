#!/usr/bin/env bash
# =============================================================================
#  run_element_learning.sh  —  end-to-end Element-Learning pipeline driver
# =============================================================================
#
#  Automates the three manual steps of an Element-Learning (EL) run so a user
#  can launch the whole thing with a single command and collect the results:
#
#    Phase 1  SAMPLE     Run the SEM solver in EL sampling mode
#                        (:lEL_Sample => true) on the 1x1 single-element mesh.
#                        → produces  input_tensor_<TAG>.csv  and
#                          output_tensor_<TAG>.csv  (TAG identifies the test so
#                          different tests / grids never clobber each other)
#
#    Phase 2  TRAIN      Hand the two tagged CSVs to the Python trainer (default:
#                        the self-contained in-repo tools/EL_training/train_CNN.py)
#                        and collect the exported ONNX model.
#                        → produces  <NNFILE>  (e.g. JX_NN_<TAG>_model.onnx)
#
#    Phase 3  INFER      Re-run the SEM solver in inference mode
#                        (:lEL_Sample => false) on the NxN (multi-element) mesh,
#                        loading the freshly trained model.
#                        → writes the solution to the case output directory
#
#  The three phases run as SEPARATE processes, exactly as they would be run by
#  hand. Phase switching (lEL_Sample, mesh, model file) is done through the
#  JEXPRESSO_EL_* environment variables that src/io/mod_inputs.jl reads, so the
#  case's user_inputs.jl is NEVER edited or committed by this script.
#
#  --------------------------------------------------------------------------
#  USAGE
#  --------------------------------------------------------------------------
#    tools/EL_training/run_element_learning.sh [options]
#
#  Options:
#    -c, --config FILE   Load settings from FILE (default:
#                        tools/EL_training/element_learning.config if present).
#    --sample-only       Run only Phase 1 (sampling).
#    --train-only        Run only Phase 2 (training).
#    --infer-only        Run only Phase 3 (inference).
#    --skip-sample       Skip Phase 1 (reuse existing *_tensor.csv).
#    --skip-train        Skip Phase 2 (reuse existing model file).
#    --skip-infer        Skip Phase 3.
#    -h, --help          Show this help and exit.
#
#  Any configuration variable (see the CONFIGURATION block below) may also be
#  overridden from the environment, e.g.:
#
#    CASE=elementLearning_hole NSAMP=2000 \
#        tools/EL_training/run_element_learning.sh
#
#  --------------------------------------------------------------------------
#  BENCHMARK (known-good)
#  --------------------------------------------------------------------------
#    tools/EL_training/run_element_learning.sh
#  runs problems/Elliptic/elementLearning_hole end-to-end with the defaults
#  below.
# =============================================================================

set -euo pipefail

# ─────────────────────────────────────────────────────────────────────────────
# Resolve paths. This script lives in <repo>/tools/EL_training, so the repo
# root (the directory Julia is launched from — where the *_tensor.csv files and
# the NNFILE live) is two levels up.
# ─────────────────────────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/../.." && pwd)"

# ─────────────────────────────────────────────────────────────────────────────
# Locate and source the config file FIRST, so precedence is:
#     environment  >  config file  >  built-in defaults below.
# The config file should assign with the `: "${VAR:=...}"` form so that a value
# already set in the environment is never clobbered.
# ─────────────────────────────────────────────────────────────────────────────
CONFIG_FILE="${SCRIPT_DIR}/element_learning.config"
_args=("$@")
for ((i=0; i<${#_args[@]}; i++)); do
    case "${_args[$i]}" in
        -c|--config) CONFIG_FILE="${_args[$((i+1))]:-}";;
    esac
done
if [[ -n "${CONFIG_FILE}" && -f "${CONFIG_FILE}" ]]; then
    echo " # [EL pipeline] loading config: ${CONFIG_FILE}"
    # shellcheck disable=SC1090
    source "${CONFIG_FILE}"
fi

# ═════════════════════════════════════════════════════════════════════════════
# CONFIGURATION  (defaults target the elementLearning_hole benchmark)
# Every value can be overridden by a --config file or by the environment.
# ═════════════════════════════════════════════════════════════════════════════

# --- Case ------------------------------------------------------------------
: "${EQ:=Elliptic}"                       # problems/<EQ>/<CASE>
: "${CASE:=elementLearning_hole}"

# --- Per-test tag ----------------------------------------------------------
# Tags the sampled tensors and the trained model so that different tests / grids
# never overwrite each other's files. Defaults to a sanitised "<EQ>_<CASE>";
# override EL_TAG (e.g. to add the polynomial order or grid) if you run the same
# case at different resolutions.
if [[ -z "${EL_TAG:-}" ]]; then
    _tag_raw="${EQ}_${CASE}"
    EL_TAG="${_tag_raw//[^A-Za-z0-9]/_}"
fi

# --- Meshes ----------------------------------------------------------------
# Sampling learns the reference-element operator from a single element (1x1).
: "${EL_SAMPLE_MESH:=./meshes/gmsh_grids/square_dirichletT_1x1.msh}"
# Inference runs on the real, multi-element mesh. Leave EMPTY to use whatever
# :gmsh_filename the case's user_inputs.jl already specifies (for the hole
# benchmark that is plate_hole_circle_unit.msh).
: "${EL_INFER_MESH:=}"

# --- Number of sampling draws (empty → use the case's :Nsamp) --------------
: "${NSAMP:=}"

# --- Tagged tensor + model filenames ---------------------------------------
# The sampler (Julia) writes these names; the trainer (Python) reads them via
# EL_INPUT_TENSOR / EL_OUTPUT_TENSOR. All live at the repo root.
: "${EL_INPUT_TENSOR:=input_tensor_${EL_TAG}.csv}"
: "${EL_OUTPUT_TENSOR:=output_tensor_${EL_TAG}.csv}"
# Model basename the trainer builds its outputs from (dataname); the exported
# ONNX is "<EL_DATANAME>_model.onnx".
: "${EL_DATANAME:=JX_NN_${EL_TAG}}"
# Model file (relative to REPO_ROOT) that Julia loads as :NNfile at inference.
: "${NNFILE:=${EL_DATANAME}_model.onnx}"

# --- Python trainer --------------------------------------------------------
# TRAIN_DIR is the directory the trainer is launched from (its CWD). The CSVs
# are staged there and the trainer reads them from its CWD; the exported ONNX is
# copied back to the repo root afterwards.
#
# Default is the SELF-CONTAINED in-repo trainer: tools/EL_training already ships
# train_CNN.py plus all of its dependencies (train_common_EL.py, IO_EL.py,
# NN_EL.py, SLmodel_EL.py), so the pipeline runs out of the box. To use your own
# external trainer instead, set TRAIN_DIR (e.g. to an EL_Jexpresso directory).
: "${TRAIN_DIR:=${REPO_ROOT}/tools/EL_training}"
: "${TRAIN_CMD:=python train_CNN.py}"     # command run inside TRAIN_DIR
: "${TRAIN_ONNX:=${NNFILE}}"              # model file the trainer writes in TRAIN_DIR
# Fallback trainer directory, used only if TRAIN_DIR does not exist.
: "${TRAIN_DIR_FALLBACK:=${REPO_ROOT}/tools/EL_training}"
: "${TRAIN_CMD_FALLBACK:=python train_CNN.py}"

# --- Julia launch ----------------------------------------------------------
: "${JULIA:=julia}"
: "${JULIA_PROJECT:=.}"
: "${NPROCS:=1}"                          # >1 → mpiexec -n NPROCS
: "${MPIEXEC:=mpiexec}"

# ─────────────────────────────────────────────────────────────────────────────
# Parse command line (a --config is applied BEFORE the phase flags so config
# never clobbers an explicit --*-only choice).
# ─────────────────────────────────────────────────────────────────────────────
DO_SAMPLE=1; DO_TRAIN=1; DO_INFER=1

# phase-selection flags (--config was already consumed above)
i=0
while [[ $i -lt ${#_args[@]} ]]; do
    case "${_args[$i]}" in
        -c|--config)  i=$((i+1));;                       # value already consumed
        --sample-only) DO_SAMPLE=1; DO_TRAIN=0; DO_INFER=0; _explicit_only=1;;
        --train-only)  DO_SAMPLE=0; DO_TRAIN=1; DO_INFER=0; _explicit_only=1;;
        --infer-only)  DO_SAMPLE=0; DO_TRAIN=0; DO_INFER=1; _explicit_only=1;;
        --skip-sample) DO_SAMPLE=0;;
        --skip-train)  DO_TRAIN=0;;
        --skip-infer)  DO_INFER=0;;
        -h|--help)     sed -n '2,58p' "${BASH_SOURCE[0]}"; exit 0;;
        *) echo "ERROR: unknown option '${_args[$i]}'" >&2; exit 2;;
    esac
    i=$((i+1))
done

# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────
_hr()   { printf '=%.0s' {1..79}; printf '\n'; }
_step() { _hr; echo " # $*"; _hr; }

# Launch Jexpresso.run_case for the current case, serial or under MPI.
run_case_julia() {
    cd "${REPO_ROOT}"
    local jl="using Jexpresso; Jexpresso.run_case(\"${EQ}\", \"${CASE}\")"
    if [[ "${NPROCS}" -gt 1 ]]; then
        "${MPIEXEC}" -n "${NPROCS}" "${JULIA}" --project="${JULIA_PROJECT}" -e "${jl}"
    else
        "${JULIA}" --project="${JULIA_PROJECT}" -e "${jl}"
    fi
}

# ─────────────────────────────────────────────────────────────────────────────
# Phase 1 — SAMPLE
# ─────────────────────────────────────────────────────────────────────────────
phase_sample() {
    _step "Phase 1/3  SAMPLE   ${EQ}/${CASE}   mesh=${EL_SAMPLE_MESH}   tag=${EL_TAG}"
    (
        export JEXPRESSO_EL_SAMPLE=true
        export JEXPRESSO_EL_MESH="${EL_SAMPLE_MESH}"
        export JEXPRESSO_EL_TAG="${EL_TAG}"
        [[ -n "${NSAMP}" ]] && export JEXPRESSO_EL_NSAMP="${NSAMP}"
        run_case_julia
    )
    for f in "${EL_INPUT_TENSOR}" "${EL_OUTPUT_TENSOR}"; do
        if [[ ! -s "${REPO_ROOT}/${f}" ]]; then
            echo "ERROR: sampling did not produce ${f}" >&2; exit 1
        fi
    done
    echo " # [EL pipeline] sampling wrote ${EL_INPUT_TENSOR} / ${EL_OUTPUT_TENSOR}"
}

# ─────────────────────────────────────────────────────────────────────────────
# Phase 2 — TRAIN  (external EL_Jexpresso/train_CNN.py, or the in-repo trainer)
# ─────────────────────────────────────────────────────────────────────────────
phase_train() {
    local tdir tcmd tonnx
    if [[ -d "${TRAIN_DIR}" ]]; then
        tdir="${TRAIN_DIR}"; tcmd="${TRAIN_CMD}"; tonnx="${TRAIN_ONNX}"
    elif [[ -d "${TRAIN_DIR_FALLBACK}" ]]; then
        echo " # [EL pipeline] TRAIN_DIR '${TRAIN_DIR}' not found — using in-repo trainer"
        tdir="${TRAIN_DIR_FALLBACK}"; tcmd="${TRAIN_CMD_FALLBACK}"; tonnx="${TRAIN_ONNX}"
    else
        echo "ERROR: no training directory found (TRAIN_DIR='${TRAIN_DIR}')." >&2
        echo "       Set TRAIN_DIR to your EL_Jexpresso directory in the config." >&2
        exit 1
    fi

    _step "Phase 2/3  TRAIN    ${tcmd}   (cwd=${tdir})   tag=${EL_TAG}"

    for f in "${EL_INPUT_TENSOR}" "${EL_OUTPUT_TENSOR}"; do
        if [[ ! -s "${REPO_ROOT}/${f}" ]]; then
            echo "ERROR: ${f} missing — run the sampling phase first." >&2; exit 1
        fi
        # Stage the CSVs into the trainer's working directory. Skip the copy if
        # the trainer already reads from the repo root (same dir).
        if [[ "$(cd "${tdir}" && pwd)" != "${REPO_ROOT}" ]]; then
            cp -f "${REPO_ROOT}/${f}" "${tdir}/${f}"
        fi
    done

    # train_CNN.py opens an interactive matplotlib window (plt.ion/plt.show);
    # force a non-interactive backend so the pipeline never blocks or fails on a
    # headless machine. Users can override by exporting MPLBACKEND themselves.
    # EL_INPUT_TENSOR / EL_OUTPUT_TENSOR / EL_DATANAME tell the trainer which
    # tagged files to read and what to name its outputs.
    (
        cd "${tdir}"
        export MPLBACKEND="${MPLBACKEND:-Agg}"
        export EL_INPUT_TENSOR="${EL_INPUT_TENSOR}"
        export EL_OUTPUT_TENSOR="${EL_OUTPUT_TENSOR}"
        export EL_DATANAME="${EL_DATANAME}"
        eval "${tcmd}"
    )

    if [[ ! -s "${tdir}/${tonnx}" ]]; then
        echo "ERROR: trainer did not produce '${tonnx}' in ${tdir}" >&2; exit 1
    fi
    # Publish the model where the inference phase (Julia CWD = REPO_ROOT) reads
    # :NNfile from.
    if [[ "${tdir}/${tonnx}" != "${REPO_ROOT}/${NNFILE}" ]]; then
        cp -f "${tdir}/${tonnx}" "${REPO_ROOT}/${NNFILE}"
    fi
    echo " # [EL pipeline] trained model ready: ${REPO_ROOT}/${NNFILE}"
}

# ─────────────────────────────────────────────────────────────────────────────
# Phase 3 — INFER
# ─────────────────────────────────────────────────────────────────────────────
phase_infer() {
    if [[ ! -s "${REPO_ROOT}/${NNFILE}" ]]; then
        echo "ERROR: model file '${NNFILE}' not found — run the training phase first." >&2
        exit 1
    fi
    _step "Phase 3/3  INFER    ${EQ}/${CASE}   model=${NNFILE}   mesh=${EL_INFER_MESH:-<case default>}"
    (
        export JEXPRESSO_EL_SAMPLE=false
        export JEXPRESSO_EL_NNFILE="${NNFILE}"
        # Inference must run on the MULTI-element mesh, never the 1x1 sampling
        # mesh. If EL_INFER_MESH is set, use it; otherwise fall back to the
        # case's user_inputs.jl :gmsh_filename by making sure no stray
        # JEXPRESSO_EL_MESH (e.g. exported in the parent shell) leaks in.
        if [[ -n "${EL_INFER_MESH}" ]]; then
            export JEXPRESSO_EL_MESH="${EL_INFER_MESH}"
        else
            unset JEXPRESSO_EL_MESH
        fi
        run_case_julia
    )
    echo " # [EL pipeline] inference complete — see the case output directory."
}

# ─────────────────────────────────────────────────────────────────────────────
# Drive the pipeline
# ─────────────────────────────────────────────────────────────────────────────
echo " # [EL pipeline] repo root : ${REPO_ROOT}"
echo " # [EL pipeline] case      : ${EQ}/${CASE}"
echo " # [EL pipeline] tag       : ${EL_TAG}"
echo " # [EL pipeline] tensors   : ${EL_INPUT_TENSOR} / ${EL_OUTPUT_TENSOR}"
echo " # [EL pipeline] model     : ${NNFILE}"
echo " # [EL pipeline] phases    : sample=${DO_SAMPLE} train=${DO_TRAIN} infer=${DO_INFER}"

[[ "${DO_SAMPLE}" -eq 1 ]] && phase_sample
[[ "${DO_TRAIN}"  -eq 1 ]] && phase_train
[[ "${DO_INFER}"  -eq 1 ]] && phase_infer

_step "Element-Learning pipeline DONE"
