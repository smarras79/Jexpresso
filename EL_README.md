# Element Learning (EL)

Element Learning trains a small neural-network surrogate for the element-level
Laplace/Poisson operator and then uses it in place of the direct SEM solve. A
full run has three steps:

1. **Sample** — run the SEM solver in sampling mode (`:lEL_Sample => true`) on a
   single-element `1x1` mesh. This writes the training data as
   `input_tensor_<TAG>.csv` and `output_tensor_<TAG>.csv`, where `<TAG>`
   identifies the test/grid so runs on different cases never overwrite each
   other's data.
2. **Train** — a Python trainer turns those two CSVs into an ONNX model
   (`JX_NN_<TAG>_model.onnx`). A self-contained copy of the trainer
   (`train_CNN.py` and its dependencies) ships in `tools/EL_training`, so no
   external directory is required; you can still point at your own trainer if
   you prefer.
3. **Infer** — re-run the SEM solver in inference mode (`:lEL_Sample => false`)
   on the real multi-element (`NxN`) mesh, loading the trained model as
   `:NNfile`.

---

## Automated pipeline (recommended)

`tools/EL_training/run_element_learning.sh` runs all three steps as separate
processes with a single command. Phase switching (`lEL_Sample`, mesh, model
file) is done through `JEXPRESSO_EL_*` environment variables read by
`src/io/mod_inputs.jl`, so **your case's `user_inputs.jl` is never edited**.

### Run the benchmark

The defaults target `problems/Elliptic/elementLearning_hole`, a known-good case:

```bash
# from the Jexpresso repo root
tools/EL_training/run_element_learning.sh
```

That runs, in order:

* **SAMPLE** — `:lEL_Sample=>true` on `square_dirichletT_1x1.msh`
  → `input_tensor.csv`, `output_tensor.csv`
* **TRAIN**  — stages the CSVs into the trainer directory (`tools/EL_training`
  by default), runs `python train_CNN.py`, copies the resulting
  `JX_NN_model.onnx` back to the repo root
* **INFER**  — `:lEL_Sample=>false` on the case's real mesh, loading
  `JX_NN_model.onnx`; the solution is written to the case output directory

### Configure it

Copy the template and edit what you need (the driver loads
`element_learning.config` automatically):

```bash
cp tools/EL_training/element_learning.config.example \
   tools/EL_training/element_learning.config
```

Key settings (see the template for the full list):

| Variable         | Meaning                                           | Default |
|------------------|---------------------------------------------------|---------|
| `EQ` / `CASE`    | `problems/<EQ>/<CASE>` to run                      | `Elliptic` / `elementLearning_hole` |
| `EL_SAMPLE_MESH` | single-element mesh used for sampling              | `.../square_dirichletT_1x1.msh` |
| `EL_INFER_MESH`  | multi-element mesh for inference (empty → case default) | *(empty)* |
| `NSAMP`          | number of sampling draws (empty → case `:Nsamp`)   | *(empty)* |
| `NNFILE`         | trained model file (`:NNfile`)                     | `JX_NN_model.onnx` |
| `TRAIN_DIR`      | directory the trainer runs in                      | `tools/EL_training` (in-repo) |
| `TRAIN_CMD`      | training command run inside `TRAIN_DIR`            | `python train_CNN.py` |
| `NPROCS`         | `>1` runs the SEM phases under `mpiexec -n NPROCS` | `1` |

Any variable can also be set inline for a one-off run:

```bash
# quick smoke test with fewer samples
NSAMP=2000 tools/EL_training/run_element_learning.sh
```

### Run individual phases

```bash
tools/EL_training/run_element_learning.sh --sample-only   # step 1 only
tools/EL_training/run_element_learning.sh --train-only     # step 2 only
tools/EL_training/run_element_learning.sh --infer-only     # step 3 only
tools/EL_training/run_element_learning.sh --skip-sample    # reuse existing CSVs
```

`run_element_learning.sh --help` lists every option.

> **Trainer.** `TRAIN_DIR` defaults to the in-repo `tools/EL_training`, which is
> self-contained: it ships `train_CNN.py` plus its dependencies
> (`train_common_EL.py`, `IO_EL.py`, `NN_EL.py`, `SLmodel_EL.py`). Install the
> Python requirements once with
> `pip install torch scikit-learn scipy matplotlib onnx`. To use an external
> trainer (e.g. an `EL_Jexpresso` directory) set `TRAIN_DIR` in the config.

---

## Diagnostics (printed at the end of inference)

After the inference run, a consolidated block reports the solve time and the
accuracy of the inferred solution:

```
 # ================== ELEMENT-LEARNING DIAGNOSTICS ==================
 #   time      : inference = 0.0123 s | direct SEM = 0.456 s | speedup = 37.1×
 #   accuracy  : inference vs numerical (direct SEM)  →  ‖e‖_L2 = … , rel = … , ‖e‖_∞ = …
 #   accuracy  : inference vs exact (manufactured)    →  ‖e‖_L2 = … , rel = … , ‖e‖_∞ = …
 #   accuracy  : direct SEM vs exact (manufactured)   →  ‖e‖_L2 = … , rel = … , ‖e‖_∞ = …
 # ==================================================================
```

* The **direct SEM** solve (`A \ RHS`) is the always-available numerical
  reference; the inferred solution is compared against it for both time
  (speedup) and accuracy.
* The **exact (manufactured)** rows appear only when the case stores an exact
  field `qe` (e.g. an MMS test); both the inference and the direct SEM errors
  against it are printed so their accuracy can be compared directly.
* Norms are mass-matrix-weighted L2 (absolute + relative) and L∞.

Set `:lEL_diagnostics => false` in the case's `user_inputs.jl` to skip the extra
direct solve (useful on very large meshes where the direct solve is exactly what
element learning is avoiding).

---

## Manual steps (equivalent to what the script automates)

Set `JEXPRESSO_EL_TAG` so the sampler tags its output (here `hole`):

```bash
# 1) Sample: :lEL_Sample => true, 1x1 mesh
JEXPRESSO_EL_TAG=hole JEXPRESSO_EL_SAMPLE=true \
JEXPRESSO_EL_MESH=./meshes/gmsh_grids/square_dirichletT_1x1.msh \
    julia --project=. -e 'using Jexpresso; Jexpresso.run_case("Elliptic","elementLearning_hole")'
#      → input_tensor_hole.csv, output_tensor_hole.csv

# 2) Train
cd tools/EL_training               # (or your own external trainer directory)
cp ../../input_tensor_hole.csv ../../output_tensor_hole.csv .
EL_INPUT_TENSOR=input_tensor_hole.csv EL_OUTPUT_TENSOR=output_tensor_hole.csv \
EL_DATANAME=JX_NN_hole python train_CNN.py     # → JX_NN_hole_model.onnx
cp JX_NN_hole_model.onnx ../../    # place it where :NNfile is read from

# 3) Infer: :lEL_Sample => false, NxN mesh
JEXPRESSO_EL_SAMPLE=false JEXPRESSO_EL_NNFILE=JX_NN_hole_model.onnx \
    julia --project=. -e 'using Jexpresso; Jexpresso.run_case("Elliptic","elementLearning_hole")'
```

The `JEXPRESSO_EL_*` overrides used by the script (all optional; unset values
leave `user_inputs.jl` untouched):

| Variable                  | Overrides         | Effect                                   |
|---------------------------|-------------------|------------------------------------------|
| `JEXPRESSO_EL_SAMPLE`     | `:lEL_Sample`     | sampling (`true`) vs inference (`false`) |
| `JEXPRESSO_EL_MESH`       | `:gmsh_filename`  | mesh for this phase                      |
| `JEXPRESSO_EL_NNFILE`     | `:NNfile`         | trained model to load at inference       |
| `JEXPRESSO_EL_OUTPUT_DIR` | `:output_dir`     | output directory                         |
| `JEXPRESSO_EL_NSAMP`      | `:Nsamp`          | number of sampling draws                 |
| `JEXPRESSO_EL_TAG`        | `:EL_tensor_tag`  | tags the tensor CSV filenames per test   |

The Python trainer reads `EL_INPUT_TENSOR` / `EL_OUTPUT_TENSOR` (which CSVs to
load) and `EL_DATANAME` (model output basename); the driver sets all three so
the trainer picks up the tagged files automatically.
