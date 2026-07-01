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

## Diagnostics

### Accuracy + solve timing (printed at the end of inference)

After the inference run, a consolidated block reports the solve time and the
accuracy of the inferred solution:

```
 # ================== ELEMENT-LEARNING DIAGNOSTICS ==================
 #   @btime minima (compilation excluded). The global matrix A is assembled
 #   once beforehand and shared by both methods (counted for neither).
 #   THIS SOLVE (single, time-independent):  EL (block-extract+infer) = 0.116 s | direct SEM (factorize+solve) = 0.0547 s | speedup = 0.47×
 #   amortized per solve (if reusing A):     EL inference (surrogate)  = 0.0227 s | direct SEM (back-solve) = 0.0031 s | speedup = 0.14×
 #   accuracy  : inference vs numerical (direct SEM)  →  ‖e‖_L2 = … , rel = … , ‖e‖_∞ = …
 #   accuracy  : inference vs exact (manufactured)    →  ‖e‖_L2 = … , rel = … , ‖e‖_∞ = …
 #   accuracy  : direct SEM vs exact (manufactured)   →  ‖e‖_L2 = … , rel = … , ‖e‖_∞ = …
 # ==================================================================
```
(illustrative numbers)

**What the two timing rows mean.** The **global matrix `A`** is assembled once by
the SEM setup, before either solve, and is shared — it is counted for neither
method. What differs is what each method does *with* `A`:

| Row | Element learning | Direct SEM |
|-----|------------------|------------|
| **THIS SOLVE** (the one solve that actually runs) | extract per-element blocks from `A` + surrogate = the full `elementLearning_Axb!` | factorize `A` + solve = `A \ RHS` |
| **amortized per solve** (only if many solves reuse the same `A`) | surrogate only (`elementLearning_infer!`) | back-solve with a pre-computed `lu(A)` factor |

Because this is a **time-independent** problem the solver runs **once**, so the
top row (**THIS SOLVE**) is the honest comparison: each method includes its own
`A`-dependent setup — EL's per-element block extraction (`elementLearning_Axb!`
Sections 1–2) and SEM's sparse LU factorization. The bottom row only applies if
you solve repeatedly with the *same* operator `A` (e.g. time-stepping or many
right-hand sides): then each method's `A`-dependent setup is done once and reused,
so per solve you pay only the EL surrogate / the SEM back-substitution. Each row's
speedup is `SEM / EL` (>1 means EL is faster).

> An earlier version timed the EL **surrogate** (block extraction excluded)
> against the SEM **full `A\RHS`** (factorization included) — not symmetric, and
> it flattered EL. The two rows above fix that by splitting both methods the same
> way.

* The **exact (manufactured)** accuracy rows appear only when the case stores an
  exact field `qe` (e.g. an MMS test); both the inference and the direct SEM
  errors against it are printed so their accuracy can be compared directly.
* Norms are mass-matrix-weighted L2 (absolute + relative) and L∞.

**Robust timing (`@btime`).** A single `@time`/`time_ns` shot is noisy — it
captures whatever GC pause, JIT compilation, or OS hiccup happened on that one
run, so repeated measurements scatter widely. The diagnostics instead time with
**BenchmarkTools' `@belapsed`** (the same engine as `@btime`): it runs each
operation many times over a time budget, excludes compilation, and reports the
**minimum**, which is far more repeatable. So a single hand-run (one cold solve)
will read noticeably higher and jump around; the diagnostics numbers are the
steady-state minima. BenchmarkTools is loaded lazily, only when the diagnostics
actually time something, so it is not pulled into the baseline of every run.

Set `:lEL_diagnostics => false` in the case's `user_inputs.jl` to skip the
timing/comparison entirely (useful on very large meshes where the direct solve
is exactly what element learning is avoiding); `:EL_timing_seconds => S`
(default 2.0) tunes the per-solve BenchmarkTools time budget.

The **same robust `@btime` timing applies to all the Laplace comparison solvers**
— direct SEM (`A\RHS` in `standard_linsolve!`), FFT, and Chebyshev — not just
element learning, so `SOLVER TIMING [...]` lines are stable run-to-run for every
solver. Because BenchmarkTools excludes compilation internally, even the first
run in a fresh session reports the steady-state minimum (no manual warm-up run
needed). Set `:lbenchmark_solve => false` to fall back to a single quick solve
(one `time_ns` shot) on very large problems where the repeated benchmark solves
would be too costly.

### Reference solution + difference in the VTU

The inference run also writes, into the same `iter_1.vtu`, extra nodal fields for
side-by-side visualisation in ParaView:

| Field                     | Meaning                                             |
|---------------------------|-----------------------------------------------------|
| `u_SEM_numerical`         | direct SEM (numerical) solution                     |
| `diff_SEM_minus_infer`    | SEM − inferred                                      |
| `u_exact_manufactured`    | exact/manufactured solution (only if the case has one) |
| `diff_exact_minus_infer`  | exact − inferred (only if the case has one)         |

### Timing hierarchy (printed by the pipeline)

When launched through `run_element_learning.sh`, a hierarchy of wall-clock and
solver timings is printed at the very end:

```
 # EL PIPELINE TIMING HIERARCHY
 #   (a) full run (phases that ran)      : 10m 32s
 #        ├─ sampling phase (SEM+startup) : 2m 00s
 #        ├─ (b) training only            : 8m 00s
 #        └─ inference phase (wall clock) : 0m 32s
 #   Solver-level (pure compute, JIT-excluded; this single solve):
 #        (c) direct SEM (factorize+solve)    : 0.0547 s
 #        (d) element learning (extract+infer): 0.116 s
 #        speedup (c)/(d)                     : 0.47×
```

(a)/(b) and the phase wall clocks come from the launcher; (c)/(d) are the pure
solver kernel times scraped from the inference run's diagnostics — the honest
like-for-like single-solve comparison (each includes its own `A`-dependent
setup: EL's per-element block extraction, SEM's LU factorization). The phase wall
clocks additionally include Julia startup, mesh read and IO. (The amortized
"if reusing `A`" numbers are in the `ELEMENT-LEARNING DIAGNOSTICS` block.)

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
