# CI Benchmark Testing

This document explains how the classical benchmark CI pipeline works, how to run
it, and how to add a new benchmark to the suite.

---

## Overview

Two GitHub Actions workflows manage benchmark testing:

| Workflow | Trigger | Purpose |
|---|---|---|
| `benchmarks.yml` | Push/PR to `master`, weekly, manual | Run all benchmarks and compare output against reference |
| `generate-ci-ref.yml` | Manual only | Generate and commit reference `.h5` files |

Each benchmark is run in **CI mode**: it reads problem inputs from
`test/CI-runs/<EQ>/<CASE>/user_inputs.jl` (reduced simulation time,
`hdf5` output) instead of the full `problems/` directory.  Output lands in
`test/CI-runs/<EQ>/<CASE>/output/`.

Pass/fail is determined in two layers:

1. **Execution** — did the Julia process exit without error?
2. **Numerical correctness** — do all HDF5 fields in the output match the
   committed reference files in `test/CI-ref/<EQ>/<CASE>/output/` within
   `atol = 1e-5`?

A markdown summary table is written to the GitHub Actions job summary page
after every run.

---

## Running the benchmark CI

### Automatic runs

`benchmarks.yml` triggers automatically on every push or pull request to
`master`, and on a weekly schedule (Sunday 02:00 UTC).  No manual action is
required.

### Manual run

1. Go to **Actions → Classical Benchmarks → Run workflow**.
2. Select the branch you want to test.
3. Click **Run workflow**.

The job summary (visible on the workflow run page under the **Summary** tab)
shows a table like:

```
### Simulation (did the code run?)
| Benchmark          | Result  |
|--------------------|---------|
| CompEuler/theta    | ✅ Pass |
| CompEuler/3d       | ✅ Pass |
| ...                | ...     |

### Numerical comparison (do results match reference?)
| Status                                                    |
|-----------------------------------------------------------|
| ✅ Pass (see Compare step log for per-benchmark detail)  |
```

If a benchmark shows **skipped** in the comparison column it means no
reference file exists for it yet — see [Bootstrapping reference files](#bootstrapping-reference-files) below.

---

## Bootstrapping reference files

Reference files must exist in `test/CI-ref/` before numerical comparison
can run.  Generate them with the dedicated workflow:

1. Go to **Actions → Generate CI Reference Solutions → Run workflow**.
2. Select the branch that holds your intended baseline code.
3. Optionally set a commit message (default: `Update CI reference solutions`).
4. Click **Run workflow**.

The workflow:
- Runs every benchmark with CI mode inputs.
- Copies the resulting `.h5` files from `test/CI-runs/<EQ>/<CASE>/output/`
  to `test/CI-ref/<EQ>/<CASE>/output/`.
- Commits and pushes those files with `[skip ci]` in the message so it
  does not trigger another benchmark run.

After this, every subsequent `benchmarks.yml` run compares fresh output
against those committed files.

### Re-generating references after an intentional change

If a physics fix, scheme change, or parameter update intentionally alters the
expected solution, run **Generate CI Reference Solutions** again on the
updated branch.  The new `.h5` files replace the old ones and become the new
golden reference.

---

## Adding a new benchmark

Follow these steps to add `<MyEquations>/<mycase>` to the CI suite.

### 1 — Create the problem files

The full problem definition lives under `problems/`:

```
problems/<MyEquations>/<mycase>/
    initialize.jl
    user_bc.jl
    user_flux.jl
    user_inputs.jl
    user_primitives.jl
    user_source.jl
```

### 2 — Create the CI-runs directory

Copy the problem directory into `test/CI-runs/` and modify
`user_inputs.jl` for fast CI execution:

```bash
cp -r problems/<MyEquations>/<mycase> test/CI-runs/<MyEquations>/<mycase>
```

Then edit `test/CI-runs/<MyEquations>/<mycase>/user_inputs.jl` and change:

| Key | CI value | Reason |
|---|---|---|
| `:outformat` | `"hdf5"` | Required for comparison |
| `:output_dir` | `"none"` | Output goes to `test/CI-runs/<EQ>/<CASE>/output/` |
| `:tend` | A small value (a few time steps) | Keep CI wall-time short |
| `:diagnostics_at_times` | Match the new `:tend` | Ensure at least one output file is written |
| `:loverwrite_output` | `true` | Avoid timestamped subdirectories |

Example — reduce a case with `Δt = 0.01` to just 5 steps:

```julia
:tend                 => 0.05,
:diagnostics_at_times => (0.05,),
:outformat            => "hdf5",
:output_dir           => "none",
:loverwrite_output    => true,
```

### 3 — Register the benchmark in the comparison script

Open `test/compare_benchmarks.jl` and add the new case to the `BENCHMARKS`
constant:

```julia
const BENCHMARKS = [
    ...
    ("MyEquations", "mycase"),   # ← add this line
]
```

### 4 — Add simulation steps to both workflows

**`benchmarks.yml`** — add a simulation step and include the new step id
in the report and fail-gate:

```yaml
- name: "Simulate: MyEquations/mycase"
  id: sim_myequations_mycase
  continue-on-error: true
  timeout-minutes: 20
  run: |
    julia --project=. -e '
      push!(empty!(ARGS), "MyEquations", "mycase", "true")
      include("src/Jexpresso.jl")
    '
```

Then add the step id to the report table:

```yaml
echo "| MyEquations/mycase | $(icon '${{ steps.sim_myequations_mycase.outcome }}') |"
```

And to the fail-gate outcomes array:

```yaml
"${{ steps.sim_myequations_mycase.outcome }}"
```

**`generate-ci-ref.yml`** — add an identical block (replacing `sim_` with
`gen_`) plus a `copy_ref MyEquations mycase` line in the
"Collect reference files" step, and add the step to the report table.

### 5 — Generate the first reference solution

Commit everything, push, then run **Generate CI Reference Solutions** on
your branch (see [Bootstrapping reference files](#bootstrapping-reference-files)).

### 6 — Verify

Trigger **Classical Benchmarks** manually on the same branch and confirm:
- The new simulation step shows ✅.
- The comparison step shows ✅ (not skipped).
