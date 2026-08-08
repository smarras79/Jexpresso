# The Jexpresso CI test suite

Everything CI runs is declared in **one file**: [`test/ci_cases.jl`](ci_cases.jl).
Adding a test to CI is adding one line to it; no workflow file has to change.

The suite currently contains a single case:

| Case | Timeout | Tolerance |
|---|---|---|
| `CompEuler/theta` | 40 min | `atol = 1e-5` |

The other cases that used to be in the suite are listed, commented out, at the
bottom of `CI_CASES` in `test/ci_cases.jl`, and their input decks are still in
`test/CI-runs/`. Re-enabling one is a matter of uncommenting its line and
regenerating its reference solution.

---

## How it fits together

```
test/ci_cases.jl          ← the registry: which cases run, timeout, tolerance
      │
      ├── test/runtests.jl           `Pkg.test()` → run + compare every case
      ├── test/compare_benchmarks.jl compare only (no solver load)
      ├── test/generate_ci_ref.jl    (re)generate the reference solutions
      ├── test/ci_compare.jl         shared run/compare helpers
      │
      ├── .github/workflows/CI.yml               merge gate (Pkg.test)
      ├── .github/workflows/benchmarks.yml       one job per case + report
      └── .github/workflows/generate-ci-ref.yml  (re)generate references
```

A case is run in **CI mode**, i.e. the solver reads its inputs from
`test/CI-runs/<EQS>/<CASE>/` (a reduced version of the deck in
`problems/<EQS>/<CASE>/`: short `:tend`, `hdf5` output) and writes to
`test/CI-runs/<EQS>/<CASE>/output/`. That output is compared field by field
against the committed reference in `test/CI-ref/<EQS>/<CASE>/output/`.

Outcome of a case:

| Result | Meaning |
|---|---|
| ✅ pass | the run finished and every field matches the reference within `atol` |
| ⏭️ skipped | no reference solution committed yet — run **Generate CI Reference Solutions** |
| ❌ fail | the run crashed, produced no HDF5 output, or a field drifted beyond `atol` |

### Which workflow runs when

| Workflow | Trigger | What it does |
|---|---|---|
| **CI** | push/PR to `master`, tags, manual | validates the registry, then `Pkg.test()` — this is the merge gate |
| **Benchmarks** | push to `master`, weekly (Sun 02:00 UTC), manual | one job per case, per-case artifacts and a summary table |
| **Generate CI Reference Solutions** | manual only | runs the cases and commits their output as the new references |

Benchmarks does not run on pull requests, because CI already runs the same
cases there. Uncomment the `pull_request:` trigger in
`.github/workflows/benchmarks.yml` if you want the per-case table on PRs too.

---

## Add a new test to CI

Assuming `problems/<EQS>/<CASE>/` already exists and runs (see
[ADD_A_NEW_TEST.md](../ADD_A_NEW_TEST.md)):

### 1. Create the CI version of the deck

```bash
mkdir -p test/CI-runs/<EQS>
cp -r problems/<EQS>/<CASE> test/CI-runs/<EQS>/<CASE>
```

Edit `test/CI-runs/<EQS>/<CASE>/user_inputs.jl` so the run is short and
writes HDF5:

| Key | CI value | Why |
|---|---|---|
| `:outformat` | `"hdf5"` | the comparison reads `.h5` files |
| `:output_dir` | `"none"` | output goes to `test/CI-runs/<EQS>/<CASE>/output/` |
| `:loverwrite_output` | `true` | no timestamped subdirectory |
| `:tend` | a few time steps' worth | keep CI wall-time short |
| `:ndiagnostics_outputs` / `:diagnostics_at_times` | consistent with `:tend` | at least one output file must be written |

All six `user_*.jl` / `initialize.jl` files must be present — the solver
includes them unconditionally.

### 2. Register the case — the one line

In `test/ci_cases.jl`:

```julia
const CI_CASES = CICase[
    CICase(eqs = "CompEuler", case = "theta", timeout = 40, atol = 1e-5),
    CICase(eqs = "<EQS>",     case = "<CASE>", timeout = 20),   # ← new
]
```

`timeout` is the per-case limit in minutes; `atol` defaults to `1e-5`.

Check the registry locally:

```bash
julia test/ci_cases.jl validate    # dirs and files present?
julia test/ci_cases.jl list        # what will run
julia test/ci_cases.jl matrix      # the JSON the workflows consume
```

### 3. Generate the reference solution

Until a reference exists the case is reported as *skipped*, not passing.
Create it either locally or on GitHub — both produce the same files.

**Locally** (recommended: you see the run, and you review the diff before it
becomes the golden answer):

```bash
julia --project=. test/generate_ci_ref.jl <EQS>/<CASE>

git status --short test/CI-ref     # look at what appeared
git add test/CI-ref
git commit -m "Add CI reference solution for <EQS>/<CASE>"
git push
```

That is all there is to it: the script runs the case with `CI_MODE=true` and
copies the resulting `.h5` files (plus the `user_inputs.jl` they came from)
into `test/CI-ref/<EQS>/<CASE>/output/`. It does not touch git — committing
and pushing is yours to do, deliberately, because a reference solution is the
definition of "correct" from then on.

**On GitHub**: commit and push first, then **Actions → Generate CI Reference
Solutions → Run workflow** on your branch (leave `cases` as `all`, or pass
just `<EQS>/<CASE>`). The workflow runs the same script and pushes one commit
with the references back to your branch, marked `[skip ci]`.

### 4. Verify

```bash
julia --project=. test/runtests.jl <EQS>/<CASE>   # run + compare, locally
```

or trigger **Actions → Benchmarks → Run workflow** and confirm the new row is
green.

---

## Regenerating references after an intentional change

If a physics or numerics change moves the expected answer, regenerate the
affected cases — locally:

```bash
julia --project=. test/generate_ci_ref.jl CompEuler/theta   # or several, or all
git add test/CI-ref && git commit -m "Update CI reference solutions" && git push
```

or with **Actions → Generate CI Reference Solutions** on the branch carrying
the change (pass a comma-separated list in `cases` to regenerate only some).

Either way the new `.h5` files become the golden reference. Review that diff
carefully — it is the record of what your change did to the solution.

### `test/generate_ci_ref.jl` in full

```
julia --project=. test/generate_ci_ref.jl [options] [<EQS>/<CASE> ...]

  (no case given)         every case registered in test/ci_cases.jl
  --copy-only, --no-run   skip the simulation and publish the output already
                          in test/CI-runs/<EQS>/<CASE>/output/ — use it after
                          having run a case by hand. Needs no packages, so
                          plain `julia` (no --project) is enough.
  --dest DIR              write somewhere other than test/CI-ref/ (this is how
                          the workflow stages references before committing
                          them from a single job)
  -h, --help              usage
```

Per case it runs the case, clears the stale `.h5` files from the destination
(a reference set must correspond to exactly one run — a leftover file from an
older run would be reported as missing output forever), then copies the fresh
ones in. It exits non-zero if a case crashed or produced no HDF5 output, and
prints a summary table of what was published.

The `generate-ci-ref` workflow calls this same script, so the local and the
CI paths cannot drift apart.

---

## Running things locally

```bash
# whole suite (run + compare), exactly what CI does
julia --project=. -e 'using Pkg; Pkg.test()'

# a single registered case
julia --project=. test/runtests.jl CompEuler/theta

# comparison only, against whatever is already in test/CI-runs/**/output
julia --project=. test/compare_benchmarks.jl CompEuler/theta

# run a case by hand
julia --project=. -e 'using Jexpresso; Jexpresso.run_case("CompEuler", "theta"; CI_MODE = true)'
```
