# CI benchmark testing

This document explains how the Jexpresso CI pipeline works, how to run it, and
how to add a case to it.

The short version: **every case CI runs is one line in
[`test/ci_cases.jl`](../test/ci_cases.jl)**. The GitHub Actions job matrices
are generated from that file at run time, so adding or removing a case never
touches a workflow file.

---

## Overview

| Workflow | Trigger | Purpose |
|---|---|---|
| `CI.yml` | push/PR to `master`, tags, manual | Validate the registry, then `Pkg.test()` — run and compare every registered case. This is the merge gate. |
| `benchmarks.yml` | push to `master`, weekly (Sun 02:00 UTC), manual | One job per case: own log, own timeout, own output artifact, plus an aggregated summary table |
| `generate-ci-ref.yml` | manual only | Run the cases and commit their output as the reference (“golden”) solutions |
| `Documentation.yml` | push/PR to `master`, tags, manual | Build and deploy the documentation |
| `cleanup.yml` | PR closed | Delete the documentation preview of that PR |

`benchmarks.yml` deliberately does not run on pull requests — `CI.yml` already
runs the same cases there, and these are real simulations. The
`pull_request:` trigger is present but commented out in the workflow if you
want the per-case table on PRs as well.

### What a case run means

Each case runs in **CI mode**: the solver reads
`test/CI-runs/<EQS>/<CASE>/user_inputs.jl` (a shortened deck, `hdf5` output)
instead of `problems/<EQS>/<CASE>/`, and writes to
`test/CI-runs/<EQS>/<CASE>/output/`.

Pass/fail has two layers:

1. **Execution** — did the run finish without error?
2. **Numerical correctness** — does every HDF5 field match the committed
   reference in `test/CI-ref/<EQS>/<CASE>/output/` within the case's `atol`
   (default `1e-5`)?

A case with no committed reference is reported as **skipped**, not failed.

---

## The registry

```julia
# test/ci_cases.jl
const CI_CASES = CICase[
    CICase(eqs = "CompEuler", case = "theta", timeout = 40, atol = 1e-5),
]
```

| Field | Meaning |
|---|---|
| `eqs`, `case` | directory names under `test/CI-runs/` |
| `timeout` | per-case limit in minutes, used for the GitHub Actions job |
| `atol` | absolute tolerance of the field-by-field comparison |

The registry has a small command line interface, used by the workflows and
useful locally:

```bash
julia test/ci_cases.jl list        # cases that will run
julia test/ci_cases.jl validate    # registry vs. what is actually on disk
julia test/ci_cases.jl matrix      # the JSON the workflows turn into jobs
julia test/ci_cases.jl matrix "CompEuler/theta"   # a subset
```

`validate` is the first thing both CI and the benchmarks workflow do, so a
typo in the registry fails in seconds instead of after a 40-minute run.

### Current suite

Only **`CompEuler/theta`** is enabled. The cases that were in the suite
previously are listed, commented out, at the bottom of `CI_CASES`; their input
decks are still in `test/CI-runs/` and several still have reference solutions
in `test/CI-ref/`. Re-enabling one is: uncomment the line, run
**Generate CI Reference Solutions** for it, done.

---

## Running the pipeline

### Locally

```bash
# everything CI does
julia --project=. -e 'using Pkg; Pkg.test()'

# one case
julia --project=. test/runtests.jl CompEuler/theta

# comparison only, against output already in test/CI-runs/**/output
julia --project=. test/compare_benchmarks.jl CompEuler/theta
```

### On GitHub

**Actions → Benchmarks → Run workflow**, pick a branch, and optionally narrow
the run with the `cases` input (`all`, or e.g.
`CompEuler/theta,Burgers/case1`). The run page shows one job per case and a
summary table:

```
## Benchmark results

| Case                | Simulation | Comparison vs reference |
|---------------------|------------|-------------------------|
| `CompEuler/theta`   | ✅ pass    | ✅ pass                 |
```

Each job also uploads its HDF5 output as an artifact (`output-<EQS>-<CASE>`),
which is what you want when a comparison fails and you need to look at the
fields.

---

## Reference solutions

References live in `test/CI-ref/<EQS>/<CASE>/output/` and are committed to the
repository.

To create or update them:

1. **Actions → Generate CI Reference Solutions → Run workflow**.
2. Choose the branch holding the intended baseline code.
3. `cases`: `all`, or a comma-separated subset.
4. Optionally set a commit message.

The workflow runs one job per case, stages each case's `.h5` files (plus the
`user_inputs.jl` they were produced with), and a single final job merges them
into `test/CI-ref/` and pushes **one** commit to the branch, marked
`[skip ci]`.

Re-run it whenever a physics or numerics change intentionally moves the
expected answer — and review the resulting diff, since it is the record of
what changed.

---

## Adding a case

Full checklist in [`test/CIdescription.md`](../test/CIdescription.md). In brief:

1. Copy the deck: `cp -r problems/<EQS>/<CASE> test/CI-runs/<EQS>/<CASE>`, then
   shorten it — `:outformat => "hdf5"`, `:output_dir => "none"`,
   `:loverwrite_output => true`, a small `:tend`, and diagnostics settings that
   still write at least one output file.
2. Add one line to `CI_CASES` in `test/ci_cases.jl`.
3. Push, then run **Generate CI Reference Solutions** for the new case.
4. Confirm the case is green in **Benchmarks** (and no longer skipped).

Nothing else — no workflow edit, no list to keep in sync in the comparison
script.

---

## Layout

```
test/ci_cases.jl              registry + matrix generator + validator
test/ci_compare.jl            run/compare helpers shared by everything below
test/runtests.jl              Pkg.test() entry point: run + compare each case
test/compare_benchmarks.jl    comparison only (HDF5 only, no solver load)
test/CI-runs/<EQS>/<CASE>/    shortened input decks; output/ is written here
test/CI-ref/<EQS>/<CASE>/     committed reference solutions
```
