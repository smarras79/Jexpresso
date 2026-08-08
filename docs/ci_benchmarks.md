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
`test/CI-runs/<EQS>/<CASE>/user_inputs.jl` (a shortened deck) instead of
`problems/<EQS>/<CASE>/`, and writes to `test/CI-runs/<EQS>/<CASE>/output/`.

CI mode also forces the three output settings the comparison depends on —
`:outformat => "hdf5"`, `:output_dir => "none"`, `:loverwrite_output => true`
— announcing each override it actually applies. A deck copied out of
`problems/` therefore produces comparable output even though it asks for VTK
somewhere else. `JEXPRESSO_CI_OUTPUT=0` disables the overrides.

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

### Locally

```bash
# run the case(s) and copy the output into test/CI-ref/
julia --project=. test/generate_ci_ref.jl CompEuler/theta   # or no argument: all cases

# review, then commit — the script never touches git
git status --short test/CI-ref
git add test/CI-ref
git commit -m "Update CI reference solutions"
git push
```

Useful variants:

| Command | Effect |
|---|---|
| `julia --project=. test/generate_ci_ref.jl` | regenerate every registered case |
| `julia test/generate_ci_ref.jl --copy-only CompEuler/theta` | publish the output of a run you already did by hand (no packages needed) |
| `julia --project=. test/generate_ci_ref.jl --dest /tmp/refs CompEuler/theta` | write elsewhere, e.g. to diff against the committed set before overwriting it |

The script clears stale `.h5` files from the destination before copying, so a
reference set always corresponds to exactly one run, and exits non-zero if a
case crashed or produced no HDF5 output.

### With the workflow

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

### Is my local reference good enough for Linux CI?

References are compared, in CI, against a run on `ubuntu-latest` x64. A
reference generated on macOS can drift past `atol = 1e-5` over a few thousand
time steps. To check before spending a CI run: dispatch **Generate CI
Reference Solutions** on your branch, download the `ciref-<EQS>-<CASE>`
artifact from the run page (it is uploaded before anything is committed), and
compare it against your local output at the CI tolerance:

```bash
unzip ~/Downloads/ciref-CompEuler-theta.zip -d /tmp/linux-ref
julia --project=. test/compare_benchmarks.jl --ref /tmp/linux-ref CompEuler/theta
```

Pass: keep your reference. Fail: `git pull` — the workflow already committed
the runner's version to your branch — and use that one. Step-by-step in
[`test/CIdescription.md`](../test/CIdescription.md#checking-a-reference-across-platforms).

---

## Adding a case

One command, from a case that already exists and runs under `problems/`:

```bash
julia --project=. test/generate_ci_ref.jl <EQS>/<CASE>

git add test/CI-ref test/CI-runs test/ci_cases.jl
git commit -m "Add <EQS>/<CASE> to CI"
git push
```

It copies `problems/<EQS>/<CASE>/` into `test/CI-runs/` (leaving out any
`output/`), runs it in CI mode — which forces HDF5 output, next to the case
inputs, one write at `:tend` — publishes the result as the reference, and
appends the `CICase(...)` line to `test/ci_cases.jl` with a timeout derived
from the measured run time. From the next push, CI runs the case.

What it cannot decide for you: how long the case should run. The copied deck
keeps the `:tend` from `problems/`; shorten it in
`test/CI-runs/<EQS>/<CASE>/user_inputs.jl` if the case is slow, then run the
script again.

### Meshes

`meshes/` is the developer's link to
[`smarras79/JexpressoMeshes`](https://github.com/smarras79/JexpressoMeshes)
and does not exist on a CI runner. A case that CI runs must keep its mesh
committed **in this repository, next to its deck**:

```julia
:gmsh_filename => "./problems/<EQS>/<CASE>/<mesh>.msh",
```

`generate_ci_ref.jl` symlinks that mesh into `test/CI-runs/<EQS>/<CASE>/`
(relative link, stored by git as mode 120000 — no second copy of the bytes)
and retargets `:gmsh_filename` at the link. Paths inside a deck are resolved
from the repository root, never from the deck. `julia test/ci_cases.jl
validate` fails in seconds when the mesh is missing, instead of letting the
job burn 15 minutes and die in `sem_setup`.

Full checklist and flags in [`test/CIdescription.md`](../test/CIdescription.md).

---

## Layout

```
test/ci_cases.jl              registry + matrix generator + validator
test/ci_compare.jl            run/compare helpers shared by everything below
test/runtests.jl              Pkg.test() entry point: run + compare each case
test/compare_benchmarks.jl    comparison only (HDF5 only, no solver load)
test/generate_ci_ref.jl       (re)generate the reference solutions
test/CI-runs/<EQS>/<CASE>/    shortened input decks; output/ is written here
test/CI-ref/<EQS>/<CASE>/     committed reference solutions
```
