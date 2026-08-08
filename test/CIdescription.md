# The Jexpresso CI test suite

Everything CI runs is declared in **one file**: [`test/ci_cases.jl`](ci_cases.jl).
No workflow file ever has to change — the GitHub Actions job matrices are
generated from that registry. And you do not edit it by hand either: adding a
test to CI is one command,
`julia --project=. test/generate_ci_ref.jl <EQS>/<CASE>`
([details below](#add-a-new-test-to-ci)).

The suite currently contains a single case:

| Case | Timeout | Tolerance |
|---|---|---|
| `CompEuler/theta` | 40 min | `atol = 1e-5` |

The other cases that used to be in the suite are listed, commented out, at the
bottom of `CI_CASES` in `test/ci_cases.jl`. Re-enabling one is uncommenting its
line and running `test/generate_ci_ref.jl` for it — that recreates the CI deck
from `problems/` if it is no longer under `test/CI-runs/`, and generates the
reference solution.

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

A CI case is self-contained inside this repository — deck, mesh and reference
solution are all committed, because a runner has nothing else.

A case is run in **CI mode**, i.e. the solver reads its inputs from
`test/CI-runs/<EQS>/<CASE>/` (a reduced version of the deck in
`problems/<EQS>/<CASE>/`: short `:tend`, `hdf5` output) and writes to
`test/CI-runs/<EQS>/<CASE>/output/`. That output is compared field by field
against the committed reference in `test/CI-ref/<EQS>/<CASE>/output/`.

Outcome of a case:

| Result | Meaning |
|---|---|
| ✅ pass | the run finished and every field matches the reference within `atol` |
| ⏭️ skipped | no reference solution committed yet — run `test/generate_ci_ref.jl` |
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

One command. Assuming `problems/<EQS>/<CASE>/` exists and runs:

```bash
julia --project=. test/generate_ci_ref.jl <EQS>/<CASE>

git add test/CI-ref test/CI-runs test/ci_cases.jl
git commit -m "Add <EQS>/<CASE> to CI"
git push
```

That script does every step of the old checklist for you:

| Step | Done by |
|---|---|
| copy `problems/<EQS>/<CASE>/` → `test/CI-runs/<EQS>/<CASE>/` (skipping any `output/`, symlinking any `.msh`) | `generate_ci_ref.jl` |
| write HDF5, next to the case inputs, no timestamped directory | CI mode (`src/run.jl`) |
| write **one** output, the final state, whatever cadence the deck asked for | CI mode (`src/run.jl`) |
| print a step heartbeat — first 5 steps, then every 100th — so a long run shows progress instead of looking hung | CI mode (`src/run.jl`) |
| run the case and publish `test/CI-ref/<EQS>/<CASE>/output/` | `generate_ci_ref.jl` |
| add the `CICase(...)` line to `test/ci_cases.jl`, with a timeout from the measured run time | `generate_ci_ref.jl` |
| commit | **you** — the diff is the new definition of "correct", so look at it |

The one judgement call left to you is **how long the case should run**. The
copied deck keeps whatever `:tend` `problems/` had; if that is minutes of
wall time, shorten `:tend` in
`test/CI-runs/<EQS>/<CASE>/user_inputs.jl` and run the script again. Nothing
else in the deck needs editing.

Useful flags:

| Flag | Effect |
|---|---|
| `--refresh-deck` | re-copy the deck from `problems/`, discarding edits to the CI copy |
| `--copy-only` | publish output from a run you already did by hand (needs no packages) |
| `--no-register` | leave `test/ci_cases.jl` alone |
| `--dest DIR` | write the references elsewhere, e.g. to diff before overwriting |

### If the case reads a GMSH mesh

`meshes/` is each developer's link to
[`smarras79/JexpressoMeshes`](https://github.com/smarras79/JexpressoMeshes).
A CI runner clones Jexpresso and nothing else, so that directory does not
exist there and a case pointing into it fails with `Msh file not found`
about fifteen minutes into the job.

A case that CI runs keeps its mesh **in this repository, next to its deck**:

```
problems/<EQS>/<CASE>/
├── user_inputs.jl
├── initialize.jl
├── …
└── <mesh>.msh          ← committed here
```

```julia
:lread_gmsh    => true,
:gmsh_filename => "./problems/<EQS>/<CASE>/<mesh>.msh",
```

`generate_ci_ref.jl` **symlinks** that mesh into the CI deck and retargets
`:gmsh_filename` at the link, so the CI case has the mesh at its own path
without a second copy of the bytes in git:

```
test/CI-runs/<EQS>/<CASE>/
├── user_inputs.jl      :gmsh_filename => "./test/CI-runs/<EQS>/<CASE>/<mesh>.msh"
└── <mesh>.msh          → ../../../../problems/<EQS>/<CASE>/<mesh>.msh
```

The link is relative, so it resolves in every clone and in the archive a
runner unpacks; git stores it as a link (mode 120000), the same trick
`test/meshes` uses. Paths inside a deck are resolved from the **repository
root**, never from the deck itself.

Because it is a link and not a copy, refining the mesh under `problems/`
changes what CI runs — the comparison against the old reference will fail,
which is the point: regenerate the reference deliberately rather than have
the two drift apart silently.

Prefer the smallest mesh that still exercises the case: it is committed
history and it is CI wall-time.

`julia test/ci_cases.jl validate` checks the mesh exists, so a mesh left
behind in `meshes/` fails the registry job in seconds instead of a quarter of
an hour into the run.

### Verify

```bash
julia --project=. test/runtests.jl <EQS>/<CASE>   # run + compare, locally
julia test/ci_cases.jl list                       # what CI will run
```

or trigger **Actions → Benchmarks → Run workflow** and confirm the new row is
green.

### Doing it on GitHub instead

Push the deck first, then **Actions → Generate CI Reference Solutions → Run
workflow** on your branch. It runs the same script and pushes one commit with
the references. Note that the workflow only knows about *registered* cases —
a brand new case is easier to add locally, where the script registers it for
you.

---

## Checking a reference across platforms

A reference generated on your machine is compared, in CI, against a run on
`ubuntu-latest` x64. Different libm, different BLAS, different vectorisation
— after a few thousand nonlinear time steps the two can disagree by more than
`atol = 1e-5`, and a reference generated on macOS (especially arm64) is the
usual suspect when CI fails on tolerance while everything passes locally.

Rather than discovering that from a 40-minute red CI run, have the runner
produce its own copy and compare the two **before** committing yours.

### 1. Get the runner's version

Push your branch (the case must be registered in `test/ci_cases.jl`; if you
generated locally, `test/generate_ci_ref.jl` already added it), then:

**Actions → Generate CI Reference Solutions → Run workflow**, branch = yours,
`cases` = `<EQS>/<CASE>`.

### 2. Download it without adopting it

Each case job uploads its result as an artifact **before** the commit job
touches the repository. On the run page, under **Artifacts**, download
`ciref-<EQS>-<CASE>` and unpack it:

```bash
unzip ~/Downloads/ciref-CompEuler-theta.zip -d /tmp/linux-ref
ls /tmp/linux-ref/CompEuler/theta/output/     # t.h5, var_1_0.h5, …
```

### 3. Diff it against your run

Your own output is still in `test/CI-runs/<EQS>/<CASE>/output/` from when you
generated the reference. Compare the Linux reference against it, at the
tolerance CI uses:

```bash
julia --project=. test/compare_benchmarks.jl --ref /tmp/linux-ref CompEuler/theta
```

| Result | Meaning | What to do |
|---|---|---|
| all tests pass | the platforms agree within `atol` | keep your reference, commit it, open the PR |
| `field 'q' differs … max │Δ│ = …` | genuine cross-platform drift | use the Linux reference instead (below) |
| `… is missing from the run output` | the two runs produced different *files*, not different numbers — usually one side is stale | regenerate locally and re-check |

A byte-level look, when you want one:

```bash
h5diff -v /tmp/linux-ref/CompEuler/theta/output/var_1_0.h5 \
          test/CI-runs/CompEuler/theta/output/var_1_0.h5
```

### 4. If they disagree, adopt the Linux reference

The workflow's commit job already pushed it to your branch, so:

```bash
git pull                                   # brings in the runner's references
julia --project=. test/runtests.jl CompEuler/theta   # now expected to FAIL locally
```

That local failure is normal and is the point: the reference now matches the
machine that gates the merge, not yours. If you would rather have a single
reference both platforms accept, widen the tolerance for that case instead —
`CICase(eqs = …, case = …, timeout = …, atol = 1e-4)` — and regenerate.

A shortcut worth knowing: the commit job only commits when the files actually
changed. If the workflow run reports *"No reference files changed — nothing to
commit"*, the runner reproduced your files byte for byte and there is nothing
left to check.

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
  --refresh-deck          re-copy the deck from problems/ even when the CI
                          copy exists (discards edits to that copy)
  --no-register           leave test/ci_cases.jl alone
  --dest DIR              write somewhere other than test/CI-ref/ (this is how
                          the workflow stages references before committing
                          them from a single job)
  -h, --help              usage
```

Per case it copies the deck from `problems/` if `test/CI-runs/` has none, runs
the case, clears the stale `.h5` files from the destination (a reference set
must correspond to exactly one run — a leftover file from an older run would
be reported as missing output forever), copies the fresh ones in, and
registers the case if it is new. It exits non-zero if a case crashed or
produced no HDF5 output, and prints a summary table of what was published.

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
