# The Jexpresso CI test suite

Everything CI runs is declared in **one file**: [`test/ci_cases.jl`](ci_cases.jl).
No workflow file ever has to change — the GitHub Actions job matrices are
generated from that registry. And you do not edit it by hand either: adding a
test to CI is one command,
`julia --project=. test/generate_ci_ref.jl <EQS>/<CASE>`
([details below](#add-a-new-test-to-ci)).

The suite currently contains:

| Case | Dim | Timeout | Tolerance | VTK check |
|---|---|---|---|---|
| `CompEuler/theta` | 2D | 40 min | `atol = 1e-5` | off |
| `CompEuler/sod1d` | 1D | 10 min | `atol = 1e-5` | n/a (1D has no VTK writer) |
| `AdvDiff/Wave_Train` | 1D | 10 min | `atol = 1e-5` | n/a (1D has no VTK writer) |

The other cases that used to be in the suite are listed, commented out, at the
bottom of `CI_CASES` in `test/ci_cases.jl`. Re-enabling one is uncommenting its
line and running `test/generate_ci_ref.jl` for it — that recreates the CI deck
from `problems/` if it is no longer under `test/CI-runs/`, and generates the
reference solution.

### A reference nobody runs against is an error

`test/runtests.jl` iterates `CI_CASES`. A case that is not in it is not run —
no failure, no warning, no line in the test summary — and committing its
reference solution does not change that. So the dangerous state is not a
missing reference, it is a **present reference with no registry entry**: the
deck, the golden files and this table all say the case is covered, and CI
reports "all passed" without ever having evaluated it. That is exactly what
happened to `CompEuler/sod1d`, disabled in place ("commented out sod1d for
now") while its reference stayed committed, so a deliberately broken 1D RHS
still produced a green suite.

`validate()` now checks that direction too, and `julia test/ci_cases.jl
validate` — the `registry` job that gates the whole CI workflow — fails on it:

```
ERROR: CompEuler/sod1d: test/CI-ref/CompEuler/sod1d/output/ holds a reference
solution, but the case is not in CI_CASES, so nothing ever compares against it
and the suite stays green however far the solution drifts. The registry has it
commented out — uncomment that line:
      CICase(eqs = "CompEuler", case = "sod1d", timeout = 10),
      If the reference is obsolete instead, delete test/CI-ref/CompEuler/sod1d/.
```

Disabling a case therefore means deciding what happens to its reference:
uncomment the line to keep the coverage, or delete `test/CI-ref/<EQS>/<CASE>/`
to say plainly that the case is not tested. `test/generate_ci_ref.jl`
re-enables a commented-out line in place rather than adding a second one for
the same case, so the registry never carries a live entry and a dead one
together.

The reverse case is only a note: a **registered case with no reference yet**
prints one from `validate`, and the comparison reports itself skipped.

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
| ✅ pass | the run finished and every field matches the reference within tolerance ([what that means](#what-the-tolerance-actually-tests)) |
| ⏭️ skipped | no reference solution committed yet — run `test/generate_ci_ref.jl` |
| ❌ fail | the run crashed, produced no HDF5 output, or a field drifted beyond tolerance |

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
| run the case and publish `test/CI-ref/<EQS>/<CASE>/output/` | `generate_ci_ref.jl` |
| add the `CICase(...)` line to `test/ci_cases.jl`, with a timeout from the measured run time | `generate_ci_ref.jl` |
| commit | **you** — the diff is the new definition of "correct", so look at it |

Nothing is printed during the solve (one write, at the end). To watch
progress or measure the step rate, run with `JEXPRESSO_STEP_HEARTBEAT=1`,
which prints the first 5 steps and then every 100th.

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

And on `test/runtests.jl`:

| Flag | Effect |
|---|---|
| `--vtk` | also run the [VTK writer check](#checking-the-vtk-writer) for the selected cases |

### Checking the VTK writer

The comparison reads HDF5, because that is the format that holds the solution
as `Float64` arrays with no geometry wrapped around it. Production runs,
though, write **VTK** — so nothing in CI would notice `write_vtk` breaking.
The VTK smoke test closes that gap. It runs the case a second time with VTK
output and checks the writer produced non-empty `.vtu`/`.pvtu` files. It
asserts nothing about the numbers: no reference, no tolerance.

Only for **2D and 3D** cases: `write_vtk` has methods for `NSD_2D` and
`NSD_3D` only, so a 1D case (which writes PNG/ASCII, plus HDF5 in CI) has no
VTK output to check and the test reports that rather than a `MethodError`.

**Run it for a case, one command:**

```bash
julia --project=. test/runtests.jl --vtk CompEuler/theta
```

That does the normal run-and-compare *and* the VTK check. Drop the case name
to do it for every registered case:

```bash
julia --project=. test/runtests.jl --vtk
```

**Turn it on permanently for a case** (it is off by default — it doubles that
case's CI wall time) by setting the flag in `test/ci_cases.jl`:

```julia
CICase(eqs = "CompEuler", case = "theta", timeout = 40, vtk_smoke = true),
```

with `timeout` raised to cover two runs. Output goes to the case's usual
`test/CI-runs/<EQS>/<CASE>/output/`, which is gitignored, and the Benchmarks
workflow uploads it as an artifact you can open in ParaView.

Under the hood the smoke run sets `JEXPRESSO_CI_OUTFORMAT=vtk`, which is the
one CI output convention that is *not* hardcoded; everything else (output
next to the case, no timestamped directory, a single write at `:tend`) still
applies. You can use it directly for a one-off:

```bash
JEXPRESSO_CI_OUTFORMAT=vtk julia --project=. \
    -e 'using Jexpresso; Jexpresso.run_case("CompEuler", "theta"; CI_MODE = true)'
```

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

## What the tolerance actually tests

Each floating-point field of each output file is compared with
`isapprox(ref, gen; atol, rtol)`, which for arrays means

    ‖gen − ref‖₂ ≤ max(atol, rtol · max(‖gen‖₂, ‖ref‖₂))

with the 2-norm taken over the **whole field** (every node in the file). Two
consequences:

* `atol` is not a pointwise bound. With N nodes, `atol = 1e-5` allows an RMS
  difference of about `1e-5/√N` per node — for the 8405 nodes of
  `CompEuler/3d`, about 1e-7. On a case run for thousands of steps that is
  close to a bit-for-bit test.
* The output holds **total** fields, not perturbations: ρθ ≈ 300·ρ. A
  relative tolerance sized for the momenta would let ρθ drift by kelvins.
  That is why `rtol` can be given per variable,
  `rtol = [ρ, ρu, ρv, ρw, ρθ, …]` (`rtol[i]` applies to `var_i_*.h5`).

The default is `rtol = 0` (absolute only). The `@error` line of a failing
field prints `‖Δ‖₂`, `‖Δ‖₂/‖ref‖₂` and `max |Δ|`, which is what you need to
tell round-off drift from a real change.

## Round-off-sensitive cases: `CompEuler/3d` and DSGS

`CompEuler/3d` (3D rising thermal, `:visc_model => DSGS()`, 2000 steps of
Δt = 0.5 s to t = 1000 s) does not reproduce across machines. The
`ubuntu-latest` pool mixes CPU types — seen in November 2026: Intel Xeon
6975P (Granite Rapids), AMD EPYC 7763 (Zen 3), AMD EPYC 9V74 (Zen 4). Each
type reproduces its own result bit for bit, but different types give
different momenta. The solutions are physically the same — same structures,
same amplitudes — but μ_dsgs differs slightly between them, and over 2000
steps that difference grows to ~1% of the momentum field.

So with `atol` alone, whether the case passes is decided by the CPU the job
lands on, not by the PR. The committed reference is, bit for bit, the Intel
output; the failing CI runs of #175–#177 (18 Nov 2026, none of which touches
DSGS) are, bit for bit, the AMD EPYC 7763 output. The reference had already
been regenerated once for the same reason (`dc8d92d9a`).

Measured differences at t = 1000 s, `‖Δ‖₂/‖ref‖₂` per field. The runner
outputs come from the Benchmarks workflow run on the same commit with the
deck unchanged (DSGS) or with only `:visc_model => SMAG()`:

| Pair compared | ρ | ρu | ρw | ρθ |
|---|---|---|---|---|
| **DSGS** Intel 6975P (= reference) vs AMD 7763 (= #175 CI) | 1.3e-5 | 1.1e-2 | 5.7e-3 | 7.7e-6 |
| **DSGS** Intel 6975P vs macOS arm64 | 1.5e-5 | 5.9e-3 | 4.1e-3 | 7.8e-6 |
| **DSGS** AMD 7763 vs macOS arm64 | 1.6e-5 | 1.1e-2 | 6.7e-3 | 6.1e-6 |
| **SMAG** AMD 7763 vs AMD 9V74 | 0 | 0 | 0 | 0 |
| **SMAG** AMD (either) vs macOS arm64 | 2.2e-14 | 1.2e-11 | 9.5e-12 | 3.0e-15 |
| **DSGS** arm64 vs arm64 with θc·(1 + 1e-12) | 7.6e-6 | 6.2e-3 | 2.8e-3 | 5.9e-6 |
| **SMAG** arm64 vs arm64 with θc·(1 + 1e-12) | 3.1e-14 | 1.4e-11 | 1.2e-11 | 3.1e-15 |
| SMAG vs DSGS (a real model change, for scale) | 4.2e-4 | 3.6e-1 | 2.2e-1 | 7.3e-5 |

(ρv is ~1e-11 everywhere — the case is a one-element-thick slab — and passes
on `atol`. No SMAG job landed on an Intel runner.)

Across platforms SMAG agrees to ~1e-11, DSGS to ~1e-2. The two perturbation
rows isolate the cause without involving hardware: changing the warm-bubble
amplitude θc by one part in 10¹² gives ~1e-11 in the SMAG momenta (no
amplification) and 6e-3 in the DSGS momenta — amplified by ~10¹⁰, to the
same level as the platform-to-platform differences. A better reference
cannot avoid this: any round-off-level change (another CPU, compiler, libm or
BLAS, a reordered sum in an unrelated part of the code) saturates at this
level under DSGS.

A plausible mechanism: the DSGS coefficient is built from the element
residual normalised by the spread of the solution and then clipped by
`Cmin`/`Cmax`. Where the flow is nearly at rest the residual itself is at
round-off level, so μ_dsgs there carries round-off noise into the dissipation
of the next step; the bubble's nonlinearity then amplifies it.

**What CI does about it.** `CompEuler/3d` is registered with
`rtol = [1e-4, 5e-2, 5e-2, 5e-2, 5e-5]` (ρ, ρu, ρv, ρw, ρθ): about 5× the
largest drift observed on each field. Checked against the outputs above, every
DSGS pair passes, and the SMAG-vs-DSGS change fails on all four of ρ, ρu, ρw
and ρθ. Note how small a real change is on the total ρ and ρθ (4e-4, 7e-5): a
loose relative tolerance on those fields would never fail.

**For a new case**: if it uses DSGS (or any closure whose coefficient depends
on a residual or on a max/min of the solution) and runs more than a few
hundred steps, expect this. Measure it before choosing a tolerance — run the
case twice locally with a 1e-12 change to the initial condition and compare
the outputs: copy the first run's `test/CI-runs/<EQS>/<CASE>/output/` to
`<DIR>/<EQS>/<CASE>/output/`, make the second run, then
`julia --project=. test/compare_benchmarks.jl --ref <DIR> <EQS>/<CASE>`. If the difference is far above round-off, set a
per-variable `rtol` a few times above it rather than regenerating the
reference whenever the runner changes.

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
