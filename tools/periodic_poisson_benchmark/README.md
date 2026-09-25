# Periodic Poisson solver benchmark — running it at high resolution

This directory benchmarks six solvers of the doubly periodic Poisson problem
`Elliptic/poisson_periodic_sem` (−∇²u = f on [0,2π]², exact solution known):

| key | solver |
|---|---|
| `sem` | SEM, full system, sparse direct (Cholesky) |
| `sem_amg` | SEM, full system, AMG-preconditioned CG |
| `sc_direct` | SEM, static condensation (element learning), skeleton system by sparse direct |
| `sc_amg` | same condensation, skeleton system by AMG-preconditioned CG |
| `ps` | pseudo-spectral (Fourier collocation, Kopriva) |
| `fft` | FFT (FFTW) |

The files:

| file | what it does |
|---|---|
| `pipeline.jl` | runs the sweep and writes `results.csv`, `results.md` and the figures |
| `plot.py` | draws the SVG figures from `results.csv` |
| `latex.py` | writes the LaTeX paper section from `results.csv` |

## Resolution: h-refinement with `:linitial_refine`

The deck's mesh has 16×16 elements. `levels` refines it uniformly
`L` times through Jexpresso's `:linitial_refine` / `:init_refine_lvl` (p4est),
and periodicity is kept. Level `L` therefore has (16·2^L)² elements and
(16·2^L·N)² unknowns at order N. The pseudo-spectral and FFT grids follow the
same count, 16·2^L·N points per direction, so all six solvers are always
compared at the same number of unknowns.

| level L | elements | unknowns at N = 4 | unknowns at N = 8 |
|---:|---:|---:|---:|
| 0 | 16² | 4 096 | 16 384 |
| 1 | 32² | 16 384 | 65 536 |
| 2 | 64² | 65 536 | 262 144 |
| 3 | 128² | 262 144 | 1 048 576 |
| 4 | 256² | 1 048 576 | 4 194 304 |
| 5 | 512² | 4 194 304 | 16 777 216 |

Refining the mesh at fixed N is the first experiment of the paper section's
recommendations: it shows whether, and where, AMG overtakes the sparse
direct solve.

## Running it on a laptop

One-time setup, from the repository root with Julia 1.11.9:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
```

Python 3 is needed only for the figures and needs no packages.

### From the REPL (the intended way)

```julia
julia --project=.
julia> using Jexpresso
julia> include("tools/periodic_poisson_benchmark/pipeline.jl")
julia> rows = run_periodic_poisson_benchmark(; levels = 0:3, nops = 2:8,
                                             outdir = "ppb_highres", warmup = :small,
                                             max_unknowns = Dict(:all => 810_000, :sem_amg => 300_000),
                                             resume = true)
```

### As a script

This is still one Julia session and follows the same protocol:

```bash
julia --project=. tools/periodic_poisson_benchmark/pipeline.jl \
      --levels 0:3 --nops 2:8 --outdir ppb_highres --warmup small \
      --max-unknowns all=810000,sem_amg=300000 --resume
```

The pipeline changes into the repository root by itself, because the deck
names its mesh relative to it. A relative `outdir` is taken from where you
started.

### Options

| keyword (REPL) | flag (script) | meaning |
|---|---|---|
| `levels = 0:3` | `--levels 0:3` | refinement levels (default `0:0`, the 16×16 mesh) |
| `nops = 2:8` | `--nops 2:8` or `--nops 4,8` | SEM orders N |
| `solvers = (:sem, :sc_amg)` | `--solvers sem,sc_amg` | subset of the six solvers |
| `outdir = "dir"` | `--outdir dir` | where results and figures go |
| `max_unknowns = Dict(:sem_amg => 300_000, :all => 4_200_000)` | `--max-unknowns sem_amg=300000,all=4200000` | skip configurations above this many unknowns, per solver; `all` caps every solver |
| `resume = true` | `--resume` | keep what is already in `outdir/results.csv` and run only the missing configurations |
| `warmup = :small` | `--warmup small` | cheaper warm-up (see the timing protocol below) |
| `amg_itmax = 10_000` | `--amg-itmax 10000` | CG iteration cap of the AMG solves |
| `plot = false` | `--noplot` | skip the figures |

`results.csv` is rewritten after every configuration. If a long sweep is
interrupted, or a configuration runs out of memory, run the same command
again with `resume` and only the missing configurations run. A configuration
that throws an error is reported and skipped, and the sweep continues.

### Timing protocol

The protocol is unchanged: every configuration runs twice in the same Julia
session and only the second run is recorded. The garbage collector runs before
that second run, the mesh and SEM caches are off, and output files are off.

With `warmup = :same` (the default), the discarded first run is the identical
configuration. At level 3 or 4 that doubles the cost of every configuration.
`warmup = :small` makes the first run use the same solver and N at level
min(L, 1). Every code path is then compiled, including the p4est refinement,
but the large problem is solved only once. The recorded run is the second
run either way.

## Resources

These were measured in the development container (Intel Xeon 2.1 GHz, 4 cores,
one Julia thread), with each solver in its own Julia process. The peak RSS
therefore covers that solver alone, and includes about 2.5 GiB that the Julia
process holds (packages and compiled code) at any resolution.

| level, N | unknowns | solver | peak RSS | time-to-solution | of which the solve step |
|---|---:|---|---:|---:|---:|
| L2, N = 8 | 262 144 | SEM direct | 5.4 GiB | 31 s | 0.11 s |
| L2, N = 8 | 262 144 | SC direct | 4.9 GiB | 31 s | 0.35 s |
| L2, N = 8 | 262 144 | SC AMG (55 CG its) | 4.9 GiB | 34 s | 3.4 s |
| L2, N = 8 | 262 144 | SEM AMG (149 CG its) | 5.4 GiB | 92 s | 62 s |
| L2, N = 8 | 262 144 | pseudo-spectral | 2.5 GiB | 0.47 s | 0.013 s |
| L2, N = 8 | 262 144 | FFT | 2.5 GiB | 0.60 s | 0.006 s |
| L3, N = 6 | 589 824 | SEM direct | 6.5 GiB | 47 s | 0.20 s |
| L3, N = 8 | 1 048 576 | pseudo-spectral | 2.5 GiB | 2.1 s | 0.10 s |

For the SEM solvers, memory is dominated by the element matrices and the
assembly. That part grows like elements × (N+1)⁴, which is 4× per level at
fixed N. From the rows above:

| level, N | unknowns | SEM, estimated peak RSS |
|---|---:|---:|
| L3, N = 7 | 802 816 | ≈ 9 GiB |
| L3, N = 8 | 1 048 576 | ≈ 14 GiB |
| L4, N = 4 | 1 048 576 | ≈ 7 GiB |
| L4, N = 6 | 2 359 296 | ≈ 18 GiB |
| L4, N = 8 | 4 194 304 | ≈ 45 GiB |

The pseudo-spectral and FFT solvers read no mesh. The driver dispatches to
them before `sem_setup`, so they build no SEM infrastructure at any
resolution, and their memory is their own grid and operators.

Full-system AMG (`sem_amg`) is the slow one. Its CG iteration count grows with
N and under refinement (149 iterations at L2, N = 8), so cap it first.
Suggested sweeps, with `warmup = :small` so no configuration is solved twice:

| laptop RAM | suggested sweep |
|---|---|
| 16 GB | `levels = 0:3, max_unknowns = Dict(:all => 810_000, :sem_amg => 300_000)`, which reaches L3 up to N = 7 |
| 32 GB | `levels = 0:4, max_unknowns = Dict(:all => 1_100_000, :sem_amg => 600_000)`, which reaches L3 N = 8 and L4 N ≤ 4 |
| 64 GB | `levels = 0:4, max_unknowns = Dict(:all => 2_400_000, :sem_amg => 1_100_000)`, which reaches L4 N ≤ 6 |

For example, on a 16 GB laptop:

```julia
julia> rows = run_periodic_poisson_benchmark(; levels = 0:3, outdir = "ppb_highres", warmup = :small,
                                             max_unknowns = Dict(:all => 810_000, :sem_amg => 300_000),
                                             resume = true)
```

Start with a level you know fits, check `results.md`, then extend the sweep
with `resume = true`. If a configuration runs out of memory, it is reported
and skipped. If the operating system kills Julia instead, run again with
`resume = true` and a lower `max_unknowns`.

## Output

The figures are SVGs, each in a light and a dark variant (`-dark`), in
`<outdir>/assets`:

- **Per level:** `ppb_error_vs_{dofs,order,solve_time,total_time}_L<L>`.
- **Across levels:** at every N run on two or more levels, the h-refinement
  figures `ppb_hrefine_{solve,cost,total}_N<N>`. They plot wall-clock against
  unknowns with one curve per solver, which is where a crossover of AMG and
  the direct solves shows up. `cost` is the solver's own setup plus solve,
  without the SEM infrastructure.

With the default `outdir` (this directory), a single-level run writes the
README figures in the repository's `assets/` instead.

To show a run's figures and table in the main README.md (each run in its own
section, figures copied to `assets/<outdir name>/`; rerunning replaces that
section):

```bash
python3 tools/periodic_poisson_benchmark/readme_figures.py ppb_1level
python3 tools/periodic_poisson_benchmark/readme_figures.py ppb_highres
git add assets/ppb_1level assets/ppb_highres README.md
```

The LaTeX section for the paper covers one mesh at a time:

```bash
python3 tools/periodic_poisson_benchmark/latex.py ppb_highres/results.csv ppb_highres/section.tex --level 3
```

Without `--level`, it takes the finest level at which all six solvers ran at
two or more orders. Orders missing a solver (skipped by `max_unknowns`) are
left out of the tables.

Report the machine you ran on: edit the `TODO(author)` line in the generated
section. For timings that are comparable across runs, keep the machine
otherwise idle and plugged in.
