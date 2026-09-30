# 3D periodic Poisson solver comparison (for Wulver)

The 2D benchmark in `tools/periodic_poisson_benchmark` is too small for iterative solvers to overtake the sparse direct solve. In 2D, nested-dissection Cholesky costs O(n^1.5), against O(n) for a good iterative method. In 3D it costs O(n²) in time and O(n^(4/3)) in memory, so the crossover moves to sizes a cluster can reach.

This directory runs the same comparison in 3D, up to about 1.7·10⁷ unknowns, as one SLURM job per configuration on NJIT's Wulver.

## The problem

$$-\nabla^2 u = f \quad\text{on } [0,2\pi]^3,\qquad u \text{ periodic in } x, y, z,$$

with the exact solution

$$u = A\big(p(x)p(y)p(z) - (c^2-1)^{-3/2}\big),\qquad p(s) = \frac{1}{c - \cos s},\qquad c = \frac{r + r^{-1}}{2}.$$

This is the 3D version of the 2D test problem. u has zero mean, A scales its peak to 1, and its Fourier coefficients decay like r^(|kx|+|ky|+|kz|), so it is not band-limited. The 3D default is r = 0.5. With the 2D value r = 0.8, the product of three sharp peaks is not resolved at feasible sizes.

## The discretization, and why it is assembled here

Jexpresso's linear-solve path is 2D-only: the Laplacian assembly, the periodic reduction and the static condensation all assume two dimensions, and they store every element matrix densely. On a Cartesian mesh with GLL quadrature (Jexpresso's default "inexact" quadrature, Q = N), the SEM stiffness matrix is exactly a sum of Kronecker products of assembled 1D matrices:

$$K = M_z\otimes M_y\otimes K_x + M_z\otimes K_y\otimes M_x + K_z\otimes M_y\otimes M_x,\qquad M = M_z\otimes M_y\otimes M_x \ \text{(lumped)}.$$

That is the same operator Jexpresso builds, and the operator of the CEED bake-off problem BP5. `poisson3d.jl` assembles it in seconds for 10⁷ unknowns.

Everything else comes from Jexpresso:
- the LGL nodes and weights (Kopriva's algorithm);
- the Lagrange basis derivatives;
- AMG (`jx_amg_setup`/`jx_amg_solve`);
- the skeleton solve of static condensation (`el_skeleton_solve`);
- the pseudo-spectral axis operators (`_collocation_axis` with Kopriva's `FourierDerivativeMatrix`);
- the FFT solver (`FFTPoissonSolver`).

**Verified** (the precompile job reruns both checks on Wulver):
- **`verify_2d.jl`**, the same code in 2D against Jexpresso's own periodic SEM solve on 16×16 elements: the errors agree to within 4·10⁻¹² at every N = 2…8. The discretizations are the same.
- **`verify_3d.jl`**:
  - the error falls exponentially in N: 0.45 at N = 2 down to 1.8·10⁻⁴ at N = 8, on 4³ elements;
  - the h-convergence rate is 4.03 at N = 3, as expected (N + 1 = 4);
  - the five SEM solvers agree to 3·10⁻¹³, and pseudo-spectral and FFT agree to 3·10⁻¹⁵.

## The seven solvers

| key | solver |
|---|---|
| `sem` | sparse Cholesky (CHOLMOD) of the full SEM system, with METIS nested-dissection ordering |
| `sem_amg` | smoothed-aggregation AMG + CG, full system |
| `sem_jacobi` | Jacobi-preconditioned CG, full system (the CEED BP5 solver) |
| `sc_direct` | element-level static condensation; skeleton system by Cholesky (METIS ordering) |
| `sc_amg` | the same condensation; skeleton system by AMG + CG |
| `ps` | pseudo-spectral Fourier collocation, matrix diagonalization, O(N_g⁴) in 3D |
| `fft` | FFT (FFTW, planned with `ESTIMATE`) |

Solver details:
- **Singular systems:** the SEM systems pin one unknown and are shifted to zero M-weighted mean afterwards.
- **CG tolerance:** CG stops at 10⁻¹² relative preconditioned residual.
- **Static condensation:** it is the algorithm of `elementLearning_Axb!` (S^e = A_bb − A_bo A_oo⁻¹ A_ob per element, then assembly of the skeleton system, then recovery of the interiors), done at element level with a dense Cholesky factorization of A_oo per element.

**Why METIS.** CHOLMOD's own default here is AMD ordering. In 3D, AMD's fill grows like n^1.6. At n = 46 656, N = 3, AMD gives a factor with 88 M nonzeros in 6.7 s, against 31 M nonzeros in 1.5 s with METIS. METIS nested dissection gives the near-optimal n^(4/3). A fair direct baseline needs it, so it is the default (`ORDERING=amd` restores CHOLMOD's choice). METIS is loaded from Jexpresso's existing dependencies.

## Running it on Wulver

1. Set up Julia 1.11.9 (for example with juliaup), check out `sm/elementLearning`, and edit `slurm/wulver3d.env`:
   - `ACCOUNT` (your PI's UCID) and `QOS`;
   - `JULIA` and `REPO`;
   - `SWEEPS`, the solvers and the thread counts.
2. Look at the plan without submitting anything:
   ```bash
   DRY_RUN=1 bash tools/poisson3d_benchmark/slurm/submit_wulver3d.sh
   ```
   It prints every job: unknowns, threads, cores, memory, time and partition. Configurations that cannot fit are listed as skipped, with the reason (memory above `BIGMEM_MAX_GB`, or time above `MAX_HOURS`).
3. Submit:
   ```bash
   bash tools/poisson3d_benchmark/slurm/submit_wulver3d.sh
   ```
   This submits one precompile + verification job, then one job per configuration (after the first succeeds), then a merge job that draws the figures after all of them end.
4. Rerun the same command any time. It resubmits only the configurations without an `ok` result: failed, timed out, preempted, or new ones.

**Parallelism.**
- *Across configurations:* every configuration is its own job, so they run in parallel across Wulver's nodes, and each job's peak memory is its own.
- *Within a configuration:* each job runs `THREADS_*` threads, chosen by problem size and the same for every solver at a given n.
  - Julia threads parallelize the element loops of static condensation.
  - BLAS threads parallelize CHOLMOD's supernodal factorization and the dense kernels.
  - AlgebraicMultigrid.jl and the CG iterations run on one thread, so threads favour the direct solvers. Set all `THREADS_*` to 1 for a strictly single-core comparison of the algorithms.
- *Not included:* distributed-memory (MPI) solves, which would need PETSc or MUMPS.

**Resource model** (in `submit_wulver3d.sh`), calibrated here in 3D for N = 2…6 up to 1.1·10⁵ unknowns:
- nnz(K) = (3N+4)·n;
- METIS factor nnz ≈ 25·n^(4/3);
- skeleton n_s = n·(1 − ((N−1)/N)³), with nnz(B) ≈ 0.7·(N+1)³·n_s and a factor of about 5·(N+1)^1.5·n_s^(4/3).

Memory requests carry `MEM_SAFETY` (×1.5). Wulver charges MAX(cores, memory/4 GB), so every job also gets memory/4 GB cores at no extra cost. Jobs above 500 GB go to `bigmem` (2 TB). With the default sweeps (N = 2 and 4, n from 4·10³ to 1.7·10⁷), there are 122 jobs:
- the direct solvers need bigmem at about 7·10⁶ unknowns (about 800 GB and about 8 h);
- they are skipped at 1.7·10⁷ unknowns (about 2.5 TB), where AMG needs 56 GB and less than an hour.

## Outputs (in `OUTDIR`)

- `results.csv` and `results.md`: one row per configuration, containing:
  - errors;
  - assembly, rhs, setup, solve and total seconds;
  - CG iterations;
  - nnz of K, of the Cholesky factor and of the skeleton matrix;
  - peak memory;
  - threads;
  - status.
- `assets/`: per SEM order N, light and dark SVG figures:
  - `ppb3d_cost_N<N>`: setup + solve against n, with dashed slopes n and n²;
  - `ppb3d_total_N<N>`, `ppb3d_memory_N<N>`, `ppb3d_iters_N<N>`, `ppb3d_error_N<N>`.
- `parts/<solver>_ne<ne>_N<N>/`: each job's own results.
- `logs/`: one log per job, plus the settings used.

## Running locally

```bash
julia --project=. -t 4 tools/poisson3d_benchmark/verify_3d.jl
julia --project=. -t 4 tools/poisson3d_benchmark/bench3d.jl \
      --solvers sem,sem_amg,sem_jacobi,sc_direct,sc_amg,ps,fft --nes 4,8,12 --nops 2,4 \
      --outdir ppb3d_local --resume
python3 tools/poisson3d_benchmark/plot3d.py ppb3d_local
```

The timing protocol is the one from 2D: a discarded warm-up run (`--warmup small`: 2³ elements, same solver and N), garbage collection, then the recorded second run. The pseudo-spectral solver needs an even grid size ne·N.
