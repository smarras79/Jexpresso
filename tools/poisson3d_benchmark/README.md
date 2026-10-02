# 3D periodic Poisson solver comparison (for Wulver)

The 2D benchmark in `tools/periodic_poisson_benchmark` is too small for iterative solvers to overtake the sparse direct solve. In 2D, nested-dissection Cholesky costs O(n^1.5), against O(n) for a good iterative method. In 3D it costs O(n²) in time and O(n^(4/3)) in memory, so the crossover moves to sizes a cluster can reach.

This directory runs the same comparison in 3D, up to about 1.7·10⁷ unknowns, in one SLURM job on one node of NJIT's Wulver. The threaded run goes through Jexpresso itself (`run_case` on a 3D deck). A separate MPI run uses distributed-memory solvers: MUMPS and hypre BoomerAMG.

## The problem

$$-\nabla^2 u = f \quad\text{on } [0,2\pi]^3,\qquad u \text{ periodic in } x, y, z,$$

with the exact solution

$$u = A\big(p(x)p(y)p(z) - (c^2-1)^{-3/2}\big),\qquad p(s) = \frac{1}{c - \cos s},\qquad c = \frac{r + r^{-1}}{2}.$$

This is the 3D version of the 2D test problem. u has zero mean, A scales its peak to 1, and its Fourier coefficients decay like r^(|kx|+|ky|+|kz|), so it is not band-limited. The 3D default is r = 0.5. With the 2D value r = 0.8, the product of three sharp peaks is not resolved at feasible sizes.

## The runs go through Jexpresso

The threaded 3D benchmark runs Jexpresso's own solvers on the 3D deck `problems/Elliptic/poisson_periodic_sem_3d`. It calls `Jexpresso.run_case("Elliptic", "poisson_periodic_sem_3d"; inputs = overrides)` once per configuration, as the 2D study does with `poisson_periodic_sem` (`jexpresso3d.jl`; `bench3d.jl --engine jexpresso`, the default). The deck uses the built-in Cartesian grid: `:nelx = :nely = :nelz = ne` hexahedra on [0,2π]³, periodic in x, y and z, with the exact solution above (r = 0.5).

| key | Jexpresso path (deck flags as in 2D) |
|---|---|
| `sem` | `sem_setup` (3D Laplacian: `DSS_laplace_sparse_3D`) → `periodic_sem_system` → `periodic_sem_factorize` (sparse Cholesky, METIS ordering: `:sparse_ordering => "metis"`) → `periodic_sem_direct_solve` |
| `sem_amg` | … → `periodic_sem_amg_solve` (`jx_amg_setup` / `jx_amg_solve`: smoothed-aggregation AMG + CG) |
| `sc_direct` | … → `periodic_sem_sc_solve` → `elementLearning_Axb!` → `el_sc_solve!` (threaded condensation into the skeleton matrix B) → Cholesky of B (METIS ordering) |
| `sc_amg` | … → `elementLearning_Axb!` → `el_sc_schur_cg!` (`:EL_skeleton_solver => "amg"`): CG on the skeleton Schur complement, matrix-free, with tensor-product interior solves, preconditioned by the full-system AMG restricted to the skeleton |
| `pmg_amg` | … → `periodic_sem_pmg_solve` (`:linsolve_pmg => "amg"`): CG on the full SEM system, preconditioned by a p-multigrid V-cycle over the SEM orders N, N/2, …, 1, with AMG on the p = 1 level |
| `pmg_gmg` | the same with `:linsolve_pmg => "gmg"`: geometric h-multigrid on the p = 1 grid in place of AMG |
| `ps` | `pseudospectral_linsolve!` → `FourierCollocationPoissonSolver3D` (O(N_g⁴)) |
| `fft` | `fft_linsolve!` → `FFTPoissonSolver` (FFTW, `ESTIMATE`) |

The timings are Jexpresso's per-phase timers of that run (`JX_TIMINGS`): `assembly` = `sem_setup` (high-order mesh, metrics, mass and Laplacian; SEM only), then `rhs`, `setup` and `solve` as the solvers record them. CG stops at 10⁻¹² relative preconditioned residual. Jexpresso's metric type carries the array sizes, so a new mesh size recompiles: every configuration is warmed up on itself, and only the second run is recorded.

**What was added to Jexpresso for 3D** (2D results are bit-for-bit unchanged, checked for all four SEM solvers):
- **`DSS_laplace_sparse_3D`** (`src/kernel/infrastructure/element_matrices.jl`): the 3D stiffness matrix, assembled straight into sparse triplets with no element matrices. With collocation (Q = N), a basis function's gradient at a quadrature point is non-zero only along the three grid lines through it, the structure the radiative-transfer operator uses (`sparse_lhs_assembly_3Dby2D`). It handles general (skewed, curved) hexahedra, with an optional diffusivity a(x,y,z). It is exactly symmetric, so `factorize` picks Cholesky. It is type-stable, sizes its triplet arrays exactly once, and fills them in parallel. Against a brute-force dense quadrature on sheared elements it agrees to 3·10⁻¹⁴.
- **The periodic solve in 3D** (`periodic_sem.jl`): periodic boundary faces are detected, and the right-hand side takes z.
- **Static condensation in 3D** (`elementLearningStructs.jl`, `periodic_sem_sc_solve`):
  - (ngl−2)³ interior nodes per element, with a boundary-first connectivity built from `connijk`.
  - The per-element blocks that held the same values in pairs now share storage (4 arrays instead of 9), in 2D too.
  - A build with no model allocates no inference buffers, and no skeleton-matrix copy for them.
- **The 3D pseudo-spectral solver** (`fourier_collocation.jl`): a solve allocates nothing; the FFT and pseudo-spectral drivers have a 3D grid for a 3D deck.
- **3D lumped mass without element matrices** (`DSS_mass_collocation_3D!`): with collocation (Q = N) the element mass matrix is diagonal, but `build_mass_matrix!(::NSD_3D, ::Inexact)` built it dense, (N+1)³ × (N+1)³ per element with (N+1)⁹ operations (17 GB at N = 8 on 16³ elements). The mass vector is now assembled straight from the weights, bit for bit the same, on the CPU without AMR. `sem_setup` at 6³ elements, N = 6 went from 17.1 to 2.3 s; at 8³, N = 6 from 46.6 to 9.7 s.
- **Static-condensation recovery reuses the condensation's inverse**: `elementLearning_Axb!` inverted each element's interior block twice, once to condense and once to recover the interiors. The condensation now keeps the inverse in place of the block (no extra memory) and the recovery uses it, bit for bit the same. The recovery at 6³ elements, N = 6 went from 0.099 to 0.013 s.
- **Static condensation rewritten for speed and memory** (`src/kernel/solvers/static_condensation.jl`), used whenever the condensation is a solver (`:lstatic_condensation`), 2D and 3D. The element-learning sampling still runs the original loop, as does `:EL_sc_kernel => "legacy"`.
  - **`el_sc_solve!`** (`sc_direct`, and `sc_amg` with `:EL_sc_amg => "assembled"`):
    - **Per-element algebra:** each element's interior block is factored with Cholesky (LAPACK `potrf`), and its Schur-complement contribution is formed with BLAS (`trsm`, `syrk`). This replaces a dense `inv` and a scalar triple loop.
    - **Only active boundary nodes enter:** with collocation the interior couples only to the face-interior nodes. Edge and vertex pairs were stored zeros in B, which made AMG-CG on B twice as expensive.
    - **Assembly into B's own pattern:** B is assembled directly into its sparsity pattern, with the elements coloured so the threaded assembly has no write conflicts. The old route went through triplets and `sparse()`, 14 GB of triplets at 16³ elements, N = 8.
    - **Exactly symmetric B:** no symmetrization copy and no pinned-submatrix copy.
    - **Nothing stored per element:** the recovery refactors each interior block instead.
    - **Speed:** the condensation is 12–20× faster than the legacy loop (4³ elements at N = 5, 8³ at N = 6).
  - **`el_sc_schur_cg!`** (`sc_amg`, the default `:EL_sc_amg => "schur"`) uses the SEM tensor product:
    - **Fast-diagonalization interior solves:** on an affine hexahedron the interior block is c_x ω̂⊗ω̂⊗K̂ + c_y ω̂⊗K̂⊗ω̂ + c_z K̂⊗ω̂⊗ω̂. Fast diagonalization inverts it in O(N⁴) from three numbers per element, read off the matrix and checked entry by entry.
    - **Matrix-free Schur complement:** S p = (Â [p; −Â_oo⁻¹ Â_ob p])_b, one sparse matvec plus the interior solves per iteration. B is never formed.
    - **Preconditioner:** the full-system AMG V-cycle restricted to the skeleton. Since (A⁻¹)_bb = S⁻¹, it is as good for S as the V-cycle is for A.
    - **Fallback:** elements that are not affine fall back to `el_sc_solve!`.
    - **Result:** SC AMG costs one full-system AMG setup and about as many iterations as SEM AMG-CG (57 against 63 at 8³ elements, N = 6), with no skeleton matrix. The assembled-B route needed 34 iterations, but each was much more expensive, and its AMG setup was on a far denser matrix.
- **p-multigrid** (`src/kernel/solvers/pmultigrid.jl`, `jx_pmg_setup` / `jx_pmg_cg`): a preconditioner built for high-order spectral elements.
  - **Levels:** the SEM orders N → ⌊N/2⌋ → … → 1, all on the same elements.
    - The order-N level is the matrix K of the solve itself (`sem_setup` → `DSS_laplace_sparse_3D` → periodic reduction).
    - Each coarser order is rediscretized with Jexpresso's LGL basis of that order (`basis_structs_ξ_ω!`, `build_Interpolation_basis!`) and the element geometry from the metrics: on a box element, A_e = c_x ω̂⊗ω̂⊗K̂ + c_y ω̂⊗K̂⊗ω̂ + c_z K̂⊗ω̂⊗ω̂.
    - The tensor form is checked against K's diagonal, and Jexpresso's per-element axis orientation is honored.
  - **Transfers between orders:** the tensor-product Lagrange interpolation ℓ_a^(p)(ξ_k^(q)), applied element by element in three 1D passes. The work is threaded, with elements coloured for the restriction.
  - **Smoother:** Chebyshev–Jacobi of degree 3 on [0.25, 1.1]·λ_max(D⁻¹A), with λ_max from power iteration. It needs only threaded matvecs, no Gauss–Seidel sweep, and keeps the V-cycle symmetric.
  - **Coarse level** (p = 1, the element vertices):
    - `amg`: smoothed-aggregation AMG.
    - `gmg`: geometric h-multigrid with trilinear interpolation and Galerkin coarse operators, coarsening by 2 down to at most 512 nodes, then a sparse Cholesky solve.
  - **Singular system:** CG runs on the singular periodic system, with the mean projected out in the preconditioner.
  - **Applicability:** uniform axis-aligned box elements on a fully periodic box; anything else stops with an error.
  - **Cost:** 9–13 CG iterations at N = 4–6, against 29–63 for AMG-CG on the SEM matrix.
  - **Tuning:** `:pmg_degree` and `:pmg_lower`.
- **AMG-CG matrix-vector products are threaded** (`JXSymCSC`, `amg.jl`): for an exactly symmetric matrix, (A x)_j is the dot product of column j with x, so the product is threaded with no write conflicts. The AMG V-cycle (AlgebraicMultigrid.jl) stays serial.
- **`:sparse_ordering => "metis"`**: opt-in METIS nested-dissection ordering for both sparse Cholesky factorizations (full system and skeleton). The default stays CHOLMOD's own (AMD). In 3D, AMD's fill grows like n^1.6; at 13 824 unknowns (N = 4) Jexpresso's SEM direct setup takes 1.6 s with AMD and 0.53 s with METIS.

**The Kronecker engine** (`poisson3d.jl`, `bench3d.jl --engine kronecker`) is now an independent cross-check. On a Cartesian mesh with GLL quadrature, the SEM stiffness matrix is exactly K = M_z⊗M_y⊗K_x + M_z⊗K_y⊗M_x + K_z⊗M_y⊗M_x, built from Jexpresso's 1D LGL and Lagrange routines. It also has Jacobi-CG (`sem_jacobi`, the CEED BP5 solver) and 2D.

**Verified** (the SLURM job reruns all of these before the benchmark):
- **`verify_3d_jexpresso.jl`**: Jexpresso's six 3D solvers against the Kronecker engine, on 4³ elements at N = 2, 3 and 6³ at N = 4. The SEM errors agree to within 6.8·10⁻¹³; FFT and pseudo-spectral agree bit for bit or to 2·10⁻¹⁶.
- **`verify_2d.jl`**: the Kronecker engine in 2D against `Jexpresso.run_case` on the 2D deck (16×16 elements, N = 2…8), run live. The errors agree to within 7·10⁻¹².
- **`verify_3d.jl`**: the Kronecker engine converges as it should. The error falls exponentially in N, the h-convergence rate is 4.03 at N = 3, and the five SEM solvers agree to 3·10⁻¹³.

## Running it on Wulver

One script, `slurm/run_wulver3d.sbatch`, on one `general` node (128 cores, 128 × 4000 MB ≈ 500 GB, up to 72 h):
```bash
cd /project/smarras/smarras/Jexpresso      # checkout of sm/elementLearning
sbatch tools/poisson3d_benchmark/slurm/run_wulver3d.sbatch
```
It follows the Jexpresso job script:
1. `module load Julia/1.11.9` and `module load GCC MPICH`, then `MPIPreferences.use_system_binary()`;
2. `Pkg.instantiate(); Pkg.precompile()`, one serial process;
3. a serial warm-up (`using MPI; using Jexpresso`), then `verify_2d.jl`, `verify_3d.jl` and `verify_3d_jexpresso.jl`; the job stops if any of these fails;
4. six of Jexpresso's solvers side by side (`SOLVERS`: the iterative SEM solvers `sem_amg`, `sc_amg`, `pmg_amg`, `pmg_gmg`, and `ps`, `fft`; the direct solvers `sem` and `sc_direct` are left out, since they are never used for large problems; run them by hand with `bench3d.jl --solver sem` if wanted), one Julia process each with `THREADS = 21` Julia and BLAS threads (6 × 21 = 126 cores). Each runs the sweep below through `run_case`, one mesh level after the other with the orders in increasing order, into `OUTDIR/parts/<solver>/results.csv`, with a log in `OUTDIR/logs/<solver>.log`;
5. when all have finished, `plot3d.py` merges the results and draws the figures.

The settings are the few variables at the top of the script: `OUTDIR` (default `ppb3d_wulver`), `SOLVERS`, `THREADS`, the sweep and its size limits:
- **The sweep is the 2D benchmark's:** mesh levels of `LEVELS` = 8³, 16³, 32³ and 64³ elements, each with the SEM orders `NOPS` = 2…8. The Fourier solvers run on the (ne·N)³ grid, which has the same number of unknowns n = (ne·N)³.
- **Size limits:** configurations above these caps are not run.
  - The iterative SEM solvers (`sem_amg`, `sc_amg`, `pmg_amg`, `pmg_gmg`) stop at 7.1·10⁶ (`MAXN_SEM`). Jexpresso's 3D SEM infrastructure (high-order mesh, metrics, matrices) took 4.1 GB at 2.6·10⁵ unknowns. That projects to about 110 GB at 7.1·10⁶ and about 260 GB at 1.7·10⁷, too much with two of them on one node.
  - The Fourier solvers go on to 1.7·10⁷ (`MAXN`).
- **Figures:** `plot3d.py` leaves the direct solvers out of the figures; their rows stay in `results.md`. Pass `--direct` to draw them.

**Resubmitting** the same script resumes: configurations with an `ok` row are skipped, and failed or unfinished ones are rerun (for example after the time limit).

**What the parallelism means:**
- Julia threads parallelize Jexpresso's 3D Laplacian assembly.
- BLAS threads parallelize CHOLMOD's supernodal factorization and the dense kernels.
- AlgebraicMultigrid.jl and the CG iterations run on one thread, so threads favour the direct solvers.
- The six processes share the node's memory bandwidth, so timings are slightly pessimistic for all of them. For cleaner timings, set `SOLVERS` to one solver and `THREADS=128`, and submit once per solver.
- For distributed-memory (MPI) solvers, see the next section.

## The MPI version: MUMPS and BoomerAMG

The solvers above run in one process, with threads. `mpi/` solves the same problem with distributed-memory solvers on MPI ranks:

| key | solver |
|---|---|
| `mumps` | MUMPS sparse Cholesky (SPD LDLᵀ), distributed matrix input (ICNTL(18) = 3), METIS nested-dissection ordering, no low-rank compression: exact |
| `boomeramg` | hypre PCG, preconditioned by one BoomerAMG V-cycle |
| `jacobi` | hypre PCG with diagonal scaling (the CEED BP5 solver) |

How it works (`mpi/poisson3d_mpi.jl`):
- **Same system:** the same SEM operator, right-hand side, exact solution and pinned unknown as the serial benchmark.
- **Assembly without communication:** every row of K = Mz⊗My⊗Kx + Mz⊗Ky⊗Mx + Kz⊗My⊗Mx follows from the periodic 1D matrices, so each rank builds its own rows.
- **Partition:** the (ne·N)³ nodes are split into boxes over an MPI process grid (`MPI.Dims_create`, e.g. 8×4×4 for 128 ranks). Each rank's nodes are numbered contiguously, which is the row range hypre and MUMPS expect.
- **BoomerAMG settings:** HMIS coarsening, extended+i interpolation (at most 4 entries per row), l1-scaled symmetric Gauss-Seidel smoothing, strong threshold 0.5. With 0.5, a 2.6·10⁵-unknown case takes 14 iterations at N = 4 and 13 at N = 2; 0.25 converges the same but sets up slower, and 0.7 needs 20 iterations.
- **CG stopping test:** 10⁻¹² relative preconditioned residual, as in the serial benchmark.
- **MUMPS ordering:** sequential METIS on the host. The MUMPS_jll binary has no parallel ordering: PT-SCOTCH returns INFOG(1) = −38 and ParMETIS fails.
- **MUMPS right-hand side and solution:** both are gathered on the host. That costs two n-vectors on rank 0 and is negligible next to the factorization.
- **Timing:** every phase starts and ends at a barrier, so a time is the slowest rank's.
- **MPICH_jll is pinned to 4.3** (`mpi/Project.toml`). With the system MPI, the MPI-dependent binaries (HYPRE_jll, MUMPS_jll, SCALAPACK32_jll) still load the Fortran MPI library of MPICH_jll on top of the system `libmpi`, so the two must be the same MPICH series. Wulver has MPICH 4.3.0; MPICH_jll 5 fails there with `libmpifort.so: undefined symbol: MPIR_fortran_false`. On a cluster with another MPICH version, change that compat entry to match.
- **No Jexpresso on the ranks:** the ranks do not load Jexpresso (several hundred MB per rank). The LGL nodes and weights are computed with Kopriva's algorithm, the one Jexpresso implements, in a separate small environment (`mpi/Project.toml`: MPI, HYPRE.jl, MUMPS.jl).

**Verified** (`mpi/verify_mpi.jl`, rerun by the SLURM job):
- **Against the serial benchmark:** the errors equal those of the serial, Jexpresso-based `sem` solve at the same (ne, N), for ne = 4, 6 and N = 2, 3, 4, to within 4.5·10⁻¹².
- **Solver agreement:** the three MPI solvers agree to 2·10⁻¹³.
- **h-convergence:** the rate is 4.02 at N = 3, over 4³ → 16³ elements.
- **Rank count:** the results are the same on 1, 3 and 4 ranks, including uneven splits, to 3·10⁻¹³.

**Running it on Wulver:**
```bash
cd /project/smarras/smarras/Jexpresso      # checkout of sm/elementLearning
sbatch tools/poisson3d_benchmark/slurm/run_wulver3d_mpi.sbatch
```
This is the same Jexpresso job template as `run_wulver3d.sbatch` (one node, 128 tasks):
1. MPIPreferences `use_system_binary()` for both environments;
2. resolve and precompile both environments;
3. a serial warm-up, then a small serial reference run (`bench3d.jl --solvers sem`), then `verify_mpi.jl` on 8 ranks against that reference;
4. `mpirun -bind-to core -np 128` for each solver in turn, so timings are clean, one `mpirun` per mesh level (a killed run only loses the rest of that level), into `ppb3d_wulver_mpi/parts/<solver>/`, with logs in `ppb3d_wulver_mpi/logs/`. The sweep is the same as in the threaded run: `LEVELS` = 8³, 16³, 32³, 64³ elements × `NOPS` = 2…8;
5. `plot3d.py`, which draws the same figures as for the threaded run (memory summed over ranks).

**Sizes:**
- The iterative solvers go to 6·10⁷ unknowns (`MAXN`).
- MUMPS stops at 2.2·10⁶ unknowns (`MAXN_DIRECT`), which fits one node's 500 GB.
- With more nodes (`--nodes`, and `NP`), raise `MAXN_DIRECT`: MUMPS spreads the factor over all ranks.
- Resubmitting resumes.

**Local run:**
```bash
julia --project=tools/poisson3d_benchmark/mpi -e 'using Pkg; Pkg.instantiate()'
julia --project=tools/poisson3d_benchmark/mpi -e 'using MPI; run(`$(MPI.mpiexec()) -n 4 julia --project=tools/poisson3d_benchmark/mpi tools/poisson3d_benchmark/mpi/verify_mpi.jl`)'
julia --project=tools/poisson3d_benchmark/mpi -e 'using MPI; run(`$(MPI.mpiexec()) -n 4 julia --project=tools/poisson3d_benchmark/mpi tools/poisson3d_benchmark/mpi/bench3d_mpi.jl --solvers mumps,boomeramg,jacobi --nop 2 --nes 4,8,16 --outdir ppb3d_mpi_local`)'
python3 tools/poisson3d_benchmark/plot3d.py ppb3d_mpi_local
```

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
  - **the figures of the 2D benchmark:**
    - per mesh level, `ppb3d_error_vs_solve_time_ne<ne>`, `ppb3d_error_vs_cost_ne<ne>` (setup + solve), `ppb3d_error_vs_total_time_ne<ne>` and `ppb3d_error_vs_order_ne<ne>`;
    - for all levels, `ppb3d_error_vs_Ng`: error against unknowns per direction ∛n = ne·N, one SEM curve per mesh, the Fourier solvers and the predicted 0.5^(∛n/2);
  - `ppb3d_total_N<N>`, `ppb3d_memory_N<N>`, `ppb3d_iters_N<N>`, `ppb3d_error_N<N>`.
- `parts/<solver>/`: each solver's own results.
- `logs/<solver>.log`: each solver's output (setup and checks go to `ppb3d.<jobid>.out` in the submit directory).

## Running locally

```bash
julia --project=. -t 4 tools/poisson3d_benchmark/verify_3d_jexpresso.jl
julia --project=. -t 4 tools/poisson3d_benchmark/bench3d.jl \
      --solvers sem,sem_amg,sc_direct,sc_amg,ps,fft --nes 4,8 --nops 2,3,4 \
      --outdir ppb3d_local --resume
python3 tools/poisson3d_benchmark/plot3d.py ppb3d_local
```

To run one configuration and look at it in ParaView:

```bash
julia --project=. -t 4 tools/poisson3d_benchmark/run_one3d.jl --solver sem --ne 8 --nop 4 --outdir vis3d
```

The SEM solvers write `vis3d/Elliptic/poisson_periodic_sem_3d/output/iter_1.pvtu` (field `u`). `ps` and `fft` write `pseudospectral_laplace.vtk` or `fft_laplace.vtk` there, on the uniform (ne·N)³ grid, with the fields `u`, `u_exact` and `error`. The script prints the path and the L∞ error.

Or one configuration straight through Jexpresso:

```julia
using Jexpresso
Jexpresso.run_case("Elliptic", "poisson_periodic_sem_3d";
                   inputs = Dict(:nop => 4, :nelx => 8, :nely => 8, :nelz => 8,
                                 :linsolve_amg => true, :sparse_ordering => "metis"))
```

The timing protocol is the one from 2D: a discarded first run (with `--engine jexpresso`, the identical configuration; with `--engine kronecker`, `--warmup small` uses 2³ elements), garbage collection, then the recorded second run. The pseudo-spectral solver needs an even grid size ne·N. `--engine kronecker` adds `sem_jacobi` and `--d 2`.
