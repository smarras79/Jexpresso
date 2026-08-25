# CompEuler/rtb3d_schur

A 3D rising thermal bubble, built to measure the **scalar Schur stage solve**
(`:imex_schur`) against the standard five-field one under IMEX3D.

## The physical problem

A warm bubble in a closed, neutrally-stratified box of dry air at rest. It is
buoyant, rises, deforms into a mushroom cap under its own shear, and is held
together by an SGS viscosity. Standard benchmark physics — the bubble is
`CompEuler/rtb_imex`'s, unchanged, so nothing here is novel except the geometry.

| | |
|---|---|
| Domain | x, y ∈ [−5, 5] km, z ∈ [0, 10] km — a 10 km cube |
| Background | θ_ref = 300 K, hydrostatic, at rest |
| Perturbation | Δθ = 2 K at the centre, linear taper θ_c(1 − r/r₀) |
| Bubble | centre (0, 0, 2500) m, radius r₀ = 2000 m |
| Boundaries | free-slip (no-flux) on **all six** faces |
| Viscosity | AV, μ = 125 m²/s on the four momentum/θ equations |
| Integrator | `IMEX_ARK(:ARS343)`, Δt = 0.6 s, t_end = 1000 s |

**The bubble is a sphere, not a cylinder**, and that is deliberate. `rtb_imex`
measures r in the x–z plane only, which leaves ρu and ρv identically zero. The
Schur reduction eliminates ρu and ρv along with everything else, so measuring it
where two of the eliminated fields are exactly zero would report a saving a real
3D flow need not repeat.

## The grid

| | |
|---|---|
| Elements | 20 × 20 × 80, nop = 4 (LGL) |
| Gridpoints | 81 × 81 × 321 = **2,106,081** |
| Element size | h_x = h_y = 500 m, h_z = 125 m |
| Smallest LGL gaps | 86.4 m horizontal, 21.6 m vertical |
| **Anisotropy** | **h_x/h_z = 4:1** |
| Element columns | 400 |

The aspect ratio is the point. The column preconditioner is exact for the
vertical acoustic operator and does nothing for the horizontal one, so what it
leaves for the Krylov iteration scales with h_x/h_z. 4:1 is where the mock sweep
found the reduction most favourable.

`rtb3d_20x20x80.msh` is committed. Rebuild it, or a variant, with
`julia --project=. problems/CompEuler/rtb3d_schur/generate_mesh.jl [nx ny nz]`
— it drives the gmsh SDK that GridapGmsh already ships, so no `gmsh` on PATH is
needed. `DBG_MESH=10x10x40` selects the smaller variant (same 4:1, 270,681
points), which fits on one rank.

## Running it

```bash
sbatch submit_Jexpresso_profile.sh full     # five fields, 5·Np unknowns
sbatch submit_Jexpresso_profile.sh schur    # scalar Schur, Np unknowns
```

The positional argument is preferred to `DBG_SCHUR=1 sbatch …`: a site that
defaults `sbatch` to `--export=NONE` drops the variable and runs the baseline
under a banner saying otherwise. A typo in the argument exits 2 rather than
guessing.

There is a third, *diagnostic* configuration, reached by adding one line to
`user_inputs.jl` rather than by an argument:

```julia
:imex_schur_kernel    => false,     # reference matvec, ~6× slower
```

This runs the Schur reduction with two full five-field operator applications
per matvec instead of the bespoke sweeps. It is not a configuration to run for
results — it is an independently written statement of the same operator,
agreeing with the fast one to 1.9e-16, and it is what the first cluster profile
of this path measured. Use it only if a Schur run looks *wrong* rather than
slow: a difference there is a kernel bug, no difference means look at the
reduction.

**25 ranks, one node** (`pick_nranks.jl`: a 5×5 rank grid over the 20×20
columns, 16 columns and 84,243 points each). Keep both arms at the same count:
the MPI reduce is the one part of the stage solve the reduction cannot shrink,
so thin ranks inflate the term that cannot improve. Memory: ~2.2 GB/rank
expected; take 8 GB (`--cpus-per-task=2`) on the first run because the *mesh
broadcast*, not the solve, is the peak.

**One rank will not work** — the full mesh peaks at 13.4 GB during setup.

## What to expect

Measured on the 10 × 10 × 40 variant, one rank, Δt = 0.6, rtol 1e-8:

| | s/step (steady) | cold iterations | total wall |
|---|---|---|---|
| five-field | 19.46 | 20 | 225.8 s |
| scalar Schur | **16.15** | **61** | 196.9 s |
| | **1.21× faster** | 3× more | |

**Expect the iteration count to go the wrong way.** The mock predicted 2.67×
*fewer* iterations at this anisotropy; on a real mesh the scalar system needs
3× *more*. The reduction still wins, entirely on per-iteration cost — one
implicit field instead of five cuts the matvec, the gather/scatter and the
orthogonalisation, which on the 64×64×60 profile are 73% of the stage solve
against 6.8% for the non-scaling MPI reduce. A rise in iterations/solve under
`DBG_SCHUR=1` is the expected result, not a fault. Read the step time.

### What 25 ranks on the full mesh actually gave

The reduction did everything it was designed to do except the one term that
dominates:

| per step, 25 ranks, 20×20×80, CFL_h 1.05 | five-field | Schur, reference H | |
|---|---|---|---|
| matvec | 4.191 | 6.104 | **0.69× — worse** |
| preconditioner | 1.688 | 0.342 | 4.9× better |
| banded solve | 1.495 | 0.217 | 6.9× better |
| orthogonalise | 2.110 | 0.171 | 12.3× better |
| MPI reduce | 0.507 | 0.057 | 8.9× better |
| warm iterations/solve | 23.1 | 12.1 | 1.92× fewer |
| **step** | **9.581** | **9.152** | 1.047× |

Note the iteration count went the *right* way here, 1.92× fewer — the opposite
of the one-rank result above, and the reason that table is now only of
historical interest. Warm-started production and a cold-start self-check
disagree about this path; production is the one that counts.

The matvec was the problem. `schur_H!` was built out of two full five-field
`hevi_apply_A!` calls — deliberately, to reuse verified code — so the *scalar*
matvec cost 2.05× **one** five-field application. Halving the iterations and
then doubling the cost of each is a wash, and it ate a 12.3× win everywhere
else.

`schur_kernel.jl` replaces those two calls with the two sweeps H actually
needs. Measured against the reference form it replaces, and agreeing with it to
1.9e-16:

| | |
|---|---|
| H, reference form | 2.09× one `hevi_apply_A!` |
| H, kernel | 0.36× one `hevi_apply_A!` |
| **whole stage solve** | **4.06× faster** (`tools/schur_kernel_e2e.jl`) |

The 4.06× is end-to-end over a complete `imex3d_solve_schur!` at identical
iteration counts, not the matvec alone.

### The measured result: 3.56×

Both arms run back to back on the same nodes, 25 ranks, 20×20×80, CFL_h 1.05,
Δt 0.6, rtol 1e-8, 50 profiled steps:

| s/step | five-field | Schur + kernel | |
|---|---|---|---|
| rhs! | 0.799 | 0.742 | untouched |
| f_imp | 0.297 | 0.351 | 0.85× — see below |
| **stage solve** | **9.542** | **1.880** | **5.08×** |
| — matvec | 4.955 | 1.092 | 4.54× |
| — preconditioner | 2.042 | 0.250 | 8.18× |
| — of which banded | 1.670 | 0.193 | 8.64× |
| — orthogonalise | 1.922 | 0.199 | 9.67× |
| — MPI reduce | 0.357 | 0.053 | 6.69× |
| — Schur setup + recover | — | 0.199 | per solve, not per iteration |
| iterations/step | 69.4 | 36.2 | 1.92× fewer |
| **step** | **10.651** | **2.990** | **3.56×** |
| stage solve share of step | 89.6% | 62.9% | |

The matvec now wins on both counts: 1.92× fewer of them, and each one 2.37×
cheaper *despite* the five-field figure being an average over a cheaper mix.

`f_imp` going the wrong way by 18% is expected, not a regression: `:imex_schur`
forces the advective Θ row, so the two arms compute a genuinely different f_imp.

**Run-to-run variation on this cluster is around 13%**, measured: two `schur`
runs of the identical configuration gave 3.398 and 2.990 s/step. Quote the
per-term ratios, which are 4–10×, rather than the third digit of the headline.

### Where the time is now

The stage solve is down from 89.6% of the step to 62.9%, and `rhs!` + `f_imp`
are now 36.5% — the same order as the matvec. Amdahl has started to bite, so
further work on the stage solve returns less than it did.

The cheapest remaining lever is the iteration count, not per-iteration cost:
12.1 iterations/solve at `rtol 1e-8`, which is tighter than a third-order
tableau needs. `DBG_RTOL=1e-6` is one job and would say whether ~2
iterations/solve come off for free. After that it is horizontal
preconditioning, which is real work — the column preconditioner is exact
vertically and does nothing for the horizontal acoustic coupling that
CFL_h = 1.05 puts in.

**Treat 1.21× as a lower bound from one rank and nothing more.** Two effects
pull opposite ways: one rank has no MPI reduce at all, which *flatters* Schur;
and that variant's CFL_h = 0.526 is half the full mesh's, keeping the stage
solve a smaller share of the step, which *understates* it. The full mesh on 25
ranks is the number that settles it.

What to compare, in order of value: **(1)** `measured step` from the profile
block; **(2)** the profile *split*, which shows where the saving comes from
rather than that there is one; **(3)** iterations/solve.

## Two things to check before believing a result

**Which arm actually ran.** The submit banner prints `DBG_SCHUR` from the
submitting shell, so under a launcher that does not propagate environment
(OpenMPI's `mpirun`) it would say `DBG_SCHUR=1` while every rank ran the
five-field solve — a dead heat with nothing explaining it. The authoritative
line comes from inside the ranks, in the setup report:

```
Stage solve: preconditioned GMRES on the SCALAR SCHUR system
    ...
    matvec: bespoke scalar sweeps (~0.36x one 5-field apply)
```

versus `... on all 5 fields`. The `matvec:` line says which form of H is in
use, which a profile cannot be read without.

**The two arms are different splittings**, not two solvers for one problem.
`:imex_schur` forces the advective Θ row, without which the elimination does not
close on a single scalar; the two rows differ by 0.06% of the flux form. That is
the right comparison for wall clock and the wrong one for reading a state
difference between `output_imex_full/` and `output_imex_schur/` as an error.
What *is* pinned: given the same operator, the two stage solves return the same
five fields to ~1e-12 (`test/imex3d/test_schur_stage.jl`).
