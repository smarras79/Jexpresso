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
DBG_SCHUR=0 sbatch submit_Jexpresso_profile.sh    # five fields, 5·Np unknowns
DBG_SCHUR=1 sbatch submit_Jexpresso_profile.sh    # scalar Schur, Np unknowns
```

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
```

versus `... on all 5 fields`.

**The two arms are different splittings**, not two solvers for one problem.
`:imex_schur` forces the advective Θ row, without which the elimination does not
close on a single scalar; the two rows differ by 0.06% of the flux form. That is
the right comparison for wall clock and the wrong one for reading a state
difference between `output_imex_full/` and `output_imex_schur/` as an error.
What *is* pinned: given the same operator, the two stage solves return the same
five fields to ~1e-12 (`test/imex3d/test_schur_stage.jl`).
