# CompEuler/rtb2d_schur — the desktop rising thermal bubble

The 2D analogue of `CompEuler/rtb3d_schur`, at 1/16 of the points: the same
bubble, the same mesh spacing, the same Δt, the same Schur A/B, sized to be
re-run on a laptop or a workstation instead of a queue slot.

## The physical problem

A warm bubble at rest in a neutrally stratified, dry, hydrostatic atmosphere,
released at t = 0 and left to rise under its own buoyancy while sound waves run
freely through the box.

| | |
|---|---|
| domain | 10 km (x) × 10 km (z), −5000 ≤ x ≤ 5000, 0 ≤ z ≤ 10000 |
| background | θ = 300 K uniform, hydrostatic p, dry ideal gas |
| bubble | Δθ = 2 K·(1 − r/r₀) for r < r₀, **cylindrical**: r = √((x−x_c)² + (z−z_c)²) |
| | centre z_c = 2500 m, radius r₀ = 2000 m, x_c the domain centre |
| initial velocity | zero everywhere |
| boundaries | free slip, no flux, on every face |
| viscosity | constant, μ = 125 m²/s on the four non-density equations |
| Δt | 0.6 s (IMEX), tend 1000 s for the physics, 60 s for a timing run |

The bubble deforms into the usual mushroom over ~1000 s. **Plot `dθ`, not `θ`** —
`θ` is the total, 300 K of background with 2 K on top, so it autoscales to the
stratification and the bubble disappears.

## It is a slab, and that is not a shortcut

HEVI and IMEX3D refuse a 2D mesh outright (`params_setup.jl`: *"IMEX3D is 3D
only"*), and the refusal is protective rather than lazy:

- `build_column_topology` reads `coords[3,:]` and builds its column catalogue
  from distinct (x, y) pairs. A 2D mesh shapes `coords` as `(2, npoin)` and puts
  the vertical along **y**, so that is both a `BoundsError` and the wrong axis.
- every element kernel is a triple `i,j,k` loop over `ngl³` nodes of
  `connijk[iel,i,j,k]`, which in 2D is shaped `(nelem, ngl, ngl, 1)`.
- the operator is contracted with `dζd{x,y,z}`, which a 2D mesh **allocates and
  never fills**. Without the guard the acoustic operator would assemble to zero
  — silently wrong, not dead. That is the one that matters.

So this case uses the idiom `rtb_hevi` and `rtb_imex` already use: a 3D mesh
**one element thick in y**, free slip front and back, and a **cylindrical**
bubble. v is zero at t = 0, free slip holds it at zero on both y faces, and
nothing in a y-invariant problem generates it — the answer is exactly the 2D
one while the solver runs its ordinary 3D path. Nothing is special-cased.

## The grid

```
20 × 1 × 80 elements, nop = 4   ->   81 × 5 × 321  =  130,005 gridpoints
h_x = 500 m      h_z = 125 m    ->   4:1
```

Held fixed from the 3D case on purpose. The column preconditioner is exact for
the vertical acoustic operator and does nothing for the horizontal one, so what
it leaves for the Krylov iteration — and therefore how much the Schur reduction
can save — scales with the **acoustic anisotropy** h_x/h_z. Same spacing, same
anisotropy, same nop, same Δt, same CFL_h ≈ 1.05 as `rtb3d_schur`. The only
change is 20 elements of y collapsing to 1.

y is one element **500 m** wide, isotropic with h_x. A thinner slab would raise
h_x/h_y and hand the horizontal operator a stiffness the 3D case does not have
— the one thing this benchmark exists to hold fixed.

`DBG_MESH=10x1x40` gives the same 4:1 at 33,005 points, for a smoke test.
Rebuild either, or a variant, with

```bash
julia --project=. problems/CompEuler/rtb2d_schur/generate_mesh.jl [nx ny nz]
```

which drives the gmsh SDK GridapGmsh already ships — no `gmsh` on PATH needed.
Keep `ny = 1` and `nz = 4·nx`.

## Running it

```bash
./problems/CompEuler/rtb2d_schur/run_rtb2d.sh full     # five fields, 5·Np
./problems/CompEuler/rtb2d_schur/run_rtb2d.sh schur    # scalar Schur, Np
./problems/CompEuler/rtb2d_schur/run_rtb2d.sh all      # both, in order
```

A third, diagnostic configuration is reached by adding `:imex_schur_kernel =>
false` to `user_inputs.jl` — the Schur reduction with the reference matvec, ~6×
slower. It is not for results; it is an independently written statement of the
same operator, kept for debugging.

`NR=4` by default, and **not your core count**. `:lxy_partition` never cuts z,
and the slab is one element thick in y, so there are 20 element columns and that
is the hard ceiling on ranks. `tools/pick_nranks.jl` lists 1, 2, 4, 5 and 10 as
the counts that leave no rank empty:

| ranks | columns/rank | points/rank | halo/elem | |
|---|---|---|---|---|
| 1 | 20 | 130,005 | 2.1 | |
| 2 | 10 | 65,002 | 2.2 | |
| 4 | 5 | 32,501 | 2.4 | default |
| 5 | 4 | 26,001 | 2.5 | |
| 10 | 2 | 13,000 | 3.0 | communication-bound |

At 10 ranks the halo per element is up 43% on 2 columns of work: the run stops
being about the stage solve, which is the thing being measured. Ten cores does
not mean ten ranks here. That ceiling is a property of a 2D problem, not of this
deck — spend the rest on threads, or on the 3D case.

Other knobs: `DBG_MESH`, `DBG_TEND`, `DBG_DT`, `DBG_RTOL`, and everything the
3D deck takes.

## Schur IMEX vs no-Schur IMEX: what to expect

The three arms differ like this:

| arm | stage solve | unknowns | H is |
|---|---|---|---|
| `full` | five-field | 5·Np | — |
| `schur` | scalar Schur | Np | bespoke scalar sweeps |

**What should carry over from 3D.** The per-iteration savings are a property of
solving for one field instead of five, and the anisotropy that sets them is
identical here. The 3D case, both arms back to back on the same nodes, 25 ranks:

| s/step | five-field | Schur + kernel | |
|---|---|---|---|
| stage solve | 9.542 | 1.880 | 5.08× |
| — matvec | 4.955 | 1.092 | 4.54× |
| — preconditioner | 2.042 | 0.250 | 8.18× |
| — orthogonalise | 1.922 | 0.199 | 9.67× |
| — MPI reduce | 0.357 | 0.053 | 6.69× |
| iterations/step | 69.4 | 36.2 | 1.92× fewer |
| **step** | **10.651** | **2.990** | **3.56×** |

### Measured here: 4.37×

Both arms back to back, 20×1×80, 2 ranks, Δt 0.6, rtol 1e-8, 50 profiled steps:

| s/step | five-field | Schur | |
|---|---|---|---|
| rhs! | 0.326 | 0.323 | untouched |
| f_imp | 0.108 | 0.119 | different splitting |
| **stage solve** | **4.585** | **0.700** | **6.55×** |
| — matvec | 1.897 | 0.384 | 4.94× |
| — preconditioner | 1.330 | 0.154 | 8.62× |
| — orthogonalise | 1.104 | 0.057 | 19.5× |
| — MPI reduce | 0.0164 | 0.0025 | 6.56× |
| iterations/step | 72.6 | 38.1 | **1.91× fewer** |
| **step** | **5.028** | **1.151** | **4.37×** |

**1.91× fewer iterations, against the 3D case's 1.92×.** That is the design
goal met: same h_x/h_z, same CFL_h, so the Krylov behaviour is the same problem
at a sixteenth of the size. This case is a faithful proxy.

**It is FASTER than the 3D case's 3.56×, and an earlier draft of this file
predicted the opposite.** That prediction was wrong, and the reasoning behind it
was backwards: it argued the MPI reduce would be a larger share of a smaller
problem. The reduce does not scale with problem size, it scales with RANK COUNT
— 0.3% of the step on 2 ranks against 3.4% on 25. Fewer ranks means *less*
communication, not relatively more. Expect the gap to narrow, not widen, as you
add ranks here.

What genuinely does not transfer: absolute s/step (130,005 points against
2,106,081), and these numbers were taken on a 4-core shared box at 2 ranks, so
they are contended. Run-to-run variation on the 3D cluster measured ~13% between
two identical runs; expect the same or worse on a desktop, and run each arm
twice before believing a step time. The per-term ratios are 5–20× and survive
it.

`rhs!` and `f_imp` are untouched by any of this and are 38.4% of the Schur step
here (36.5% in 3D). That is the ceiling on anything further.

**Read the profile split, not just the step time.** That is the part that
transfers. `JEXPRESSO_HEVI_PROFILE=1` (which `run_rtb2d.sh` sets) prints the
per-term breakdown; compare `matvec`, `precond`, `orthogonalise` and
`MPI reduce` between arms rather than one headline number.

If a Schur run ever looks wrong rather than slow, `:imex_schur_kernel => false`
runs the same reduction through the reference matvec. A difference there is a
kernel bug; no difference means the reduction itself is what to look at.

## Two things to check before believing a result

**Which arm actually ran.** The setup report prints it from inside the ranks:

```
Stage solve: preconditioned GMRES on the SCALAR SCHUR system
    ...
    matvec: bespoke scalar sweeps (~0.36x one 5-field apply)
```

versus `... on all 5 fields`, and versus `reference form -- two 5-field applies`.

**The two arms are different splittings**, not two solvers for one problem.
`:imex_schur` forces the advective Θ row, without which the elimination does not
close on a single scalar; the two rows differ by 0.06% of the flux form
(`test/hevi/test_theta_advective.jl`). That is the right comparison for wall
clock and the wrong one for reading a state difference between
`output_imex_full/` and `output_imex_schur/` as an error. What *is* pinned: given
the same operator, the two stage solves return the same five fields to ~1e-12
(`test/imex3d/test_schur_stage.jl`), and the two forms of H agree to 1.9e-16
(`test/hevi/test_schur_kernel.jl`).
