# MPI Parallel Scaling Audit

Static reading of the MPI parallelization on branch
`claude/sharp-thompson-OmWOs` after the wholesale infra swap from
`origin/hw/giga_les`. **This is a code review, not a benchmark.** It
identifies anti-patterns and synchronization points that would prevent
good strong/weak scaling; it does *not* prove scaling. Real scaling
claims require a measured strong/weak sweep (1, 2, 4, … ranks at fixed
problem size, wall-time per timestep, plot efficiency).

## Punch list

| Severity | Item | Location |
|---|---|---|
| **BROKEN** | xy-bin partitioner: prime `nparts` (7, 11, 13, …) collapses to 1-D strip; non-column-extruded 3-D meshes get layer-multiplied imbalance; uniform bins on graded meshes pile cells into a few bins. No Metis/ParMETIS fallback. | `src/kernel/mesh/mesh.jl:21-41` |
| **BROKEN** | NetCDF writer `MPI.gather`s full per-rank fields to rank 0. | `src/io/write_output.jl:708-712` |
| **BROKEN** | `find_gip_owner`: `Gatherv` all `ip2gip` to rank 0, dict-build on rank 0, `Scatterv` back. Init/AMR only, but O(global mesh) on rank 0 → memory wall on ≥10 M-cell meshes. | `src/kernel/mesh/mesh.jl:1829-1879` |
| **BROKEN** | Gmsh read on rank 0 holds the full mesh in memory before partitioning. Mesh-size ceiling regardless of rank count. | `src/kernel/mesh/mesh.jl:65,70` |
| BOTTLENECK > ~512 ranks | `computeCFL` does 4–7 sequential `Allreduce`s per print (umax, vmax, wmax, μmax, plus 1–2 in `soundSpeed`). Each a separate latency round-trip — trivially batchable into one Allreduce of a small vector. | `src/kernel/physics/soundSpeed.jl:111-172, 183` |
| BOTTLENECK > ~1024 ranks (init only) | `restructure_for_periodicity` has ~12 root-funneled Bcast/Gather/Scatterv rounds. Heavy but one-shot. | `src/kernel/mesh/restructure_for_periodicity.jl` |
| **OK** | DSS halo exchange: pre-cached `Isend`/`Irecv!` on neighbour-only `active_send/recv_ranks`, `Waitall`. No collectives in the hot path. | `src/kernel/mpi/mpi_communications.jl:237-399`, `src/kernel/operators/rhs.jl:473, 609, 630, 649` |
| **OK** | `rhs!` itself: **zero** Allreduce / Barrier / Bcast / Gather per timestep. | `src/kernel/operators/rhs.jl:106` |
| **OK** | Parallel VTK writer (`pvtk_grid`). | `src/io/write_output.jl:296-314` |
| **OK** | Boundary conditions: no MPI calls in BC code; cross-rank consistency handled by the subsequent DSS. | `src/kernel/boundaryconditions/*` |

## Detail by section

### 1. Partition quality — `src/kernel/mesh/mesh.jl:21-41`

**BROKEN at moderate scale (worst-case factor ~2× imbalance; total
failure on prime `nparts`).**

- Algorithm: bin cell centroids into a rectangular `nx × ny` grid in the
  xy plane only.
- `divisors = [d for d in 1:nparts if nparts % d == 0]` → for prime
  `nparts` (7, 11, 13, 17, 23, …) the only divisors are `1` and
  `nparts`, forcing a 1×N or N×1 strip decomposition regardless of mesh
  aspect ratio.
- 3-D meshes that are *not* column-extruded collapse multiple z-layers
  into one bin: cell count per bin = (#xy cells in bin) × (#layers).
  xy distribution variations are multiplied by layer count. There is no
  z partitioning available.
- Bins are uniform Cartesian over the bounding box (`mesh.jl:37-38`),
  not equipopulated. A graded mesh (refined near boundaries, sponge
  layers, terrain) puts most cells in a few bins → severe imbalance.
- No fallback to Metis/ParMETIS anywhere (`grep Metis` empty).
- A diagnostic exists (`measure_elements_per_rank`, `mesh.jl:1775`) but
  it only *prints*; `needs_redistribution` is returned but unused on the
  `lxy_partition` path. Watch its min/avg/max ratio at runtime to
  quantify the imbalance.

### 2. DSS / halo exchange — `src/kernel/mpi/mpi_communications.jl`, `src/kernel/operators/rhs.jl`

**OK for thousands of ranks.**

- `DSS_global_RHS!` (`element_matrices.jl:1119`) →
  `assemble_mpi!` (`mpi_communications.jl:237-315`) uses cached
  non-blocking pairs: `MPI.Isend` / `MPI.Irecv!` on
  `active_send_ranks` / `active_recv_ranks` only, followed by
  `MPI.Waitall`. No collectives in the hot path.
- Buffers, request handles, and the active-rank lists are pre-allocated
  in `setup_assembler` (called once during init). `MPI.Alltoall` and
  two `MPI.Barrier`s at setup (`mpi_communications.jl:116-117, 150`)
  are init-only.
- Both forward-gather and scatter-back phases are point-to-point —
  the correct neighbour-only DSS pattern.
- Per RHS evaluation: up to **3 DSS rounds** (`rhs.jl:473, 609, 630,
  649`) plus an `inputs[:ladapt]` gather/scatter pair (`rhs.jl:600,
  662`). For an RK method with `s` stages, that is `~3s` neighbour
  exchanges per timestep.
- **Cannot determine from static reading**: actual neighbour count per
  rank for `assemble_mpi!` — depends on the xy partitioner output for
  the specific mesh. If the strip-decomposition path is hit (prime
  `nparts`, see §1), neighbour count is O(1) which is good, but DSS
  volume per rank scales as the long-side ghost layer and won't shrink
  with `nparts`. Recommend instrumenting
  `length(cache.active_send_ranks)` and `sum(cache.send_data_sizes)`
  at runtime.

### 3. Time-loop synchronization — `src/kernel/solvers/TimeIntegrators.jl`

**Per step (always): zero collectives.** `rhs!` (`rhs.jl:106`) is
collective-free.

**Per output (`dosetimes` cadence, `affect!`
`TimeIntegrators.jl:136-173`):**

- 1 `Allreduce(ad_lvl_max)` (`TimeIntegrators.jl:145`).
- Inside `computeCFL` (`physics/soundSpeed.jl:111, 114, 117, 120,
  157-172, 183`): **4–7 sequential `MPI.Allreduce` calls per CFL
  print** (umax, vmax, [wmax], μmax, plus 1–2 inside `soundSpeed`).
  Each is a separate latency-bound round trip.
  **BOTTLENECK at > ~512 ranks** if outputs are frequent. Trivially
  batchable into one `Allreduce` of a small vector.
- `write_output` (VTK path → `pvtk_grid`, `write_output.jl:296-314`)
  is properly parallel per-rank → **OK**.
- NetCDF path (`write_output.jl:702-787`) uses **`MPI.gather` of full
  per-rank arrays to rank 0**, then rank 0 writes.
  **BROKEN beyond ~32 ranks for large meshes** — both memory
  (rank 0 holds the full global field) and serial write time. HDF5 path
  (`write_output.jl:579-602`) writes one file per rank per variable —
  scales, but produces `N × nvar` files.
- Restart (`do_restart!` `TimeIntegrators.jl:58`): 2 `MPI.Barrier`,
  `cp` and `rm` only on rank 0 (lines 79-82) — serial filesystem ops.
  OK if rare, expensive if `restart_time` is small.
- AMR refinement (`TimeIntegrators.jl:243-258`): `Allreduce` + `solve`
  restart per AMR step plus 2 `Barrier`s. Per AMR event, not per step.

### 4. Serial chokepoints

- **`src/io/write_output.jl:708-712`** — `MPI.gather` of `ip2gip`,
  `x`, `y`, subelem, qout to rank 0 in the NetCDF writer.
  **BROKEN at scale.**
- **`src/kernel/mesh/mesh.jl:1829-1879` `find_gip_owner`** —
  `Gatherv` of all `ip2gip` to rank 0, dictionary build on rank 0,
  `Scatterv` back. Called from `mesh.jl:1090, 1315, 1545-46, 1639-40`
  and `Projection.jl:1094`. **Init/AMR only**, but on very large
  meshes rank 0 memory and dict ops will be a wall (> ~1024 ranks).
- **`src/io/les_statistics.jl:71, 307`** — `MPI.gather`s on rank 0
  during cache build (init) and final write. Init-only, fine.
- **`src/kernel/mesh/restructure_for_periodicity.jl`** — 91 MPI calls;
  12+ `MPI.Bcast` / `Gather` / `Scatterv` rounds (lines 1304, 1376-79,
  1422, 1570-72, 1615, 1685-87, 1725, 1771-72, 1991-2016, 2589-94,
  2716-21) all root-funnelled. Init-only but heavy.
- **`conformity4ncf_q!`** (`Projection.jl:3001+`, called from
  `params_setup.jl:257, 267`): uses the same neighbour-only DSS path
  → **parallelizes correctly**. Init-only.

### 5. Initialization

- Gmsh read: `mesh.jl:65, 70` — `@outputrootonly` reads on rank 0 then
  partitions (`DiscreteModel(parts, smodel, cell_to_part)` line 68).
  Rank 0 holds the full mesh in memory → **BROKEN for meshes > ~10 M
  cells regardless of rank count.**
- `find_gip_owner` and `restructure_for_periodicity` are O(global mesh)
  on rank 0. Expect minutes at > 1024 ranks for > 10 M-DOF meshes.
- `Axb_rad_mpi.jl:29, 32, 214, 234, 263` and
  `build_rad_3d.jl:285, 921-930` use `Allreduce!` over global vectors
  (`y_local`, `b_global`). Radiation-only; if `radiation_time_step`
  fires often this becomes **per-timestep BROKEN**.

### 6. Boundary conditions

**OK.** `grep MPI src/kernel/boundaryconditions/` returns nothing.
`apply_boundary_conditions_*!` (`rhs.jl:543, 633`) operate on local
`poin_in_bdy_*` arrays only; cross-rank consistency is handled by the
subsequent `DSS_global_RHS!` (`rhs.jl:649`). No serial reduction in BC
code.

## Bottom line

The per-timestep core is genuinely good — pure neighbour-only
point-to-point DSS, zero collectives in `rhs!`. That is the right
architecture for thousands of ranks. The failure modes are at the
edges:

1. **Partitioner quality** is the most likely thing to bite first. The
   xy-bin partitioner is "good enough for benchmarks", not a serious
   partitioner. Plug in Metis/ParMETIS before claiming scaling
   results on irregular meshes.
2. **Output writers**: the NetCDF path will OOM rank 0 at moderate
   scale. Use VTK (parallel) or per-rank HDF5 for scaling runs.
3. **Mesh size** is bounded by rank-0 memory because the gmsh read is
   serial. Below ~10 M cells you are fine.
4. **CFL prints** are easy to fix (batch the Allreduces). Worth doing
   if you print every step at > 512 ranks.

## Suggested low-effort fixes

- **CFL Allreduce batching** — DONE in commit `23d7c82d`. Collapsed
  2D `[u, v]` → 1 Allreduce and 3D `[u, v, w, μmax]` → 1 Allreduce;
  also folded local `cmax` from `soundSpeed` into the same reduction
  to fix a pre-existing per-rank-local acoustic-CFL bug.
- **Runtime imbalance diagnostic** — already wired:
  `measure_elements_per_rank` (`mesh.jl:1775`) is called from
  `mod_mesh_read_gmsh!` (`mesh.jl:342`) and prints min/avg/max/imbalance
  ratio on rank 0 during init. What is missing is the *follow-through* —
  the returned `lneed_redistribute` flag is set but no code path acts
  on it. Acting on it is a real feature, not a small fix.
- **NetCDF parallel I/O** — give the NetCDF path the same per-rank-write
  treatment the VTK path already has. Bigger change; parallel
  HDF5/NetCDF on macOS arm64 has its own portability issues. Defer
  until someone actually needs the NetCDF format at scale.

## Follow-up: ParMETIS integration plan

The xy-bin partitioner (`_compute_xy_partition`, `mesh.jl:21-41`) is
the most fragile piece of the parallel pipeline (see §1 of the punch
list above). It works for column-extruded LES boxes on composite
`nparts`; it falls over on prime `nparts`, on graded meshes, on
arbitrary 3-D geometries, and on anything where the xy projection is
not representative of cell density. Plugging in METIS / ParMETIS would
fix all of these.

This is **not a small edit** and the right way to make it happen is
to land it as its own PR after the design questions below are
answered. Filing the plan here so it does not get lost.

### Open design questions (answer before coding)

1. **Which Julia binding?** Options as of writing:
   - `Metis.jl` (mature, well-maintained) — wraps the **serial** METIS
     library. Adequate if rank 0 partitions and broadcasts; not adequate
     at the very scale where partitioning quality matters most.
   - `ParMETIS_jll` — Yggdrasil ships the C library, but I am not aware
     of a high-level Julia wrapper. May need to write the ccalls
     directly, which is real work.
   - Native partitioning hooks in `GridapDistributed` /
     `PartitionedArrays` (already deps in `Project.toml`). Worth a
     careful read before reinventing.
   - **Decision needed**: do we accept the serial-METIS limitation
     (partition on rank 0, broadcast assignment), or commit to
     ParMETIS_jll and the wrapper work?

2. **Required vs optional dependency?**
   - As **required**: simpler, but inherits the binding's platform
     support. `ParMETIS_jll` may not build cleanly on macOS arm64
     (this needs verification on the target dev machine before
     committing).
   - As **optional** (Project.toml `[extras]` / package extension):
     keeps the build portable but doubles the partitioner code paths.
   - **Recommendation**: optional, with the existing xy-bin partitioner
     as the universal fallback. New input flag selects:
     ```julia
     :partitioner => :parmetis | :metis | :xy | :gmsh   # default :xy
     ```

3. **Adjacency-graph source.** ParMETIS needs a CSR-format cell graph
   (each cell + neighbours via shared face). Gridap should expose this
   via `get_cell_faces` / `get_face_cells`, but the exact wiring needs
   verification. For 3-D high-order spectral elements, neighbour
   relationship is via faces (not vertices), so use the **3-D-face**
   topology, not vertex topology.

4. **Interaction with p4est / AMR.** When `inputs[:lamr] == true` the
   mesh is an `OctreeDistributedDiscreteModel` and p4est handles
   partitioning internally on refine/coarsen. A METIS partition would
   only apply to the initial coarse mesh; p4est takes over afterwards.
   Confirm this does not double-partition or invalidate p4est's
   internal connectivity.

5. **Interaction with periodic restructuring.**
   `restructure_for_periodicity.jl` already does heavy rank-0 work
   (~12 root-funnelled rounds, audit §4). It assumes a specific
   ownership pattern. Need to check whether METIS-assigned ownership
   breaks any of its invariants.

### Concrete sequencing once the questions are answered

1. **Spike** (~half day): in a scratch script, partition a 10k-cell
   3-D theta mesh with `Metis.jl` on rank 0, broadcast the assignment,
   construct `DiscreteModel(parts, smodel, cell_to_part)` exactly as
   `_compute_xy_partition` does today. Compare min/avg/max/imbalance
   against the xy-bin partition. If imbalance is comparable or worse,
   stop — the win is too small to justify the dependency.

2. **Hook up the input flag** (~1 day): add `:partitioner` default
   in `mod_inputs.jl`, branch in `mod_mesh_read_gmsh!` between
   `_compute_xy_partition` / `_compute_metis_partition` /
   `_compute_parmetis_partition` / gmsh-default based on the flag.
   Keep `lxy_partition` for backwards compatibility (`:xy` when set).

3. **Verify partition quality** (~1 day): run the existing
   `measure_elements_per_rank` diagnostic on a representative sweep
   (theta 2-D, theta 3-D, BoMex, a graded-mesh case if available) at
   2/4/8/16/32 ranks. Record imbalance ratios before / after. Reject
   any partitioner that does not strictly improve over xy-bin on
   non-column-extruded meshes.

4. **Measure DSS communication volume** (~1 day): instrument
   `assemble_mpi!` to log `length(cache.active_send_ranks)` and
   `sum(cache.send_data_sizes)` per rank at startup. Compare the
   partitioners. The ratio of communication-volume reduction to
   imbalance-reduction is what actually matters for scaling, not
   imbalance alone.

5. **ParMETIS only if METIS-on-rank-0 hits a memory wall** (~3 days
   if needed): if step 3 shows METIS wins on quality but step 4
   shows the rank-0 partition phase is too expensive at the target
   mesh size, then escalate to ParMETIS_jll and write the ccall
   wrapper. Until that data exists, **do not** start on ParMETIS —
   it is the higher-cost option and may not be needed.

### Risks / unknowns

- `ParMETIS_jll` build status on macOS arm64 not verified.
- Gridap's cell-adjacency exposure may not be CSR-shaped; conversion
  cost is O(nelem) but adds complexity.
- p4est-managed AMR paths may render the partition irrelevant for
  AMR runs (METIS only helps the initial coarse mesh).
- `restructure_for_periodicity` rank-0 invariants need verification
  against arbitrary partitions.

### Rough effort estimate

- Steps 1–4 (METIS-on-rank-0 path with quality + comm measurements):
  **3–4 days** of focused work.
- Step 5 (ParMETIS escalation if needed): **+3 days**, with platform
  risk on macOS arm64.

Total: **3–7 days** depending on whether ParMETIS is needed.
