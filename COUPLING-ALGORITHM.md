# Coupling Algorithm: Jexpresso -- Alya

This document describes the protocol between Jexpresso (Julia SEM solver) and
Alya (Fortran CFD code), as implemented by the Alya proxy in
[`AlyaProxy/alya_all2all_time_loop.f90`](AlyaProxy/alya_all2all_time_loop.f90)
(built as `AlyaProxy/Alya.x`) and by Jexpresso in
[`src/kernel/coupling/couplingStructs.jl`](src/kernel/coupling/couplingStructs.jl).
It covers only the coupling logic and the data that flows between the two
codes. For how to build and launch a coupled run, see
[RUN-COUPLED.md](RUN-COUPLED.md).

In one sentence: at every time step Jexpresso interpolates its velocity onto
Alya's structured grid and sends each Alya rank the values at the points that
rank owns. Data flows one way, Jexpresso → Alya; Alya sends nothing back.

References below name functions rather than line numbers. Unless another file
is given, Julia functions are in `couplingStructs.jl`.

---

## 1. Launch and handshake

Both codes are launched as a single MPMD MPI job sharing one `MPI_COMM_WORLD`,
**Alya first**:

```bash
export JEXPRESSO_COUPLED=1
mpirun -np 2 ./AlyaProxy/Alya.x : -np 2 julia --project=. ./src/Jexpresso.jl CompEuler 3dAlya
```

(`mpirun` here is the launcher of the MPI both codes use. With MPI.jl on the
MPICH bundled with Julia, `MPICH_jll`, that is the `mpiexec` shipped with it;
`run_coupled.sh` picks the right one either way — see RUN-COUPLED.md §6.)

> **Prerequisite — one MPI for both codes.** Because Alya (Fortran) and
> Jexpresso (Julia) share a single `MPI_COMM_WORLD`, they must be built and run
> against the **same MPI implementation, version, and ABI** — Alya via `mpif90`,
> Jexpresso via `MPI.jl`/`MPIPreferences`. A mismatch makes the handshake below
> hang in `MPI_Init` or corrupt the exchanged data. See
> [RUN-COUPLED.md](RUN-COUPLED.md) for how to identify the MPI on both sides,
> **reconcile them if they differ** (rebind MPI.jl or recompile Alya), and verify
> the match before launching.

Alya ranks occupy world ranks `0 .. NA-1` and Jexpresso ranks
`NA .. NA+NJ-1`. Two facts follow from that layout, and both codes rely on them:

- **Alya's world rank 0 is the root** of every setup broadcast and of the name
  gather. Put Alya first on the `mpirun` line.
- **Alya's rank 0 is a master that owns no grid points.** The grid is split over
  its ranks `1 .. NA-1`, so `NA ≥ 2`. With one Alya rank no point has an owner;
  Jexpresso stops with an error rather than run uncoupled
  (`extract_local_alya_coordinates`).

`JEXPRESSO_COUPLED=1` is what makes Jexpresso split its communicator
(`src/run.jl`). Without it Jexpresso runs standalone on the whole world,
Alya's ranks included, and deadlocks in its first collective. Export it for the
whole job; `-x` after the `:` is not honoured by every launcher.

### 1.1 Communicator split

Each code splits `MPI_COMM_WORLD` by its own color:

- **Alya:** `MPI_Comm_split(MPI_COMM_WORLD, 1, rank, PAR_COMM_FINAL)`
- **Jexpresso:** `MPI.Comm_split(world, APPID, wrank)` in
  `je_init_mpi_and_split_comm`, with `APPID` from the environment (default 2).

Jexpresso then knows it is coupled because its communicator is smaller than the
world (`lsize < wsize`). From here on, every Jexpresso-internal collective must
use that local communicator, `get_mpi_comm()`; only the operations in this
document use the world (§7).

### 1.2 Application names

Every rank sends a 128-character name to world rank 0 with `MPI_Gather`:
`"ALYA"` from Alya, `"JEXPRESSO"` from `je_perform_coupling_handshake`. Alya's
rank 0 prints them under `=== Coupling labels (world size= N ) ===`.

---

## 2. Setup: what is exchanged, in order

Every row is an operation that **both** codes perform, in exactly this order.
One side skipping or reordering a row deadlocks the job.

| # | Operation on `MPI_COMM_WORLD` | Alya (Fortran) | Jexpresso (Julia) | Payload |
|---|---|---|---|---|
| 1 | `Comm_split` | color 1 | color `APPID` = 2 | — |
| 2 | `Gather` to world rank 0 | `app_name` | `je_perform_coupling_handshake` | 128 chars |
| 3 | `Bcast` from world rank 0 | `ndime` | `je_receive_alya_data` | 1 × `Int32` |
| 4 | `Bcast` × 3, for `d = 1:3` | `rem_min(d)`, `rem_max(d)`, `rem_nx(d)` | `je_receive_alya_data` | `Float64`, `Float64`, `Int32` |
| 5 | `Allreduce(SUM)` | `alya_to_world` | `je_receive_alya_data` (zeros) | `NA` × `Int32` |
| 6 | `Barrier` | before the `Alltoall` | `setup_coupling_and_mesh` or `je_early_coupling_sync!` | — |
| 7 | `Alltoall` | sends zeros | points it will send to each world rank | 1 × `Int32` per rank |
| 8 | point-to-point, tag 0 | `MPI_Recv` from each Jexpresso rank with points | `je_send_node_list` | the global IDs of those points, `Int32` |

Nothing else is exchanged before the time loop. In particular **neither side
sends `neqs`, the number of steps, `Δt` or `tend`**; see §4.

- `ndime` is the spatial dimension (the proxy sends 3).
- `rem_min`, `rem_max`, `rem_nx` describe Alya's structured grid: bounding box
  and number of nodes per direction. All three components are always sent; a 2D
  grid has `rem_nx(3) = 1`.
- `alya_to_world[a]` is the world rank of Alya rank `a`. Jexpresso uses it to
  address messages and to tell Alya's workers from its master (world rank 0).
- After row 7, Alya's `npoin_recv(j)` is the number of points world rank `j`
  will send it every step.
- The ID list of row 8 tells Alya where each value it will receive belongs.
  Alya receives it once and uses it at every step.

Proxy values: `rem_min = [-5000, -3000, 0]`, `rem_max = [5000, 1500, 10000]`,
`rem_nx = [10, 10, 10]`. These are the bounds of the mesh the `3dAlya` case
reads, `hexa_TFI_10x1x10.msh`, so all 1000 points lie inside Jexpresso's domain.

---

## 3. Which points each Jexpresso rank sends

`extract_local_alya_coordinates` decides, for every Jexpresso rank, which Alya
grid points it is responsible for. Alya sends no coordinates: Jexpresso rebuilds
the grid from the metadata of §2.

**Coordinates and IDs.** For 0-based structured indices `(i1, i2, i3)`:

```
x = rem_min[1] + i1*dx,   dx = (rem_max[1] - rem_min[1]) / (rem_nx[1] - 1)   (same for y, z)
id = i1 + rem_nx[1]*(i2 + rem_nx[2]*i3) + 1                                   (1-based)
```

**Selection.** Each rank keeps the points that lie in its part of the mesh: in
its bounding box in 3D, in one of its elements in 2D (after index-space cropping
and block-wise tests). Points on a shared face are claimed by several ranks; an
`Allreduce(MAX)` on the Jexpresso communicator keeps exactly one claimant. A
second pass gives the points that no rank claimed to the lowest rank whose
bounding box contains them. `verify_coupling_communication_pattern` prints the
result: `[VERIFY] … CHECK 2` must report `All Alya points accounted for`, or
some Alya points lie outside Jexpresso's domain and will receive nothing.

**Owner.** Alya assigns points to its workers `1 .. NA-1` in contiguous chunks
of the ID range, the first `mod(nmax, NA-1)` workers taking one extra point.
Jexpresso applies the same rule to find the Alya rank that receives each point.

**Order.** Points are sorted by `(owner, id)`. The ID list of §2 row 8 and every
later field message use this same order, which is how Alya matches values to
points.

**When.** Normally inside `setup_coupling_and_mesh`, after `sem_setup`. If a
mesh cache from an earlier run exists, it is done before `with_mpi` instead
(`je_prefetch_caches!` → `_je_prefetch_geometry!`), and rows 6–8 then run early
too (`je_early_coupling_sync!`), so that Alya is released from its barrier while
Jexpresso is still compiling. A case that rescales its grid (`:xscale`,
`:xdisp`, `:yscale`, `:ydisp`) never takes the early path, because the cache
holds the grid before rescaling.

---

## 4. Time loop: one exchange per step

### 4.1 Alya

The proxy computes its step count locally from hard-coded values:
`t0 = 0`, `dt = 0.5`, `tend = 1000`, so `nsteps = int((tend - t0)/dt) = 2000`.
Its rank 0 prints it at startup (`Steps: 2000`). At every step it:

1. posts one `MPI_Irecv` per Jexpresso rank with `npoin_recv > 0`:
   `npoin_recv * nfields` doubles, tag 0, where `nfields = ndime`;
2. waits for all of them (`MPI_Waitall`);
3. scatters the values into its own chunk of the grid using the ID list;
4. every `out_dt = 100` time units, gathers the grid on its rank 0 and writes
   `alya_grid_NNNNNN.vts`, with the received fields as `var1`, `var2`, `var3`.

### 4.2 Jexpresso

`setup_coupling_callback` registers a `DiscreteCallback` that fires after every
accepted step with `t > tinit + :couple_time_tol`. Each call runs
`je_perform_coupling_exchange_3d` (or `je_perform_coupling_exchange` in 2D):

1. **Output variables.** The state is converted with the case's `user_uout!`.
   For `3dAlya` that is `(ρ, u, v, w, θ)`.
2. **Interpolation.** Each local Alya point is located in the SEM mesh (element
   bins, bounding-box test, Newton solve for the reference coordinates) and the
   tensor-product Lagrange interpolant is evaluated there. Locating the points
   depends only on the two grids, so on a static mesh it is done once, at the
   first exchange (`[coupling] interpolation cache built: …`); every later step
   is one dot product per point and variable. The cache is off under `:lamr` or
   `:ladapt`, or with `:lcouple_cache_interp => false`. A point that no element
   contains takes its nearest node's value.
3. **Packing.** Columns `2 : neqs-1`, the velocity, are packed point by point
   (`[u1, v1, w1, u2, v2, w2, …]`) into one buffer per destination Alya rank
   (`pack_velocity_data!`). `neqs - 2 == ndime` is asserted at setup, because
   that is the `nfields` Alya expects.
4. **Send.** One `MPI.Isend` per destination, tag 0, then `MPI.Waitall`
   (`coupling_exchange_data!`).

With `SEND_COORDS = true` (top of `couplingStructs.jl`) Jexpresso sends the
interpolated mesh coordinates instead of the velocity: same buffer shape, and
Alya's VTS files should then reproduce its own grid coordinates. It is a check
of the point location and interpolation.

### 4.3 The step contract

The only thing that keeps the two loops in step is that **Jexpresso performs
exactly `nsteps` exchanges**. Neither side tells the other how many steps it
takes:

- Fewer (a larger `:Δt`, a shorter `:tend`): Alya waits forever in
  `MPI_Waitall` while Jexpresso waits in the final barrier.
- More: Alya never receives the extra messages. While they are small enough to
  be sent eagerly they are silently lost, and if the extra exchange comes first
  every step reaches Alya one step late (§6). Once a message is too large to be
  sent eagerly, Jexpresso blocks in `MPI_Waitall` while Alya waits in the final
  barrier.

So `(:tend - :tinit) / :Δt` in the case's `user_inputs.jl` must equal Alya's
`nsteps`. Change both sides together: `t0`/`dt`/`tend` in
`alya_all2all_time_loop.f90`, then rebuild `Alya.x`. Two more things add steps
on the Jexpresso side and must be avoided:

- **Output times off the step grid.** Every entry of `:diagnostics_at_times`
  must be `tinit` plus a multiple of `:Δt`; otherwise the integrator shortens a
  step to land on it, and that is one more exchange.
- **Adaptive time stepping** (`:ode_adaptive_solver => true`).

Jexpresso's integrator warm-up in `time_loop!` runs one throw-away step with the
real callback set. The exchange is suspended for that step
(`CouplingData.exchange_enabled`), so only the real steps exchange. At the end
of the solve Jexpresso prints

```
 # Coupling: 2000 exchanges sent to Alya
```

which must equal Alya's `Steps:` line.

---

## 5. Shutdown

After its loop Alya prints `Alya: time loop complete, syncing with Julia...` and
waits in `MPI_Barrier(MPI_COMM_WORLD)`. Jexpresso issues the matching
`MPI.Barrier(world)` at the end of its `with_mpi` block in `src/run.jl`, after
the solve and its output. Then both call `MPI_Finalize`.

---

## 6. How messages are matched

Every point-to-point message uses tag 0 on `MPI_COMM_WORLD`. They are told
apart by order. MPI never lets a message overtake an earlier one from the same
sender with the same tag and communicator. From each Jexpresso rank, Alya
therefore receives the ID list first (the blocking `MPI_Recv` of §2 row 8), then
exactly one field message per step, in step order. That is also why an extra or
missing exchange shifts every later message rather than failing.

---

## 7. Rules for code that runs in a coupled job

- **Never run a Jexpresso-internal collective on `MPI.COMM_WORLD`.** Use
  `get_mpi_comm()`: it is Jexpresso's own communicator under coupling and
  `COMM_WORLD` standalone. A collective on `COMM_WORLD` waits for Alya's ranks,
  which never arrive. This typically shows up only at the first diagnostic
  output, and never in a standalone test.
- **Take all of the Jexpresso ranks, or none, into a world operation.** The
  prefetch and early-sync paths first agree across ranks (`_je_all_ranks`) on
  whether to enter one.
- **Anything that adds or removes steps changes the exchange count** (§4.3).

---

## 8. Code map

| Stage | Jexpresso | Alya proxy (`alya_all2all_time_loop.f90`) |
|---|---|---|
| Split, coupled or not | `src/run.jl` → `je_init_mpi_and_split_comm` | `MPI_Comm_split` |
| Names | `je_perform_coupling_handshake` | STEP 0: HANDSHAKE |
| Grid metadata, rank map | `je_receive_alya_data` | STEP 2: GRID METADATA; Alya -> World rank map |
| Points and owners | `extract_local_alya_coordinates` | point range of each worker (`i_start`, `i_end`) |
| Counts, ID list | `setup_coupling_and_mesh` or `je_early_coupling_sync!`; `je_send_node_list` | STEP 3: COUNT EXCHANGE; STEP 3b: RECEIVE … NODE LIST |
| Mesh + coupling object | `problems/drivers.jl` → `setup_coupling_and_mesh` | — |
| Per-step exchange | `setup_coupling_callback` → `je_perform_coupling_exchange_3d` → `coupling_exchange_data!` | TIME LOOP: `MPI_Irecv` + `MPI_Waitall` |
| Callback registration, warm-up | `time_loop!` in `src/kernel/solvers/TimeIntegrators.jl` | — |
| Output | Jexpresso's own VTK | `write_alya_grid_vts` → `alya_grid_NNNNNN.vts` |
| Shutdown | `MPI.Barrier(world)` in `src/run.jl` | `MPI_Barrier(MPI_COMM_WORLD)` |
