# Jexpresso-Alya Coupling: Data Exchange Reference

This document describes the MPI coupling protocol between **Jexpresso** (Julia, Spectral Element Method solver) and **Alya** (Fortran, CFD proxy/unit-test program). It identifies every quantity exchanged, the direction of communication, the MPI calls used, and the exact source-code locations.

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [MPI Communicators](#2-mpi-communicators)
3. [Initialization Phase: Alya Broadcasts to All Ranks](#3-initialization-phase-alya-broadcasts-to-all-ranks)
4. [Rank Mapping: Alya-to-World](#4-rank-mapping-alya-to-world)
5. [Communication Pattern Exchange (Alltoall)](#5-communication-pattern-exchange-alltoall)
6. [Time Loop: Jexpresso Sends to Alya](#6-time-loop-jexpresso-sends-to-alya)
7. [Time Loop: Alya Sends to Jexpresso](#7-time-loop-alya-sends-to-jexpresso)
8. [VTS Output: Alya Gatherv and File Write](#8-vts-output-alya-gatherv-and-file-write)
9. [Complete Data Flow Summary](#9-complete-data-flow-summary)
10. [File Reference](#10-file-reference)

---

## 1. Architecture Overview

The coupled simulation runs as a single `mpirun` job with two SPMD programs sharing `MPI_COMM_WORLD`:

```
mpirun -np <N_alya> ./AlyaProxy/Alya_enhanced.x : -np <N_jexpresso> julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
```

Alya ranks occupy world ranks `0 .. N_alya-1`, Jexpresso ranks occupy `N_alya .. N_alya+N_jexpresso-1`. Each code splits `MPI_COMM_WORLD` into a local communicator (`PAR_COMM_FINAL` in Fortran, `local_comm` in Julia) using `MPI_Comm_split` with distinct color values.

The coupling has three phases:
1. **Initialization** -- Alya broadcasts grid metadata; all ranks exchange rank mappings and point counts.
2. **Time loop** -- At every timestep, Jexpresso interpolates its SEM solution onto Alya's structured grid coordinates, and sends both the interpolated field values **and** the corresponding grid coordinates to the appropriate Alya ranks. Alya sends a return field back to Jexpresso.
3. **Output** -- Alya gathers the received data across its own ranks and writes VTS files at scheduled intervals.

---

## 2. MPI Communicators

| Communicator | Fortran name | Julia name | Scope |
|---|---|---|---|
| World | `MPI_COMM_WORLD` | `world` / `cpg.comm_world` | All ranks from both codes |
| Alya local | `PAR_COMM_FINAL` | -- | Alya ranks only |
| Jexpresso local | -- | `local_comm` / `get_mpi_comm()` | Jexpresso ranks only |

**Code references:**
- Fortran split: `alya_all2all_time_loop.f90:78`
- Julia split: `Jexpresso-mini-coupled.jl:31`

---

## 3. Initialization Phase: Alya Broadcasts to All Ranks

During initialization, Alya rank 0 broadcasts structured grid metadata to **all** ranks in `MPI_COMM_WORLD`. These broadcasts must occur in the exact same order on both sides.

### 3.1 Handshake (MPI_Gather)

Each rank sends a 128-character application name to world rank 0.

| Quantity | Value | Direction | MPI call |
|---|---|---|---|
| `app_name` (128 chars) | `"ALYA"` or `"JEXPRESSO"` | All ranks -> rank 0 | `MPI_Gather` |

**Code references:**
- Fortran: `alya_all2all_time_loop.f90:85-87`
- Julia: `Jexpresso-mini-coupled.jl:42`

### 3.2 Grid Metadata Broadcasts

The following quantities are broadcast sequentially from Alya rank 0 (world rank 0) to all ranks. The order is critical -- both sides must match exactly.

| # | Quantity | Fortran variable | Julia variable | Type | Count | Description |
|---|---|---|---|---|---|---|
| 1 | Spatial dimension | `ndime` | `ndime_buf` | `MPI_INTEGER` / `Int32` | 1 | 2 or 3 |
| 2a | Grid minimum (x) | `rem_min(1)` | `rem_min[1]` | `MPI_REAL` / `Float32` | 1 | Min x-coordinate of Alya grid |
| 2b | Grid maximum (x) | `rem_max(1)` | `rem_max[1]` | `MPI_REAL` / `Float32` | 1 | Max x-coordinate of Alya grid |
| 2c | Grid points (x) | `rem_nx(1)` | `rem_nx[1]` | `MPI_INTEGER` / `Int32` | 1 | Number of grid points in x |
| 3a | Grid minimum (y) | `rem_min(2)` | `rem_min[2]` | `MPI_REAL` / `Float32` | 1 | Min y-coordinate |
| 3b | Grid maximum (y) | `rem_max(2)` | `rem_max[2]` | `MPI_REAL` / `Float32` | 1 | Max y-coordinate |
| 3c | Grid points (y) | `rem_nx(2)` | `rem_nx[2]` | `MPI_INTEGER` / `Int32` | 1 | Number of grid points in y |
| 4a | Grid minimum (z) | `rem_min(3)` | `rem_min[3]` | `MPI_REAL` / `Float32` | 1 | Min z-coordinate |
| 4b | Grid maximum (z) | `rem_max(3)` | `rem_max[3]` | `MPI_REAL` / `Float32` | 1 | Max z-coordinate |
| 4c | Grid points (z) | `rem_nx(3)` | `rem_nx[3]` | `MPI_INTEGER` / `Int32` | 1 | Number of grid points in z |
| 5 | Number of equations | `neqs` | `neqs_buf` | `MPI_INTEGER` / `Int32` | 1 | Equations per grid point (e.g. 4) |
| 6 | Number of timesteps | `nsteps` | `nsteps_buf` | `MPI_INTEGER` / `Int32` | 1 | Total number of time steps |

Items 2-4 are broadcast inside a loop over dimensions (idime = 1, 2, 3), with one `MPI_Bcast` per scalar per dimension.

**Code references:**
- Fortran broadcasts: `alya_all2all_time_loop.f90:124-135`
  - `ndime`: line 124
  - `rem_min/rem_max/rem_nx` loop: lines 126-130
  - `neqs`: line 133
  - `nsteps`: line 135
- Julia receives: `couplingStructs.jl:857-879` (`je_receive_alya_data`)
  - Also duplicated in `Jexpresso-mini-coupled.jl:68-80` for the entry-point path

### 3.3 Quantities Defined by Alya (Hardcoded in Fortran)

These values are set in the Fortran program and define the Alya grid:

```fortran
! alya_all2all_time_loop.f90, lines 102-113
t0     = 0.0d0
dt     = 0.25d0
tend   = 1000.0d0
nsteps = int((tend - t0) / dt)  ! = 4000

rem_min = [-5000.0,     0.0, 0.0]
rem_max = [ 5000.0, 10000.0, 0.0]
rem_nx  = [100,      100,      1]
ndime   = 2
neqs    = 4
```

The four equation variables correspond to:
1. `rho` -- density
2. `u` -- x-velocity
3. `v` -- y-velocity
4. `theta` -- potential temperature

---

## 4. Rank Mapping: Alya-to-World

After the grid metadata broadcasts, all ranks collectively build a mapping from Alya's local rank indices to world rank indices using `MPI_Allreduce` with `MPI_SUM`.

| Quantity | Type | Direction | MPI call |
|---|---|---|---|
| `alya_to_world[0:asize-1]` | `Int32` array | All ranks (collective) | `MPI_Allreduce(MPI_SUM)` |

Each Alya rank sets its own entry to its world rank; all other entries are zero. The sum produces the complete mapping.

**Code references:**
- Fortran: `alya_all2all_time_loop.f90:140-145`
- Julia: `Jexpresso-mini-coupled.jl:86-87` and `couplingStructs.jl:882-887`

---

## 5. Communication Pattern Exchange (Alltoall)

Before the time loop, all ranks exchange point-count arrays via `MPI_Alltoall` to establish how many points each rank will send to / receive from every other rank.

| Quantity | Type | Direction | MPI call |
|---|---|---|---|
| `npoin_send[0:size-1]` / `npoin_recv[0:size-1]` | `Int32` array | All ranks (collective) | `MPI_Alltoall` |

- **Jexpresso side**: Each Jexpresso rank fills `npoin_send[alya_world_rank]` with the count of Alya grid points it will service (computed by `extract_local_alya_coordinates`).
- **Alya side**: Initializes `npoin_send` to all zeros (Alya does not initiate point distribution), then receives counts via `Alltoall` into `npoin_recv`.

After this exchange, each rank knows exactly how many points it communicates with every other rank.

**Code references:**
- Fortran: `alya_all2all_time_loop.f90:158-166`
- Julia: `couplingStructs.jl:120-123` (inside `setup_coupling_and_mesh`)

---

## 6. Time Loop: Jexpresso Sends to Alya

At every timestep, Jexpresso interpolates its SEM solution onto Alya's structured grid points and sends **two** arrays to each Alya rank:

### 6.1 Interpolated Field Values (TAG_DATA = 2000)

| Quantity | Description | Size per Alya rank | Data type | MPI tag |
|---|---|---|---|---|
| Interpolated solution | `neqs` values per grid point (rho, u, v, theta) | `npoin_send[dest] * neqs` doubles | `Float64` | `TAG_DATA + dest_world_rank` |

The interpolation pipeline:
1. Reshape Jexpresso solution vector `u` into matrix `[npoin x neqs]`
2. Call `interpolate_solution_to_alya_coords()` -- for each Alya grid point in this rank's domain, find the containing SEM element, map to reference coordinates via Newton iteration, evaluate Lagrange basis functions, and interpolate all `neqs` solution components.
3. Pack into per-destination send buffers via `pack_interpolated_data!()`

**Code references:**
- Interpolation: `couplingAuxiliaryFunctions.jl:619-704` (`interpolate_solution_to_alya_coords`)
- Packing: `couplingAuxiliaryFunctions.jl:465-508` (`pack_interpolated_data!`)
- MPI Isend (data): `couplingAuxiliaryFunctions.jl:377-380`

### 6.2 Grid Coordinates (TAG_COORD = 3000)

| Quantity | Description | Size per Alya rank | Data type | MPI tag |
|---|---|---|---|---|
| Alya grid coordinates | `ndime` coordinate components per point (x, y [, z]) | `npoin_send[dest] * ndime` doubles | `Float64` | `TAG_COORD + dest_world_rank` |

The coordinates are sent alongside the data so that Alya can correctly place each interpolated value at its grid position, regardless of the order in which points arrive from different Jexpresso ranks.

For each point `i` sent to a given Alya rank:
- `send_coord_bufs[dest][(i-1)*ndime + 1]` = x-coordinate
- `send_coord_bufs[dest][(i-1)*ndime + 2]` = y-coordinate
- `send_coord_bufs[dest][(i-1)*ndime + 3]` = z-coordinate (if 3D)

**Code references:**
- Coordinate packing: `couplingAuxiliaryFunctions.jl:497-503` (inside `pack_interpolated_data!`)
- MPI Isend (coords): `couplingAuxiliaryFunctions.jl:383-386`

### 6.3 MPI Communication Pattern (Julia -> Alya)

```
Julia rank J                                Alya rank A (world rank W)
  |                                              |
  |--- MPI_Isend(data,  tag=TAG_DATA+W)  ------->|  MPI_Irecv(data,  tag=TAG_DATA+W)
  |--- MPI_Isend(coord, tag=TAG_COORD+W) ------->|  MPI_Irecv(coord, tag=TAG_COORD+W)
  |                                              |
```

All communication is **non-blocking** (`MPI_Isend` / `MPI_Irecv`), followed by `MPI_Waitall`.

**Code references:**
- Julia send side: `couplingAuxiliaryFunctions.jl:348-401` (`coupling_exchange_data!`)
- Fortran receive side: `alya_all2all_time_loop.f90:281-321`

---

## 7. Time Loop: Alya Sends to Jexpresso

At the same timestep, Alya sends a return field back to each Jexpresso rank. In the current proxy implementation, this is a dummy field (`42.0 + arank + 0.01 * step`), but the protocol supports arbitrary field data.

### 7.1 Return Field Values (TAG_DATA = 2000)

| Quantity | Description | Size per Julia rank | Data type | MPI tag |
|---|---|---|---|---|
| Alya field data | `neqs` values per point | `npoin_recv[dest] * neqs` doubles | `Float64` | `TAG_DATA + dest_world_rank` |

```
Alya rank A (world rank W)                  Julia rank J (world rank V)
  |                                              |
  |--- MPI_Isend(data, tag=TAG_DATA+V)  -------->|  MPI_Irecv(data, tag=TAG_DATA+V)
  |                                              |
```

**Code references:**
- Fortran send side: `alya_all2all_time_loop.f90:323-340`
  - Fill buffer: line 275 (`sendbuf_all = dummy_field`)
  - Tag: `TAG_DATA + i` where `i` is the destination world rank (line 331/335)
- Julia receive side: `couplingAuxiliaryFunctions.jl:364-371`
  - Tag: `TAG_DATA + wrank` (line 367)

### 7.2 Unpacking Received Data in Jexpresso

The received Alya data is unpacked via `unpack_received_data!()`. This function is currently a placeholder -- the physical coupling (how Alya's field modifies Jexpresso's solution, e.g., as forcing or boundary conditions) is application-dependent.

**Code reference:** `couplingAuxiliaryFunctions.jl:527-577`

---

## 8. VTS Output: Alya Gatherv and File Write

At scheduled output times, Alya gathers all received data to Alya rank 0 and writes a VTK StructuredGrid (`.vts`) file.

### 8.1 Output Scheduling

Output is controlled by the Fortran-side variables:

```fortran
! alya_all2all_time_loop.f90, lines 245-247
out_tstart = t0         ! first eligible output time
out_dt     = 10.0d0     ! write every 10 seconds
out_tend   = tend        ! last eligible output time
```

A time-bucket algorithm determines whether the current step crosses an output time (lines 264-272).

### 8.2 MPI_Gatherv Within Alya

Two `MPI_Gatherv` calls collect data from all Alya ranks to Alya rank 0:

| Quantity | Source buffer | Destination buffer | Size | Communicator |
|---|---|---|---|---|
| Field values | `recvbuf_all` (each Alya rank's portion) | `full_field` (rank 0 only) | `nmax * neqs` doubles | `PAR_COMM_FINAL` |
| Coordinates | `recvcoord_all` (each Alya rank's portion) | `full_coord` (rank 0 only) | `nmax * ndime` doubles | `PAR_COMM_FINAL` |

**Code references:**
- Data Gatherv: `alya_all2all_time_loop.f90:554-562`
- Coordinate Gatherv: `alya_all2all_time_loop.f90:569-577`

### 8.3 Coordinate-Based Grid Placement

Alya rank 0 uses the received coordinates to compute grid indices and place values at the correct structured grid positions:

```fortran
! alya_all2all_time_loop.f90, lines 667-696
do k = 0, nmax - 1
   x = full_coord(k * ndime + 1)
   y = full_coord(k * ndime + 2)
   ix = nint((x - rem_min(1)) / dx)
   iy = nint((y - rem_min(2)) / dy)
   iz = 0
   if (ndime >= 3) iz = nint((z - rem_min(3)) / dz)
   ! clamp to valid range
   ix = max(0, min(rem_nx(1)-1, ix))
   iy = max(0, min(rem_nx(2)-1, iy))
   iz = max(0, min(rem_nx(3)-1, iz))
   grid_var(ix, iy, iz) = full_field(k * neqs + ieq)
end do
```

This coordinate-based placement is necessary because points from different Jexpresso ranks arrive in arbitrary order (determined by spatial domain decomposition and binning), not in flat-index order.

### 8.4 VTS File Structure

Each output file (`alya_grid_NNNNNN.vts`) contains:
- **Points**: Structured grid coordinates computed from `rem_min`, `rem_max`, `rem_nx`
- **PointData**: One `DataArray` per equation variable:
  - `rho` (density)
  - `u` (x-velocity)
  - `v` (y-velocity)
  - `theta` (potential temperature)

Written in VTS order: ix varies fastest, then iy, then iz.

---

## 9. Complete Data Flow Summary

### Phase 1: Initialization (one-time)

```
Direction: Alya -> All ranks (MPI_Bcast on MPI_COMM_WORLD)
============================================================
 ndime        Int32     spatial dimension (2 or 3)
 rem_min[1]   Float32   grid minimum x
 rem_max[1]   Float32   grid maximum x
 rem_nx[1]    Int32     grid points in x
 rem_min[2]   Float32   grid minimum y
 rem_max[2]   Float32   grid maximum y
 rem_nx[2]    Int32     grid points in y
 rem_min[3]   Float32   grid minimum z
 rem_max[3]   Float32   grid maximum z
 rem_nx[3]    Int32     grid points in z
 neqs         Int32     equations per point (4)
 nsteps       Int32     total number of timesteps

Direction: All ranks (collective MPI_Allreduce on MPI_COMM_WORLD)
============================================================
 alya_to_world[0:asize-1]   Int32[]   Alya local rank -> world rank map

Direction: All ranks (collective MPI_Alltoall on MPI_COMM_WORLD)
============================================================
 npoin_send / npoin_recv     Int32[]   points exchanged per rank pair
```

### Phase 2: Time Loop (every timestep)

```
Direction: Jexpresso -> Alya (non-blocking, MPI_COMM_WORLD)
============================================================
 Interpolated field data     Float64[]   neqs values per point    tag = 2000 + dest_rank
 Grid coordinates            Float64[]   ndime values per point   tag = 3000 + dest_rank

Direction: Alya -> Jexpresso (non-blocking, MPI_COMM_WORLD)
============================================================
 Return field data           Float64[]   neqs values per point    tag = 2000 + dest_rank
```

### Phase 3: VTS Output (at scheduled intervals, Alya-internal)

```
Direction: All Alya ranks -> Alya rank 0 (MPI_Gatherv on PAR_COMM_FINAL)
============================================================
 full_field    Float64[]   all field values     (nmax * neqs elements)
 full_coord    Float64[]   all coordinates      (nmax * ndime elements)
```

---

## 10. File Reference

### Julia Coupling Files

| File | Key functions / structures | Purpose |
|---|---|---|
| `src/Jexpresso-mini-coupled.jl` | Entry point, `MPI_Comm_split`, receives Alya broadcasts | Coupled driver: sets up communicators, receives grid metadata, launches Jexpresso |
| `src/kernel/coupling/couplingStructs.jl` | `CouplingData` struct, `setup_coupling_and_mesh()`, `extract_local_alya_coordinates()`, `allocate_coupling_buffers()`, `je_receive_alya_data()` | Data structures, Alya grid coordinate extraction, buffer allocation |
| `src/kernel/coupling/couplingAuxiliaryFunctions.jl` | `coupling_exchange_data!()`, `pack_interpolated_data!()`, `interpolate_solution_to_alya_coords()`, `je_perform_coupling_exchange()`, `unpack_received_data!()` | Interpolation, data packing, MPI exchange, unpacking |
| `src/mpi/JexpressoCoupling.jl` | `CouplingContext`, `CouplingBuffer`, `initialize_coupling()` | Generic coupling module (communicator management, buffered exchanges) |
| `src/kernel/solvers/TimeIntegrators.jl` | `time_loop!()`, `cb_coupling` callback | ODE solver with coupling callback at every timestep |
| `problems/CompEuler/thetaAlya/user_inputs.jl` | `user_inputs()` | Problem configuration (dt, tend, diagnostics_at_times, mesh file, etc.) |

### Fortran Coupling Files

| File | Key sections | Purpose |
|---|---|---|
| `AlyaProxy/alya_all2all_time_loop.f90` | Main program, `write_alya_grid_vts` subroutine | Alya proxy: broadcasts metadata, exchanges data at each step, writes VTS output |

### MPI Tags

| Tag constant | Value | Used for |
|---|---|---|
| `TAG_DATA` | `2000` | Field data exchange (both directions) |
| `TAG_COORD` | `3000` | Coordinate data (Jexpresso -> Alya only) |

Actual tag used in messages: `TAG + destination_world_rank` (ensures uniqueness when multiple senders communicate with the same receiver).

### Key Data Structures

**`CouplingData`** (Julia, `couplingStructs.jl:42-73`):
```
npoin_recv       Int32[]                  Points to receive from each world rank
npoin_send       Int32[]                  Points to send to each world rank
recv_from_ranks  Int32[]                  World ranks we receive from
send_to_ranks    Int32[]                  World ranks we send to
comm_world       MPI.Comm                 MPI_COMM_WORLD
lrank            Int32                    Local rank in Jexpresso communicator
neqs             Int                      Number of equations per point
ndime            Int                      Spatial dimension (2 or 3)
send_bufs        Vector{Vector{Float64}}  Per-destination data send buffers
recv_bufs        Vector{Vector{Float64}}  Per-source data receive buffers
send_coord_bufs  Vector{Vector{Float64}}  Per-destination coordinate send buffers
alya_local_coords Matrix{Float64}         Local Alya grid coordinates [n_local x ndime]
alya_local_ids    Vector{Int32}           Global Alya point IDs (1-based)
alya_owner_ranks  Vector{Int32}           Owning Alya world rank for each local point (0-based)
```
