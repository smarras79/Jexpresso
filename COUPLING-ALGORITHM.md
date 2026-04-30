
# Coupling Algorithm: Jexpresso -- Alya

This document describes the coupling algorithm between Jexpresso (Julia SEM solver) and Alya (Fortran CFD code), from the initial handshake through the time-loop exchange of interpolated fields. It focuses exclusively on the coupling logic and the data that flows between the two codes.

---

## 1. Handshake

Both codes are launched as a single MPI job sharing `MPI_COMM_WORLD`:

```
mpirun -np <N_alya> ./AlyaProxy/Alya.x \
     : -np <N_jexpresso> julia --project=. ./src/Jexpresso.jl CompEuler thetaAlya
```

Alya ranks occupy world ranks `0 .. N_alya-1`; Jexpresso occupies `N_alya .. N_alya+N_jexpresso-1`.

### 1.1 Communicator split

Each code calls `MPI_Comm_split` with its own color to create a local communicator:

- **Alya:** `MPI_Comm_split(MPI_COMM_WORLD, 1, rank, PAR_COMM_FINAL)` -- `alya_all2all_time_loop.f90:78`
- **Jexpresso:** `MPI.Comm_split(world, appid, wrank)` -- `Jexpresso-mini-coupled.jl:31`

### 1.2 Application name exchange

Every rank sends a 128-character application name to world rank 0 via `MPI_Gather`:

```
rank -> "ALYA      ..."   (128 bytes)     alya_all2all_time_loop.f90:85-87
rank -> "JEXPRESSO ..."   (128 bytes)     Jexpresso-mini-coupled.jl:42
```

This confirms how many ranks belong to each code.

### 1.3 Coupling detection

In `drivers.jl:20`, Jexpresso calls `je_perform_coupling_handshake(world, lsize)` (`couplingStructs.jl:826-845`). The logic is simple: if `wsize > nparts` (more world ranks than Jexpresso ranks), coupling is active.

**Code path:** `drivers.jl:20` -> `couplingStructs.jl:826`

---

## 2. Receiving Alya's Grid Metadata

Once the handshake confirms coupled mode, the driver calls `setup_coupling_and_mesh()` (`drivers.jl:25`), which immediately calls `je_receive_alya_data(world, lsize)` (`couplingStructs.jl:847`).

Alya rank 0 broadcasts the following scalars to all ranks in `MPI_COMM_WORLD`, using `MPI_Bcast`. The broadcasts must happen in exactly this order on both sides:

| # | Quantity | Fortran source | Julia receiver | Type |
|---|---|---|---|---|
| 1 | `ndime` | `alya_all2all_time_loop.f90:124` | `couplingStructs.jl:857-859` | `Int32` |
| 2 | `rem_min(1), rem_max(1), rem_nx(1)` | `:127-129` (idime=1) | `:865-869` | `Float32, Float32, Int32` |
| 3 | `rem_min(2), rem_max(2), rem_nx(2)` | `:127-129` (idime=2) | `:865-869` | `Float32, Float32, Int32` |
| 4 | `rem_min(3), rem_max(3), rem_nx(3)` | `:127-129` (idime=3) | `:865-869` | `Float32, Float32, Int32` |
| 5 | `neqs` | `:133` | `:872-874` | `Int32` |
| 6 | `nsteps` | `:135` | `:877-879` | `Int32` |

**Physical meaning:**

- `ndime`: spatial dimension (2 or 3)
- `rem_min[1:3]`, `rem_max[1:3]`: bounding box of the Alya structured grid
- `rem_nx[1:3]`: number of grid nodes in each direction
- `neqs`: number of field variables per grid node (4 for the compressible Euler system: rho, rho*u, rho*v, rho*theta)
- `nsteps`: total number of time steps Alya will perform

Jexpresso stores these in the global dictionary `JEXPRESSO_COUPLING_DATA` (`couplingStructs.jl:889-897`).

---

## 3. Alya-to-World Rank Mapping

Immediately after the grid metadata broadcasts, all ranks collectively build a map from Alya local rank indices to world rank indices:

```
alya_to_world[arank] = world_rank
```

Each Alya rank sets its own entry; all other entries are zero. An `MPI_Allreduce(MPI_SUM)` combines them into the complete map.

- **Fortran:** `alya_all2all_time_loop.f90:140-145`
- **Julia:** `couplingStructs.jl:882-887`

This map is essential because Jexpresso needs to know which world rank to send data to for a given Alya rank.

---

## 4. Extracting Alya Grid Coordinates on Each Jexpresso Rank

After receiving the grid metadata, Jexpresso sets up its own SEM mesh (`sem_setup`) and then determines which Alya grid points fall within each Jexpresso rank's spatial domain. This is done by `extract_local_alya_coordinates()` (`couplingStructs.jl:534-801`).

### 4.1 Computing Alya grid coordinates

Alya does **not** send actual coordinate arrays. Instead, Jexpresso reconstructs every Alya grid point from the metadata:

```
dx = (rem_max[d] - rem_min[d]) / (rem_nx[d] - 1)    for each dimension d

For a point with structured indices (i1, i2, i3) (all 0-based):
    x = rem_min[1] + i1 * dx
    y = rem_min[2] + i2 * dy
    z = rem_min[3] + i3 * dz
```

The flat (1-based) index is `ipoin = i1 + nx1*(i2 + nx2*i3) + 1`.

### 4.2 Spatial filtering with index-space cropping and binning

Rather than looping over all `nmax = nx1 * nx2 * nx3` Alya points, the function uses two levels of pruning:

1. **Index-space cropping** (`couplingStructs.jl:614-657`): For each dimension, compute the range of Alya grid indices whose physical coordinates could overlap the local Jexpresso bounding box. Only scan indices in `[i_lo, i_hi]` per dimension.

2. **Block-based spatial bins** (`couplingStructs.jl:713-763`): The cropped index range is divided into blocks of size `block_size` (default 64x64x64). For each block:
   - Compute the block's physical bounding box.
   - If the block is entirely outside the local domain: skip it.
   - If the block is entirely inside: accept all points without per-point checks.
   - Otherwise: check each point individually against the local bounding box.

### 4.3 Determining the Alya owner rank

For each accepted point, Jexpresso computes which Alya rank "owns" it. Alya distributes points in contiguous chunks across its ranks:

```julia
# couplingStructs.jl:604-609
r     = mod(nmax, nranks_alya)
npoin = div(nmax, nranks_alya)
# First r ranks get (npoin+1) points; remaining get npoin points
if ipoin <= r * (npoin + 1)
    alya_rank = div(ipoin - 1, npoin + 1) + 1
else
    alya_rank = r + div(ipoin - r*(npoin+1) - 1, npoin) + 1
end
alya_world_rank = alya2world[alya_rank]
```

### 4.4 Sorting for consistent ordering

After extraction, the local point list is sorted by `(owner_rank, global_id)` (`couplingStructs.jl:786-791`). This ensures that points sent to each Alya rank arrive in ascending flat-index order, consistent with the `MPI_Gatherv` displacement layout on the Alya side.

### 4.5 Output

The function returns three arrays:

- `alya_local_coords[n_local, ndime]`: physical coordinates of each accepted point
- `alya_local_ids[n_local]`: 1-based global Alya point ID
- `alya_owner_ranks[n_local]`: 0-based world rank of the owning Alya rank

---

## 5. Building `npoin_recv` and `npoin_send`

With the extracted point list in hand, `setup_coupling_and_mesh()` (`couplingStructs.jl:94-129`) builds the communication-count arrays that Alya needs.

### 5.1 Building `npoin_recv` from the extracted points

```julia
# couplingStructs.jl:97-102
npoin_recv = zeros(Int32, wsize)   # one entry per world rank
for owner_wrank in alya_owner_ranks
    npoin_recv[owner_wrank + 1] += 1
end
```

`npoin_recv[r+1]` tells Jexpresso how many points it needs to exchange with world rank `r`. This is computed directly from the extracted point ownership -- no separate counting pass over the full grid is needed.

### 5.2 Exchanging counts with Alya via `MPI_Alltoall`

Both codes participate in a collective `MPI_Alltoall` on `MPI_COMM_WORLD`:

```julia
# Jexpresso (couplingStructs.jl:120-123)
send_counts_to_alya = Vector{Int32}(npoin_recv)  # what I will send
recv_counts_from_alya = zeros(Int32, wsize)
MPI.Alltoall!(send_counts_to_alya, recv_counts_from_alya, 1, world)
```

```fortran
! Alya (alya_all2all_time_loop.f90:158-166)
npoin_send = 0   ! Alya initializes to 0 (never initiates)
call MPI_Alltoall(npoin_send, 1, MPI_INTEGER4, &
                  npoin_recv, 1, MPI_INTEGER4, MPI_COMM_WORLD, ierr)
```

After `Alltoall`:
- **Alya** has `npoin_recv[j]` = number of points it will receive from Jexpresso world rank `j`.
- **Jexpresso** has `recv_counts_from_alya[j]` = 0 for all `j` (since Alya sent zeros), but this is expected -- the return exchange uses the same point counts as the forward exchange.

### 5.3 Setting `npoin_send`

Because the coupling is symmetric (Alya sends back exactly the same points it receives), Jexpresso sets:

```julia
# couplingStructs.jl:128-129
npoin_send   = copy(npoin_recv)
send_to_ranks = copy(recv_from_ranks)
```

### 5.4 Buffer allocation

With the final counts known, `allocate_coupling_buffers()` (`couplingStructs.jl:803-824`) creates per-destination send/receive buffers:

- `send_bufs[r]`: `npoin_send[r] * neqs` doubles (field values to send to world rank `r`)
- `recv_bufs[r]`: `npoin_recv[r] * neqs` doubles (field values to receive from world rank `r`)
- `send_coord_bufs[r]`: `npoin_send[r] * ndime` doubles (coordinates to send to world rank `r`)

---

## 6. Time-Loop: Interpolation and Exchange

The coupling exchange happens at every time step via a `DiscreteCallback` registered in the ODE solver (`TimeIntegrators.jl:124-149`).

### 6.1 Coupling callback trigger

```julia
# TimeIntegrators.jl:136
coupling_condition(u_state, t, integrator) = t > t0 + tol0
```

The callback fires at every step after the initial time. When it fires, it calls `je_perform_coupling_exchange()` (`couplingAuxiliaryFunctions.jl:841-869`).

### 6.2 Solution preparation

Before interpolation, the conservative state vector `u` is converted to output (primitive) variables:

```julia
# couplingAuxiliaryFunctions.jl:847-849
u2uaux!(u_mat, u, neqs, npoin)                        # reshape flat vector to [npoin x neqs]
call_user_uout(qout, u_mat, u_mat, 0, SOL_VARS_TYPE, npoin, neqs, neqs)
```

The `user_uout!` function (defined in `user_primitives.jl:24-30`) converts conservative to primitive variables:

```
qout[ip, 1] = rho          = u[1]
qout[ip, 2] = velocity_x   = u[2] / u[1]    (= rho*u / rho)
qout[ip, 3] = velocity_y   = u[3] / u[1]    (= rho*v / rho)
qout[ip, 4] = theta        = u[4] / u[1]    (= rho*theta / rho)
```

### 6.3 Interpolation onto Alya coordinates

`interpolate_solution_to_alya_coords()` (`couplingAuxiliaryFunctions.jl:619-704`) interpolates the SEM solution to each local Alya grid point.

**Algorithm for each Alya point `(px, py)`:**

1. **Find the containing SEM element:**
   - Build spatial bins over element bounding boxes (`couplingStructs.jl:457-501`).
   - Look up candidate elements from the bin containing `(px, py)`.
   - For each candidate, do a quick bounding-box check, then map `(px, py)` to reference coordinates.

2. **Map to reference coordinates via Newton iteration** (`couplingAuxiliaryFunctions.jl:735-787`):
   - Starting from `(xi, eta) = (0, 0)`, iterate:
     ```
     Evaluate basis and derivatives at (xi, eta)
     Compute physical position:  x_curr = sum( psi_i * x_elem_i )
     Compute Jacobian:           J = [dx/dxi  dx/deta; dy/dxi  dy/deta]
     Residual:                   r = (px - x_curr, py - y_curr)
     Newton update:              (xi, eta) += J^{-1} * r
     ```
   - Accept if `|r| < 1e-12` and `|xi|, |eta| <= 1 + 1e-10`.

3. **Evaluate Lagrange basis and interpolate** (`couplingAuxiliaryFunctions.jl:674-687`):
   - Compute 1D basis values at `xi_ref` and `eta_ref` using the barycentric formula.
   - The 2D basis is the tensor product: `psi_{ij}(xi, eta) = psi_i(xi) * psi_j(eta)`.
   - Interpolate each equation:
     ```
     u_interp[pt, q] = sum_{i,j} psi_i(xi_ref) * psi_j(eta_ref) * qout[node_{ij}, q]
     ```

4. **Fallback:** If no element contains the point (should not happen if domains overlap correctly), use nearest-neighbor (`couplingAuxiliaryFunctions.jl:693-700`).

**Result:** `u_interp[n_local, neqs]` -- interpolated primitive variables at every local Alya grid point.

### 6.4 Packing data and coordinates into send buffers

`pack_interpolated_data!()` (`couplingAuxiliaryFunctions.jl:465-508`) packs the interpolated values **and** the Alya grid coordinates into per-destination buffers.

For each local Alya point `i`:

```julia
owner_rank = alya_owner_ranks[i]     # 0-based world rank

# Pack neqs field values contiguously
send_bufs[owner_rank+1][offset+1 : offset+neqs] = u_interp[i, 1:neqs]

# Pack ndime coordinate values contiguously
send_coord_bufs[owner_rank+1][coffset+1 : coffset+ndime] = alya_local_coords[i, 1:ndime]
```

**Why send coordinates?** Different Jexpresso ranks extract Alya points in an order determined by spatial binning, not by flat grid index. When Alya receives data from multiple Jexpresso ranks and concatenates the buffers, the values are not in flat-index order. By also receiving the coordinates, Alya can compute the correct grid position for each value: `ix = nint((x - rem_min) / dx)`.

### 6.5 MPI exchange

`coupling_exchange_data!()` (`couplingAuxiliaryFunctions.jl:348-401`) performs the non-blocking MPI exchange.

**Jexpresso posts receives** for the return field from Alya:

```julia
# For each Alya rank that sends to us:
tag = TAG_DATA + my_world_rank          # = 2000 + wrank
MPI.Irecv!(recv_bufs[src+1], src, tag, comm_world)
```

**Jexpresso posts sends** of interpolated data and coordinates:

```julia
# For each Alya rank we send to:
tag_data  = TAG_DATA  + dest_world_rank   # = 2000 + dest
tag_coord = TAG_COORD + dest_world_rank   # = 3000 + dest
MPI.Isend(send_bufs[dest+1],       dest, tag_data,  comm_world)
MPI.Isend(send_coord_bufs[dest+1], dest, tag_coord, comm_world)
```

**Alya posts matching receives** (`alya_all2all_time_loop.f90:286-321`):

```fortran
! Data receive from each sender:
tag = TAG_DATA + rank                    ! 2000 + my_world_rank
MPI_Irecv(recvbuf_all(offset+1), npoin_recv(i)*neqs, ..., i, tag, ...)

! Coordinate receive from each sender:
tag = TAG_COORD + rank                   ! 3000 + my_world_rank
MPI_Irecv(recvcoord_all(offset+1), npoin_recv(i)*ndime, ..., i, tag, ...)
```

**Alya posts sends** of its return field (`alya_all2all_time_loop.f90:326-340`):

```fortran
tag = TAG_DATA + i                       ! 2000 + dest_world_rank
MPI_Isend(sendbuf_all(offset+1), npoin_recv(i)*neqs, ..., i, tag, ...)
```

All operations complete with `MPI_Waitall`.

### 6.6 Tag matching summary

| Message | Sender | Receiver | Tag formula |
|---|---|---|---|
| Interpolated field data | Jexpresso rank J | Alya rank A (world W) | `2000 + W` |
| Grid coordinates | Jexpresso rank J | Alya rank A (world W) | `3000 + W` |
| Return field from Alya | Alya rank A (world W) | Jexpresso rank J (world V) | `2000 + V` |

The tag always contains the **receiver's** world rank, ensuring uniqueness when multiple senders communicate with the same receiver.

### 6.7 Receiving Alya's return field

After `MPI_Waitall`, `unpack_received_data!()` (`couplingAuxiliaryFunctions.jl:527-577`) processes the field received from Alya. In the current implementation this is a placeholder -- the physical coupling (applying Alya's field as forcing or boundary conditions on Jexpresso) is application-dependent and left for future development.

---

## 7. Summary: Step-by-Step Coupling Flow

```
INITIALIZATION (one-time)
=========================

1. Handshake
   drivers.jl:20  ->  couplingStructs.jl:826
   All ranks: MPI_Gather(app_name)
   Jexpresso detects coupling: wsize > nparts

2. Receive Alya's grid metadata
   drivers.jl:25  ->  couplingStructs.jl:75  ->  couplingStructs.jl:847
   Alya rank 0: MPI_Bcast(ndime, rem_min, rem_max, rem_nx, neqs, nsteps)
   All ranks:   MPI_Allreduce(alya_to_world map)

3. Setup SEM mesh
   couplingStructs.jl:84
   Standard Jexpresso mesh partitioning (independent of coupling)

4. Extract local Alya coordinates
   couplingStructs.jl:90-92  ->  couplingStructs.jl:534
   Each Jexpresso rank: identify Alya grid points inside local domain
   Result: coordinates, global IDs, owner Alya ranks

5. Build communication arrays
   couplingStructs.jl:94-129
   Build npoin_recv from extracted ownership
   MPI_Alltoall to inform Alya of incoming point counts
   Set npoin_send = npoin_recv (symmetric exchange)
   Allocate send/recv/coord buffers

6. Store coupling object
   couplingStructs.jl:163-181
   CouplingData struct holds all arrays and buffers
   Attached to params for use in time loop


TIME LOOP (every step)
======================

7. Coupling callback fires
   TimeIntegrators.jl:136-148
   Condition: t > t0 + tol

8. Prepare solution
   couplingAuxiliaryFunctions.jl:847-849
   Convert conservative -> primitive variables
   [rho, rho*u, rho*v, rho*theta] -> [rho, u, v, theta]

9. Interpolate onto Alya coordinates
   couplingAuxiliaryFunctions.jl:852-856  ->  :619
   For each local Alya point:
     Find containing SEM element (spatial bins)
     Newton iteration: (px,py) -> (xi,eta) in reference space
     Evaluate tensor-product Lagrange basis
     u_interp = sum( psi_ij * qout[node_ij] )

10. Pack into send buffers
    couplingAuxiliaryFunctions.jl:859  ->  :465
    Pack neqs field values + ndime coordinates per point
    Grouped by destination Alya rank

11. MPI exchange
    couplingAuxiliaryFunctions.jl:862  ->  :348
    Jexpresso -> Alya: field data (tag 2000+dest) + coordinates (tag 3000+dest)
    Alya -> Jexpresso: return field (tag 2000+dest)
    All non-blocking (Irecv/Isend + Waitall)

12. Unpack received Alya data
    couplingAuxiliaryFunctions.jl:865  ->  :527
    (Placeholder for applying Alya's field to Jexpresso solution)
```
