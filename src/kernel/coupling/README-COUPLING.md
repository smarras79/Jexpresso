# Jexpresso-Alya Coupling Algorithm Documentation

**Version:** 1.0  
**Date:** 2024  
**Authors:** Jexpresso Development Team

---

## Table of Contents

1. [Overview](#overview)
2. [Algorithm Breakdown](#algorithm-breakdown)
3. [Where to Add Interpolation](#where-to-add-interpolation)
4. [Required New Functions](#required-new-functions)
5. [Summary Table](#summary-table)
6. [Next Steps](#next-steps)

---

## Overview

This document describes the coupling initialization workflow between Jexpresso and Alya, focusing on the coordinate exchange and where interpolation will be implemented.

The `setup_coupling_and_mesh()` function orchestrates the complete coupling initialization between Jexpresso and Alya. It establishes the communication pattern, identifies which Jexpresso ranks own which Alya grid points, and prepares buffers for data exchange.

### Function Signature
```julia
coupling, sem, partitioned_model, qp = setup_coupling_and_mesh(
    world,        # MPI.COMM_WORLD
    lsize,        # Number of Jexpresso ranks
    inputs,       # Input parameters dictionary
    nranks,       # Number of ranks for mesh setup
    distribute,   # Distribution function
    rank,         # Current rank
    OUTPUT_DIR,   # Output directory path
    TFloat        # Floating point type (Float32/Float64)
)
```

### Returns

- `coupling`: CouplingData structure with communication pattern and buffers
- `sem`: SEM mesh and solver structures  
- `partitioned_model`: Partitioned mesh model
- `qp`: Initialized solution arrays

---

## Algorithm Breakdown

### Step 1: Receive Alya's Domain Metadata
```julia
je_receive_alya_data(world, lsize)
```

#### What it does

- Receives structured grid information from Alya via `MPI.Bcast` on `COMM_WORLD`
- Alya rank 0 broadcasts:
  - `ndime`: Spatial dimension (2 or 3)
  - `rem_min`: Grid bounding box minimum `[xmin, ymin, zmin]`
  - `rem_max`: Grid bounding box maximum `[xmax, ymax, zmax]`
  - `rem_nx`: Grid resolution `[nx, ny, nz]`
- Builds `alya2world` mapping (Alya local rank → world rank)
- Stores this metadata in global `JEXPRESSO_COUPLING_DATA`

#### Key Points

> ⚠️ **Important:** Alya does **NOT** send actual coordinate arrays. It only sends the bounding box and resolution for its structured grid. Jexpresso will compute Alya coordinates on-the-fly from this metadata.

---

### Step 2: Setup Jexpresso's SEM Mesh
```julia
sem, partitioned_model = sem_setup(inputs, nranks, distribute, rank)
```

#### What it does

- Creates the Spectral Element Method (SEM) mesh for Jexpresso
- Partitions the domain across Jexpresso ranks
- Each rank knows its local mesh coordinates: `sem.mesh.x`, `sem.mesh.y`, `sem.mesh.z`
- Converts mesh arrays to GPU if using GPU backend

#### Key Points

- After this step, each Jexpresso rank knows its **local domain bounding box**

---

### Step 3: 🔴 Build Coupling Communication Pattern

**THIS IS WHERE ALYA COORDINATES ARE RECEIVED**
```julia
npoin_recv, npoin_send, recv_from_ranks, send_to = 
    build_coupling_communication_arrays(sem.mesh, coupling_data, 
                                        local_comm, world)
```

#### What it does

This is the **critical function** where Alya's grid coordinates are computed and ownership is determined.

#### Detailed Workflow

1. **Computes Alya grid spacing** from metadata:
```julia
   rem_dx[idim] = (rem_max[idim] - rem_min[idim]) / (rem_nx[idim] - 1)
```

2. **Loops over ALL Alya grid points** (e.g., 100 points in a 10×10 grid):
```julia
   for ipoin in 1:nmax
       # Convert flat index to 3D grid indices
       i0 = ipoin - 1
       ri3 = div(i0, rem_nx[1] * rem_nx[2])
       ri2 = div(i0 - ri3 * rem_nx[1] * rem_nx[2], rem_nx[1])
       ri1 = mod(i0 - ri3 * rem_nx[1] * rem_nx[2] - ri2 * rem_nx[1], rem_nx[1])
       
       # 🔴 COMPUTE ALYA COORDINATE from structured grid
       x[1] = rem_min[1] + ri1 * rem_dx[1]
       x[2] = rem_min[2] + ri2 * rem_dx[2]
       x[3] = rem_min[3] + ri3 * rem_dx[3]
       
       # Check if this coordinate is in my Jexpresso domain
       in_my_rank = point_in_my_domain(x, mesh.bounds)
       
       if in_my_rank
           # Determine which Alya rank owns this point
           alya_rank = compute_alya_rank_owner(ipoin, nmax, nranks_alya)
           alya_world_rank = alya2world[alya_rank]
           
           # This Alya point belongs to me!
           npoin_recv[alya_world_rank + 1] += 1
       end
   end
```

3. **Determines which Alya rank owns each point** using the distribution formula

4. **Builds receive counts**: `npoin_recv[i]` = how many points to receive from world rank `i`

#### Key Points

> ⚠️ **THIS IS WHERE ALYA COORDINATES ARE COMPUTED**  
> - Each Jexpresso rank identifies which Alya grid points fall within its domain
> - However, **the actual coordinate values are not yet stored** for later use
> - We only count how many points we receive from each rank

---

### Step 4: Exchange Communication Pattern with Alya
```julia
MPI.Barrier(world)
MPI.Alltoall!(send_counts_to_alya, recv_counts_from_alya, 1, world)
```

#### What it does

- Synchronizes all processes with a barrier
- Uses `MPI_Alltoall` to transpose the communication pattern:
  - **Jexpresso → Alya**: "I need X points from you"
  - **Alya → Jexpresso**: "You want Y points from me"
- After this, both codes know:
  - `npoin_recv[i]`: How many points I receive from world rank i
  - `npoin_send[i]`: How many points I send to world rank i

#### Key Points

- This establishes the **bidirectional communication pattern**
- Both codes now know the data exchange topology

---

### Step 5: Print Communication Pattern (Verification)
```julia
print_coupling_pattern(lrank, npoin_recv, npoin_send, recv_from_ranks, 
                      send_to_ranks, coupling_data, world)
```

#### What it does

- Prints detailed information about communication pattern
- Shows which ranks send/receive to/from which ranks
- Displays world rank breakdown (Alya vs Jexpresso ranks)
- Useful for debugging and verification

#### Example Output
```
=========================================
Julia local rank 0, world rank 2
=========================================
Total points to RECV from Alya: 50
Total points to SEND to Alya: 0
-----------------------------------------
World rank breakdown (size=4):
  World rank 0 = Alya rank 0 (ALYA)
  World rank 1 = Alya rank 1 (ALYA)
  World rank 2 = Julia rank 0 (JULIA)
  World rank 3 = Julia rank 1 (JULIA)
-----------------------------------------
RECV pattern (npoin_recv array):
  <- World rank 0 (Alya rank 0): 25 points
  <- World rank 1 (Alya rank 1): 25 points
=========================================
```

---

### Step 6: Build Alya Point Ownership Map (Verification)
```julia
ownership = build_alya_point_ownership_map(sem.mesh, coupling_data, 
                                           local_comm, world)
```

#### What it does

- **Re-loops over ALL Alya grid points** (similar to Step 3)
- For each point:
  - 🔴 **Computes the Alya coordinate** from grid metadata
  - Determines which Jexpresso rank owns it
  - **Stores the coordinate** in `ownership.alya_coords`
  - Stores ownership information
- Uses `MPI.Allreduce` to gather ownership information across all ranks

#### Key Points

- This creates a **complete map** of Alya coordinates and their owners
- Used for verification and debugging
- **This is currently only for printing, not used in actual coupling yet**

---

### Step 7: Print Ownership Information (Verification)
```julia
print_alya_point_ownership(ownership, coupling_data; 
                          max_points_to_print=20, only_owned=true)
verify_alya_point_distribution(ownership, coupling_data)
```

#### What it does

- Prints which Jexpresso rank owns each Alya coordinate
- Shows distribution statistics
- Verifies that all Alya points are owned by some Jexpresso rank

#### Example Output
```
======================================================================
Alya Point Ownership Map
======================================================================
Total Alya points: 100
Spatial dimension: 2
----------------------------------------------------------------------
Points owned by Jexpresso: 100
Points not owned: 0
----------------------------------------------------------------------
Point ID | Coordinates (x, y)          | Owner
----------------------------------------------------------------------
       1 | (  -5000.00,       0.00) | Jexpresso rank 1 (world rank 3)
       2 | (  -3888.89,       0.00) | Jexpresso rank 1 (world rank 3)
       6 | (    555.56,       0.00) | Jexpresso rank 0 (world rank 2)
       7 | (   1666.67,       0.00) | Jexpresso rank 0 (world rank 2)
...

======================================================================
Alya Point Distribution Statistics
======================================================================
Points per Jexpresso rank:
  Jexpresso rank 0: 50 points (50.0%)
  Jexpresso rank 1: 50 points (50.0%)
----------------------------------------------------------------------
✓ All Alya points are owned by Jexpresso ranks
======================================================================
```

---

### Step 8: Initialize Solution Arrays
```julia
qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)
```

#### What it does

- Initializes the solution arrays for the PDE solver
- Creates:
  - `qp.qn`: Numerical solution vector
  - `qp.qe`: Exact solution (if available)
  - Other solution-related arrays
- **After this step, we have the solution vector `u[:]` or `qp.qn[:]`**

#### Key Points

> ✅ **Solution arrays are now available for interpolation**

---

### Step 9: Create Coupling Data Structure
```julia
coupling = CouplingData(
    npoin_recv = npoin_recv,
    npoin_send = npoin_send,
    recv_from_ranks = recv_from_ranks,
    send_to_ranks = send_to_ranks,
    comm_world = world,
    lrank = lrank,
    neqs = qp.neqs
)
```

#### What it does

- Packages all coupling information into a single structure
- Stores:
  - Communication pattern arrays
  - MPI communicator and rank information
  - Number of equations per point
- Will hold send/receive buffers for time-stepping

---

### Step 10: Allocate Communication Buffers
```julia
coupling.send_bufs, coupling.recv_bufs = 
    allocate_coupling_buffers(npoin_recv, npoin_send, coupling.neqs)
```

#### What it does

- Allocates `Vector{Vector{Float64}}` buffers for each world rank
- `recv_bufs[i]`: buffer to receive data from world rank `i`
- `send_bufs[i]`: buffer to send data to world rank `i`
- Each buffer sized for `npoin[i] × neqs` values

#### Key Points

- These buffers will hold **interpolated velocities** to send to Alya
- Currently empty, will be filled during time-stepping

---

## Where to Add Interpolation

### Current Gap in the Algorithm

Currently, the code:

- ✅ Computes Alya coordinates (in Step 3 and Step 6)
- ✅ Knows which Alya points belong to each Jexpresso rank
- ✅ Has solution arrays `u[:]` or `qp.qn[:]` available
- ❌ **Does NOT store the Alya coordinates that belong to each rank for coupling**
- ❌ **Does NOT interpolate the solution onto these coordinates**
- ❌ **Does NOT pack interpolated values into send buffers**

---

### Where Interpolation Should Be Added

#### Option A: Add to `setup_coupling_and_mesh` (Initialization)

**Location:** Between Steps 8 and 9
```julia
# Step 8: Initialize solution arrays
qp = initialize(sem.mesh.SD, 0, sem.mesh, inputs, OUTPUT_DIR, TFloat)

# ⭐ NEW STEP 8.5: Store Alya coordinates that belong to this rank
# and interpolate initial solution onto them
alya_local_coords, alya_local_indices = extract_local_alya_coordinates(
    sem.mesh, coupling_data, local_comm, world
)

# ⭐ NEW STEP 8.6: Interpolate Jexpresso solution onto Alya coordinates
alya_interpolated_values = interpolate_solution_to_alya_coords(
    alya_local_coords, sem.mesh, qp.qn, sem.basis, qp.neqs
)

# Step 9: Create coupling structure (add the new data)
coupling = CouplingData(
    npoin_recv = npoin_recv,
    npoin_send = npoin_send,
    recv_from_ranks = recv_from_ranks,
    send_to_ranks = send_to_ranks,
    comm_world = world,
    lrank = lrank,
    neqs = qp.neqs,
    alya_local_coords = alya_local_coords,           # NEW
    alya_local_indices = alya_local_indices,         # NEW
    alya_interp_values = alya_interpolated_values    # NEW
)
```

**Advantages:**
- Stores Alya coordinates once during initialization
- Can verify interpolation before time-stepping begins
- Useful for debugging and validation

---

#### Option B: Add to Time Loop (Every Time Step)

**Location:** In `time_loop!` function, within the coupling callback
```julia
function do_coupling_exchange!(integrator)
    # integrator.u contains current solution
    
    # ⭐ Interpolate current solution onto Alya coordinates
    for ip in 1:length(coupling.alya_local_coords)
        alya_x = coupling.alya_local_coords[ip, 1]
        alya_y = coupling.alya_local_coords[ip, 2]
        
        # Interpolate all equations at this point
        for ieq in 1:coupling.neqs
            u_interp = interpolate_at_point(
                alya_x, alya_y, params.mesh, 
                view(integrator.u, :, ieq), params.basis
            )
            
            # Pack into send buffer for appropriate Alya rank
            alya_rank = ... # determine from stored info
            idx = (ip - 1) * coupling.neqs + ieq
            coupling.send_bufs[alya_rank][idx] = u_interp
        end
    end
    
    # Exchange data with Alya via MPI
    exchange_coupling_data!(coupling)
end
```

**Advantages:**
- Interpolates with current solution at each time step
- Necessary for time-dependent problems
- Reflects actual physical coupling

---

### Recommended Approach

**Use BOTH options:**

1. **In `setup_coupling_and_mesh` (Initialization):**
   - Extract and **store** the Alya coordinates that belong to this rank
   - Perform **initial interpolation** for verification
   - Add coordinates to `CouplingData` structure
   - Verify that interpolation works correctly

2. **In `time_loop!` coupling callback:**
   - **Re-interpolate** at each time step using stored coordinates
   - Pack interpolated values into send buffers
   - Exchange with Alya via `MPI_Alltoallv` or similar

This two-stage approach allows for:
- Early detection of interpolation issues
- Efficient time-stepping (coordinates stored, not recomputed)
- Clear separation of initialization and runtime coupling

---

## Required New Functions

### 1. Extract Local Alya Coordinates
```julia
"""
    extract_local_alya_coordinates(mesh, coupling_data, local_comm, world_comm)

Extract the subset of Alya grid coordinates that belong to this Jexpresso rank.

# Arguments
- `mesh`: Jexpresso mesh structure
- `coupling_data`: Dict with Alya grid metadata
- `local_comm`: Jexpresso local communicator
- `world_comm`: MPI.COMM_WORLD

# Returns
- `alya_local_coords`: Matrix [n_local_points × ndime] of Alya coordinates
- `alya_local_ids`: Vector [n_local_points] of global Alya point IDs
- `alya_owner_ranks`: Vector [n_local_points] of Alya world ranks that own each point

# Description
Similar to `build_coupling_communication_arrays` but stores coordinates instead 
of just counting. This is necessary for interpolation during coupling.
"""
function extract_local_alya_coordinates(mesh, coupling_data, local_comm, world_comm)
    # TODO: Implement
    # Loop over Alya grid points (like Step 3)
    # Store coordinates that fall in my domain
    # Return coordinate matrix and associated metadata
end
```

---

### 2. Interpolate Solution to Alya Coordinates
```julia
"""
    interpolate_solution_to_alya_coords(alya_coords, mesh, u, basis, neqs)

Interpolate Jexpresso solution onto Alya grid coordinates using SEM basis functions.

# Arguments
- `alya_coords`: Matrix [n_points × ndime] of Alya coordinates
- `mesh`: Jexpresso SEM mesh structure
- `u`: Solution vector [npoin × neqs]
- `basis`: SEM basis function information
- `neqs`: Number of equations per point

# Returns
- `u_interp`: Matrix [n_points × neqs] of interpolated values

# Description
For each Alya coordinate:
1. Find which SEM element contains it
2. Map physical coordinates to reference coordinates (ξ, η) ∈ [-1,1]²
3. Evaluate Lagrange basis functions at (ξ, η)
4. Interpolate solution: u_interp = Σᵢ uᵢ φᵢ(ξ, η)

# Algorithm Details
- Use element bounding boxes for fast point location
- Handle edge cases (points on element boundaries)
- Use tolerance for floating-point comparisons
- Raise error if point not found in any element
"""
function interpolate_solution_to_alya_coords(alya_coords, mesh, u, basis, neqs)
    # TODO: Implement
    # For each Alya coordinate:
    #   1. Find containing element
    #   2. Map to reference space
    #   3. Evaluate basis functions
    #   4. Interpolate solution
    # Return interpolated values
end
```

---

### 3. Helper Function: Find Containing Element
```julia
"""
    find_containing_element(x, y, mesh)

Find which SEM element contains the given physical coordinate.

# Arguments
- `x, y`: Physical coordinates
- `mesh`: SEM mesh structure

# Returns
- `iel`: Element index (1-based), or 0 if not found
- `ξ, η`: Reference coordinates within element [-1,1]²

# Description
Uses bounding box checks followed by inverse mapping for accuracy.
"""
function find_containing_element(x, y, mesh)
    # TODO: Implement
    # 1. Quick bounding box rejection
    # 2. For candidate elements, compute inverse mapping
    # 3. Check if (ξ, η) ∈ [-1-tol, 1+tol]²
    # 4. Return element and reference coordinates
end
```

---

### 4. Helper Function: Evaluate Basis Functions
```julia
"""
    evaluate_lagrange_basis(ξ, η, basis)

Evaluate all Lagrange basis functions at reference coordinates (ξ, η).

# Arguments
- `ξ, η`: Reference coordinates in [-1, 1]²
- `basis`: Basis function information (nodes, degrees)

# Returns
- `φ`: Vector of basis function values at (ξ, η)

# Description
For tensor-product Lagrange polynomials:
φᵢⱼ(ξ, η) = Lᵢ(ξ) × Lⱼ(η)

where Lᵢ are 1D Lagrange polynomials.
"""
function evaluate_lagrange_basis(ξ, η, basis)
    # TODO: Implement
    # Evaluate tensor product of 1D Lagrange polynomials
    # Return flattened array of basis values
end
```

---

## Summary Table

| Step | Function | Alya Coords | Solution Available | Action |
|------|----------|-------------|-------------------|---------|
| 1 | `je_receive_alya_data` | Metadata only | ❌ | Receive grid info |
| 2 | `sem_setup` | - | ❌ | Build mesh |
| 3 | `build_coupling_communication_arrays` | 🔴 **Computed here** | ❌ | Count ownership |
| 4 | `MPI_Alltoall` | - | ❌ | Exchange pattern |
| 5 | `print_coupling_pattern` | - | ❌ | Verification |
| 6 | `build_alya_point_ownership_map` | 🔴 **Computed & stored** | ❌ | Verification only |
| 7 | `print_alya_point_ownership` | - | ❌ | Verification |
| 8 | `initialize` | - | ✅ **Available** | Initialize u[:] |
| **8.5 (NEW)** | `extract_local_alya_coordinates` | 🔴 **Store for coupling** | ✅ | **⭐ STORE COORDS** |
| **8.6 (NEW)** | `interpolate_solution_to_alya_coords` | ✅ | ✅ | **⭐ INTERPOLATE** |
| 9 | `CouplingData` constructor | ✅ | ✅ | Package data |
| 10 | `allocate_coupling_buffers` | ✅ | ✅ | Allocate buffers |

---

## Next Steps

### Immediate Tasks

1. **Implement `extract_local_alya_coordinates`**
   - Extract the coordinate loop from `build_coupling_communication_arrays`
   - Store coordinates instead of just counting
   - Return matrix of local Alya coordinates

2. **Implement `interpolate_solution_to_alya_coords`**
   - Implement element location (bounding box + inverse mapping)
   - Implement basis function evaluation
   - Implement interpolation formula
   - Handle edge cases and error conditions

3. **Add coordinates to `CouplingData` structure**
```julia
   mutable struct CouplingData
       # ... existing fields ...
       
       # NEW: Alya coordinate information
       alya_local_coords::Union{Nothing, Matrix{Float64}}
       alya_local_ids::Union{Nothing, Vector{Int32}}
       alya_owner_ranks::Union{Nothing, Vector{Int32}}
   end
```

4. **Integrate into `setup_coupling_and_mesh`**
   - Call new functions between Steps 8 and 9
   - Store coordinates in coupling structure
   - Verify interpolation works

5. **Implement time-stepping coupling callback**
   - Re-use stored coordinates
   - Interpolate at each time step
   - Pack into send buffers
   - Exchange via MPI

### Testing Strategy

1. **Unit Tests**
   - Test `find_containing_element` with known points
   - Test `evaluate_lagrange_basis` against analytical values
   - Test `interpolate_solution_to_alya_coords` with polynomial solutions

2. **Integration Tests**
   - Run coupled simulation with simple test case
   - Verify conservation properties
   - Check interpolation accuracy

3. **Performance Tests**
   - Profile interpolation overhead
   - Optimize element search if needed
   - Consider caching strategies

---

## References

- Jexpresso SEM solver documentation
- Alya coupling interface specification
- MPI documentation for collective operations

---

**Document Version History:**

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024 | Initial documentation |

---

**End of Document**