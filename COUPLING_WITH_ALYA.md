# Coupling Jexpresso with Alya

This guide explains how to couple Jexpresso with Alya using MPI.

## Prerequisites

- Both Jexpresso and Alya must use the same MPI installation
- Verify MPI compatibility:
  ```bash
  # Check Jexpresso's MPI
  julia --project=. -e 'using MPI; MPI.versioninfo()'

  # Check Alya's MPI
  ldd /path/to/alya/bin/alya.x | grep mpi
  ```

## Launch Methods

### Method 1: MPMD Mode (Recommended)

Launch both codes in one `mpirun` command:

```bash
mpirun -np 4 julia --project=. run_jexpresso.jl \
    --coupling \
    --code-id 1 \
    --n-codes 2 \
    --code-name "Jexpresso" \
    CompEuler theta : \
    -np 4 /path/to/alya/bin/alya.x alya_input.dat
```

**Requirements:**
- Alya must support coupling via `MPI_COMM_WORLD` splitting
- Alya must be able to identify its code_id and n_codes (either via command line or config file)

**Advantages:**
- Clean separation of executables
- Each code runs independently
- Easy to debug

### Method 2: Custom Launcher Script

Use the provided `launch_jexpresso_alya.jl` script:

```bash
# Edit the script first to configure paths and ranks
vim launch_jexpresso_alya.jl

# Then launch
mpirun -np 8 julia --project=. launch_jexpresso_alya.jl
```

**Edit these sections in the launcher:**

1. **Rank distribution** (lines ~20-21):
   ```julia
   JEXPRESSO_RANKS = 4  # Number of ranks for Jexpresso
   ALYA_RANKS = 4       # Number of ranks for Alya
   ```

2. **Jexpresso problem** (lines ~24-25):
   ```julia
   JEXPRESSO_EQS = "CompEuler"   # Your equation type
   JEXPRESSO_CASE = "theta"      # Your test case
   ```

3. **Alya configuration** (lines ~28-29):
   ```julia
   ALYA_INPUT_FILE = "alya.dat"
   ALYA_EXECUTABLE = "/path/to/alya/bin/alya.x"
   ```

4. **Alya launch implementation** (lines ~90-130):
   - Currently has placeholder code
   - You need to implement how Alya is launched
   - See options below

## Implementing Alya Launch

### Option A: Alya has Julia Interface

If Alya provides a Julia wrapper or interface:

```julia
# In launch_jexpresso_alya.jl, in the Alya section:
include("/path/to/alya/julia_interface.jl")

run_alya_coupled(
    input_file = ALYA_INPUT_FILE,
    code_id = 2,
    n_codes = 2,
    comm = MPI.COMM_WORLD
)
```

### Option B: Call Alya via System Command

If Alya is a standalone executable:

```julia
# This is tricky because we're already in MPI context
# Alya needs to use the existing MPI ranks, not spawn new ones

# Check if Alya supports being called as a library
# Otherwise, you may need to use the MPMD method instead
```

### Option C: Alya Fortran Interface via ccall

If you have access to Alya's source and can create a wrapper:

```julia
# Create a Fortran wrapper in Alya:
# subroutine alya_main_coupled(code_id, n_codes)
#     integer, intent(in) :: code_id, n_codes
#     ! Initialize Alya in coupling mode
#     ! Run Alya solver
# end subroutine

# Then call from Julia:
ccall((:alya_main_coupled_, "libalya.so"), Cvoid,
      (Ref{Int32}, Ref{Int32}),
      Int32(2), Int32(2))
```

## Coupling Interface in Jexpresso

### Step 1: Identify Coupling Points

Where does Jexpresso need to exchange data with Alya?
- Boundary conditions? (most common)
- Initial conditions?
- Every timestep?
- Every N timesteps?

### Step 2: Add Coupling Code

**In your boundary condition file** (e.g., `problems/CompEuler/theta/user_bc.jl`):

```julia
function user_bc!(...)
    # ... existing BC code ...

    # Add coupled BC if coupling is active
    if !isnothing(Jexpresso.coupling_ctx)
        apply_alya_coupling!(mesh, state, t)
    end
end

function apply_alya_coupling!(mesh, state, t)
    using Jexpresso.JexpressoCoupling

    ctx = Jexpresso.coupling_ctx
    alya_code_id = 2

    # Extract boundary data
    local_velocity = extract_boundary_velocity(mesh, state)

    # Gather to coupling root
    if is_coupling_root(ctx)
        global_velocity = zeros(total_boundary_points)
    else
        global_velocity = zeros(0)
    end
    gather_to_coupling_root!(ctx, local_velocity, global_velocity)

    # Exchange with Alya
    if is_coupling_root(ctx)
        global_pressure = zeros(total_boundary_points)

        # Send velocity (tag 200), receive pressure (tag 201)
        coupling_exchange_workflow!(ctx, global_velocity, global_pressure, alya_code_id;
                                    tag_send=200, tag_recv=201)
    else
        global_pressure = zeros(0)
    end

    # Distribute to all Jexpresso ranks
    local_pressure = zeros(length(local_velocity))
    scatter_from_coupling_root!(ctx, global_pressure, local_pressure)

    # Apply pressure BC
    apply_pressure_boundary!(mesh, state, local_pressure)
end
```

## Coordination with Alya Team

You need to agree on these parameters with whoever maintains Alya:

### 1. Code IDs
- Jexpresso: `code_id = 1`
- Alya: `code_id = 2`

### 2. MPI Tags
Document what tags are used for each data type:

| Data Type | Direction | MPI Tag |
|-----------|-----------|---------|
| Velocity  | Jexpresso → Alya | 200 |
| Pressure  | Alya → Jexpresso | 201 |
| Temperature | Jexpresso → Alya | 202 |
| ... | ... | ... |

### 3. Data Format
- Array sizes (number of interface points)
- Data type (Float64? Float32?)
- Array layout (contiguous? strided?)
- Units (SI? CGS? Non-dimensional?)

### 4. Communication Pattern
- Who sends first?
- Blocking or non-blocking?
- Synchronization points?

### 5. Initialization
- Exchange application names?
- Exchange mesh information?
- Verify interface compatibility?

## Testing

### Phase 1: Verify Both Codes Start

```bash
# Test that launcher works and both codes initialize
mpirun -np 8 julia --project=. launch_jexpresso_alya.jl
```

Look for:
- Both codes print startup messages
- No MPI errors
- Both codes initialize coupling correctly

### Phase 2: Test Name Exchange

Add to Jexpresso's initialization:

```julia
if !isnothing(Jexpresso.coupling_ctx)
    ctx = Jexpresso.coupling_ctx
    partner_name = exchange_app_names!(ctx, "Jexpresso", 2)
    @info "Coupled with: $partner_name"
end
```

Verify:
- Jexpresso receives "Alya" (or Alya's name)
- Alya receives "Jexpresso"

### Phase 3: Test Dummy Data Exchange

Before implementing real data extraction, test with dummy arrays:

```julia
if is_coupling_root(ctx)
    dummy_send = rand(100)
    dummy_recv = zeros(100)
    coupling_exchange_workflow!(ctx, dummy_send, dummy_recv, 2)
    @info "Exchange successful: recv mean=$(mean(dummy_recv))"
end
```

### Phase 4: Real Coupling

Once dummy exchange works, integrate with actual simulation data.

## Troubleshooting

### Problem: "Unable to find entry point"
**Cause:** Alya executable not found or wrong path
**Solution:** Check `ALYA_EXECUTABLE` path in launcher

### Problem: Deadlock (program hangs)
**Cause:** Mismatched send/recv order between Jexpresso and Alya
**Solution:**
- Ensure both codes call exchange in same order
- Use `coupling_exchange_workflow!` which handles ordering internally
- Check MPI tags match

### Problem: "Incorrect number of codes"
**Cause:** Alya didn't initialize coupling or has wrong n_codes
**Solution:** Verify Alya initializes with `n_codes=2`

### Problem: Wrong data received
**Cause:** Buffer size mismatch
**Solution:**
- Print array sizes on both sides
- Ensure both allocate same-sized buffers
- Check data types match (Float64 vs Float32)

### Problem: File I/O conflicts
**Cause:** Both codes writing to same output directory
**Solution:** Set separate output directories:
```julia
# In Jexpresso's user_inputs.jl
if !isnothing(Jexpresso.coupling_ctx)
    OUTPUT_DIR = "output_jexpresso"
else
    OUTPUT_DIR = "output"
end
```

## Example: Minimal Alya Interface

If Alya needs a coupling interface, here's a minimal Fortran example:

```fortran
module alya_coupling
    use mpi
    implicit none

    integer :: code_id, n_codes
    integer :: comm_world, comm_local, comm_inter
    integer :: world_rank, local_rank

contains

    subroutine alya_coupling_init(cid, nc)
        integer, intent(in) :: cid, nc
        integer :: ierr

        code_id = cid
        n_codes = nc
        comm_world = MPI_COMM_WORLD

        ! Get world rank
        call MPI_Comm_rank(comm_world, world_rank, ierr)

        ! Split communicator (similar to JexpressoCoupling)
        call MPI_Comm_split(comm_world, code_id, world_rank, comm_local, ierr)
        call MPI_Comm_rank(comm_local, local_rank, ierr)

        print *, '[Alya] Coupling initialized: code_id=', code_id, &
                 'local_rank=', local_rank
    end subroutine

    subroutine send_to_jexpresso(data, n, tag)
        real(8), intent(in) :: data(n)
        integer, intent(in) :: n, tag
        integer :: jexpresso_root, ierr

        ! Only rank 0 sends (coupling root)
        if (local_rank == 0) then
            jexpresso_root = 0  ! World rank of Jexpresso's root
            call MPI_Send(data, n, MPI_DOUBLE_PRECISION, &
                         jexpresso_root, tag, comm_world, ierr)
        endif
    end subroutine

    subroutine recv_from_jexpresso(data, n, tag)
        real(8), intent(out) :: data(n)
        integer, intent(in) :: n, tag
        integer :: jexpresso_root, ierr

        if (local_rank == 0) then
            jexpresso_root = 0
            call MPI_Recv(data, n, MPI_DOUBLE_PRECISION, &
                         jexpresso_root, tag, comm_world, &
                         MPI_STATUS_IGNORE, ierr)
        endif

        ! Broadcast to all Alya ranks
        call MPI_Bcast(data, n, MPI_DOUBLE_PRECISION, 0, comm_local, ierr)
    end subroutine

end module alya_coupling
```

## Getting Help

See also:
- `COUPLING_WITH_EXTERNAL_CODE.md` - General external code coupling
- `COUPLING_USAGE.md` - API reference
- `COUPLING_TESTS.md` - Testing procedures
- `src/mpi/coupling_helpers.jl` - Helper functions

For Alya-specific questions, coordinate with the Alya development team.
