# Intercommunicator Coupling for Jexpresso

This guide explains how to use Jexpresso with intercommunicator-based coupling for truly independent execution with external codes like Alya.

## Overview

**Intercommunicator coupling** provides:
- ✅ **True independence**: Each code has its own communicator for internal MPI operations
- ✅ **No blocking**: Codes run independently without synchronizing every step
- ✅ **Optional communication**: Codes communicate only when needed via intercommunicator
- ✅ **Clean separation**: Internal collective operations don't interfere between codes

## When to Use This Mode

Use `--intercomm-coupling` when:
- ✅ You need truly independent execution (different time steps, different termination criteria)
- ✅ You want to communicate occasionally without tight synchronization
- ✅ Your external code runs its own simulation with its own MPI operations
- ✅ You need maximum flexibility in coupling

**Comparison with other modes:**

| Mode | Independence | Synchronization | Use Case |
|------|-------------|----------------|----------|
| `--intercomm-coupling` | Full | Optional (point-to-point) | Production coupling |
| `--coupling` (JexpressoCoupling) | Moderate | Manual (comm_split) | Jexpresso-to-Jexpresso |
| `--gather-coupling` | None | All collective ops | Testing only |

## Requirements

### Environment Variable

Both codes must set `APPID` environment variable:
- First code (e.g., Alya): `APPID=0`
- Second code (e.g., Jexpresso): `APPID=1`

### Fortran Code Structure

Your Fortran code must follow this pattern:

```fortran
program your_code
  use mpi
  implicit none

  integer :: ierr, world_rank, world_size, appid
  integer :: local_comm, local_rank, local_size
  integer :: inter_comm, remote_leader_world
  character(len=256) :: appid_str
  integer, allocatable :: world_appids(:)
  integer :: my_appid

  ! Initialize MPI
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)

  ! Read APPID from environment
  call get_environment_variable("APPID", appid_str)
  read(appid_str, *) appid

  ! Split by APPID to create local communicator
  call MPI_Comm_split(MPI_COMM_WORLD, appid, world_rank, local_comm, ierr)
  call MPI_Comm_rank(local_comm, local_rank, ierr)
  call MPI_Comm_size(local_comm, local_size, ierr)

  ! Discover remote leader
  allocate(world_appids(world_size))
  my_appid = appid
  call MPI_Allgather(my_appid, 1, MPI_INTEGER, &
                     world_appids, 1, MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

  remote_leader_world = -1
  do i = 1, world_size
    if (world_appids(i) /= appid) then
      remote_leader_world = i - 1
      exit
    end if
  end do

  ! Create intercommunicator
  call MPI_Intercomm_create(local_comm, 0, MPI_COMM_WORLD, &
                            remote_leader_world, 12345, &
                            inter_comm, ierr)

  ! Exchange names (for verification)
  if (local_rank == 0) then
    ! Fortran sends tag 100, receives tag 101
    call MPI_Sendrecv(send_name, 128, MPI_CHARACTER, 0, 100, &
                      recv_name, 128, MPI_CHARACTER, 0, 101, &
                      inter_comm, MPI_STATUS_IGNORE, ierr)
  end if

  !============================================================
  ! YOUR SIMULATION CODE HERE
  ! Use local_comm for internal operations (barriers, reductions, etc.)
  ! Use inter_comm for communicating with Jexpresso (when needed)
  !============================================================

  do time_step = 1, max_steps
    ! Your simulation code
    call MPI_Barrier(local_comm, ierr)  ! Only synchronizes within Fortran ranks

    ! Optional: Send/receive data to/from Jexpresso
    if (local_rank == 0 .and. mod(time_step, exchange_frequency) == 0) then
      call MPI_Send(data, size, MPI_DOUBLE_PRECISION, 0, tag, inter_comm, ierr)
      call MPI_Recv(data, size, MPI_DOUBLE_PRECISION, 0, tag+1, inter_comm, &
                    MPI_STATUS_IGNORE, ierr)
    end if
  end do

  ! Cleanup
  call MPI_Comm_free(inter_comm, ierr)
  call MPI_Comm_free(local_comm, ierr)
  deallocate(world_appids)
  call MPI_Finalize(ierr)

end program your_code
```

## Usage with Jexpresso

### Basic Syntax

```bash
mpirun --tag-output \
    -np N1 -x APPID=0 ./your_fortran_code : \
    -np N2 -x APPID=1 julia --project=. src/Jexpresso.jl <equations> <case> <CI_mode> \
        --intercomm-coupling --code-name "Jexpresso"
```

### Example: Coupling with Alya

```bash
# 4 Alya ranks + 4 Jexpresso ranks
mpirun --tag-output \
    -np 4 -x APPID=0 ./alya/Alya.x : \
    -np 4 -x APPID=1 julia --project=. src/Jexpresso.jl CompEuler theta false \
        --intercomm-coupling --code-name "Jexpresso"
```

### What Happens

1. **MPI initialization**: All ranks start in `MPI_COMM_WORLD` (size = 8)
2. **Comm split**:
   - Ranks 0-3 (APPID=0, Alya): Create `local_comm_alya` (size = 4)
   - Ranks 4-7 (APPID=1, Jexpresso): Create `local_comm_jexpresso` (size = 4)
3. **Intercommunicator creation**: Both codes create `inter_comm` for cross-talk
4. **Independent execution**:
   - Alya uses `local_comm_alya` for all internal MPI operations
   - Jexpresso uses `local_comm_jexpresso` for all internal MPI operations
   - No interference - codes run independently!
5. **Optional communication**: Codes can exchange data via `inter_comm` when needed

## Testing

### Step 1: Test with Minimal Versions

First test that intercommunicator coupling works:

```bash
# Compile Fortran test
mpif90 alya_mini_coupler.f90 -o alya_mini_coupler

# Run test
mpirun --tag-output \
    -np 2 -x APPID=0 ./alya_mini_coupler : \
    -np 2 -x APPID=1 julia --project=. Jexpresso-mini-coupled.jl
```

**Expected output:**
```
[Alya-mini] world_size=4 appid=0 remote_leader_world=2
[Jexpresso] world_size=4 appid=1 remote_leader_world=0
[Alya-mini root] Partner name = JULIA
[Jexpresso root] Partner name = ALYA
[Alya-mini] Step 1/5
[Jexpresso] Step 1/5
...
[Both codes complete 5 steps and exit cleanly]
```

### Step 2: Test with Full Jexpresso

Once minimal test works:

```bash
mpirun --tag-output \
    -np 2 -x APPID=0 ./alya_mini_coupler : \
    -np 2 -x APPID=1 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --intercomm-coupling --code-name "Jexpresso"
```

### Step 3: Production Coupling

With real Alya simulation:

```bash
mpirun --tag-output \
    -np 4 -x APPID=0 ./alya/Alya.x : \
    -np 4 -x APPID=1 julia --project=. src/Jexpresso.jl CompEuler theta false \
        --intercomm-coupling --code-name "Jexpresso"
```

## Accessing the Intercommunicator in Jexpresso Code

After initialization, you can access the intercommunicator context:

```julia
# In user code (e.g., user_source.jl or custom coupling routine)
if !isnothing(Jexpresso.intercomm_ctx)
    inter_comm = Jexpresso.intercomm_ctx.inter_comm
    local_rank = Jexpresso.intercomm_ctx.local_rank

    # Example: Send data to Fortran root (rank 0 in inter_comm)
    if local_rank == 0 && time_step % exchange_frequency == 0
        data = [1.0, 2.0, 3.0]  # Your coupling data
        MPI.Send(data, 0, 200, inter_comm)  # tag=200

        # Receive response
        recv_data = Vector{Float64}(undef, 3)
        MPI.Recv!(recv_data, 0, 201, inter_comm)  # tag=201
    end
end
```

## Communication Patterns

### Point-to-Point (Root-to-Root)

Most common pattern - only root ranks communicate:

**Fortran (local_rank=0):**
```fortran
if (local_rank == 0) then
  call MPI_Send(data, size, MPI_DOUBLE_PRECISION, 0, 200, inter_comm, ierr)
  call MPI_Recv(data, size, MPI_DOUBLE_PRECISION, 0, 201, inter_comm, &
                MPI_STATUS_IGNORE, ierr)
end if
```

**Julia (local_rank=0):**
```julia
if local_rank == 0
    MPI.Send(data, 0, 201, inter_comm)  # Note: opposite tag
    MPI.Recv!(data, 0, 200, inter_comm)
end if
```

### Collective on Intercommunicator

Synchronize all ranks of both codes:

```julia
# Barrier across intercommunicator (all ranks participate)
MPI.Barrier(inter_comm)
```

**Note**: This synchronizes all Alya ranks with all Jexpresso ranks. Use sparingly.

### Asynchronous Exchange

Non-blocking communication for maximum independence:

```julia
if local_rank == 0
    req = MPI.Isend(send_data, 0, 300, inter_comm)
    # Continue with computation...
    MPI.Wait(req)
end
```

## Tag Conventions

**IMPORTANT**: Tags must be opposite between codes:

| Operation | Fortran Sends | Fortran Receives | Julia Sends | Julia Receives |
|-----------|---------------|------------------|-------------|----------------|
| Name exchange | 100 | 101 | 101 | 100 |
| Velocity field | 200 | 201 | 201 | 200 |
| Pressure field | 202 | 203 | 203 | 202 |

## Troubleshooting

### ERROR: APPID not set

**Symptom**: `APPID environment variable not set`

**Solution**: Use `-x APPID=N` in mpirun:
```bash
mpirun -np 2 -x APPID=0 ./code1 : -np 2 -x APPID=1 ./code2
```

### Deadlock during name exchange

**Symptom**: Codes hang after creating intercommunicator

**Cause**: Tag mismatch in initial exchange

**Solution**: Verify tags are opposite:
- Fortran: `MPI_Sendrecv(..., 100, ..., 101, ...)`
- Julia: `MPI.Sendrecv!(..., 101, ..., 100, ...)`

### Jexpresso hangs during simulation

**Symptom**: Codes exchange names successfully but Jexpresso hangs later

**Cause**: Jexpresso is calling collective operations on wrong communicator

**Solution**: Jexpresso should use `comm` (which is set to `local_comm`) for all internal operations. The intercommunicator is stored in `intercomm_ctx.inter_comm` and should only be used for cross-code communication.

### Cannot enable multiple coupling modes

**Symptom**: Error about multiple coupling modes

**Solution**: Use only one coupling flag:
- `--intercomm-coupling` (recommended for independence)
- `--coupling` (JexpressoCoupling, for Jexpresso-to-Jexpresso)
- `--gather-coupling` (for testing only)

## Advantages Over Gather-Based Coupling

| Feature | Intercommunicator | Gather-based |
|---------|------------------|--------------|
| Independence | ✅ Full | ❌ None |
| Internal MPI ops | ✅ No interference | ❌ Must match |
| Time stepping | ✅ Independent | ❌ Must synchronize |
| Termination | ✅ Independent criteria | ❌ Must finish together |
| Complexity | Moderate | Simple |
| Best for | Production coupling | Quick tests only |

## Reference Files

- **Module**: `src/mpi/IntercommCoupling.jl`
- **Integration**: `src/run.jl` (search for `intercomm_coupling`)
- **Fortran template**: `alya_mini_coupler.f90`
- **Julia template**: `Jexpresso-mini-coupled.jl`
- **Working test**: `INTERCOMMUNICATOR_TEST.md`

## Summary

Intercommunicator coupling is the **recommended approach** for production coupling where:
1. Codes run independently with their own time steps and termination criteria
2. Communication happens only when needed (not every step)
3. Internal MPI operations don't interfere between codes

Use it with: `--intercomm-coupling` and `APPID` environment variable.
