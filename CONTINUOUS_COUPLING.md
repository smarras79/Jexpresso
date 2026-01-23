# Continuous Coupling: Running Both Codes in Parallel

This guide explains how to set up continuous execution where both Alya and Jexpresso continue running after the coupling initialization.

## Overview

For a **full coupled simulation** where both codes run together:
- Both codes participate in coupling initialization (`MPI_Comm_split` + `MPI_Gather`)
- Both codes continue with their simulation loops
- Both codes must synchronize at key points (barriers)
- Both codes must reach `MPI_Finalize` at approximately the same time

## Required Code Structure

### Fortran Side (Alya)

Your Fortran code must have this structure:

```fortran
! 1. MPI Initialization
call MPI_Init(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

! 2. Coupling Initialization
call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)
call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                app_dumm, 128, MPI_CHARACTER, &
                0, MPI_COMM_WORLD, ierr)

! 3. Simulation Loop - CRITICAL: Keep running!
do time_step = 1, max_steps
    ! Your simulation code here

    ! Synchronize with all ranks (including Julia)
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
end do

! 4. Finalization - Julia must also reach this
call MPI_Barrier(MPI_COMM_WORLD, ierr)  ! Final sync
call MPI_Finalize(ierr)
```

**Key Requirements:**
- ✅ Continue running after `MPI_Gather` (don't exit immediately)
- ✅ Use `MPI_Barrier(MPI_COMM_WORLD)` to synchronize with Julia ranks
- ✅ Both codes must reach `MPI_Finalize` together

### Julia Side (Jexpresso)

For full Jexpresso simulation, use:

```bash
mpirun -np 2 ./alya/Alya.x : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --gather-coupling --code-name "Jexpresso"
```

**Note**: This will encounter the dyld/package loading issue on some systems. See workaround below.

## Testing with Simplified Versions

### Step 1: Test with Simplified Continuous Versions

First, test that continuous coupling works using simplified versions:

```bash
# Compile Fortran continuous test
mpif90 -cpp -DUSEMPIF08 alya_gather_continuous.f90 -o alya_gather_continuous

# Run coupled test (10 seconds)
./test_continuous_coupling.sh
```

Or manually:
```bash
mpirun --tag-output \
    -np 2 ./alya_gather_continuous : \
    -np 2 julia --project=. jexpresso_gather_continuous.jl --code-name "Jexpresso"
```

**Expected Output:**
```
[Alya rank 0] Coupling initialization complete. Starting simulation...
[Jexpresso rank 0] Coupling initialization complete. Starting simulation...
[Alya rank 0] Time step 1, elapsed: 1.00s
[Jexpresso rank 0] Time step 1, elapsed: 1.00s
[Alya rank 0] Time step 2, elapsed: 2.00s
[Jexpresso rank 0] Time step 2, elapsed: 2.00s
...
[Alya rank 0] Simulation complete after 10 steps
[Jexpresso rank 0] Simulation complete after 10 steps
[Alya rank 0] Done
[Jexpresso rank 0] Done
```

### Step 2: Run with Full Jexpresso (if needed)

Once simplified versions work, try with full Jexpresso:

```bash
mpirun --tag-output \
    -np 2 ./alya_gather_continuous : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --gather-coupling --code-name "Jexpresso"
```

**If this crashes with dyld errors**, it means package loading conflicts with MPI. Solutions:
1. Use the simplified `jexpresso_gather_continuous.jl` instead (recommended for testing)
2. Fix package precompilation issues (advanced)
3. Use the intercommunicator approach instead (see `INTERCOMMUNICATOR_TEST.md`)

## Critical Synchronization Points

Both codes **must** call the same collective operations:

### 1. MPI_Comm_split
```fortran
! Fortran
call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)
```
```julia
# Julia
MPI.Comm_split(world, MPI_UNDEFINED, rank)
```

### 2. MPI_Gather
```fortran
! Fortran
call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                app_dumm, 128, MPI_CHARACTER, &
                0, MPI_COMM_WORLD, ierr)
```
```julia
# Julia
MPI.Gather!(send_buf, nothing, 0, world)
```

### 3. Barriers During Simulation (if used)
```fortran
! Fortran - in simulation loop
call MPI_Barrier(MPI_COMM_WORLD, ierr)
```
```julia
# Julia - must match each Fortran barrier
MPI.Barrier(world)
```

### 4. Final Barrier Before MPI_Finalize
```fortran
! Fortran - before finalize
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Finalize(ierr)
```
```julia
# Julia
MPI.Barrier(world)
MPI.Finalize()
```

## Modifying Your Alya Code

If your current Alya code looks like this (unit test version):

```fortran
call MPI_Gather(app_name, ...)
call MPI_Finalize(ierr)  ! ❌ Exits immediately - Julia hangs!
```

You need to change it to:

```fortran
call MPI_Gather(app_name, ...)

! Add simulation loop here
do time_step = 1, num_time_steps
    ! Your Alya simulation code

    ! Optionally synchronize with Julia
    if (mod(time_step, sync_frequency) == 0) then
        call MPI_Barrier(MPI_COMM_WORLD, ierr)
    end if
end do

! Finalize together
call MPI_Barrier(MPI_COMM_WORLD, ierr)
call MPI_Finalize(ierr)
```

## Comparison: Unit Test vs Continuous

| Feature | Unit Test Version | Continuous Version |
|---------|------------------|-------------------|
| Fortran behavior | Exits after gather | Continues running |
| Julia command | `coupling_test_minimal.jl` | `--gather-coupling` or `jexpresso_gather_continuous.jl` |
| Simulation | No simulation | Both run simulations |
| Barriers | Only during init | Multiple during simulation |
| Use case | Testing coupling setup | Production coupling |

## Troubleshooting

### Julia hangs after "Initialize fields"

**Cause**: Alya exited early (called `MPI_Finalize` too soon)

**Solution**: Ensure Alya continues running after the gather. See "Modifying Your Alya Code" above.

### Deadlock during simulation

**Symptom**: Both codes hang after some time steps

**Cause**: Mismatch in collective operations (one code calls barrier, other doesn't)

**Solution**: Ensure both codes call the same collective operations in the same order:
- Check all `MPI_Barrier` calls
- Check all `MPI_Bcast`, `MPI_Allreduce`, etc.
- Add logging to identify where they diverge

### Fortran completes but Julia still running

**Cause**: Julia simulation time longer than Fortran

**Solution**: Ensure both codes exit at the same time:
- Use the same number of time steps
- Or use time-based exit criteria
- Always call final barrier before `MPI_Finalize`

### dyld crash when loading Jexpresso

**Cause**: Package loading issues in MPI environment

**Solution**: Use simplified Julia version:
```bash
mpirun -np 2 ./alya : \
    -np 2 julia jexpresso_gather_continuous.jl
```

## Production Deployment

For production coupling with real Alya simulation:

1. **Ensure Alya continues running** after the gather
2. **Add synchronization points** where data exchange occurs
3. **Use same run duration** for both codes
4. **Test with simplified versions first** to verify coupling works
5. **Scale up gradually** (2+2 ranks → 4+4 ranks → production scale)

## Reference Files

- **Fortran continuous test**: `alya_gather_continuous.f90`
- **Julia continuous test**: `jexpresso_gather_continuous.jl`
- **Test script**: `test_continuous_coupling.sh`
- **Quick coupling test**: `coupling_test_minimal.jl` (exits after init)

## Example: Full Workflow

```bash
# 1. Test coupling only (quick test)
mpirun -np 2 ./alya/Alya.x : \
    -np 2 julia coupling_test_minimal.jl

# 2. Test continuous execution (simplified)
./test_continuous_coupling.sh

# 3. Run with real Alya (if it continues after gather)
mpirun -np 2 ./alya/Alya.x : \
    -np 2 julia jexpresso_gather_continuous.jl

# 4. Production: Full Jexpresso simulation
mpirun -np 4 ./alya/Alya.x : \
    -np 4 julia --project=. src/Jexpresso.jl CompEuler theta false \
        --gather-coupling --code-name "Jexpresso"
```
