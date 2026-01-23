# Gather-Based Coupling Usage Guide

This guide explains how to use Jexpresso's gather-based coupling mode to couple with external codes (like Alya) that use `MPI_Gather` on a shared `MPI_COMM_WORLD`.

## Overview

**Gather-based coupling** is a simpler coupling approach where:
- Both codes share the same `MPI_COMM_WORLD`
- All ranks from both codes participate in collective operations together
- No separate communicators are created (unlike intercommunicator approach)
- Suitable for codes that cannot be modified and use shared collective operations

## When to Use This Mode

Use `--gather-coupling` when:
- ✅ Your external code uses `MPI_Gather`, `MPI_Comm_split`, or other collective ops on `MPI_COMM_WORLD`
- ✅ You cannot modify the external code to use intercommunicators
- ✅ The external code expects all ranks to participate in its collective operations
- ✅ You want a simpler coupling setup without separate communicators

Use `--coupling` (intercommunicator mode) when:
- ✅ You need independent execution of both codes with occasional synchronization
- ✅ You want separate communicators for each code
- ✅ The external code supports intercommunicator-based coupling

## Command Line Usage

### Basic Syntax

```bash
julia --project=. src/Jexpresso.jl <equations> <case> <CI_mode> --gather-coupling [--code-name NAME]
```

### Parameters

- `--gather-coupling`: Enable gather-based coupling mode (required)
- `--code-name NAME`: Name to send during coupling (default: "Jexpresso")
- `--coupling-test-only`: Exit after coupling initialization without running simulation (optional, for testing)

### Example 1: Testing Coupling Only (Alya unit test)

When your external code exits immediately after the gather (e.g., unit test), use the **minimal coupling test**:

```bash
mpirun --tag-output \
    -np 2 ./alya/Alya.x : \
    -np 2 julia --project=. coupling_test_minimal.jl --code-name "Jexpresso"
```

This lightweight script:
1. Participates in `MPI_Comm_split`
2. Sends application name via `MPI_Gather`
3. Exits cleanly without loading Jexpresso simulation code

**Why use this instead of `--coupling-test-only`?**
- Avoids loading heavy Jexpresso modules (plotting, physics solvers, etc.)
- Prevents package loading issues in MPI environment
- Much faster for quick coupling tests
- Use this when testing with external code unit tests that exit quickly

### Example 2: Testing with Partial Jexpresso Initialization

If you need to test that coupling works after Jexpresso initializes but before running simulation:

```bash
mpirun --tag-output \
    -np 2 ./alya/Alya.x : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --gather-coupling --coupling-test-only --code-name "Jexpresso"
```

**Note**: This loads all Jexpresso modules (slower) but exits before running simulation. Use Example 1 (minimal test) unless you specifically need to verify module loading.

### Example 3: Full Coupled Simulation

When your external code continues running (full Alya simulation):

```bash
mpirun --tag-output \
    -np 2 ./alya/Alya.x : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false \
        --gather-coupling --code-name "Jexpresso"
```

This launches:
- Ranks 0-1: Alya (runs full simulation)
- Ranks 2-3: Jexpresso (runs full simulation)

**Important**: Both codes must call `MPI_Finalize` at approximately the same time to avoid hanging.

## What Happens During Initialization

When `--gather-coupling` is enabled, Jexpresso:

1. **Participates in `MPI_Comm_split`** with `MPI_UNDEFINED` color
   - Creates a NULL communicator (as external code may expect)
   - Ensures collective operation completes for all ranks

2. **Participates in `MPI_Gather`** to send application name
   - Sends 128-byte application name to global rank 0
   - Global rank 0 (typically external code) receives names from all ranks

3. **Continues with normal execution** on full `MPI_COMM_WORLD`
   - Uses all assigned ranks for computation
   - Participates in any subsequent collective operations

## Rank Distribution in MPMD Mode

When using MPMD launch with `:`, ranks are assigned sequentially:

```
mpirun -np 2 ./prog1 : -np 2 ./prog2 : -np 3 ./prog3

Rank distribution:
  Ranks 0-1: prog1
  Ranks 2-3: prog2
  Ranks 4-6: prog3
```

**Important**: The first program always gets rank 0 (global root).

## Example: Coupling with Alya

### Alya Code Pattern

Your Alya code does:
```fortran
! Initialize MPI
call MPI_Init(ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

! Split communicator (creates NULL for all with MPI_UNDEFINED)
call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)

! Gather application names
call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                app_dumm, 128, MPI_CHARACTER, &
                0, MPI_COMM_WORLD, ierr)

! Continue with simulation...
```

### Jexpresso Participation

When you run:
```bash
mpirun -np 2 ./Alya.x : -np 2 julia src/Jexpresso.jl ... --gather-coupling
```

Jexpresso automatically:
- Participates in the `MPI_Comm_split` (line 1 of Alya code)
- Participates in the `MPI_Gather` (line 2 of Alya code)
- Sends "Jexpresso" (or custom `--code-name`) to Alya's rank 0

### Expected Output

```
[Alya rank 0] World size: 4
[Alya rank 1] World size: 4
[Jexpresso rank 2] World size: 4
[Jexpresso rank 3] World size: 4
[ Info: Initializing gather-based MPI coupling mode
[ Info: Participated in MPI_Comm_split with MPI_UNDEFINED
[ Info: Gather-based coupling initialization complete (sent 'Jexpresso' to global root)
[Alya rank 0] Received from Rank 1: ALYA
[Alya rank 0] Received from Rank 2: JEXPRESSO
[Alya rank 0] Received from Rank 3: JEXPRESSO
```

## Testing with Standalone Program

Before coupling with your full Alya code, test with the provided test program:

```bash
# Compile test program
mpif90 -cpp -DUSEMPIF08 alya_gather_coupling.f90 -o alya_test

# Run test
mpirun --tag-output \
    -np 2 ./alya_test : \
    -np 2 julia --project=. src/Jexpresso.jl CompEuler wave1d false --gather-coupling
```

## Troubleshooting

### Hang after coupling initialization

**Symptom**: Jexpresso completes field initialization but then hangs

**Cause**: External code called `MPI_Finalize` and exited while Jexpresso continues running

**Solution**: Use `--coupling-test-only` flag:
```bash
mpirun -np 2 ./alya/Alya.x : \
    -np 2 julia src/Jexpresso.jl ... --gather-coupling --coupling-test-only
```

This makes Jexpresso exit after coupling initialization, matching the external code's behavior.

### Deadlock at MPI_Comm_split

**Symptom**: Program hangs, no output after MPI initialization

**Cause**: Jexpresso is not participating in the `MPI_Comm_split` collective operation

**Solution**: This should be automatic with `--gather-coupling`. If still hanging:
- Check that external code calls `MPI_Comm_split` with `MPI_UNDEFINED`
- Verify both codes are using the same MPI_COMM_WORLD

### Wrong data received by external code

**Symptom**: External code receives garbled application names

**Cause**: Data type or buffer size mismatch

**Solution**:
- Jexpresso sends 128 bytes of `Int8` (MPI_CHARACTER equivalent)
- External code must expect 128 bytes per rank
- Use `--code-name` to customize the sent name

### Jexpresso uses wrong number of ranks

**Symptom**: Jexpresso reports wrong number of ranks or errors about mesh partitioning

**Cause**: In gather-coupling mode, Jexpresso uses ALL ranks in MPI_COMM_WORLD

**Solution**: This is expected behavior. If you need Jexpresso to use only some ranks:
- Launch with correct rank count: `mpirun -np 2 alya : -np 4 jexpresso` gives Jexpresso 4 ranks
- Or use `--coupling` (intercommunicator mode) instead for independent communicators

### Cannot use both coupling modes

**Symptom**: Error "Cannot enable both --coupling and --gather-coupling simultaneously"

**Cause**: Both flags were specified

**Solution**: Choose one:
- `--gather-coupling` for shared COMM_WORLD approach
- `--coupling` for intercommunicator approach

## Comparison: Gather vs Intercommunicator Coupling

| Feature | `--gather-coupling` | `--coupling` |
|---------|-------------------|--------------|
| Communicator | Shared MPI_COMM_WORLD | Separate local + inter |
| Setup complexity | Simple | More complex |
| Code modification | None (external code) | May require changes |
| Collective ops | All ranks participate | Independent per code |
| Use case | Unmodifiable external codes | Flexible coupling scenarios |
| Synchronization | Tight coupling | Loose coupling |

## Additional Collective Operations

If your external code performs additional collective operations (barriers, broadcasts, reductions) on `MPI_COMM_WORLD`, Jexpresso must also participate in them.

**Example: Adding a barrier**

If Alya does:
```fortran
call MPI_Barrier(MPI_COMM_WORLD, ierr)
```

You need to add to Jexpresso's main loop:
```julia
MPI.Barrier(Jexpresso.comm)
```

This ensures both codes synchronize at the barrier.

## Reference Files

- **Implementation**: `src/run.jl` (lines ~100-160)
- **Test program**: `alya_gather_coupling.f90`
- **Standalone Julia test**: `julia_couple_with_alya_gather.jl`
- **Test script**: `test_jexpresso_alya_gather.sh`

## Support and Issues

If you encounter issues with gather-based coupling:

1. Test with the standalone programs first (`alya_gather_coupling.f90` + `julia_couple_with_alya_gather.jl`)
2. Check that your external code uses the expected collective operation pattern
3. Verify MPI ranks are distributed as expected in MPMD mode
4. Report issues with detailed MPI output using `--tag-output`
