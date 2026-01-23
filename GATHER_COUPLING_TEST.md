# MPI_Gather-based Coupling Test

This directory contains test programs for coupling Julia and Fortran codes using `MPI_Gather` on a shared `MPI_COMM_WORLD`.

## Overview

This coupling approach is different from the intercommunicator method:
- **Both codes run in the same MPI_COMM_WORLD**
- **All ranks from both codes participate in collective operations together**
- **No intercommunicators are created**
- **Simpler but requires coordination on all collective operations**

## Files

- `alya_gather_coupling.f90` - Fortran code that gathers app names from all ranks
- `julia_couple_with_alya_gather.jl` - Julia code that participates in the gather

## How It Works

When launched with MPMD:
```bash
mpirun -np 2 ./alya_gather_coupling : -np 2 julia julia_couple_with_alya_gather.jl
```

The MPI world looks like:
```
MPI_COMM_WORLD (size = 4):
  Rank 0, 1: Fortran (app_name = "OTRO5")
  Rank 2, 3: Julia (app_name = "JEXPRESSO")
```

### The Gather Operation

Both codes call `MPI_Gather` to send their application name to rank 0:

**Fortran:**
```fortran
call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                app_dumm, 128, MPI_CHARACTER, &
                0, MPI_COMM_WORLD, ierr)
```

**Julia:**
```julia
if rank == 0
    MPI.Gather!(send_buf, recv_buf, 0, world)
else
    MPI.Gather!(send_buf, nothing, 0, world)
end
```

### What Rank 0 Receives

Rank 0 (Fortran) receives an array with names from all 4 ranks:
```
app_dumm[1] = "OTRO5"       (from Fortran rank 0)
app_dumm[2] = "OTRO5"       (from Fortran rank 1)
app_dumm[3] = "JEXPRESSO"   (from Julia rank 2)
app_dumm[4] = "JEXPRESSO"   (from Julia rank 3)
```

## Compilation and Execution

### Step 1: Compile the Fortran code

```bash
mpif90 -o alya_gather_coupling alya_gather_coupling.f90
```

### Step 2: Run the coupled test

```bash
mpirun --tag-output \
    -np 2 ./alya_gather_coupling : \
    -np 2 julia --project=. julia_couple_with_alya_gather.jl
```

### Alternative: Different rank distribution

```bash
# 3 Fortran ranks + 1 Julia rank
mpirun --tag-output \
    -np 3 ./alya_gather_coupling : \
    -np 1 julia --project=. julia_couple_with_alya_gather.jl
```

## Expected Output

```
[Fortran rank 0] World size: 4
[Fortran rank 1] World size: 4
[Julia rank 2] World size: 4
[Julia rank 3] World size: 4

[Fortran rank 0] Received application names from all ranks:
  Rank 0: OTRO5
  Rank 1: OTRO5
  Rank 2: JEXPRESSO
  Rank 3: JEXPRESSO

[Julia rank 2] Received application names from all ranks:
  Rank 0: OTRO5
  Rank 1: OTRO5
  Rank 2: JEXPRESSO
  Rank 3: JEXPRESSO

[Fortran rank 0] Done
[Fortran rank 1] Done
[Julia rank 2] Done
[Julia rank 3] Done
```

## Key Differences from Intercommunicator Approach

| Feature | Intercommunicator | Shared COMM_WORLD |
|---------|------------------|-------------------|
| Communicator | Separate local + inter comms | One shared MPI_COMM_WORLD |
| Rank numbering | Independent in each code | Continuous across codes |
| Collective ops | Independent within each code | All ranks must participate together |
| Complexity | Higher (multiple communicators) | Lower (single communicator) |
| Flexibility | High (codes can operate independently) | Low (tight coupling required) |
| Use case | Loosely coupled codes | Tightly coupled codes with shared operations |

## Important Notes

1. **All collective operations must include ALL ranks from both codes**
   - If Fortran does `MPI_Barrier(MPI_COMM_WORLD)`, Julia must also call it
   - If Fortran does `MPI_Gather`, Julia must participate

2. **MPI data types must match**
   - Fortran `MPI_CHARACTER` = Julia `Int8` (signed char)
   - Both use 128-byte buffers

3. **No APPID needed**
   - Unlike intercommunicator approach, no environment variables required
   - Rank distribution is automatic based on MPMD launch order

4. **Synchronization is critical**
   - Both codes must reach collective operations in the same order
   - Deadlock will occur if one code skips a collective operation

## Adapting to Real Alya Code

If your Alya code uses this pattern:
1. Julia must participate in EVERY collective operation Alya performs on MPI_COMM_WORLD
2. Julia needs to know the sequence of collective operations Alya will perform
3. Consider using barriers between major steps to ensure synchronization
4. Data exchange can be done via:
   - `MPI_Gather` / `MPI_Scatter` (as shown here)
   - `MPI_Sendrecv` with known rank pairs
   - `MPI_Bcast` for broadcasting data

## Troubleshooting

**Deadlock during gather:**
- Check that both codes call MPI_Gather with the same root rank (0)
- Ensure both codes use MPI_COMM_WORLD (not local communicators)

**Wrong data received:**
- Verify buffer size is 128 bytes in both codes
- Check data type: Fortran MPI_CHARACTER vs Julia Int8

**Rank mismatch errors:**
- Total ranks must equal sum of ranks in MPMD launch
- Example: `-np 2 : -np 2` gives total 4 ranks (not 2)
