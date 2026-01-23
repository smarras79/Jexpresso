# Intercommunicator Coupling Test

This guide tests MPI intercommunicator coupling with minimal programs.

## Files

1. **`alya_mini_coupler.f90`** - Fortran program using MPI_Intercomm_create
2. **`Jexpresso-mini-coupled.jl`** - Julia program using MPI_Intercomm_create

Both use APPID environment variables to identify applications and dynamically discover remote leaders.

## Compile Fortran

```bash
# With MPI Fortran compiler
mpif90 -o alya_mini_coupler alya_mini_coupler.f90

# Or with gfortran + MPI
gfortran -o alya_mini_coupler alya_mini_coupler.f90 -lmpi

# If using MPI_F08 interface
mpif90 -DUSEMPIF08 -o alya_mini_coupler alya_mini_coupler.f90
```

## Run with MPMD

```bash
mpirun --tag-output \
    -np 2 -x APPID=0 ./alya_mini_coupler : \
    -np 2 -x APPID=1 julia --project=. Jexpresso-mini-coupled.jl
```

**Key points:**
- `--tag-output` prefixes output with rank numbers
- `-x APPID=0` sets environment variable for Fortran (APPID=0)
- `-x APPID=1` sets environment variable for Julia (APPID=1)
- `:` separates the two programs

## Expected Output

```
[1,0]<stdout>:[Alya-mini rank 0] World size: 4
[1,1]<stdout>:[Alya-mini rank 1] World size: 4
[1,2]<stdout>:Jexpresso-mini starting...
[1,2]<stdout>:[Jexpresso rank 2] World size: 4
[1,3]<stdout>:Jexpresso-mini starting...
[1,3]<stdout>:[Jexpresso rank 3] World size: 4
[1,0]<stdout>:[Alya-mini] world_size=4 appid=0 remote_leader_world=2
[1,2]<stdout>:[Jexpresso] world_size=4 appid=1 remote_leader_world=0
[1,0]<stdout>:[Alya-mini root] Partner name = JULIA
[1,2]<stdout>:[Jexpresso root] Partner name = ALYA
[1,0]<stdout>:[Alya-mini] Step 1/5
[1,2]<stdout>:[Jexpresso] Step 1/5
[1,0]<stdout>:[Alya-mini] Step 2/5
[1,2]<stdout>:[Jexpresso] Step 2/5
...
[1,0]<stdout>:[Alya-mini rank 0] Done
[1,2]<stdout>:[Jexpresso rank 2] Done
```

## What Was Fixed

### Original Problem: Tag Mismatch

**Both codes had:**
```
Send tag 100, Receive tag 101
```

This causes deadlock because neither receives what the other sends!

### Solution: Opposite Tags

**Fortran (alya_mini_coupler.f90):**
```fortran
call MPI_Sendrecv(snd, NAME_LEN, dt_schar, 0, 100, &  ! send tag 100
                  rcv, NAME_LEN, dt_schar, 0, 101, &  ! recv tag 101
```

**Julia (Jexpresso-mini-coupled.jl) - FIXED:**
```julia
MPI.Sendrecv!(send, 0, 101, recv, 0, 100, inter_comm)  # send tag 101, recv tag 100
```

Now:
- Fortran sends tag 100 → Julia receives tag 100 ✅
- Julia sends tag 101 → Fortran receives tag 101 ✅

## How It Works

1. **APPID Assignment**: Each application gets a unique ID via environment variable
   - Fortran: `APPID=0`
   - Julia: `APPID=1`

2. **Comm Splitting**: Each app splits MPI_COMM_WORLD by APPID
   - Fortran ranks form `local_comm` with APPID=0
   - Julia ranks form `local_comm` with APPID=1

3. **Dynamic Discovery**: All ranks do Allgather to share their APPID
   - Each app finds first rank of other app
   - This becomes `remote_leader_world`

4. **Intercommunicator**: Both create intercommunicator
   ```fortran
   MPI_Intercomm_create(local_comm, 0, world, remote_leader, TAG, inter_comm)
   ```

5. **Communication**: Use `inter_comm` for between-app communication
   - Rank 0 in each app's `local_comm` = rank 0 in `inter_comm`
   - Sendrecv, Barrier work on `inter_comm`

## Troubleshooting

### Problem: "APPID not set" error
**Solution**: Use `-x APPID=N` flag with mpirun

### Problem: "did not find remote leader"
**Cause**: Both apps have same APPID
**Solution**: Ensure Fortran has APPID=0, Julia has APPID=1

### Problem: Hangs at MPI_Intercomm_create
**Cause**: TAG mismatch or rank issues
**Solution**: Both must use same TAG (12345 in our case)

### Problem: Hangs at Sendrecv
**Cause**: Tag mismatch (original problem)
**Solution**: Use opposite tags as shown above

### Problem: "libmpi not found" in Julia
**Solution**: Check MPI.jl is using correct MPI
```bash
julia --project=. -e 'using MPI; MPI.versioninfo()'
```

## Differences from JexpressoCoupling Approach

This test uses **MPI intercommunicators** instead of the JexpressoCoupling module:

**Intercommunicator approach:**
- ✅ Standard MPI construct
- ✅ Explicit separation of applications
- ✅ Works with any MPI implementation
- ❌ More complex setup (requires APPID, manual discovery)
- ❌ Need to manage multiple communicators

**JexpressoCoupling approach (used in full Jexpresso):**
- ✅ Simpler API (initialize_coupling, exchange_field_data, etc.)
- ✅ Automatic root rank discovery
- ✅ Helper functions for common patterns
- ✅ Better for production coupling
- ❌ Requires both codes to use similar approach

For your real Alya coupling, you'll use the **JexpressoCoupling module** approach, but this intercommunicator test verifies the basic MPI infrastructure works.

## Next Steps

Once this test works:

1. ✅ Verify both codes can communicate via intercommunicators
2. ✅ Test with different rank counts
3. ✅ Move to full Jexpresso with JexpressoCoupling module
4. ✅ Integrate coupling into real Alya
5. ✅ Test actual physics data exchange

## See Also

- `MINIMAL_COUPLING_TEST.md` - JexpressoCoupling module approach
- `MPMD_COUPLING_GUIDE.md` - Complete MPMD guide
- `COUPLING_WITH_ALYA.md` - Full Alya integration guide
