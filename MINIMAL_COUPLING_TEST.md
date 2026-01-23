# Minimal Coupling Test

This guide tests MPI coupling with minimal programs before trying full Jexpresso.

## Files Created

1. **`Jexpresso-mini.jl`** - Minimal Julia program that initializes coupling
2. **`alya_mini_coupled.f90`** - Modified Alya test with coupling support

## Step 1: Compile Alya Mini

Compile the Fortran test program:

```bash
cd /path/to/Jexpresso

# With MPI Fortran compiler
mpif90 -o alya_mini_coupled alya_mini_coupled.f90

# Or with gfortran + MPI
gfortran -o alya_mini_coupled alya_mini_coupled.f90 -lmpi
```

**Verify it compiles:**
```bash
ls -lh alya_mini_coupled
```

## Step 2: Test Standalone

Test each code separately first:

### Test Jexpresso-mini alone
```bash
mpirun -np 2 julia --project=. Jexpresso-mini.jl
```

**Expected output:**
```
Jexpresso-mini starting...
[Jexpresso-mini rank 0] World size: 2
[Jexpresso-mini rank 1] World size: 2
[Jexpresso-mini rank 0] Initializing coupling...
[ Info: [Jexpresso] Coupling initialized: local_size=2, world_rank=0, ...
[Jexpresso-mini rank 0] Coupling initialized!
...
```

### Test alya_mini_coupled alone
```bash
mpirun -np 2 ./alya_mini_coupled
```

**Expected output:**
```
[Alya-mini rank 0] World size: 2
[Alya-mini rank 1] World size: 2
[Alya-mini rank 0] Initializing coupling...
[Alya-mini rank 0] Coupling initialized!
  Local rank: 0 / 2
...
```

## Step 3: Test MPMD Coupling

Now test them together:

```bash
mpirun -np 2 julia --project=. Jexpresso-mini.jl : \
       -np 2 ./alya_mini_coupled
```

**Expected output:**
```
Jexpresso-mini starting...
[Alya-mini rank 2] World size: 4
[Jexpresso-mini rank 0] World size: 4
[Alya-mini rank 3] World size: 4
[Jexpresso-mini rank 1] World size: 4
[Jexpresso-mini rank 0] Initializing coupling...
[Alya-mini rank 2] Initializing coupling...
[ Info: [Jexpresso] Coupling initialized: local_size=2, world_rank=0, root_ranks=[0, 2]
[ Info: [Alya-mini] Coupling initialized: local_size=2, world_rank=2, root_ranks=[0, 2]
[Jexpresso-mini rank 0] Coupling initialized!
  Local rank: 0 / 2
  World rank: 0
  Is root: true
[Alya-mini rank 2] Coupling initialized!
  Local rank: 0 / 2
  World rank: 2
  Is root: true
[Jexpresso-mini] Testing name exchange...
[Alya-mini] Testing name exchange...
[Jexpresso-mini] Exchanged names successfully!
  My name: Jexpresso
  Partner name: Alya-mini
[Alya-mini] Exchanged names successfully!
  My name: Alya-mini
  Partner name: Jexpresso
[Jexpresso-mini] Step 1/5
[Alya-mini] Step 1/5
...
[Jexpresso-mini] Test complete!
[Alya-mini] Test complete!
```

**Success criteria:**
- ✅ Both codes initialize coupling
- ✅ Names exchanged successfully
- ✅ Both codes synchronize and run 5 steps
- ✅ Both codes finish cleanly

## Step 4: Modify Your Real Alya

Once the minimal test works, add coupling to your real Alya code.

### Changes needed in `unitt_alya_with_another_code.f`:

```fortran
program unitt_alya_with_another_code
#ifdef USEMPIF08
  use mpi_f08
  implicit none
#define MY_MPI_COMM      type(MPI_Comm)
#else
  implicit none
  include 'mpif.h'
#define MY_MPI_COMM      integer(4)
#endif

  MY_MPI_COMM            :: comm_world, comm_local     ! ADD comm_local
  integer(4)             :: ierr, rank, size
  integer(4)             :: local_rank, local_size     ! ADD local vars
  integer(4)             :: code_id, n_codes           ! ADD coupling vars
  character(128)         :: app_name
  character(128), allocatable :: app_dumm(:)

  ! Initialize MPI
  app_name = 'OTRO5'
  call MPI_Init(ierr)

  comm_world = MPI_COMM_WORLD                          ! SAVE world comm
  call MPI_COMM_RANK(comm_world, rank, ierr)
  call MPI_COMM_SIZE(comm_world, size, ierr)

  ! =================================================================
  ! ADD COUPLING INITIALIZATION
  ! =================================================================
  code_id = 2   ! Alya is code 2
  n_codes = 2   ! Total codes: Jexpresso + Alya

  ! Split communicator
  call MPI_COMM_SPLIT(comm_world, code_id, rank, comm_local, ierr)

  call MPI_COMM_RANK(comm_local, local_rank, ierr)
  call MPI_COMM_SIZE(comm_local, local_size, ierr)

  if (rank == size/2) then  ! First rank of Alya
     write(*,*) '[Alya] Coupling initialized: local_size=', local_size
  endif
  ! =================================================================

  ! Allocate on root of THIS CODE (use local_rank, not rank!)
  if (local_rank == 0) then              ! CHANGE: use local_rank
     allocate(app_dumm(local_size))      ! CHANGE: use local_size
  end if

  ! Use comm_local for internal Alya operations
  call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                  app_dumm, 128, MPI_CHARACTER, &
                  0, comm_local, ierr)         ! CHANGE: use comm_local

  ! Print result
  if (local_rank == 0) then              ! CHANGE: use local_rank
     print *, "Alya received from Rank 0: ", app_dumm(1)
     deallocate(app_dumm)
  end if

  call MPI_Finalize(ierr)

end program unitt_alya_with_another_code
```

**Key changes:**
1. ✅ Add `comm_local` and local rank variables
2. ✅ Call `MPI_COMM_SPLIT` with code_id=2
3. ✅ Use `comm_local` instead of `MPI_COMM_WORLD` for internal operations
4. ✅ Use `local_rank` and `local_size` instead of global rank/size

## Troubleshooting

### Problem: "mpif90: command not found"
**Solution:** Install MPI or use full path:
```bash
# Find MPI compiler
which mpifort
which mpif90

# Or install
brew install open-mpi  # macOS
apt install libopenmpi-dev  # Linux
```

### Problem: Hangs at initialization
**Cause:** One code not calling `MPI_COMM_SPLIT`

**Solution:** Ensure both codes split the communicator with matching code_id

### Problem: "Invalid communicator"
**Cause:** Using `MPI_COMM_WORLD` instead of `comm_local`

**Solution:** Replace all `MPI_COMM_WORLD` with `comm_local` for internal ops

### Problem: Wrong data in gather
**Cause:** Gathering from wrong communicator

**Solution:** Use `comm_local` for gathering within Alya

## Next Steps

Once minimal test works:

1. ✅ Verify name exchange works
2. ✅ Add data exchange (send/receive arrays)
3. ✅ Test with different rank counts
4. ✅ Integrate into full Jexpresso
5. ✅ Integrate into full Alya

## See Also

- `MPMD_COUPLING_GUIDE.md` - Complete MPMD guide
- `COUPLING_WITH_ALYA.md` - Full Alya coupling documentation
- `src/mpi/coupling_helpers.jl` - Helper functions for Jexpresso
