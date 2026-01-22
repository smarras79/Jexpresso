# MPMD Mode Coupling Guide

This guide shows how to couple Jexpresso with Alya using MPMD (Multiple Program Multiple Data) mode.

## Testing Strategy

Test incrementally to ensure each step works before moving to the next.

### Step 1: Test Jexpresso Alone (n_codes=1)

Verify Jexpresso runs with coupling enabled but no partner:

```bash
mpirun -np 4 julia --project=. run_jexpresso.jl \
    --coupling \
    --code-id 1 \
    --n-codes 1 \
    CompEuler theta
```

**Expected behavior:**
- All ranks belong to code 1
- Coupling initializes successfully
- No data exchange (no partner exists)
- Simulation runs normally

**Success criteria:**
- ✅ No errors during coupling initialization
- ✅ Simulation completes successfully
- ✅ Output shows: `[ Info: [Jexpresso] Coupling initialized: local_size=4, ...`

### Step 2: Test MPMD with Alya Mimic (Julia)

Test MPMD launch with a Julia program that mimics Alya:

```bash
mpirun -np 4 julia --project=. run_jexpresso.jl \
    --coupling --code-id 1 --n-codes 2 CompEuler theta : \
    -np 4 julia --project=. test_alya_mimic.jl
```

**Expected behavior:**
- Ranks 0-3: Jexpresso
- Ranks 4-7: Alya mimic
- Both codes initialize coupling
- Both codes synchronize
- Optional data exchange test at step 3

**Success criteria:**
- ✅ Both codes start and initialize
- ✅ No deadlocks or hangs
- ✅ Both codes complete successfully
- ✅ If data exchange enabled, both codes send/receive

### Step 3: Launch with Real Alya

Once steps 1-2 work, launch with actual Alya:

```bash
mpirun -np 4 julia --project=. run_jexpresso.jl \
    --coupling --code-id 1 --n-codes 2 CompEuler theta : \
    -np 4 /path/to/alya/bin/alya.x alya_input.dat
```

**Prerequisites:**
- Alya must be compiled with same MPI library as Julia
- Alya must initialize coupling (see Fortran code below)
- Alya input file must specify coupling parameters

## What Alya Needs to Do

### Minimal Alya Coupling Interface

Add this to Alya's main program or initialization:

```fortran
module alya_jexpresso_coupling
    use mpi
    implicit none

    ! Coupling context
    integer :: code_id = 2
    integer :: n_codes = 2
    integer :: comm_world
    integer :: comm_local
    integer :: world_rank, world_size
    integer :: local_rank, local_size
    logical :: is_coupling_root
    integer :: jexpresso_root_rank = 0  ! World rank of Jexpresso's root

contains

    subroutine initialize_alya_coupling()
        integer :: ierr

        comm_world = MPI_COMM_WORLD

        ! Get world rank/size
        call MPI_Comm_rank(comm_world, world_rank, ierr)
        call MPI_Comm_size(comm_world, world_size, ierr)

        ! Split communicator by code_id
        call MPI_Comm_split(comm_world, code_id, world_rank, comm_local, ierr)

        ! Get local rank/size
        call MPI_Comm_rank(comm_local, local_rank, ierr)
        call MPI_Comm_size(comm_local, local_size, ierr)

        ! Only rank 0 of each code is a coupling root
        is_coupling_root = (local_rank == 0)

        if (world_rank == world_size/2) then
            write(*,*) '[Alya] Coupling initialized:'
            write(*,*) '  Local size:', local_size
            write(*,*) '  World rank:', world_rank
            write(*,*) '  Is root:', is_coupling_root
        endif

    end subroutine initialize_alya_coupling

    subroutine send_to_jexpresso(data, n, tag)
        real(8), intent(in) :: data(n)
        integer, intent(in) :: n, tag
        integer :: ierr

        if (is_coupling_root) then
            call MPI_Send(data, n, MPI_DOUBLE_PRECISION, &
                         jexpresso_root_rank, tag, comm_world, ierr)
        endif
    end subroutine send_to_jexpresso

    subroutine recv_from_jexpresso(data, n, tag)
        real(8), intent(out) :: data(n)
        integer, intent(in) :: n, tag
        integer :: ierr

        if (is_coupling_root) then
            call MPI_Recv(data, n, MPI_DOUBLE_PRECISION, &
                         jexpresso_root_rank, tag, comm_world, &
                         MPI_STATUS_IGNORE, ierr)
        endif

        ! Broadcast to all Alya ranks
        call MPI_Bcast(data, n, MPI_DOUBLE_PRECISION, 0, comm_local, ierr)
    end subroutine recv_from_jexpresso

    subroutine exchange_with_jexpresso(send_data, recv_data, n, send_tag, recv_tag)
        real(8), intent(in) :: send_data(n)
        real(8), intent(out) :: recv_data(n)
        integer, intent(in) :: n, send_tag, recv_tag
        integer :: ierr, req(2), status(MPI_STATUS_SIZE, 2)

        if (is_coupling_root) then
            ! Non-blocking send/recv to avoid deadlock
            call MPI_Isend(send_data, n, MPI_DOUBLE_PRECISION, &
                          jexpresso_root_rank, send_tag, comm_world, req(1), ierr)
            call MPI_Irecv(recv_data, n, MPI_DOUBLE_PRECISION, &
                          jexpresso_root_rank, recv_tag, comm_world, req(2), ierr)
            call MPI_Waitall(2, req, status, ierr)
        endif

        ! Broadcast received data to all Alya ranks
        call MPI_Bcast(recv_data, n, MPI_DOUBLE_PRECISION, 0, comm_local, ierr)
    end subroutine exchange_with_jexpresso

    subroutine synchronize_with_jexpresso()
        integer :: ierr
        call MPI_Barrier(comm_world, ierr)
    end subroutine synchronize_with_jexpresso

    subroutine finalize_alya_coupling()
        integer :: ierr

        ! Final synchronization
        call MPI_Barrier(comm_world, ierr)

        if (world_rank == world_size/2) then
            write(*,*) '[Alya] Coupling finalized'
        endif
    end subroutine finalize_alya_coupling

end module alya_jexpresso_coupling
```

### Using the Coupling Module in Alya

```fortran
program alya_main
    use mpi
    use alya_jexpresso_coupling
    implicit none

    integer :: ierr, timestep
    real(8), allocatable :: velocity(:), pressure(:)

    ! Initialize MPI
    call MPI_Init(ierr)

    ! Initialize coupling
    call initialize_alya_coupling()

    ! Allocate interface arrays
    allocate(velocity(100), pressure(100))

    ! Initialize Alya solver (use comm_local, not MPI_COMM_WORLD!)
    call alya_init(comm_local)

    ! Time stepping loop
    do timestep = 1, 100

        ! Solve Alya timestep (use comm_local)
        call alya_solve_timestep(comm_local)

        ! Couple with Jexpresso every N steps
        if (mod(timestep, 10) == 0) then

            ! Extract pressure at coupling interface
            call alya_get_interface_pressure(pressure, 100)

            ! Exchange: send pressure, receive velocity
            call exchange_with_jexpresso(pressure, velocity, 100, &
                                        send_tag=201, recv_tag=200)

            ! Apply velocity BC
            call alya_apply_velocity_bc(velocity, 100)

        endif

    enddo

    ! Finalize
    call finalize_alya_coupling()
    call MPI_Finalize(ierr)

end program alya_main
```

## Key Points

### CRITICAL: Use comm_local in Alya

Alya must use `comm_local` for all internal MPI operations, NOT `MPI_COMM_WORLD`:

```fortran
! WRONG - will include Jexpresso ranks!
call MPI_Bcast(data, n, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)

! CORRECT - only Alya ranks
call MPI_Bcast(data, n, MPI_DOUBLE, 0, comm_local, ierr)
```

### MPI Tag Convention

Coordinate tags between Jexpresso and Alya:

| Data Type | Direction | MPI Tag |
|-----------|-----------|---------|
| Velocity | Jexpresso → Alya | 200 |
| Pressure | Alya → Jexpresso | 201 |
| Temperature | (as needed) | 202 |
| Custom data | (as needed) | 300+ |

### Jexpresso Side Coupling

In your Jexpresso user code (e.g., `user_bc.jl`):

```julia
function couple_with_alya!(mesh, state, t)
    ctx = Jexpresso.coupling_ctx

    if isnothing(ctx) || ctx.n_codes == 1
        return  # Not coupling or single code mode
    end

    # Extract velocity at interface
    local_velocity = extract_boundary_velocity(mesh, state)

    # Gather to root
    if is_coupling_root(ctx)
        global_velocity = zeros(100)
    else
        global_velocity = zeros(0)
    end
    gather_to_coupling_root!(ctx, local_velocity, global_velocity)

    # Exchange with Alya (code_id=2)
    if is_coupling_root(ctx)
        global_pressure = zeros(100)
        coupling_exchange_workflow!(ctx, global_velocity, global_pressure, 2;
                                    tag_send=200, tag_recv=201)
    else
        global_pressure = zeros(0)
    end

    # Distribute to all ranks
    local_pressure = zeros(length(local_velocity))
    scatter_from_coupling_root!(ctx, global_pressure, local_pressure)

    # Apply BC
    apply_pressure_bc!(mesh, state, local_pressure)
end
```

## Troubleshooting

### Problem: Hangs at initialization
**Cause:** Alya not calling `MPI_Comm_split`
**Solution:** Ensure Alya initializes coupling before any other MPI operations

### Problem: Wrong data size
**Cause:** Buffer size mismatch between codes
**Solution:** Print array sizes on both sides, ensure they match

### Problem: Segmentation fault in Alya
**Cause:** Using `MPI_COMM_WORLD` instead of `comm_local`
**Solution:** Search Alya code for `MPI_COMM_WORLD` and replace with `comm_local`

### Problem: "Invalid communicator" error
**Cause:** Communicator not properly initialized
**Solution:** Check that `MPI_Comm_split` completed successfully

### Problem: Different MPI versions
**Cause:** Julia and Alya compiled with different MPI libraries
**Solution:**
```bash
# Check Julia's MPI
julia --project=. -e 'using MPI; MPI.versioninfo()'

# Check Alya's MPI
ldd /path/to/alya/bin/alya.x | grep mpi

# They should match (both OpenMPI or both MPICH, etc.)
```

## Testing Checklist

- [ ] Step 1: Jexpresso alone with n_codes=1 works
- [ ] Step 2: MPMD with test_alya_mimic.jl works
- [ ] Step 3: Real Alya coupling interface implemented
- [ ] Step 4: Both codes initialize without errors
- [ ] Step 5: Data exchange completes without deadlock
- [ ] Step 6: Received data is correct
- [ ] Step 7: Both codes finalize successfully

## Next Steps

1. **Run Step 1 test** to verify Jexpresso works with coupling
2. **Run Step 2 test** to verify MPMD launch works
3. **Implement Alya coupling** using the Fortran code above
4. **Coordinate with Alya team** on data formats and timing
5. **Test incrementally** with increasing complexity

## See Also

- `COUPLING_WITH_ALYA.md` - General Alya coupling guide
- `COUPLING_WITH_EXTERNAL_CODE.md` - External code patterns
- `src/mpi/coupling_helpers.jl` - Helper functions
- `test_alya_mimic.jl` - Test program
