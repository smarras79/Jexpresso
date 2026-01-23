program alya_gather_continuous
  ! Fortran code that continues running after gather
  ! Simulates a full Alya simulation with time steps

#ifdef USEMPIF08
  use mpi_f08
  implicit none
#define MY_MPI_COMM      type(MPI_Comm)
#else
  implicit none
  include 'mpif.h'
#define MY_MPI_COMM      integer(4)
#endif

  MY_MPI_COMM            :: PAR_COMM_FINAL
  integer(4)             :: ierr, rank, size
  character(128)         :: app_name
  character(128), allocatable :: app_dumm(:)
  integer                :: i, time_step
  real(8)                :: start_time, current_time
  real(8), parameter     :: max_time = 10.0  ! Run for 10 seconds

  ! Initialize MPI
  app_name = 'ALYA'
  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

  if (rank == 0) then
    print '(A,I0)', "[Alya rank 0] Starting with world size: ", size
    call flush(6)
  end if

  ! Allocate memory on rank 0
  if (rank == 0) then
     allocate(app_dumm(size))
  end if

  ! Participate in MPI_Comm_split with MPI_UNDEFINED
  if (rank == 0) then
    print *, "[Alya rank 0] Participating in MPI_Comm_split"
    call flush(6)
  end if

  call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)

  if (rank == 0) then
    print *, "[Alya rank 0] MPI_Comm_split completed"
    call flush(6)
  end if

  ! Gather application names
  if (rank == 0) then
    print *, "[Alya rank 0] Gathering application names"
    call flush(6)
  end if

  call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                  app_dumm, 128, MPI_CHARACTER, &
                  0, MPI_COMM_WORLD, ierr)

  ! Print gathered names
  if (rank == 0) then
     print *, ""
     print *, "[Alya rank 0] Received application names:"
     do i = 1, size
        print '(A,I0,A,A)', "  Rank ", i-1, ": ", trim(app_dumm(i))
        call flush(6)
     end do
     print *, ""
     print *, "[Alya rank 0] Coupling initialization complete. Starting simulation..."
     call flush(6)
  end if

  ! ============================================================
  ! SIMULATION LOOP - Alya runs independently from Jexpresso
  ! NO barriers during simulation - each code runs at its own pace
  ! ============================================================

  start_time = MPI_Wtime()
  time_step = 0

  do while (.true.)
    current_time = MPI_Wtime() - start_time
    if (current_time >= max_time) exit

    time_step = time_step + 1

    ! Simulate work (Alya's own computation)
    call sleep(1)

    ! NO MPI_Barrier here - let Jexpresso run independently

    if (rank == 0) then
      print '(A,I0,A,F6.2,A)', "[Alya rank 0] Time step ", time_step, &
                               ", elapsed: ", current_time, "s"
      call flush(6)
    end if

    ! Optional: Send/receive data to/from Jexpresso ranks using point-to-point
    ! Example: MPI_Send to rank 2 (Jexpresso rank 0)
    ! Example: MPI_Recv from rank 2 (Jexpresso rank 0)
    ! This allows asynchronous data exchange without blocking
  end do

  ! ============================================================
  ! FINALIZATION - Both codes must reach this together
  ! ============================================================

  if (rank == 0) then
    print *, ""
    print '(A,I0,A)', "[Alya rank 0] Simulation complete after ", time_step, " steps"
    print *, "[Alya rank 0] Waiting for all processes before finalize"
    call flush(6)
    deallocate(app_dumm)
  end if

  ! Final barrier to ensure all ranks finish together
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  call MPI_Finalize(ierr)

  if (rank == 0) then
    print *, "[Alya rank 0] Done"
    call flush(6)
  end if

end program alya_gather_continuous
