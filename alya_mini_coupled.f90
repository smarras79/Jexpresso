program alya_mini_coupled
! Minimal Alya coupling test - matches Jexpresso-mini.jl
!
! Usage: mpirun -np 2 julia --project=. Jexpresso-mini.jl : -np 2 ./alya_mini_coupled
!
! Compile:
!   mpif90 -o alya_mini_coupled alya_mini_coupled.f90
!

#ifdef USEMPIF08
  use mpi_f08
  implicit none
#define MY_MPI_COMM      type(MPI_Comm)
#else
  implicit none
  include 'mpif.h'
#define MY_MPI_COMM      integer(4)
#endif

  MY_MPI_COMM            :: comm_world, comm_local, comm_inter
  integer(4)             :: ierr, world_rank, world_size
  integer(4)             :: local_rank, local_size
  integer(4)             :: code_id, n_codes, color, key
  integer(4)             :: i
  character(128)         :: app_name
  character(128), allocatable :: app_dumm(:)
  logical                :: is_root

  ! Initialize MPI
  app_name = 'Alya-mini'
  call MPI_Init(ierr)

  comm_world = MPI_COMM_WORLD
  call MPI_COMM_RANK(comm_world, world_rank, ierr)
  call MPI_COMM_SIZE(comm_world, world_size, ierr)

  write(*,*) '[Alya-mini rank', world_rank, '] World size:', world_size

  ! ==========================================================================
  ! CRITICAL: Initialize coupling
  ! ==========================================================================

  code_id = 2
  n_codes = 2

  write(*,*) '[Alya-mini rank', world_rank, '] Initializing coupling...'

  ! Split communicator by code_id to create local communicator
  color = code_id
  key = world_rank
  call MPI_COMM_SPLIT(comm_world, color, key, comm_local, ierr)

  ! Get local rank info
  call MPI_COMM_RANK(comm_local, local_rank, ierr)
  call MPI_COMM_SIZE(comm_local, local_size, ierr)

  is_root = (local_rank == 0)

  ! Create inter-communicator (only for coupling roots)
  if (is_root) then
     color = 1  ! Participate in inter-comm
  else
     color = 0  ! Don't participate
  endif
  call MPI_COMM_SPLIT(comm_world, color, key, comm_inter, ierr)

  write(*,*) '[Alya-mini rank', world_rank, '] Coupling initialized!'
  write(*,*) '  Local rank:', local_rank, '/', local_size
  write(*,*) '  World rank:', world_rank
  write(*,*) '  Is root:', is_root

  ! Barrier to ensure all ranks initialized
  call MPI_BARRIER(comm_world, ierr)

  ! ==========================================================================
  ! Test: Exchange application names
  ! ==========================================================================

  if (is_root) then
     write(*,*) ''
     write(*,*) '[Alya-mini] Testing name exchange...'

     ! Allocate buffer for exchange (only root needs this)
     allocate(app_dumm(1))

     ! Exchange with Jexpresso (rank 0 in world)
     ! Use Sendrecv to avoid deadlock
     call MPI_Sendrecv(app_name, 128, MPI_CHARACTER, 0, 100, &
                       app_dumm(1), 128, MPI_CHARACTER, 0, 100, &
                       comm_world, MPI_STATUS_IGNORE, ierr)

     write(*,*) '[Alya-mini] Exchanged names successfully!'
     write(*,*) '  My name: ', trim(app_name)
     write(*,*) '  Partner name: ', trim(app_dumm(1))

     deallocate(app_dumm)
  endif

  ! ==========================================================================
  ! Simulate some work with synchronization
  ! ==========================================================================

  do i = 1, 5
     if (is_root) then
        write(*,*) ''
        write(*,*) '[Alya-mini] Step', i, '/ 5'
     endif

     ! Synchronize all codes
     call MPI_BARRIER(comm_world, ierr)

     ! Simulate work
     call sleep(1)
  end do

  write(*,*) ''
  write(*,*) '[Alya-mini rank', world_rank, '] Test complete!'

  ! ==========================================================================
  ! Cleanup
  ! ==========================================================================

  call MPI_BARRIER(comm_world, ierr)
  call MPI_Finalize(ierr)

  write(*,*) '[Alya-mini rank', world_rank, '] Done'

end program alya_mini_coupled
