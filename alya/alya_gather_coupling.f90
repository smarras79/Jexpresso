program unitt_alya_with_another_code
  ! Fortran code using MPI_Gather on MPI_COMM_WORLD
  ! Couples with Julia code in the same MPI world

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
  integer                :: i

  ! Initialize MPI
  app_name = 'OTRO5'
  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)

  print '(A,I0,A,I0)', "[Fortran rank ", rank, "] World size: ", size
  call flush(6)

  ! Allocate memory on the Root process (Rank 0) only
  if (rank == 0) then
     allocate(app_dumm(size))
  end if

  ! Split (creates NULL communicator for all ranks - kept as in original)
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)

  ! Gather all application names to rank 0
  call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                  app_dumm, 128, MPI_CHARACTER, &
                  0, MPI_COMM_WORLD, ierr)

  ndime = 3
  call MPI_Bcast(ndime, 1, INTEGER, 0, MPI_COMM_WORLD) 
  !do idime = 1,3
  !   call PAR_BROADCAST(rem_min(idime),'IN THE UNIVERSE') ! real
  !   call PAR_BROADCAST(rem_max(idime),'IN THE UNIVERSE') ! real
  !   call PAR_BROADCAST(rem_nx(idime) ,'IN THE UNIVERSE') ! integer
  !end do


  ! Print result to verify
  if (rank == 0) then
     print *, ""
     print *, "[Fortran rank 0] Received application names from all ranks:"
     do i = 1, size
        print '(A,I0,A,A)', "  Rank ", i-1, ": ", trim(app_dumm(i))
        call flush(6)
     end do
     deallocate(app_dumm)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  print '(A,I0,A)', "[Fortran rank ", rank, "] Done"
  call flush(6)

  call MPI_Finalize(ierr)

end program unitt_alya_with_another_code
