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

  MY_MPI_COMM            :: PAR_COMM_FINAL
  integer(4)             :: ierr, rank, size
  character(128)         :: app_name
  ! FIX: Make app_dumm an allocatable array of strings
  character(128), allocatable :: app_dumm(:) 
  
  ! Initialize MPI
  app_name = 'Jexpresso'
  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr) ! Get total number of procs

  ! 1. Allocate memory on the Root process (Rank 0) only
  if (rank == 0) then
     allocate(app_dumm(size)) 
  end if

  ! 2. Fix the Split (Optional context note)
  ! Passing MPI_UNDEFINED to ALL ranks makes PAR_COMM_FINAL = MPI_COMM_NULL for everyone.
  ! This is valid but useless. Assuming you keep it for now:
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)  

  ! 3. Correct MPI_Gather Call
  ! Note: We send 128 characters (1 string)
  ! We receive 128 characters (1 string) FROM EACH PROCESS
  call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                  app_dumm, 128, MPI_CHARACTER, &
                  0, MPI_COMM_WORLD, ierr)
                                                    
  ! Print result to verify
  if (rank == 0) then
     print *, "Received from Rank 1: ", app_dumm(1)
     deallocate(app_dumm)
  end if

!!!!!  call MPI_Finalize(ierr)

end program unitt_alya_with_another_code
