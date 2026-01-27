
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

  !--------------------------
  ! Declarations (specification part)
  !--------------------------
  MY_MPI_COMM                        :: PAR_COMM_FINAL
  integer(4)                         :: ierr, rank, size
  character(len=128)                 :: app_name
  character(len=128), allocatable    :: app_dumm(:)   ! receive buffer on root only
  integer                            :: i, ndime      ! <-- moved here
  character(len=128)                 :: s             ! <-- moved here

  !===================== execution part =====================
#ifdef USEMPIF08
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
#else
  call MPI_Init(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
#endif

  ! This program's label (Julia ranks will send their own label)
  app_name = 'ALYA'

  ! Allocate receive buffer only on root for GATHER
  if (rank == 0) then
     allocate(app_dumm(size))
  end if

  ! Your split (kept as-is; MPI_UNDEFINED -> MPI_COMM_NULL for everyone)
#ifdef USEMPIF08
  call MPI_Comm_split(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)
#else
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, MPI_UNDEFINED, rank, PAR_COMM_FINAL, ierr)
#endif

  ! Gather: collect 128 chars from each rank to root
#ifdef USEMPIF08
  call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                  app_dumm,  128, MPI_CHARACTER, &
                  0, MPI_COMM_WORLD, ierr)
#else
  call MPI_GATHER(app_name, 128, MPI_CHARACTER, &
                  app_dumm,  128, MPI_CHARACTER, &
                  0, MPI_COMM_WORLD, ierr)
#endif

  ! Root prints all contributions (from world ranks 0..size-1)
  if (rank == 0) then
     print *, '=== Coupling labels gathered on root (world size=', size, ') ==='
     do i = 1, size
        s = cstr_trim(app_dumm(i))           ! strip NULs and spaces
        write(*,'(A,I0,A,1X,A)') 'From world rank ', i-1, ':', s
     end do
     deallocate(app_dumm)
  end if

  ndime = 3
  call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD)
  
  call MPI_Finalize(ierr)

contains
  pure function cstr_trim(str) result(out)
    ! Trim both trailing spaces and any trailing NUL bytes (CHAR(0))
    character(len=*), intent(in) :: str
    character(len=len(str))      :: out
    integer :: nulpos

    nulpos = index(str, char(0))
    if (nulpos > 0) then
       out = str(:nulpos-1)
    else
      out = str
    end if
    out = trim(adjustl(out))
  end function cstr_trim

end program unitt_alya_with_another_code
