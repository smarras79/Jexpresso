
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
  integer(4)                         :: arank, asize
  character(len=128)                 :: app_name
  character(len=128), allocatable    :: app_dumm(:)     ! receive buffer on root only
  integer                            :: i, idime, ndime ! <-- moved here
  character(len=128)                 :: s               ! <-- moved here

  real,    dimension(1:3)            :: rem_min, rem_max
  integer, dimension(1:3)            :: rem_nx
  
  integer(4),    contiguous, pointer :: recv_comm(:,:)
  integer(4),    contiguous, pointer :: recv_comm_snd(:,:)
  integer(4),    contiguous, pointer :: alya_to_world(:)
  integer(4),    contiguous, pointer :: alya_to_world_snd(:)

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
  nullify(recv_comm,recv_comm_snd,alya_to_world,alya_to_world_snd)

  ! This program's label (Julia ranks will send their own label)
  app_name = 'ALYA'

  ! Allocate receive buffer only on root for GATHER
  if (rank == 0) then
     allocate(app_dumm(size))
  end if

  ! Your split (kept as-is; MPI_UNDEFINED -> MPI_COMM_NULL for everyone)
#ifdef USEMPIF08
  call MPI_Comm_split(MPI_COMM_WORLD, 1_4, rank, PAR_COMM_FINAL, ierr)
#else
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, 1_4, rank, PAR_COMM_FINAL, ierr)
#endif
  
  call MPI_Comm_size(PAR_COMM_FINAL,asize,ierr)
  call MPI_Comm_rank(PAR_COMM_FINAL,arank,ierr)

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
  do idime = 1,3
     rem_min(idime) = 10.1 + idime
     rem_max(idime) = 100.1 + 10.0*idime
     rem_nx(idime)  = idime

     call MPI_Bcast(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end do

  allocate(alya_to_world_snd(0:asize-1))
  allocate(alya_to_world    (0:asize-1))
  alya_to_world_snd = 0
  alya_to_world_snd(arank) = rank
  call MPI_AllReduce(alya_to_world_snd,alya_to_world,asize,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)

  allocate(recv_comm_snd(0:size-1,0:size-1))
  allocate(recv_comm    (0:size-1,0:size-1))
  recv_comm_snd = 0
  call MPI_AllReduce(recv_comm_snd,recv_comm,size*size,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)

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
