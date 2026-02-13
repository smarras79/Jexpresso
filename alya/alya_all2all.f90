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
  ! Declarations
  !--------------------------
  MY_MPI_COMM                        :: PAR_COMM_FINAL
  integer(4)                         :: ierr, rank, size
  integer(4)                         :: arank, asize
  character(len=128)                 :: app_name
  character(len=128), allocatable    :: app_dumm(:)
  integer                            :: i, idime, ndime
  character(len=128)                 :: s

  real,    dimension(1:3)            :: rem_min, rem_max
  integer, dimension(1:3)            :: rem_nx
  
  integer(4),    contiguous, pointer :: alya_to_world(:)
  integer(4),    contiguous, pointer :: alya_to_world_snd(:)
  
  ! Coupling communication arrays
  integer(4),    allocatable         :: npoin_send(:)     ! Points I send to each world rank
  integer(4),    allocatable         :: npoin_recv(:)     ! Points I receive from each world rank
  integer(4)                         :: total_send, total_recv

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
  nullify(alya_to_world,alya_to_world_snd)

  app_name = 'ALYA'

  if (rank == 0) then
     allocate(app_dumm(size))
  end if

  ! Split communicator: Alya uses appid=1
#ifdef USEMPIF08
  call MPI_Comm_split(MPI_COMM_WORLD, 1_4, rank, PAR_COMM_FINAL, ierr)
#else
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, 1_4, rank, PAR_COMM_FINAL, ierr)
#endif
  
  call MPI_Comm_size(PAR_COMM_FINAL, asize, ierr)
  call MPI_Comm_rank(PAR_COMM_FINAL, arank, ierr)

  !---------------------------------------------------------------------------
  ! STEP 1: HANDSHAKE - Gather app identities
  !---------------------------------------------------------------------------
#ifdef USEMPIF08
  call MPI_Gather(app_name, 128, MPI_CHARACTER, &
                  app_dumm,  128, MPI_CHARACTER, &
                  0, MPI_COMM_WORLD, ierr)
#else
  call MPI_GATHER(app_name, 128, MPI_CHARACTER, &
                  app_dumm,  128, MPI_CHARACTER, &
                  0, MPI_COMM_WORLD, ierr)
#endif

  if (rank == 0) then
     print *, '=== Coupling labels gathered on root (world size=', size, ') ==='
     do i = 1, size
        s = cstr_trim(app_dumm(i))
        write(*,'(A,I0,A,1X,A)') 'From world rank ', i-1, ':', s
     end do
     deallocate(app_dumm)
  end if

  !---------------------------------------------------------------------------
  ! STEP 2: SEND GRID METADATA to Jexpresso
  !---------------------------------------------------------------------------
  rem_min = [-5000.0,     0.0, 0.0]
  rem_max = [ 5000.0, 10000.0, 0.0]
  rem_nx  = [10,      10,        1]
  ndime = 2
  
  if (rank == 0) then
     print *, '=== Sending grid metadata to Jexpresso ==='
     write(*,'(A,I0)') 'ndime = ', ndime
     write(*,'(A,3F12.2)') 'rem_min = ', rem_min(1:ndime)
     write(*,'(A,3F12.2)') 'rem_max = ', rem_max(1:ndime)
     write(*,'(A,3I6)')    'rem_nx  = ', rem_nx(1:ndime)
  end if
  
  call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  do idime = 1, 3
     call MPI_Bcast(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end do

  !---------------------------------------------------------------------------
  ! Build Alya to World rank mapping (for reference)
  !---------------------------------------------------------------------------
  allocate(alya_to_world_snd(0:asize-1))
  allocate(alya_to_world    (0:asize-1))
  alya_to_world_snd = 0
  alya_to_world_snd(arank) = rank
  call MPI_AllReduce(alya_to_world_snd, alya_to_world, asize, &
                     MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)
  
  if (rank == 0) then
     print *, '=== Alya to World rank mapping ==='
     do i = 0, asize-1
        write(*,'(A,I0,A,I0)') 'Alya rank ', i, ' -> World rank ', alya_to_world(i)
     end do
  end if

  !---------------------------------------------------------------------------
  ! STEP 3, 4, 5: RECEIVE COUPLING PATTERN from Jexpresso via MPI_Alltoall
  ! Jexpresso has computed which Alya grid points fall in which Julia ranks
  ! and now tells us via Alltoall
  !---------------------------------------------------------------------------
  allocate(npoin_send(0:size-1))
  allocate(npoin_recv(0:size-1))
  npoin_send = 0
  npoin_recv = 0

  if (rank == 0) then
     print *, '=== Waiting for coupling pattern from Jexpresso ==='
  end if

  ! Jexpresso sends how many points it needs from each Alya rank
  ! After Alltoall:
  !   npoin_recv(i) = how many points world rank i needs from me
  !   npoin_send(i) = how many points I need from world rank i (should be 0 for Alya->Alya)
  call MPI_Alltoall(npoin_send, 1, MPI_INTEGER4, &
                    npoin_recv, 1, MPI_INTEGER4, &
                    MPI_COMM_WORLD, ierr)

  ! Transpose: what Jexpresso needs from me = what I send to Jexpresso
  ! For Alya, npoin_send should contain counts for sending to Julia ranks
  ! and npoin_recv contains what we receive from Julia (likely 0 for now)
  
  total_send = sum(npoin_recv)  ! This is what Julia needs from me
  total_recv = sum(npoin_send)  ! This is what I need from Julia (should be 0)
  
  if (rank == 0 .or. total_send > 0) then
     write(*,'(A,I0,A,I0,A,I0)') '[Alya rank ', arank, ', world ', rank, &
                                 '] Will send to Julia: ', total_send, ' points'
     
     ! Print non-zero send counts
     do i = 0, size-1
        if (npoin_recv(i) > 0) then
           write(*,'(A,I0,A,I0,A)') '  Will send ', npoin_recv(i), &
                                    ' points to world rank ', i
        end if
     end do
  end if

  if (rank == 0) then
     print *, '=== Coupling pattern received successfully ==='
  end if

  deallocate(npoin_send, npoin_recv)
  deallocate(alya_to_world, alya_to_world_snd)
  
  call MPI_Finalize(ierr)

contains

  !---------------------------------------------------------------------------
  pure function cstr_trim(str) result(out)
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
