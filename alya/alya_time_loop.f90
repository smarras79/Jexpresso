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
  
  integer(4),    contiguous, pointer :: recv_comm(:,:)
  integer(4),    contiguous, pointer :: recv_comm_snd(:,:)
  integer(4),    contiguous, pointer :: alya_to_world(:)
  integer(4),    contiguous, pointer :: alya_to_world_snd(:)

  !--------------------------
  ! Per-step coupling variables
  !--------------------------
  integer, parameter                 :: TAG_DATA = 1001
  integer, parameter                 :: TAG_BUFSIZE = 999
  integer(4)                         :: nranks_alya, partner
  integer(4)                         :: nvals, step, nsteps
  real(kind=8)                       :: t0, dt, tend, t
  real(kind=8), allocatable          :: sendbuf(:), recvbuf(:)
  real(kind=8)                       :: dummy_field
#ifdef USEMPIF08
  type(MPI_Status)                   :: status
#else
  integer                            :: status(MPI_STATUS_SIZE)
#endif

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

  app_name = 'ALYA'

  if (rank == 0) then
     allocate(app_dumm(size))
  end if

#ifdef USEMPIF08
  call MPI_Comm_split(MPI_COMM_WORLD, 1_4, rank, PAR_COMM_FINAL, ierr)
#else
  call MPI_COMM_SPLIT(MPI_COMM_WORLD, 1_4, rank, PAR_COMM_FINAL, ierr)
#endif
  
  call MPI_Comm_size(PAR_COMM_FINAL,asize,ierr)
  call MPI_Comm_rank(PAR_COMM_FINAL,arank,ierr)

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

     ! Count Alya ranks and identify partner
     nranks_alya = 0
     do i = 1, size
        s = cstr_trim(app_dumm(i))
        if (trim(s) == 'ALYA') nranks_alya = nranks_alya + 1
     end do
     partner = nranks_alya
     deallocate(app_dumm)
  end if

#ifdef USEMPIF08
  call MPI_Bcast(nranks_alya, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
  call MPI_Bcast(partner,     1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
#else
  call MPI_BCAST(nranks_alya, 1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
  call MPI_BCAST(partner,     1, MPI_INTEGER4, 0, MPI_COMM_WORLD, ierr)
#endif

  !======================================================================
  ! INITIAL METADATA BROADCAST
  !======================================================================
  rem_min = [-5000.0,     0.0, 0.0]
  rem_max = [ 5000.0, 10000.0, 0.0]
  rem_nx  = [10,      10,        1]
  ndime = 2
  
#ifdef USEMPIF08
  call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
  call MPI_BCAST(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
  do idime = 1,3
#ifdef USEMPIF08
     call MPI_Bcast(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
     call MPI_BCAST(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_BCAST(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
  end do

  allocate(alya_to_world_snd(0:asize-1))
  allocate(alya_to_world    (0:asize-1))
  alya_to_world_snd = 0
  alya_to_world_snd(arank) = rank
#ifdef USEMPIF08
  call MPI_AllReduce(alya_to_world_snd,alya_to_world,asize,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
  call MPI_ALLREDUCE(alya_to_world_snd,alya_to_world,asize,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

  allocate(recv_comm_snd(0:size-1,0:size-1))
  allocate(recv_comm    (0:size-1,0:size-1))
  recv_comm_snd = 0
#ifdef USEMPIF08
  call MPI_AllReduce(recv_comm_snd,recv_comm,size*size,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)
#else
  call MPI_ALLREDUCE(recv_comm_snd,recv_comm,size*size,MPI_INTEGER4,MPI_SUM,MPI_COMM_WORLD,ierr)
#endif

  !======================================================================
  ! TIME INTEGRATION SETUP
  !======================================================================
  t0     = 0.0d0
  dt     = 0.25d0     ! Must match Julia's dt
  tend   = 2.0d0
  nsteps = int( (tend - t0) / dt )

  if (rank == 0) then
     write(*,'(A,I0,A,F6.3,A,F6.3)') 'Alya: will run ', nsteps, &
        ' timesteps (dt=', dt, ', tend=', tend, ')'
  end if

  !======================================================================
  ! BARRIER 1: Wait for Julia to finish ALL initialization
  !======================================================================
  if (rank == 0) then
     write(*,'(A)') 'Alya: initialization complete, waiting for Julia...'
  end if
#ifdef USEMPIF08
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#else
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  if (rank == 0) then
     write(*,'(A)') 'Alya: both codes ready for handshake'
  end if

  !======================================================================
  ! HANDSHAKE: Receive buffer size from Julia (AFTER barrier)
  !======================================================================
  if (rank == 0) then
#ifdef USEMPIF08
     call MPI_Recv(nvals, 1, MPI_INTEGER4, partner, TAG_BUFSIZE, MPI_COMM_WORLD, status, ierr)
#else
     call MPI_RECV(nvals, 1, MPI_INTEGER4, partner, TAG_BUFSIZE, MPI_COMM_WORLD, status, ierr)
#endif
     write(*,'(A,I0)') 'Alya: Received buffer size from Jexpresso: nvals=', nvals
     allocate(sendbuf(nvals), recvbuf(nvals))
     sendbuf = 0.0d0
     recvbuf = 0.0d0
  end if

  !======================================================================
  ! BARRIER 2: Ensure handshake complete before time loop
  !======================================================================
#ifdef USEMPIF08
  call MPI_Barrier(MPI_COMM_WORLD, ierr)
#else
  call MPI_BARRIER(MPI_COMM_WORLD, ierr)
#endif
  if (rank == 0) then
     write(*,'(A)') 'Alya: starting time loop'
  end if

  !======================================================================
  ! TIME LOOP WITH PER-STEP COUPLING
  !======================================================================
  dummy_field = 42.0d0

  do step = 1, nsteps
     t = t0 + step*dt

     !-------------------------------------------------------------------
     ! ALYA PHYSICS: Advance solution by one timestep
     !-------------------------------------------------------------------
     if (rank == 0) then
        dummy_field = dummy_field + 0.01d0 * dble(step)
     end if

     !-------------------------------------------------------------------
     ! STEP 1: BROADCAST METADATA (ndime, rem_min, rem_max, rem_nx)
     !-------------------------------------------------------------------
     if (rank == 0) then
        write(*,'(A,F6.3,A)') 'Alya: t=', t, ' broadcasting metadata...'
     end if

#ifdef USEMPIF08
     call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
     call MPI_BCAST(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
     do idime = 1,3
#ifdef USEMPIF08
        call MPI_Bcast(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#else
        call MPI_BCAST(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
        call MPI_BCAST(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
#endif
     end do

     !-------------------------------------------------------------------
     ! STEP 2: EXCHANGE SOLUTION DATA (only Alya root with Julia leader)
     !-------------------------------------------------------------------
     if (rank == 0) then
        sendbuf(1:nvals) = dummy_field

        write(*,'(A,F6.3,A)') 'Alya: t=', t, ' exchanging solution data...'
#ifdef USEMPIF08
        call MPI_Sendrecv(sendbuf, nvals, MPI_DOUBLE_PRECISION, partner, TAG_DATA, &
                          recvbuf, nvals, MPI_DOUBLE_PRECISION, partner, TAG_DATA, &
                          MPI_COMM_WORLD, status, ierr)
#else
        call MPI_SENDRECV(sendbuf, nvals, MPI_DOUBLE_PRECISION, partner, TAG_DATA, &
                          recvbuf, nvals, MPI_DOUBLE_PRECISION, partner, TAG_DATA, &
                          MPI_COMM_WORLD, status, ierr)
#endif
        write(*,'(A,F6.3,A,F12.5)') 'Alya: t=', t, &
           ' coupling complete; recv[1]=', recvbuf(1)
     end if
  end do

  if (rank == 0) then
     write(*,'(A)') 'Alya: time loop complete'
     if (allocated(sendbuf)) deallocate(sendbuf, recvbuf)
  end if

  call MPI_Finalize(ierr)

contains
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
