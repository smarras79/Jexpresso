program alya_mini_coupler
  ! alya_mini_coupler.f90 - Fortran side with dynamic remote-leader discovery
  ! Pairs with Jexpresso-mini-coupled.jl
  ! Tags: Fortran sends 100, receives 101 (opposite of Julia)

  use mpi
  implicit none

  integer :: ierr, world_rank, world_size, appid
  integer :: local_comm, local_rank, local_size
  integer :: inter_comm, remote_leader_world
  integer, parameter :: NAME_LEN = 128
  integer, parameter :: TAG_CREATE = 12345
  character(len=NAME_LEN) :: snd, rcv, partner_name
  character(len=256) :: appid_str
  integer, allocatable :: world_appids(:)
  integer :: i, my_appid
  integer :: dt_schar
  integer :: status(MPI_STATUS_SIZE)

  print *, "[Alya-mini] Starting..."
  call flush(6)

  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, world_rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, world_size, ierr)

  ! Read APPID from environment
  call get_environment_variable("APPID", appid_str)
  if (len_trim(appid_str) == 0) then
    if (world_rank == 0) then
      print *, "[Alya-mini] ERROR: APPID not set. Launch with -x APPID=0 (Fortran) and -x APPID=1 (Julia)."
    end if
    call MPI_Abort(MPI_COMM_WORLD, 1, ierr)
  end if
  read(appid_str, *) appid

  ! Split COMM_WORLD by APPID to create local communicator
  call MPI_Comm_split(MPI_COMM_WORLD, appid, world_rank, local_comm, ierr)
  call MPI_Comm_rank(local_comm, local_rank, ierr)
  call MPI_Comm_size(local_comm, local_size, ierr)

  ! --- Dynamic remote leader discovery via Allgather on WORLD ---
  allocate(world_appids(world_size))
  my_appid = appid
  call MPI_Allgather(my_appid, 1, MPI_INTEGER, &
                     world_appids, 1, MPI_INTEGER, &
                     MPI_COMM_WORLD, ierr)

  remote_leader_world = -1
  do i = 1, world_size
    if (world_appids(i) /= appid) then
      remote_leader_world = i - 1  ! Convert to 0-based rank
      exit
    end if
  end do

  if (remote_leader_world < 0) then
    if (world_rank == 0) then
      print *, "[Alya-mini] ERROR: did not find remote leader"
    end if
    call MPI_Abort(MPI_COMM_WORLD, 2, ierr)
  end if

  if (world_rank == 0) then
    print '(A,I0,A,I0,A,I0)', "[Alya-mini] world_size=", world_size, &
                              " appid=", appid, &
                              " remote_leader_world=", remote_leader_world
    call flush(6)
  end if

  ! Create intercommunicator
  call MPI_Intercomm_create(local_comm, 0, MPI_COMM_WORLD, &
                            remote_leader_world, TAG_CREATE, &
                            inter_comm, ierr)
  if (ierr /= MPI_SUCCESS) then
    print *, "[Alya-mini rank", world_rank, "] ERROR: MPI_Intercomm_create failed"
    call MPI_Abort(MPI_COMM_WORLD, 3, ierr)
  end if

  ! Exchange application names using MPI_SIGNED_CHAR equivalent
  ! Fortran CHARACTER is similar to C's char, use MPI_CHARACTER
  snd = "ALYA"  ! Space-padded automatically
  rcv = ""

  if (local_rank == 0) then
    ! Fortran sends tag 100, receives tag 101 (opposite of Julia which sends 101, receives 100)
    call MPI_Sendrecv(snd, NAME_LEN, MPI_CHARACTER, 0, 100, &
                      rcv, NAME_LEN, MPI_CHARACTER, 0, 101, &
                      inter_comm, status, ierr)
    if (ierr /= MPI_SUCCESS) then
      print *, "[Alya-mini] ERROR: MPI_Sendrecv failed"
      call MPI_Abort(MPI_COMM_WORLD, 4, ierr)
    end if

    partner_name = trim(adjustl(rcv))
    print '(A,A)', "[Alya-mini root] Partner name = ", trim(partner_name)
    call flush(6)
  end if

  ! Synchronization loop
  do i = 1, 5
    call MPI_Barrier(inter_comm, ierr)
    if (local_rank == 0) then
      print '(A,I0,A)', "[Alya-mini] Step ", i, "/5"
      call flush(6)
    end if
  end do

  ! Cleanup
  call MPI_Comm_free(inter_comm, ierr)
  call MPI_Comm_free(local_comm, ierr)
  deallocate(world_appids)
  call MPI_Finalize(ierr)

  print '(A,I0,A)', "[Alya-mini rank ", world_rank, "] Done"
  call flush(6)

end program alya_mini_coupler
