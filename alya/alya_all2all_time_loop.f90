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

  !--------------------------
  ! Time loop variables
  !--------------------------
  integer, parameter                 :: TAG_DATA = 2000  ! Base tag for data exchange
  integer(4)                         :: nranks_julia, partner_base
  integer(4)                         :: step, nsteps, ipartner
  real(kind=8)                       :: t0, dt, tend, t
  real(kind=8), allocatable          :: sendbuf(:), recvbuf(:)
  real(kind=8), allocatable          :: sendbuf_all(:), recvbuf_all(:)
  integer(4)                         :: nsend, nrecv
  real(kind=8)                       :: dummy_field
#ifdef USEMPIF08
  type(MPI_Status), allocatable      :: send_status(:), recv_status(:)
  type(MPI_Request), allocatable     :: send_requests(:), recv_requests(:)
#else
  integer, allocatable               :: send_status(:,:), recv_status(:,:)
  integer, allocatable               :: send_requests(:), recv_requests(:)
#endif
  integer                            :: nactive_send, nactive_recv, ireq

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
  if (rank == 0) then
     print *, '=== Alya: Before MPI_Gather ==='
     flush(6)
  end if

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
     print *, '=== Alya: After MPI_Gather ==='
     flush(6)
  end if

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
  rem_nx  = [4,      4,        1]
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
     print *, '=== Alya: Before barrier before Alltoall ==='
     flush(6)
  end if

  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  if (rank == 0) then
     print *, '=== Alya: After barrier, about to call MPI_Alltoall ==='
     flush(6)
  end if

  call MPI_Alltoall(npoin_send, 1, MPI_INTEGER4, &
       npoin_recv, 1, MPI_INTEGER4, &
       MPI_COMM_WORLD, ierr)

  if (rank == 0) then
     print *, '=== Alya: MPI_Alltoall completed ==='
     flush(6)
  end if

  !---------------------------------------------------------------------------
  ! VERIFICATION: Print detailed coupling pattern
  !---------------------------------------------------------------------------
  total_send = sum(npoin_recv)  ! What Julia needs from me (I will send)
  total_recv = sum(npoin_send)  ! What I need from Julia (should be 0)

  ! Print summary for each Alya rank
  write(*,'(A)') '========================================='
  write(*,'(A,I0,A,I0)') 'Alya local rank ', arank, ', world rank ', rank
  write(*,'(A)') '========================================='
  write(*,'(A,I0)') 'Total points to SEND to Julia: ', total_send
  write(*,'(A,I0)') 'Total points to RECV from Julia: ', total_recv
  write(*,'(A)') '-----------------------------------------'

  ! Print world rank breakdown
  write(*,'(A,I0,A)') 'World rank breakdown (size=', size, '):'
  do i = 0, asize-1
     write(*,'(A,I0,A,I0,A)') '  World rank ', alya_to_world(i), &
          ' = Alya rank ', i, ' (ALYA)'
  end do
  do i = 0, size - asize - 1
     write(*,'(A,I0,A,I0,A)') '  World rank ', asize + i, &
          ' = Julia rank ', i, ' (JULIA)'
  end do
  write(*,'(A)') '-----------------------------------------'

  ! Print detailed send pattern (npoin_recv tells me what to send)
  if (total_send > 0) then
     write(*,'(A)') 'SEND pattern (npoin_recv array):'
     do i = 0, size-1
        if (npoin_recv(i) > 0) then
           ! Determine if this is Alya or Julia rank
           if (i < asize) then
              write(*,'(A,I0,A,I0,A,I0,A)') '  -> World rank ', i, &
                   ' (Alya rank ', i, '): ', npoin_recv(i), ' points'
           else
              write(*,'(A,I0,A,I0,A,I0,A)') '  -> World rank ', i, &
                   ' (Julia rank ', i - asize, '): ', npoin_recv(i), ' points'
           end if
        end if
     end do
  else
     write(*,'(A)') 'SEND pattern: No points to send'
  end if

  ! Print detailed receive pattern (npoin_send tells me what I'll receive)
  if (total_recv > 0) then
     write(*,'(A)') 'RECV pattern (npoin_send array):'
     do i = 0, size-1
        if (npoin_send(i) > 0) then
           if (i < asize) then
              write(*,'(A,I0,A,I0,A,I0,A)') '  <- World rank ', i, &
                   ' (Alya rank ', i, '): ', npoin_send(i), ' points'
           else
              write(*,'(A,I0,A,I0,A,I0,A)') '  <- World rank ', i, &
                   ' (Julia rank ', i - asize, '): ', npoin_send(i), ' points'
           end if
        end if
     end do
  else
     write(*,'(A)') 'RECV pattern: No points to receive'
  end if
  write(*,'(A)') '========================================='
  flush(6)

  !===========================================================================
  ! NEW: TIME INTEGRATION SETUP AND LOOP
  !===========================================================================
  
  ! Calculate number of Julia ranks
  nranks_julia = size - asize
  partner_base = asize  ! First Julia rank in world communicator
  
  ! Time integration parameters
  t0     = 0.0d0
  dt     = 0.25d0     ! Must match Julia's dt
  tend   = 2.0d0
  nsteps = int((tend - t0) / dt)

  if (rank == 0) then
     write(*,'(A)') ' '
     write(*,'(A)') '==========================================='
     write(*,'(A)') '=== STARTING TIME INTEGRATION ==='
     write(*,'(A)') '==========================================='
     write(*,'(A,I0,A,F6.3,A,F6.3)') 'Alya: will run ', nsteps, &
          ' timesteps (dt=', dt, ', tend=', tend, ')'
     write(*,'(A,I0)') 'Number of Julia ranks: ', nranks_julia
     flush(6)
  end if

  ! Allocate send/receive buffers for ALL ranks
  ! Each rank needs buffers for all partners it communicates with
  nsend = total_send  ! Total points this rank sends
  nrecv = total_recv  ! Total points this rank receives
  
  allocate(sendbuf_all(nsend))
  allocate(recvbuf_all(nrecv))
  sendbuf_all = 0.0d0
  recvbuf_all = 0.0d0

  ! Allocate MPI request arrays
  nactive_send = count(npoin_recv > 0)
  nactive_recv = count(npoin_send > 0)
  
  if (nactive_send > 0) then
#ifdef USEMPIF08
     allocate(send_requests(nactive_send))
     allocate(send_status(nactive_send))
#else
     allocate(send_requests(nactive_send))
     allocate(send_status(MPI_STATUS_SIZE, nactive_send))
#endif
  end if
  
  if (nactive_recv > 0) then
#ifdef USEMPIF08
     allocate(recv_requests(nactive_recv))
     allocate(recv_status(nactive_recv))
#else
     allocate(recv_requests(nactive_recv))
     allocate(recv_status(MPI_STATUS_SIZE, nactive_recv))
#endif
  end if

  ! Initialize dummy field
  dummy_field = 42.0d0 + dble(arank)

  !===========================================================================
  ! TIME LOOP
  !===========================================================================
  do step = 1, nsteps
     t = t0 + step * dt

     !------------------------------------------------------------------------
     ! ALYA PHYSICS: Advance solution by one timestep
     !------------------------------------------------------------------------
     dummy_field = dummy_field + 0.01d0 * dble(step)
     
     if (rank == 0) then
        write(*,'(A,F6.3,A)') 'Alya: t=', t, ' - starting coupling exchange'
        flush(6)
     end if

     !------------------------------------------------------------------------
     ! STEP 1: PREPARE SEND DATA
     ! Fill send buffer with data that Julia needs
     !------------------------------------------------------------------------
     ! For now, just fill with dummy field
     ! In real application, this would be actual Alya solution data
     sendbuf_all(1:nsend) = dummy_field

     !------------------------------------------------------------------------
     ! STEP 2: NON-BLOCKING SEND/RECV WITH JULIA RANKS
     !------------------------------------------------------------------------
     
     ! Post all receives first
     ireq = 0
     do i = 0, size-1
        if (npoin_send(i) > 0) then  ! I receive from rank i
           ireq = ireq + 1
           ! Calculate offset in recvbuf_all (simplified - assumes single equation)
           ! In production code, you'd need proper offset calculation
#ifdef USEMPIF08
           call MPI_Irecv(recvbuf_all, nrecv, MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + rank, MPI_COMM_WORLD, recv_requests(ireq), ierr)
#else
           call MPI_IRECV(recvbuf_all, nrecv, MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + rank, MPI_COMM_WORLD, recv_requests(ireq), ierr)
#endif
        end if
     end do

     ! Post all sends
     ireq = 0
     do i = 0, size-1
        if (npoin_recv(i) > 0) then  ! I send to rank i
           ireq = ireq + 1
#ifdef USEMPIF08
           call MPI_Isend(sendbuf_all, nsend, MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + i, MPI_COMM_WORLD, send_requests(ireq), ierr)
#else
           call MPI_ISEND(sendbuf_all, nsend, MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + i, MPI_COMM_WORLD, send_requests(ireq), ierr)
#endif
        end if
     end do

     ! Wait for all receives
     if (nactive_recv > 0) then
#ifdef USEMPIF08
        call MPI_Waitall(nactive_recv, recv_requests, recv_status, ierr)
#else
        call MPI_WAITALL(nactive_recv, recv_requests, recv_status, ierr)
#endif
     end if

     ! Wait for all sends
     if (nactive_send > 0) then
#ifdef USEMPIF08
        call MPI_Waitall(nactive_send, send_requests, send_status, ierr)
#else
        call MPI_WAITALL(nactive_send, send_requests, send_status, ierr)
#endif
     end if

     if (rank == 0) then
        write(*,'(A,F6.3,A,F12.5)') 'Alya: t=', t, &
             ' coupling complete; recv[1]=', recvbuf_all(1)
        flush(6)
     end if

     !------------------------------------------------------------------------
     ! STEP 3: WRITE VTU OUTPUT (all ranks participate)
     !------------------------------------------------------------------------
     call write_alya_grid_vtu(rank, arank, asize, ndime, rem_min, rem_max, rem_nx, &
                              recvbuf_all, nrecv, step, t)

     ! Barrier before master file write
     call MPI_Barrier(PAR_COMM_FINAL, ierr)

     ! Rank 0 writes PVTU master file
     if (rank == 0) then
        call write_pvtu_master(asize, step, t)
     end if

  end do  ! End time loop

  if (rank == 0) then
     write(*,'(A)') ' '
     write(*,'(A)') '==========================================='
     write(*,'(A)') 'Alya: time loop complete'
     write(*,'(A)') '==========================================='
     flush(6)
  end if

  ! Cleanup
  if (allocated(sendbuf_all)) deallocate(sendbuf_all)
  if (allocated(recvbuf_all)) deallocate(recvbuf_all)
  if (allocated(send_requests)) deallocate(send_requests)
  if (allocated(recv_requests)) deallocate(recv_requests)
  if (allocated(send_status)) deallocate(send_status)
  if (allocated(recv_status)) deallocate(recv_status)
  
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

  !---------------------------------------------------------------------------
  subroutine write_alya_grid_vtu(rank, arank, asize, ndime, rem_min, rem_max, rem_nx, &
                                 recvbuf, nvals, timestep, time)
    implicit none

    ! Arguments
    integer(4), intent(in)     :: rank, arank, asize, ndime
    real,       intent(in)     :: rem_min(3), rem_max(3)
    integer,    intent(in)     :: rem_nx(3)
    real(8),    intent(in)     :: recvbuf(:)
    integer(4), intent(in)     :: nvals, timestep
    real(8),    intent(in)     :: time

    ! Local variables
    character(len=256) :: filename
    integer            :: iunit, ierr_local
    integer            :: i, j, k, ipoin
    integer            :: i_start, i_end
    integer            :: npoin_local
    real(8)            :: x, y, z, dx, dy, dz
    real(8)            :: field_value
    integer            :: nmax, r, npoin_per_rank

    ! Calculate total number of grid points
    nmax = rem_nx(1) * rem_nx(2) * rem_nx(3)

    ! Calculate grid spacing
    if (rem_nx(1) > 1) then
       dx = dble(rem_max(1) - rem_min(1)) / dble(rem_nx(1) - 1)
    else
       dx = 0.0d0
    end if

    if (rem_nx(2) > 1) then
       dy = dble(rem_max(2) - rem_min(2)) / dble(rem_nx(2) - 1)
    else
       dy = 0.0d0
    end if

    if (rem_nx(3) > 1) then
       dz = dble(rem_max(3) - rem_min(3)) / dble(rem_nx(3) - 1)
    else
       dz = 0.0d0
    end if

    ! Distribute grid points among Alya ranks
    r = mod(nmax, asize)
    npoin_per_rank = nmax / asize

    ! Determine range of points for this rank
    if (arank < r) then
       i_start = arank * (npoin_per_rank + 1)
       i_end   = i_start + npoin_per_rank
    else
       i_start = r * (npoin_per_rank + 1) + (arank - r) * npoin_per_rank
       i_end   = i_start + npoin_per_rank - 1
    end if

    npoin_local = i_end - i_start + 1

    ! Create filename for this rank
    write(filename, '(A,I6.6,A,I4.4,A)') 'alya_grid_', timestep, '_rank', rank, '.vtu'

    ! Open file
    iunit = 100 + rank
    open(unit=iunit, file=trim(filename), status='replace', action='write', iostat=ierr_local)
    if (ierr_local /= 0) then
       write(*,*) 'Error opening file: ', trim(filename)
       return
    end if

    ! Write VTU header
    write(iunit,'(A)') '<?xml version="1.0"?>'
    write(iunit,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iunit,'(A,E16.8,A)') '  <!-- Time = ', time, ' -->'
    write(iunit,'(A)') '  <UnstructuredGrid>'
    write(iunit,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', npoin_local, &
                                  '" NumberOfCells="', max(0, npoin_local), '">'

    ! Write points
    write(iunit,'(A)') '      <Points>'
    write(iunit,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'

    do ipoin = i_start, i_end
       ! Convert global point index to 3D structured indices (0-based)
       k = ipoin / (rem_nx(1) * rem_nx(2))
       j = (ipoin - k * rem_nx(1) * rem_nx(2)) / rem_nx(1)
       i = mod(ipoin, rem_nx(1))

       ! Calculate physical coordinates
       x = dble(rem_min(1)) + i * dx
       y = dble(rem_min(2)) + j * dy
       z = dble(rem_min(3)) + k * dz

       write(iunit,'(3(E16.8,1X))') x, y, z
    end do

    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '      </Points>'

    ! Write cells (vertices)
    write(iunit,'(A)') '      <Cells>'
    write(iunit,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do ipoin = 0, npoin_local-1
       write(iunit,'(I0)') ipoin
    end do
    write(iunit,'(A)') '        </DataArray>'

    write(iunit,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    do ipoin = 1, npoin_local
       write(iunit,'(I0)') ipoin
    end do
    write(iunit,'(A)') '        </DataArray>'

    write(iunit,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    do ipoin = 1, npoin_local
       write(iunit,'(I0)') 1  ! VTK_VERTEX = 1
    end do
    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '      </Cells>'

    ! Write point data (received from Julia)
    write(iunit,'(A)') '      <PointData Scalars="received_field">'
    write(iunit,'(A)') '        <DataArray type="Float64" Name="received_field" format="ascii">'

    do ipoin = i_start, i_end
       if (ipoin - i_start + 1 <= nvals) then
          field_value = recvbuf(ipoin - i_start + 1)
       else
          field_value = 0.0d0
       end if
       write(iunit,'(E16.8)') field_value
    end do

    write(iunit,'(A)') '        </DataArray>'

    ! Add global point ID for debugging
    write(iunit,'(A)') '        <DataArray type="Int32" Name="global_id" format="ascii">'
    do ipoin = i_start, i_end
       write(iunit,'(I0)') ipoin
    end do
    write(iunit,'(A)') '        </DataArray>'

    write(iunit,'(A)') '      </PointData>'

    ! Close Piece and UnstructuredGrid
    write(iunit,'(A)') '    </Piece>'
    write(iunit,'(A)') '  </UnstructuredGrid>'
    write(iunit,'(A)') '</VTKFile>'

    close(iunit)

  end subroutine write_alya_grid_vtu

  !---------------------------------------------------------------------------
  subroutine write_pvtu_master(nranks, timestep, time)
    implicit none
    integer(4), intent(in) :: nranks, timestep
    real(8),    intent(in) :: time

    character(len=256) :: pvtu_filename, vtu_filename
    integer :: iunit, irank, ierr_local

    ! Create PVTU master file (only rank 0)
    write(pvtu_filename, '(A,I6.6,A)') 'alya_grid_', timestep, '.pvtu'

    iunit = 99
    open(unit=iunit, file=trim(pvtu_filename), status='replace', action='write', iostat=ierr_local)
    if (ierr_local /= 0) then
       write(*,*) 'Error opening PVTU file: ', trim(pvtu_filename)
       return
    end if

    ! Write PVTU header
    write(iunit,'(A)') '<?xml version="1.0"?>'
    write(iunit,'(A)') '<VTKFile type="PUnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iunit,'(A,E16.8,A)') '  <!-- Time = ', time, ' -->'
    write(iunit,'(A)') '  <PUnstructuredGrid GhostLevel="0">'

    ! Point data structure
    write(iunit,'(A)') '    <PPointData Scalars="received_field">'
    write(iunit,'(A)') '      <PDataArray type="Float64" Name="received_field"/>'
    write(iunit,'(A)') '      <PDataArray type="Int32" Name="global_id"/>'
    write(iunit,'(A)') '    </PPointData>'

    ! Points structure
    write(iunit,'(A)') '    <PPoints>'
    write(iunit,'(A)') '      <PDataArray type="Float64" NumberOfComponents="3"/>'
    write(iunit,'(A)') '    </PPoints>'

    ! List all piece files
    do irank = 0, nranks - 1
       write(vtu_filename, '(A,I6.6,A,I4.4,A)') 'alya_grid_', timestep, '_rank', irank, '.vtu'
       write(iunit,'(A,A,A)') '    <Piece Source="', trim(vtu_filename), '"/>'
    end do

    write(iunit,'(A)') '  </PUnstructuredGrid>'
    write(iunit,'(A)') '</VTKFile>'

    close(iunit)

  end subroutine write_pvtu_master

end program unitt_alya_with_another_code
