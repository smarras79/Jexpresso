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

  MY_MPI_COMM                        :: PAR_COMM_FINAL
  integer(4)                         :: ierr, rank, size
  integer(4)                         :: arank, asize
  character(len=128)                 :: app_name
  character(len=128), allocatable    :: app_dumm(:)
  integer                            :: i, idime, ndime, itime
  integer(4)                         :: neqs
  character(len=128)                 :: s

  real(8), dimension(1:3)            :: rem_min, rem_max
  integer, dimension(1:3)            :: rem_nx
  
  integer(4),    contiguous, pointer :: alya_to_world(:)
  integer(4),    contiguous, pointer :: alya_to_world_snd(:)

  integer(4),    allocatable         :: npoin_send(:)
  integer(4),    allocatable         :: npoin_recv(:)
  integer(4)                         :: total_pts

  integer(4)                         :: nranks_julia

  ! Jexpresso SEM node list received once during setup
  integer(8), allocatable            :: je_gids_all(:)     ! global node IDs from connijk
  integer(4)                         :: je_total_nodes, je_gid_offset

  integer(4)                         :: step, nsteps
  real(kind=8)                       :: t0, dt, tend, t
  
  ! Output scheduling
  real(kind=8)                       :: out_dt, out_tend
  real(kind=8)                       :: next_t, tol, dt_step, t_plus
  logical                            :: write_now
  
  real(kind=8), allocatable          :: recvbuf_all(:)
  real(kind=8), allocatable          :: ordered_buf(:)  ! recvbuf_all scattered to [i_start..i_end] order
  integer                            :: jpt, local_pos

#ifdef USEMPIF08
  type(MPI_Status),  allocatable     :: recv_status(:)
  type(MPI_Request), allocatable     :: recv_requests(:)
#else
  integer,           allocatable     :: recv_status(:,:)
  integer,           allocatable     :: recv_requests(:)
#endif
  integer                            :: nactive, ireq
  integer                            :: recv_offset

  integer                            :: nmax, r_rem, npoin_per_rank
  integer                            :: i_start, i_end, npoin_local
  integer                            :: nworkers, iworker
  integer                            :: ipoin, ix, iy, iz
  real(kind=8)                       :: xc, yc, dx, dy, dz

  !==========================================================================
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
  call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)
  nullify(alya_to_world, alya_to_world_snd)

  app_name = 'ALYA'
  if (rank == 0) allocate(app_dumm(size))

  call MPI_Comm_split(MPI_COMM_WORLD, 1_4, rank, PAR_COMM_FINAL, ierr)
  call MPI_Comm_size(PAR_COMM_FINAL, asize, ierr)
  call MPI_Comm_rank(PAR_COMM_FINAL, arank, ierr)

  !--------------------------------------------------------------------------
  ! STEP 0: HANDSHAKE
  !--------------------------------------------------------------------------
  call MPI_Gather(app_name, 128, MPI_CHARACTER, &
       app_dumm,  128, MPI_CHARACTER, &
       0, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
     print *, '=== Coupling labels (world size=', size, ') ==='
     do i = 1, size
        s = cstr_trim(app_dumm(i))
        write(*,'(A,I0,A,1X,A)') '  world rank ', i-1, ': ', s
     end do
     deallocate(app_dumm)
     flush(6)
  end if

  !--------------------------------------------------------------------------
  ! STEP 1: Time quantities
  !--------------------------------------------------------------------------
  t0     = 0.0d0
  dt     = 0.5d0
  tend   = 50.0d0
  nsteps = int((tend - t0) / dt)
  out_dt = 10.0d0       ! <-- write every 10 seconds
  !--------------------------------------------------------------------------
  ! STEP 2: GRID METADATA  (broadcast order must match Julia exactly)
  !--------------------------------------------------------------------------
  rem_min = [-5000.0,     0.0, 0.0]
  rem_max = [ 5000.0, 10000.0, 0.0]
  rem_nx  = [50,      50,      1]
  ndime   = 2

  if (rank == 0) then
     write(*,'(A,I0)')     'ndime   = ', ndime
     write(*,'(A,3F12.2)') 'rem_min = ', rem_min(1:ndime)
     write(*,'(A,3F12.2)') 'rem_max = ', rem_max(1:ndime)
     write(*,'(A,3I6)')    'rem_nx  = ', rem_nx(1:ndime)
     flush(6)
  end if

  ! 2a. ndime
  call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  ! 2b. rem_min / rem_max / rem_nx  (full 3-element arrays, matching Julia)
  do idime = 1, 3
     call MPI_Bcast(rem_min(idime), 1, MPI_DOUBLE_PRECISION,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_max(idime), 1, MPI_DOUBLE_PRECISION,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end do
  
  !
  ! Thsi can be useful
  !
  ! 2c. neqs
  !neqs = 4
  !call MPI_Bcast(neqs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  !! 2d. nsteps
  !call MPI_Bcast(nsteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  !--------------------------------------------------------------------------
  ! Alya -> World rank map
  !--------------------------------------------------------------------------
  allocate(alya_to_world_snd(0:asize-1))
  allocate(alya_to_world    (0:asize-1))
  alya_to_world_snd        = 0
  alya_to_world_snd(arank) = rank
  call MPI_AllReduce(alya_to_world_snd, alya_to_world, asize, &
       MPI_INTEGER4, MPI_SUM, MPI_COMM_WORLD, ierr)

  if (rank == 0) then
     print *, '=== Alya -> World rank map ==='
     do i = 0, asize-1
        write(*,'(A,I0,A,I0)') '  Alya ', i, ' -> world ', alya_to_world(i)
     end do
     flush(6)
  end if

  !--------------------------------------------------------------------------
  ! STEP 3: COUPLING COUNT EXCHANGE (Alltoall)
  !--------------------------------------------------------------------------
  allocate(npoin_send(0:size-1))
  allocate(npoin_recv(0:size-1))
  npoin_send = 0
  npoin_recv = 0

  call MPI_Barrier(MPI_COMM_WORLD, ierr)
  call MPI_Alltoall(npoin_send, 1, MPI_INTEGER4, &
       npoin_recv, 1, MPI_INTEGER4, &
       MPI_COMM_WORLD, ierr)

  total_pts = sum(npoin_recv)

  write(*,'(A)') repeat('=',60)
  write(*,'(A,I0,A,I0)') 'Alya local rank ', arank, '  world rank ', rank
  write(*,'(A,I0)') '  Pts to RECV/SEND (Julia<->Alya): ', total_pts
  do i = 0, size-1
     if (npoin_recv(i) > 0) then
        if (i < asize) then
           write(*,'(A,I0,A,I0,A,I0,A)') '  <-> world ', i, &
                ' (Alya ', i, '): ', npoin_recv(i), ' pts'
        else
           write(*,'(A,I0,A,I0,A,I0,A)') '  <-> world ', i, &
                ' (Julia ', i-asize, '): ', npoin_recv(i), ' pts'
        end if
     end if
  end do
  write(*,'(A)') repeat('=',60)
  flush(6)

  !--------------------------------------------------------------------------
  ! STEP 3b: RECEIVE JEXPRESSO SEM NODE LIST (once, before time loop)
  !
  ! The Alltoall in STEP 3 already told each Alya rank how many nodes to
  ! expect from each Julia rank (npoin_recv(i) for i >= asize).
  ! Julia now sends one message per partner: the sorted unique global node
  ! IDs (integer(8)) taken from connijk.
  !--------------------------------------------------------------------------
  nranks_julia = size - asize

  je_total_nodes = sum(npoin_recv(asize:size-1))

  allocate(je_gids_all(max(1, je_total_nodes)))
  je_gids_all = 0_8

  je_gid_offset = 0
  do i = asize, size-1
     if (npoin_recv(i) > 0) then
        call MPI_Recv(je_gids_all(je_gid_offset + 1), npoin_recv(i), &
                      MPI_INTEGER8, i, 0, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        write(*,'(A,I0,A,I0,A,I0,A)') &
             '[node_list] Alya world_rank=', rank, &
             ' received global node IDs from Julia world_rank=', i, &
             ' (', npoin_recv(i), ' IDs)'
        flush(6)
        je_gid_offset = je_gid_offset + npoin_recv(i)
     end if
  end do

  write(*,'(A,I0,A,I0)') &
       '[node_list] Total Jexpresso SEM nodes received by Alya world_rank=', rank, &
       ': ', je_total_nodes
  if (je_total_nodes > 0) then
     write(*,'(A)') '[node_list]   ID list:'
     do i = 1, je_total_nodes
        write(*,'(A,I0,A,I0)') '    [', i, '] gid = ', je_gids_all(i)
     end do
  end if
  flush(6)

  !--------------------------------------------------------------------------
  ! Pre-compute this rank's point range.
  ! Alya local rank 0 is the driving/master rank: it coordinates execution
  ! but never participates in data exchange with Jexpresso, and therefore
  ! owns no remote grid points.  Points are distributed only over the
  ! nworkers = asize-1 worker ranks (Alya local ranks 1..asize-1).
  !--------------------------------------------------------------------------
  nmax        = rem_nx(1) * rem_nx(2) * rem_nx(3)
  nworkers    = max(1, asize - 1)
  r_rem          = mod(nmax, nworkers)
  npoin_per_rank = nmax / nworkers

  i_start     = -1
  i_end       = -2
  npoin_local = 0
  if (arank > 0) then
     iworker = arank - 1   ! 0-based worker index
     if (iworker < r_rem) then
        i_start = iworker * (npoin_per_rank + 1)
        i_end   = i_start + npoin_per_rank
     else
        i_start = r_rem * (npoin_per_rank + 1) + (iworker - r_rem) * npoin_per_rank
        i_end   = i_start + npoin_per_rank - 1
     end if
     npoin_local = i_end - i_start + 1
  end if

  dx = 0.0d0; dy = 0.0d0; dz = 0.0d0
  if (rem_nx(1) > 1) dx = dble(rem_max(1)-rem_min(1)) / dble(rem_nx(1)-1)
  if (rem_nx(2) > 1) dy = dble(rem_max(2)-rem_min(2)) / dble(rem_nx(2)-1)
  if (rem_nx(3) > 1) dz = dble(rem_max(3)-rem_min(3)) / dble(rem_nx(3)-1)

  !==========================================================================
  ! TIME INTEGRATION SETUP
  !==========================================================================
  nranks_julia = size - asize

  if (rank == 0) then
     write(*,'(A,I0,A,F6.3,A,F6.3)') 'Steps: ', nsteps, '  dt=', dt, '  tend=', tend
     flush(6)
  end if

  !--------------------------------------------------------------------------
  ! Receive neqs from Julia via Allreduce(MAX).
  ! Alya contributes 0; Julia contributes its qp.neqs.
  ! The collective synchronises both sides before buffer allocation.
  !--------------------------------------------------------------------------
  neqs = 0
  call MPI_Allreduce(MPI_IN_PLACE, neqs, 1, MPI_INTEGER4, MPI_MAX, &
                     MPI_COMM_WORLD, ierr)
  if (rank == 0) then
     write(*,'(A,I0)') '  neqs received from Julia: ', neqs
     flush(6)
  end if

  allocate(recvbuf_all(max(1, total_pts * neqs)))
  recvbuf_all = 0.0d0

  ! ordered_buf: values placed at the correct [i_start..i_end] position
  ! using je_gids_all as the scatter index map.
  allocate(ordered_buf(max(1, npoin_local * neqs)))
  ordered_buf = 0.0d0

  nactive = count(npoin_recv > 0)

  if (nactive > 0) then
#ifdef USEMPIF08
     allocate(recv_requests(nactive), recv_status(nactive))
#else
     allocate(recv_requests(nactive), recv_status(MPI_STATUS_SIZE, nactive))
#endif
  end if

  !==========================================================================
  ! TIME LOOP
  !==========================================================================
  ! --- VTS output schedule: write once every out_dt seconds
  out_tend   = tend         ! last time eligible for output

  ! Tolerance for floating-point comparisons
  tol = max(1.0d-10, 1.0d-6 * out_dt)

  ! next_t: the next simulation time at which a VTS file should be written.
  ! Start at t0 + out_dt so the first write is at t = out_dt (not at t = 0).
  t      = t0
  next_t = t0 + out_dt

  ! Optional: disable writes if out_dt <= 0
  if (out_dt <= 0.0d0) next_t = huge(1.0d0)

  itime = 1
  do step = 1, nsteps

     ! Use the dt of this step (constant or variable)
     dt_step = dt
     t_plus  = t + dt_step

     ! Write when t_plus reaches or passes the next scheduled output time.
     write_now = (out_dt > 0.0d0 .and. t_plus >= next_t - tol .and. next_t <= out_tend + tol)

     !------------------------------------------------------------------------
     ! MPI RECEIVE: DATA (Julia -> Alya)
     ! Each sender's data goes to a distinct offset in recvbuf_all so that
     ! multiple senders do not overwrite each other.
     !------------------------------------------------------------------------
     ireq = 0
     recv_offset = 0
     do i = 0, size-1
        if (npoin_recv(i) > 0) then
           ireq = ireq + 1
#ifdef USEMPIF08
           call MPI_Irecv(recvbuf_all(recv_offset + 1), npoin_recv(i) * neqs, &
                MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, &
                recv_requests(ireq), ierr)
#else
           call MPI_IRECV(recvbuf_all(recv_offset + 1), npoin_recv(i) * neqs, &
                MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, &
                recv_requests(ireq), ierr)
#endif
           recv_offset = recv_offset + npoin_recv(i) * neqs
        end if
     end do

     ! --- Wait for data receives ---
     if (nactive > 0) then
#ifdef USEMPIF08
        call MPI_Waitall(nactive, recv_requests, recv_status, ierr)
#else
        call MPI_WAITALL(nactive, recv_requests, recv_status, ierr)
#endif
        if (ierr /= 0) write(*,'(A,I0,A,I0)') &
             'ERROR MPI_Waitall recv rank ', rank, ' step ', step
     end if

     !------------------------------------------------------------------------
     ! PRINT: received interpolated data at coordinates
     !------------------------------------------------------------------------
!!$     
!!$     if (total_pts > 0) then
!!$        write(*,'(A,I0,A,F10.4,A,I0,A,I0)') &
!!$             '--- Recv step=', step, '  t=', t, &
!!$             '  Alya_rank=', arank, '  world_rank=', rank
!!$        write(*,'(A8,A6,A6,A18,A18,A18)') 'GlobalID','ix','iy','x','y','rho'
!!$        write(*,'(A)') repeat('-', 80)
!!$        do ipoin = i_start, i_end
!!$           iz = ipoin / (rem_nx(1)*rem_nx(2))
!!$           iy = (ipoin - iz*rem_nx(1)*rem_nx(2)) / rem_nx(1)
!!$           ix = mod(ipoin, rem_nx(1))
!!$           xc = dble(rem_min(1)) + ix*dx
!!$           yc = dble(rem_min(2)) + iy*dy
!!$           i  = ipoin - i_start + 1
!!$           write(*,'(I8,I6,I6,E18.8,E18.8,E18.8)') &
!!$                ipoin, ix, iy, xc, yc, recvbuf_all((i-1)*neqs + 1)
!!$        end do
!!$        write(*,'(A)') repeat('=', 80)
!!$        flush(6)
!!$     end if

     !------------------------------------------------------------------------
     ! VTS OUTPUT: structured grid with filled contours
     !   - all Alya ranks gather to rank 0 via MPI_Gatherv (PAR_COMM_FINAL)
     !   - rank 0 writes a single .vts file per timestep
     !------------------------------------------------------------------------
     if (write_now) then
        ! Advance the target time for the next write.
        next_t = next_t + out_dt

        !------------------------------------------------------------------------
        ! REORDER: scatter recvbuf_all → ordered_buf using je_gids_all.
        ! Done only when writing output, not at every timestep.
        ! je_gids_all(jpt) is the 1-based global Alya point ID for the jpt-th
        ! received value.  local_pos = gid - 1 - i_start gives its 0-based
        ! position within this rank's [i_start..i_end] range.
        !------------------------------------------------------------------------
        ordered_buf = 0.0d0
        do jpt = 1, total_pts
           local_pos = int(je_gids_all(jpt)) - 1 - i_start
           if (local_pos >= 0 .and. local_pos < npoin_local) then
              ordered_buf(local_pos * neqs + 1 : (local_pos + 1) * neqs) = &
                   recvbuf_all((jpt - 1) * neqs + 1 : jpt * neqs)
           end if
        end do

        call write_alya_grid_vts(arank, asize, PAR_COMM_FINAL, ndime, &
             rem_min, rem_max, rem_nx, &
             ordered_buf, npoin_local, neqs, step, itime, t_plus)
        call MPI_Barrier(PAR_COMM_FINAL, ierr)
        itime = itime + 1
     end if

     t = t_plus
     
  end do  ! End time loop

  if (rank == 0) then
     write(*,'(A)') repeat('=',60)
     write(*,'(A)') 'Alya: time loop complete, syncing with Julia...'
     write(*,'(A)') repeat('=',60)
     flush(6)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)

  ! Cleanup
  if (allocated(recvbuf_all))   deallocate(recvbuf_all)
  if (allocated(recv_requests)) deallocate(recv_requests)
  if (allocated(recv_status))   deallocate(recv_status)
  if (allocated(je_gids_all))    deallocate(je_gids_all)
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
    if (nulpos > 0) then; out = str(:nulpos-1); else; out = str; end if
    out = trim(adjustl(out))
  end function cstr_trim

  !===========================================================================
  !  write_alya_grid_vts
  !
  !  Writes a VTK StructuredGrid (.vts) file so ParaView renders filled
  !  contours instead of a point cloud.
  !
  !  Strategy:
  !    1. Each Alya rank contributes its coupling points to a MPI_Gatherv
  !       within PAR_COMM_FINAL (for both data values AND coordinates).
  !    2. Alya rank 0 assembles the full field and writes one .vts file.
  !
  !  Data ordering:
  !    After MPI_Gatherv, full_field(k*neqs + ieq) holds the value for
  !    global point k (0-based). Grid indices (ix,iy,iz) are computed
  !    directly from k using the regular grid layout (ix fastest):
  !       iz = k / (rem_nx(1)*rem_nx(2))
  !       iy = (k - iz*rem_nx(1)*rem_nx(2)) / rem_nx(1)
  !       ix = mod(k, rem_nx(1))
  !
  !  Each variable (rho, u, v, theta) is written as a separate DataArray
  !  so ParaView can visualise them independently.
  !===========================================================================
  subroutine write_alya_grid_vts(arank, asize, PAR_COMM, ndime, &
                                  rem_min, rem_max, rem_nx, &
                                  recvbuf, total_pts, neqs, &
                                  timestep, itime, time)
    implicit none
    integer(4),   intent(in) :: arank, asize, ndime, itime
    MY_MPI_COMM,  intent(in) :: PAR_COMM
    real(8),      intent(in) :: rem_min(3), rem_max(3)
    integer,      intent(in) :: rem_nx(3)
    real(8),      intent(in) :: recvbuf(:)
    integer(4),   intent(in) :: total_pts, neqs, timestep
    real(8),      intent(in) :: time

    ! Local variables
    integer          :: nmax, r_rem_loc, np_base
    integer          :: i_start_r, count_r
    integer          :: irank, iunit, ierr_loc, ierr_mpi
    integer          :: nworkers_loc, iworker_loc
    integer          :: ix, iy, iz, k, ieq
    real(8)          :: x, y, z, dx, dy, dz
    character(len=256)          :: filename
    integer,  allocatable       :: sendcounts(:), displs(:)
    real(8),  allocatable       :: full_field(:)
    real(8),  allocatable       :: grid_var(:,:,:)

    ! Variable names for the neqs fields
    character(len=16) :: varname

    !--- Grid parameters ---
    ! Alya local rank 0 is the driving rank and owns no points.
    ! Distribute nmax over nworkers = asize-1 worker ranks (ranks 1..asize-1).
    nmax          = rem_nx(1) * rem_nx(2) * rem_nx(3)
    nworkers_loc  = max(1, asize - 1)
    r_rem_loc     = mod(nmax, nworkers_loc)
    np_base       = nmax / nworkers_loc

    dx = 0.0d0; dy = 0.0d0; dz = 0.0d0
    if (rem_nx(1) > 1) dx = dble(rem_max(1)-rem_min(1)) / dble(rem_nx(1)-1)
    if (rem_nx(2) > 1) dy = dble(rem_max(2)-rem_min(2)) / dble(rem_nx(2)-1)
    if (rem_nx(3) > 1) dz = dble(rem_max(3)-rem_min(3)) / dble(rem_nx(3)-1)

    !--- Compute Gatherv arrays for data ---
    allocate(sendcounts(0:asize-1))
    allocate(displs    (0:asize-1))

    do irank = 0, asize-1
       if (irank == 0) then
          ! Master rank owns no points
          sendcounts(irank) = 0
          displs(irank)     = 0
       else
          iworker_loc = irank - 1   ! 0-based worker index
          if (iworker_loc < r_rem_loc) then
             i_start_r = iworker_loc * (np_base + 1)
             count_r   = np_base + 1
          else
             i_start_r = r_rem_loc * (np_base + 1) + (iworker_loc - r_rem_loc) * np_base
             count_r   = np_base
          end if
          sendcounts(irank) = count_r * neqs
          displs(irank)     = i_start_r * neqs
       end if
    end do

    !--- This rank's local send count ---
    count_r = sendcounts(arank)

    !--- Allocate receive buffer only on root ---
    if (arank == 0) then
       allocate(full_field(nmax * neqs))
       full_field = 0.0d0
    else
       allocate(full_field(1))
    end if

    !--- Gather DATA to Alya rank 0 ---
#ifdef USEMPIF08
    call MPI_Gatherv(recvbuf(1), count_r, MPI_DOUBLE_PRECISION, &
                     full_field(1), sendcounts, displs, MPI_DOUBLE_PRECISION, &
                     0, PAR_COMM, ierr_mpi)
#else
    call MPI_GATHERV(recvbuf(1), count_r, MPI_DOUBLE_PRECISION, &
                     full_field(1), sendcounts, displs, MPI_DOUBLE_PRECISION, &
                     0, PAR_COMM, ierr_mpi)
#endif
    if (ierr_mpi /= 0) then
       write(*,'(A,I0,A,I0)') 'ERROR MPI_Gatherv data: arank=', arank, ' step=', timestep
       flush(6)
    end if

    !--- Only Alya rank 0 writes the file ---
    if (arank /= 0) then
       deallocate(sendcounts, displs, full_field)
       return
    end if

    !--- Open file ---
    write(filename,'(A,I6.6,A)') 'alya_grid_', itime, '.vts'
    iunit = 99
    open(unit=iunit, file=trim(filename), status='replace', action='write', iostat=ierr_loc)
    if (ierr_loc /= 0) then
       write(*,*) 'ERROR opening VTS file: ', trim(filename)
       deallocate(sendcounts, displs, full_field)
       return
    end if
    !=========================================================================
    ! VTK StructuredGrid XML
    !=========================================================================
    write(iunit,'(A)') '<?xml version="1.0"?>'
    write(iunit,'(A)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iunit,'(A,E16.8,A)') '  <!-- Time = ', time, ' -->'

    write(iunit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
         '  <StructuredGrid WholeExtent="', &
         0, ' ', rem_nx(1)-1, ' ', &
         0, ' ', rem_nx(2)-1, ' ', &
         0, ' ', rem_nx(3)-1, '">'

    write(iunit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
         '    <Piece Extent="', &
         0, ' ', rem_nx(1)-1, ' ', &
         0, ' ', rem_nx(2)-1, ' ', &
         0, ' ', rem_nx(3)-1, '">'

    !--- Points: compute coordinates from grid indices in VTS order ---
    write(iunit,'(A)') '      <Points>'
    write(iunit,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'

    do iz = 0, rem_nx(3)-1
       do iy = 0, rem_nx(2)-1
          do ix = 0, rem_nx(1)-1
             x = dble(rem_min(1)) + ix * dx
             y = dble(rem_min(2)) + iy * dy
             z = dble(rem_min(3)) + iz * dz
             write(iunit,'(3(E16.8,1X))') x, y, z
          end do
       end do
    end do

    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '      </Points>'

    !--- Point data: one DataArray per equation variable ---
    !
    ! For each point k in the gathered arrays, read its coordinates from
    write(iunit,'(A)') '      <PointData>'

    allocate(grid_var(0:rem_nx(1)-1, 0:rem_nx(2)-1, 0:rem_nx(3)-1))

    do ieq = 1, neqs

       ! Variable name
       select case(ieq)
          case(1); varname = 'rho'
          case(2); varname = 'u'
          case(3); varname = 'v'
          case(4); varname = 'theta'
          case default
             write(varname,'(A,I0)') 'var', ieq
       end select

       write(iunit,'(A,A,A)') &
            '        <DataArray type="Float64" Name="', trim(varname), '" format="ascii">'

       ! Map each global point index k to its (ix,iy,iz) grid position directly.
       ! Points are numbered with ix fastest, then iy, then iz (VTS order).
       grid_var = 0.0d0
       do k = 0, nmax - 1
          iz = k / (rem_nx(1) * rem_nx(2))
          iy = (k - iz * rem_nx(1) * rem_nx(2)) / rem_nx(1)
          ix = mod(k, rem_nx(1))
          grid_var(ix, iy, iz) = full_field(k * neqs + ieq)
       end do

       ! Write in VTS order (ix fastest, then iy, then iz)
       do iz = 0, rem_nx(3)-1
          do iy = 0, rem_nx(2)-1
             do ix = 0, rem_nx(1)-1
                write(iunit,'(E16.8)') grid_var(ix, iy, iz)
             end do
          end do
       end do

       write(iunit,'(A)') '        </DataArray>'

    end do  ! ieq

    deallocate(grid_var)

    write(iunit,'(A)') '      </PointData>'
    write(iunit,'(A)') '    </Piece>'
    write(iunit,'(A)') '  </StructuredGrid>'
    write(iunit,'(A)') '</VTKFile>'

    close(iunit)

    write(*,'(A,A)') 'Written VTS file: ', trim(filename)
    flush(6)

    deallocate(sendcounts, displs, full_field)

  end subroutine write_alya_grid_vts

end program unitt_alya_with_another_code
