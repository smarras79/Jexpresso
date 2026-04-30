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
  integer(4)                         :: nfields
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
  integer(4), allocatable            :: je_gids_all(:)
  integer(4)                         :: je_total_nodes, je_gid_offset

  integer(4)                         :: step, nsteps
  real(kind=8)                       :: t0, dt, tend, t

  ! Output scheduling
  real(kind=8)                       :: out_dt, out_tstart, out_tend
  real(kind=8)                       :: next_t, tol, dt_step, t_plus
  logical                            :: write_now

  real(kind=8), allocatable          :: recvbuf_all(:)
  ! ordered_buf: field values scattered into [i_start..i_end] order
  ! via je_gids_all.  Updated EVERY timestep (matches original Alya behaviour).
  real(kind=8), allocatable          :: ordered_buf(:)
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
  real(kind=8)                       :: dx, dy, dz

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
  tend   = 1000.0d0
  nsteps = int((tend - t0) / dt)
  out_dt = 100.0d0

  !--------------------------------------------------------------------------
  ! STEP 2: GRID METADATA  (broadcast order must match Julia exactly)
  !--------------------------------------------------------------------------
  rem_min = [-5000.0,     0.0, 0.0]
  rem_max = [ 5000.0, 10000.0, 0.0]
  rem_nx  = [50,      50,      1  ]
  ndime   = 2

  if (rank == 0) then
     write(*,'(A)')        ' ALYA-2-JEXP'
     write(*,'(A,I0)')     'ndime   = ', ndime
     write(*,'(A,3F12.2)') 'rem_min = ', rem_min(1:ndime)
     write(*,'(A,3F12.2)') 'rem_max = ', rem_max(1:ndime)
     write(*,'(A,3I6)')    'rem_nx  = ', rem_nx(1:ndime)
     flush(6)
  end if

  ! 2a. ndime
  call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  ! 2b. rem_min / rem_max / rem_nx (full 3-element arrays)
  do idime = 1, 3
     call MPI_Bcast(rem_min(idime), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_max(idime), 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierr)
  end do

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
  write(*,'(A,I0)') '  Pts to RECV from Julia: ', total_pts
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
  !--------------------------------------------------------------------------
  nranks_julia   = size - asize
  je_total_nodes = sum(npoin_recv(asize:size-1))

  allocate(je_gids_all(max(1, je_total_nodes)))
  je_gids_all   = 0_4
  je_gid_offset = 0

  do i = asize, size-1
     if (npoin_recv(i) > 0) then
        call MPI_Recv(je_gids_all(je_gid_offset + 1), npoin_recv(i), &
                      MPI_INTEGER4, i, 0, &
                      MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)
        write(*,'(A,I0,A,I0,A,I0,A)') &
             '[node_list] Alya world_rank=', rank, &
             ' received node IDs from Julia world_rank=', i, &
             ' (', npoin_recv(i), ' IDs)'
        flush(6)
        je_gid_offset = je_gid_offset + npoin_recv(i)
     end if
  end do

  !--------------------------------------------------------------------------
  ! Pre-compute this rank's point range.
  ! Alya local rank 0 owns no points; work is split over ranks 1..asize-1.
  !--------------------------------------------------------------------------
  nmax           = rem_nx(1) * rem_nx(2) * rem_nx(3)
  nworkers       = max(1, asize - 1)
  r_rem          = mod(nmax, nworkers)
  npoin_per_rank = nmax / nworkers

  i_start     = -1
  i_end       = -2
  npoin_local = 0
  if (arank > 0) then
     iworker = arank - 1
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
  ! nfields: number of field components received per point from Julia.
  ! Set to ndime (one value per spatial dimension).
  ! Julia must send data in field-fastest order: [f1_p1, f2_p1, f1_p2, ...]
  ! regardless of whether f is velocity or coordinates — the layout is the
  ! same either way and Alya does not need to know the difference.
  !--------------------------------------------------------------------------
  nfields = ndime
  if (rank == 0) then
     write(*,'(A,I0)') '  nfields per point from Julia: ', nfields
     flush(6)
  end if

  allocate(recvbuf_all(max(1, total_pts * nfields)))
  recvbuf_all = 0.0d0

  allocate(ordered_buf(max(1, npoin_local * nfields)))
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
  out_tstart = t0
  out_tend   = tend
  tol        = max(1.0d-10, 1.0d-6 * out_dt)
  t          = t0
  next_t     = t0 + out_dt
  if (out_dt <= 0.0d0) next_t = huge(1.0d0)

  itime = 1
  do step = 1, nsteps

     dt_step = dt
     t_plus  = t + dt_step

     write_now = (out_dt > 0.0d0 .and. &
                  t_plus >= next_t - tol .and. &
                  next_t <= out_tend + tol)

     !------------------------------------------------------------------------
     ! MPI RECEIVE: DATA (Julia -> Alya)
     ! nfields = ndime values per point, field-fastest.
     ! Julia decides what those values represent (velocity or coordinates);
     ! the MPI calls here are identical either way.
     !------------------------------------------------------------------------
     ireq        = 0
     recv_offset = 0
     do i = 0, size-1
        if (npoin_recv(i) > 0) then
           ireq = ireq + 1
#ifdef USEMPIF08
           call MPI_Irecv(recvbuf_all(recv_offset + 1), npoin_recv(i) * nfields, &
                MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, &
                recv_requests(ireq), ierr)
#else
           call MPI_IRECV(recvbuf_all(recv_offset + 1), npoin_recv(i) * nfields, &
                MPI_DOUBLE_PRECISION, i, 0, MPI_COMM_WORLD, &
                recv_requests(ireq), ierr)
#endif
           recv_offset = recv_offset + npoin_recv(i) * nfields
        end if
     end do

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
     ! SCATTER: recvbuf_all -> ordered_buf  (every timestep, unconditionally)
     !------------------------------------------------------------------------
     ordered_buf = 0.0d0
     do jpt = 1, total_pts
        local_pos = int(je_gids_all(jpt)) - 1 - i_start
        if (local_pos >= 0 .and. local_pos < npoin_local) then
           ordered_buf(local_pos * nfields + 1 : (local_pos + 1) * nfields) = &
                recvbuf_all((jpt - 1) * nfields + 1 : jpt * nfields)
        end if
     end do

     if (write_now) then
        next_t = next_t + out_dt
        call write_alya_grid_vts(arank, asize, PAR_COMM_FINAL, ndime, &
             rem_min, rem_max, rem_nx, &
             ordered_buf, npoin_local, nfields, step, itime, t_plus)
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

  if (allocated(recvbuf_all))   deallocate(recvbuf_all)
  if (allocated(ordered_buf))   deallocate(ordered_buf)
  if (allocated(recv_requests)) deallocate(recv_requests)
  if (allocated(recv_status))   deallocate(recv_status)
  if (allocated(je_gids_all))   deallocate(je_gids_all)
  deallocate(npoin_send, npoin_recv)
  deallocate(alya_to_world, alya_to_world_snd)
  call MPI_Finalize(ierr)

contains

  pure function cstr_trim(str) result(out)
    character(len=*), intent(in) :: str
    character(len=len(str))      :: out
    integer :: nulpos
    nulpos = index(str, char(0))
    if (nulpos > 0) then; out = str(:nulpos-1); else; out = str; end if
    out = trim(adjustl(out))
  end function cstr_trim

  subroutine write_alya_grid_vts(arank, asize, PAR_COMM, ndime, &
                                  rem_min, rem_max, rem_nx, &
                                  recvbuf, total_pts, nfields, &
                                  timestep, itime, time)
    implicit none
    integer(4),   intent(in) :: arank, asize, ndime, itime
    MY_MPI_COMM,  intent(in) :: PAR_COMM
    real(8),      intent(in) :: rem_min(3), rem_max(3)
    integer,      intent(in) :: rem_nx(3)
    real(8),      intent(in) :: recvbuf(:)
    integer(4),   intent(in) :: total_pts, nfields, timestep
    real(8),      intent(in) :: time

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
    character(len=16)           :: varname

    nmax         = rem_nx(1) * rem_nx(2) * rem_nx(3)
    nworkers_loc = max(1, asize - 1)
    r_rem_loc    = mod(nmax, nworkers_loc)
    np_base      = nmax / nworkers_loc

    dx = 0.0d0; dy = 0.0d0; dz = 0.0d0
    if (rem_nx(1) > 1) dx = dble(rem_max(1)-rem_min(1)) / dble(rem_nx(1)-1)
    if (rem_nx(2) > 1) dy = dble(rem_max(2)-rem_min(2)) / dble(rem_nx(2)-1)
    if (rem_nx(3) > 1) dz = dble(rem_max(3)-rem_min(3)) / dble(rem_nx(3)-1)

    allocate(sendcounts(0:asize-1))
    allocate(displs    (0:asize-1))

    do irank = 0, asize-1
       if (irank == 0) then
          sendcounts(irank) = 0
          displs(irank)     = 0
       else
          iworker_loc = irank - 1
          if (iworker_loc < r_rem_loc) then
             i_start_r = iworker_loc * (np_base + 1)
             count_r   = np_base + 1
          else
             i_start_r = r_rem_loc * (np_base + 1) + (iworker_loc - r_rem_loc) * np_base
             count_r   = np_base
          end if
          sendcounts(irank) = count_r * nfields
          displs(irank)     = i_start_r * nfields
       end if
    end do

    count_r = sendcounts(arank)

    if (arank == 0) then
       allocate(full_field(nmax * nfields))
       full_field = 0.0d0
    else
       allocate(full_field(1))
    end if

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
       write(*,'(A,I0,A,I0)') 'ERROR MPI_Gatherv: arank=', arank, ' step=', timestep
       flush(6)
    end if

    if (arank /= 0) then
       deallocate(sendcounts, displs, full_field)
       return
    end if

    write(filename,'(A,I6.6,A)') 'alya_grid_', itime, '.vts'
    iunit = 99
    open(unit=iunit, file=trim(filename), status='replace', action='write', iostat=ierr_loc)
    if (ierr_loc /= 0) then
       write(*,*) 'ERROR opening VTS file: ', trim(filename)
       deallocate(sendcounts, displs, full_field)
       return
    end if

    write(iunit,'(A)') '<?xml version="1.0"?>'
    write(iunit,'(A)') '<VTKFile type="StructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iunit,'(A,E16.8,A)') '  <!-- Time = ', time, ' -->'

    write(iunit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
         '  <StructuredGrid WholeExtent="', &
         0, ' ', rem_nx(1)-1, ' ', 0, ' ', rem_nx(2)-1, ' ', 0, ' ', rem_nx(3)-1, '">'

    write(iunit,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
         '    <Piece Extent="', &
         0, ' ', rem_nx(1)-1, ' ', 0, ' ', rem_nx(2)-1, ' ', 0, ' ', rem_nx(3)-1, '">'

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

    write(iunit,'(A)') '      <PointData>'
    allocate(grid_var(0:rem_nx(1)-1, 0:rem_nx(2)-1, 0:rem_nx(3)-1))

    do ieq = 1, nfields
       select case(ieq)
          case(1); varname = 'var1'
          case(2); varname = 'var2'
          case(3); varname = 'var3'
          case default; write(varname,'(A,I0)') 'var', ieq
       end select

       write(iunit,'(A,A,A)') &
            '        <DataArray type="Float64" Name="', trim(varname), '" format="ascii">'

       grid_var = 0.0d0
       do k = 0, nmax - 1
          iz = k / (rem_nx(1) * rem_nx(2))
          iy = (k - iz * rem_nx(1) * rem_nx(2)) / rem_nx(1)
          ix = mod(k, rem_nx(1))
          grid_var(ix, iy, iz) = full_field(k * nfields + ieq)
       end do

       do iz = 0, rem_nx(3)-1
          do iy = 0, rem_nx(2)-1
             do ix = 0, rem_nx(1)-1
                write(iunit,'(E16.8)') grid_var(ix, iy, iz)
             end do
          end do
       end do

       write(iunit,'(A)') '        </DataArray>'
    end do

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
