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
  integer                            :: i, idime, ndime
  integer(4)                         :: neqs
  character(len=128)                 :: s

  real,    dimension(1:3)            :: rem_min, rem_max
  integer, dimension(1:3)            :: rem_nx

  integer(4),    contiguous, pointer :: alya_to_world(:)
  integer(4),    contiguous, pointer :: alya_to_world_snd(:)

  integer(4),    allocatable         :: npoin_send(:)
  integer(4),    allocatable         :: npoin_recv(:)
  integer(4)                         :: total_pts

  integer, parameter                 :: TAG_DATA = 2000

  integer(4)                         :: nranks_julia
  integer(4)                         :: step, nsteps
  real(kind=8)                       :: t0, dt, tend, t
  real(kind=8), allocatable          :: recvbuf_all(:)
  real(kind=8), allocatable          :: sendbuf_all(:)
  real(kind=8)                       :: dummy_field

#ifdef USEMPIF08
  type(MPI_Status),  allocatable     :: recv_status(:), send_status(:)
  type(MPI_Request), allocatable     :: recv_requests(:), send_requests(:)
#else
  integer,           allocatable     :: recv_status(:,:), send_status(:,:)
  integer,           allocatable     :: recv_requests(:), send_requests(:)
#endif
  integer                            :: nactive, ireq
  integer                            :: recv_offset, send_offset

  integer                            :: nmax, r_rem, npoin_per_rank
  integer                            :: i_start, i_end, npoin_local
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
  dt     = 0.25d0
  tend   = 400.0d0
  nsteps = int((tend - t0) / dt)

  !--------------------------------------------------------------------------
  ! STEP 2: GRID METADATA  (broadcast order must match Julia exactly)
  !--------------------------------------------------------------------------
  rem_min = [-5000.0,     0.0, 0.0]
  rem_max = [ 5000.0, 10000.0, 0.0]
  rem_nx  = [41,      41,        1]
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
  ! 2b. rem_min / rem_max / rem_nx  (one scalar per call, matching Julia)
  do idime = 1, 3
     call MPI_Bcast(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end do
  ! 2c. neqs
  neqs = 4
  call MPI_Bcast(neqs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  ! 2d. nsteps
  call MPI_Bcast(nsteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

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
  ! Pre-compute this rank's point range (for screen printing)
  !--------------------------------------------------------------------------
  nmax           = rem_nx(1) * rem_nx(2) * rem_nx(3)
  r_rem          = mod(nmax, asize)
  npoin_per_rank = nmax / asize

  if (arank < r_rem) then
     i_start = arank * (npoin_per_rank + 1)
     i_end   = i_start + npoin_per_rank
  else
     i_start = r_rem * (npoin_per_rank + 1) + (arank - r_rem) * npoin_per_rank
     i_end   = i_start + npoin_per_rank - 1
  end if
  npoin_local = i_end - i_start + 1

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

  allocate(recvbuf_all(max(1, total_pts * neqs)))
  allocate(sendbuf_all(max(1, total_pts * neqs)))
  recvbuf_all = 0.0d0
  sendbuf_all = 0.0d0

  nactive = count(npoin_recv > 0)

  if (nactive > 0) then
#ifdef USEMPIF08
     allocate(recv_requests(nactive), recv_status(nactive))
     allocate(send_requests(nactive), send_status(nactive))
#else
     allocate(recv_requests(nactive), recv_status(MPI_STATUS_SIZE, nactive))
     allocate(send_requests(nactive), send_status(MPI_STATUS_SIZE, nactive))
#endif
  end if

  dummy_field = 42.0d0 + dble(arank)

  !==========================================================================
  ! TIME LOOP
  !==========================================================================
  do step = 1, nsteps
     t = t0 + step * dt

     dummy_field = dummy_field + 0.01d0 * dble(step)
     sendbuf_all(1:total_pts * neqs) = dummy_field

     !------------------------------------------------------------------------
     ! BIDIRECTIONAL MPI EXCHANGE  (all Irecv before all Isend)
     !------------------------------------------------------------------------

     ! --- Post all receives (Julia -> Alya) ---
     ! Each sender's data goes to a distinct offset in recvbuf_all so that
     ! multiple senders do not overwrite each other.
     ireq = 0
     recv_offset = 0
     do i = 0, size-1
        if (npoin_recv(i) > 0) then
           ireq = ireq + 1
#ifdef USEMPIF08
           call MPI_Irecv(recvbuf_all(recv_offset + 1), npoin_recv(i) * neqs, &
                MPI_DOUBLE_PRECISION, i, TAG_DATA + rank, MPI_COMM_WORLD, &
                recv_requests(ireq), ierr)
#else
           call MPI_IRECV(recvbuf_all(recv_offset + 1), npoin_recv(i) * neqs, &
                MPI_DOUBLE_PRECISION, i, TAG_DATA + rank, MPI_COMM_WORLD, &
                recv_requests(ireq), ierr)
#endif
           recv_offset = recv_offset + npoin_recv(i) * neqs
        end if
     end do

     ! --- Post all sends (Alya -> Julia) ---
     ireq = 0
     send_offset = 0
     do i = 0, size-1
        if (npoin_recv(i) > 0) then
           ireq = ireq + 1
#ifdef USEMPIF08
           call MPI_Isend(sendbuf_all(send_offset + 1), npoin_recv(i) * neqs, &
                MPI_DOUBLE_PRECISION, i, TAG_DATA + i, MPI_COMM_WORLD, &
                send_requests(ireq), ierr)
#else
           call MPI_ISEND(sendbuf_all(send_offset + 1), npoin_recv(i) * neqs, &
                MPI_DOUBLE_PRECISION, i, TAG_DATA + i, MPI_COMM_WORLD, &
                send_requests(ireq), ierr)
#endif
           send_offset = send_offset + npoin_recv(i) * neqs
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

     if (nactive > 0) then
#ifdef USEMPIF08
        call MPI_Waitall(nactive, send_requests, send_status, ierr)
#else
        call MPI_WAITALL(nactive, send_requests, send_status, ierr)
#endif
        if (ierr /= 0) write(*,'(A,I0,A,I0)') &
             'ERROR MPI_Waitall send rank ', rank, ' step ', step
     end if

     !------------------------------------------------------------------------
     ! PRINT: received interpolated data at coordinates
     !------------------------------------------------------------------------
     if (total_pts > 0) then
        write(*,'(A,I0,A,F10.4,A,I0,A,I0)') &
             '--- Recv step=', step, '  t=', t, &
             '  Alya_rank=', arank, '  world_rank=', rank
        write(*,'(A8,A6,A6,A18,A18,A18)') 'GlobalID','ix','iy','x','y','rho'
        write(*,'(A)') repeat('-', 80)
        do ipoin = i_start, i_end
           iz = ipoin / (rem_nx(1)*rem_nx(2))
           iy = (ipoin - iz*rem_nx(1)*rem_nx(2)) / rem_nx(1)
           ix = mod(ipoin, rem_nx(1))
           xc = dble(rem_min(1)) + ix*dx
           yc = dble(rem_min(2)) + iy*dy
           i  = ipoin - i_start + 1
           write(*,'(I8,I6,I6,E18.8,E18.8,E18.8)') &
                ipoin, ix, iy, xc, yc, recvbuf_all((i-1)*neqs + 1)
        end do
        write(*,'(A)') repeat('=', 80)
        flush(6)
     end if

     !------------------------------------------------------------------------
     ! VTS OUTPUT: structured grid with filled contours
     !   - all Alya ranks gather to rank 0 via MPI_Gatherv (PAR_COMM_FINAL)
     !   - rank 0 writes a single .vts file per timestep
     !------------------------------------------------------------------------
     call write_alya_grid_vts(arank, asize, PAR_COMM_FINAL, ndime, &
          rem_min, rem_max, rem_nx, recvbuf_all, total_pts, neqs, step, t)

     call MPI_Barrier(PAR_COMM_FINAL, ierr)

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
  if (allocated(sendbuf_all))   deallocate(sendbuf_all)
  if (allocated(recv_requests)) deallocate(recv_requests)
  if (allocated(send_requests)) deallocate(send_requests)
  if (allocated(recv_status))   deallocate(recv_status)
  if (allocated(send_status))   deallocate(send_status)
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
  !    1. Each Alya rank contributes its coupling points (i_start..i_end)
  !       to a MPI_Gatherv within PAR_COMM_FINAL.
  !    2. Alya rank 0 assembles the full field and writes one .vts file.
  !
  !  Data ordering:
  !    The Julia side sorts extracted coordinates by (owner_rank, global_id)
  !    so that each Alya rank receives values in ascending flat-index order.
  !    After MPI_Gatherv, full_field is indexed as:
  !       full_field( ipoin_0 * neqs + ieq )
  !    where ipoin_0 = ix + iy*nx1 + iz*nx1*nx2  (0-based flat index)
  !    and   ieq = 1..neqs                        (1-based equation index)
  !
  !  The VTS writer uses the coordinate position (ix, iy, iz) to compute
  !  the flat index for each grid point, ensuring coordinates and variable
  !  values are written in the correct VTS order (ix fastest).
  !
  !  Each variable (rho, u, v, theta) is written as a separate DataArray
  !  so ParaView can visualise them independently.
  !===========================================================================
  subroutine write_alya_grid_vts(arank, asize, PAR_COMM, ndime, &
                                  rem_min, rem_max, rem_nx, &
                                  recvbuf, total_pts, neqs, timestep, time)
    implicit none
    integer(4),   intent(in) :: arank, asize, ndime
    MY_MPI_COMM,  intent(in) :: PAR_COMM
    real,         intent(in) :: rem_min(3), rem_max(3)
    integer,      intent(in) :: rem_nx(3)
    real(8),      intent(in) :: recvbuf(:)
    integer(4),   intent(in) :: total_pts, neqs, timestep
    real(8),      intent(in) :: time

    ! Local variables
    integer          :: nmax, r_rem_loc, np_base
    integer          :: i_start_r, count_r
    integer          :: irank, iunit, ierr_loc, ierr_mpi
    integer          :: ix, iy, iz, ipoin_0, ieq
    real(8)          :: x, y, z, dx, dy, dz
    character(len=256)          :: filename
    integer,  allocatable       :: sendcounts(:), displs(:)
    real(8),  allocatable       :: full_field(:)
    real(8),  allocatable       :: grid_var(:,:,:)

    ! Variable names for the neqs fields
    character(len=16) :: varname

    !--- Grid parameters ---
    nmax     = rem_nx(1) * rem_nx(2) * rem_nx(3)
    r_rem_loc = mod(nmax, asize)
    np_base  = nmax / asize

    dx = 0.0d0; dy = 0.0d0; dz = 0.0d0
    if (rem_nx(1) > 1) dx = dble(rem_max(1)-rem_min(1)) / dble(rem_nx(1)-1)
    if (rem_nx(2) > 1) dy = dble(rem_max(2)-rem_min(2)) / dble(rem_nx(2)-1)
    if (rem_nx(3) > 1) dz = dble(rem_max(3)-rem_min(3)) / dble(rem_nx(3)-1)

    !--- Compute Gatherv arrays (all ranks compute, Gatherv ignores on non-root) ---
    allocate(sendcounts(0:asize-1))
    allocate(displs    (0:asize-1))

    do irank = 0, asize-1
       if (irank < r_rem_loc) then
          i_start_r = irank * (np_base + 1)
          count_r   = np_base + 1
       else
          i_start_r = r_rem_loc * (np_base + 1) + (irank - r_rem_loc) * np_base
          count_r   = np_base
       end if
       sendcounts(irank) = count_r * neqs   ! doubles this rank sends
       displs(irank)     = i_start_r * neqs ! displacement in full_field
    end do

    !--- This rank's local send count ---
    count_r = sendcounts(arank)

    !--- Allocate receive buffer only on root ---
    if (arank == 0) then
       allocate(full_field(nmax * neqs))
       full_field = 0.0d0
    else
       allocate(full_field(1))  ! dummy, not used on non-root
    end if

    !--- Gather: each rank sends its coupling recv buffer to Alya rank 0 ---
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
       write(*,'(A,I0,A,I0)') 'ERROR MPI_Gatherv vts: arank=', arank, ' step=', timestep
       flush(6)
    end if

    !--- Only Alya rank 0 writes the file ---
    if (arank /= 0) then
       deallocate(sendcounts, displs, full_field)
       return
    end if

    !--- Open file ---
    write(filename,'(A,I6.6,A)') 'alya_grid_', timestep, '.vts'
    iunit = 99
    open(unit=iunit, file=trim(filename), status='replace', action='write', iostat=ierr_loc)
    if (ierr_loc /= 0) then
       write(*,*) 'ERROR opening VTS file: ', trim(filename)
       deallocate(sendcounts, displs, full_field)
       return
    end if

    !=========================================================================
    ! VTK StructuredGrid XML
    !
    ! WholeExtent / Extent: "i0 i1 j0 j1 k0 k1"  (index bounds, inclusive)
    ! For a 41x41x1 grid: "0 40 0 40 0 0"
    !
    ! Point ordering in VTS: ix varies fastest, then iy, then iz.
    ! This matches our flat index: flat = ix + iy*nx1 + iz*nx1*nx2.
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
    ! Use the coordinate position (ix, iy, iz) to compute the flat index:
    !   ipoin_0 = ix + iy * rem_nx(1) + iz * rem_nx(1) * rem_nx(2)
    ! Then look up the variable value from full_field at:
    !   full_field( ipoin_0 * neqs + ieq )
    !
    ! This explicit coordinate -> flat-index -> value mapping ensures each
    ! grid point's data is written at the position matching its coordinates.
    write(iunit,'(A)') '      <PointData>'

    ! Allocate a 3D work array to scatter interleaved full_field values
    ! into (ix, iy, iz) structured positions for each variable
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

       ! Scatter: use each point's flat index to compute its grid position (ix,iy,iz)
       ! and place the variable value at the correct structured location
       grid_var = 0.0d0
       do ipoin_0 = 0, nmax - 1
          iz = ipoin_0 / (rem_nx(1) * rem_nx(2))
          iy = (ipoin_0 - iz * rem_nx(1) * rem_nx(2)) / rem_nx(1)
          ix = mod(ipoin_0, rem_nx(1))
          grid_var(ix, iy, iz) = full_field(ipoin_0 * neqs + ieq)
       end do

       ! Write in VTS order (ix fastest, then iy, then iz)
       ! using the coordinate positions stored in grid_var
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
