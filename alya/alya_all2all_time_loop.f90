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

  !---------------------------------------------------------------------------
  ! COUPLING ARRAYS
  !
  !  npoin_send(i) : always 0 (Alya does not initiate data); used only for Alltoall.
  !  npoin_recv(i) : after Alltoall = pts that Julia world rank i sends TO Alya.
  !
  !  The exchange is BIDIRECTIONAL with the same point set:
  !    Julia  -> Alya : npoin_recv(i) interpolated values  (Irecv in Alya)
  !    Alya   -> Julia: npoin_recv(i) Alya field values    (Isend from Alya)
  !---------------------------------------------------------------------------
  integer(4),    allocatable         :: npoin_send(:)
  integer(4),    allocatable         :: npoin_recv(:)
  integer(4)                         :: total_pts   ! = sum(npoin_recv)

  !---------------------------------------------------------------------------
  ! TAG PROTOCOL (must match Julia identically):
  !
  !   Direction Julia -> Alya:
  !     Julia  posts Isend  with tag = TAG_DATA + alya_world_rank
  !     Alya   posts Irecv  with tag = TAG_DATA + rank   (my world rank)  ✓
  !
  !   Direction Alya -> Julia:
  !     Alya   posts Isend  with tag = TAG_DATA + julia_world_rank  (= i)
  !     Julia  posts Irecv  with tag = TAG_DATA + wrank (its world rank)   ✓
  !---------------------------------------------------------------------------
  integer, parameter                 :: TAG_DATA = 2000

  integer(4)                         :: nranks_julia
  integer(4)                         :: step, nsteps
  real(kind=8)                       :: t0, dt, tend, t
  real(kind=8), allocatable          :: recvbuf_all(:)   ! interpolated field from Julia -> Alya
  real(kind=8), allocatable          :: sendbuf_all(:)   ! Alya field -> Julia
  real(kind=8)                       :: dummy_field

#ifdef USEMPIF08
  type(MPI_Status),  allocatable     :: recv_status(:), send_status(:)
  type(MPI_Request), allocatable     :: recv_requests(:), send_requests(:)
#else
  integer,           allocatable     :: recv_status(:,:), send_status(:,:)
  integer,           allocatable     :: recv_requests(:), send_requests(:)
#endif
  integer                            :: nactive, ireq

  ! Coordinate helpers
  integer                            :: nmax, r_rem, npoin_per_rank
  integer                            :: i_start, i_end, npoin_local
  integer                            :: ipoin, ix, iy, iz
  real(kind=8)                       :: xc, yc, zc, dx, dy, dz

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
  ! STEP 2: GRID METADATA
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

  call MPI_Bcast(ndime, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  do idime = 1, 3
     call MPI_Bcast(rem_min(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_max(idime), 1, MPI_REAL,    0, MPI_COMM_WORLD, ierr)
     call MPI_Bcast(rem_nx(idime),  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  end do

  neqs = 4   ! CompEuler 2D: rho, u, v, theta  <-- Alya sets this, broadcasts it
  call MPI_Bcast(neqs, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  !nsteps = int((tend - t0) / dt)
  !call MPI_Bcast(nsteps, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  
  !if (rank == 0) then
  !write(*,'(A,I0)') 'neqs = ', neqs
  !flush(6)
  !end if

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

  ! Diagnostics
  write(*,'(A)') repeat('=',60)
  write(*,'(A,I0,A,I0)') 'Alya local rank ', arank, '  world rank ', rank
  write(*,'(A,I0)') '  Pts to RECV from Julia: ', total_pts
  write(*,'(A,I0)') '  Pts to SEND to   Julia: ', total_pts, ' (bidirectional, same pts)'
  do i = 0, size-1
     if (npoin_recv(i) > 0) then
        if (i < asize) then
           write(*,'(A,I0,A,I0,A,I0,A)') '  <-> world ', i, ' (Alya ', i, '): ', npoin_recv(i), ' pts'
        else
           write(*,'(A,I0,A,I0,A,I0,A)') '  <-> world ', i, ' (Julia ', i-asize, '): ', npoin_recv(i), ' pts'
        end if
     end if
  end do
  write(*,'(A)') repeat('=',60)
  flush(6)

  !--------------------------------------------------------------------------
  ! Pre-compute point range for this Alya rank
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

  ! Buffers: same size in both directions (total_pts coupling points)
  allocate(recvbuf_all(max(1, total_pts * neqs)))
  allocate(sendbuf_all(max(1, total_pts * neqs)))
  recvbuf_all = 0.0d0
  sendbuf_all = 0.0d0

  ! One request per active partner, for both recv and send
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

     !--- Alya physics (placeholder) ----------------------------------------
     dummy_field = dummy_field + 0.01d0 * dble(step)

     !--- Fill send buffer with Alya's current field -------------------------
     ! In production: copy actual Alya solution at the coupling points here.
     sendbuf_all(1:total_pts * neqs) = dummy_field

     !------------------------------------------------------------------------
     ! BIDIRECTIONAL EXCHANGE
     !   Rule: post ALL Irecv BEFORE any Isend to avoid any risk of deadlock.
     !
     !   Alya Irecv (Julia -> Alya):  tag = TAG_DATA + rank  (my world rank)
     !   Alya Isend (Alya  -> Julia): tag = TAG_DATA + i     (dest world rank)
     !------------------------------------------------------------------------

     ! --- Post all receives ---
     ireq = 0
     do i = 0, size-1
        if (npoin_recv(i) > 0) then
           ireq = ireq + 1
#ifdef USEMPIF08
           call MPI_Irecv(recvbuf_all(1), npoin_recv(i) * neqs, MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + rank, MPI_COMM_WORLD, recv_requests(ireq), ierr)
#else
           call MPI_IRECV(recvbuf_all(1), npoin_recv(i) * neqs, MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + rank, MPI_COMM_WORLD, recv_requests(ireq), ierr)
#endif
        end if
     end do

     ! --- Post all sends ---
     ireq = 0
     do i = 0, size-1
        if (npoin_recv(i) > 0) then     ! Same partners, same point count
           ireq = ireq + 1
#ifdef USEMPIF08
           call MPI_Isend(sendbuf_all(1), npoin_recv(i), MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + i, MPI_COMM_WORLD, send_requests(ireq), ierr)
#else
           call MPI_ISEND(sendbuf_all(1), npoin_recv(i), MPI_DOUBLE_PRECISION, &
                i, TAG_DATA + i, MPI_COMM_WORLD, send_requests(ireq), ierr)
#endif
        end if
     end do

     ! --- Wait for receives ---
     if (nactive > 0) then
#ifdef USEMPIF08
        call MPI_Waitall(nactive, recv_requests, recv_status, ierr)
#else
        call MPI_WAITALL(nactive, recv_requests, recv_status, ierr)
#endif
        if (ierr /= 0) write(*,'(A,I0,A,I0)') &
             'ERROR MPI_Waitall recv rank ', rank, ' step ', step
     end if

     ! --- Wait for sends ---
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
     ! PRINT: received interpolated data at grid coordinates
     !------------------------------------------------------------------------
     if (total_pts > 0) then
        write(*,'(A,I0,A,F10.4,A,I0,A,I0)') &
             '--- Recv interp data  step=', step, '  t=', t, &
             '  Alya_rank=', arank, '  world_rank=', rank
        write(*,'(A8,A6,A6,A18,A18,A18)') 'GlobalID','ix','iy','x','y','field'
        write(*,'(A)') repeat('-', 80)
        do ipoin = i_start, i_end
           iz = ipoin / (rem_nx(1)*rem_nx(2))
           iy = (ipoin - iz*rem_nx(1)*rem_nx(2)) / rem_nx(1)
           ix = mod(ipoin, rem_nx(1))
           xc = dble(rem_min(1)) + ix*dx
           yc = dble(rem_min(2)) + iy*dy
           i  = ipoin - i_start + 1    ! 1-based local offset into recvbuf_all

           write(*,'(I8,I6,I6,E18.8,E18.8,E18.8)') &
                ipoin, ix, iy, xc, yc, recvbuf_all((i-1)*neqs + 1)
        end do
        write(*,'(A)') repeat('=', 80)
        flush(6)
     end if

     !------------------------------------------------------------------------
     ! VTU OUTPUT
     !------------------------------------------------------------------------
     call write_alya_grid_vtu(rank, arank, asize, ndime, rem_min, rem_max, rem_nx, &
          recvbuf_all, total_pts, neqs, step, t)
     
     call MPI_Barrier(PAR_COMM_FINAL, ierr)

  end do  ! End time loop
  
  ! After end do (time loop)
  if (rank == 0) then
     write(*,'(A)') 'Alya: time loop complete, waiting for Julia...'
     flush(6)
  end if
  call MPI_Barrier(MPI_COMM_WORLD, ierr)   ! <-- add this
  
  ! Cleanup
  if (allocated(recvbuf_all))   deallocate(recvbuf_all)
  if (allocated(sendbuf_all))   deallocate(sendbuf_all)
  if (allocated(recv_requests)) deallocate(recv_requests)
  if (allocated(send_requests)) deallocate(send_requests)
  if (allocated(recv_status))   deallocate(recv_status)
  if (allocated(send_status))   deallocate(send_status)
  deallocate(npoin_send, npoin_recv)
  deallocate(alya_to_world, alya_to_world_snd)
 ! call MPI_Finalize(ierr)
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

  !---------------------------------------------------------------------------
  subroutine write_alya_grid_vtu(rank, arank, asize, ndime, rem_min, rem_max, rem_nx, &
                                 recvbuf, nvals, neqs, timestep, time)
    implicit none
    integer(4), intent(in) :: rank, arank, asize, ndime
    real,       intent(in) :: rem_min(3), rem_max(3)
    integer,    intent(in) :: rem_nx(3)
    real(8),    intent(in) :: recvbuf(:)
    integer(4), intent(in) :: nvals, timestep
    real(8),    intent(in) :: time
    integer(4), intent(in) :: neqs   ! add this argument

    character(len=256) :: filename
    integer :: iunit, ierr_local, i, j, k, ipoin
    integer :: i_start, i_end, npoin_local
    real(8) :: x, y, z, dx, dy, dz, fv
    integer :: nmax, r, np

    nmax = rem_nx(1)*rem_nx(2)*rem_nx(3)
    dx = 0.0d0; dy = 0.0d0; dz = 0.0d0
    if (rem_nx(1)>1) dx = dble(rem_max(1)-rem_min(1))/dble(rem_nx(1)-1)
    if (rem_nx(2)>1) dy = dble(rem_max(2)-rem_min(2))/dble(rem_nx(2)-1)
    if (rem_nx(3)>1) dz = dble(rem_max(3)-rem_min(3))/dble(rem_nx(3)-1)

    r  = mod(nmax, asize)
    np = nmax / asize
    if (arank < r) then
       i_start = arank*(np+1);  i_end = i_start + np
    else
       i_start = r*(np+1) + (arank-r)*np;  i_end = i_start + np - 1
    end if
    npoin_local = i_end - i_start + 1

    write(filename,'(A,I6.6,A,I4.4,A)') 'alya_grid_', timestep, '_rank', rank, '.vtu'
    iunit = 100 + rank
    open(unit=iunit, file=trim(filename), status='replace', action='write', iostat=ierr_local)
    if (ierr_local /= 0) then; write(*,*) 'Error opening: ', trim(filename); return; end if

    write(iunit,'(A)') '<?xml version="1.0"?>'
    write(iunit,'(A)') '<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">'
    write(iunit,'(A,E16.8,A)') '  <!-- Time = ', time, ' -->'
    write(iunit,'(A)') '  <UnstructuredGrid>'
    write(iunit,'(A,I0,A,I0,A)') '    <Piece NumberOfPoints="', npoin_local, &
         '" NumberOfCells="', max(0,npoin_local), '">'
    write(iunit,'(A)') '      <Points>'
    write(iunit,'(A)') '        <DataArray type="Float64" NumberOfComponents="3" format="ascii">'
    do ipoin = i_start, i_end
       k = ipoin / (rem_nx(1)*rem_nx(2))
       j = (ipoin - k*rem_nx(1)*rem_nx(2)) / rem_nx(1)
       i = mod(ipoin, rem_nx(1))
       x = dble(rem_min(1)) + i*dx;  y = dble(rem_min(2)) + j*dy;  z = dble(rem_min(3)) + k*dz
       write(iunit,'(3(E16.8,1X))') x, y, z
    end do
    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '      </Points>'
    write(iunit,'(A)') '      <Cells>'
    write(iunit,'(A)') '        <DataArray type="Int32" Name="connectivity" format="ascii">'
    do ipoin = 0, npoin_local-1; write(iunit,'(I0)') ipoin; end do
    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '        <DataArray type="Int32" Name="offsets" format="ascii">'
    do ipoin = 1, npoin_local; write(iunit,'(I0)') ipoin; end do
    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '        <DataArray type="UInt8" Name="types" format="ascii">'
    do ipoin = 1, npoin_local; write(iunit,'(I0)') 1; end do
    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '      </Cells>'
    write(iunit,'(A)') '      <PointData Scalars="received_field">'
    write(iunit,'(A)') '        <DataArray type="Float64" Name="received_field" format="ascii">'
    do ipoin = i_start, i_end
       if (ipoin-i_start+1 <= nvals) then; fv = recvbuf((ipoin - i_start)*neqs + 1); else; fv = 0.0d0; end if
       write(iunit,'(E16.8)') fv
    end do
    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '        <DataArray type="Int32" Name="global_id" format="ascii">'
    do ipoin = i_start, i_end; write(iunit,'(I0)') ipoin; end do
    write(iunit,'(A)') '        </DataArray>'
    write(iunit,'(A)') '      </PointData>'
    write(iunit,'(A)') '    </Piece>'
    write(iunit,'(A)') '  </UnstructuredGrid>'
    write(iunit,'(A)') '</VTKFile>'
    close(iunit)
  end subroutine write_alya_grid_vtu

end program unitt_alya_with_another_code
