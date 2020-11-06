module spectral_wave_data_shape_4_impl_1_def

use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float

use kind_values, only: knd => kind_swd_interface, wp => kind_swd_internal

use open_swd_file_def, only: open_swd_file, swd_validate_binary_convention, &
                             swd_magic_number
use spectral_wave_data_def, only: spectral_wave_data
use spectral_interpolation_def, only: spectral_interpolation
use swd_version, only: version
use swd_fft_def, only: swd_fft

implicit none
private

! This module provides an implemention of the shape 4 class of the 
! spectral-wave-data API.
!
! Written by Jens Bloch Helmers, August, 19. 2019
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: spectral_wave_data_shape_4_impl_1
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type, extends(spectral_wave_data) :: spectral_wave_data_shape_4_impl_1
    integer            :: nx    ! Number of spectral components in x-direction
    integer            :: ny    ! Number of spectral components in y-direction
    integer            :: nsumx ! Number of spectral components in x-summation (<=nx)
    integer            :: nsumy ! Number of spectral components in y-summation (<=ny)
    real(wp)           :: dkx   ! Constant spacing of wave numbers in x-direction
    real(wp)           :: dky   ! Constant spacing of wave numbers in y-direction
    integer            :: icur     ! 'pointer' to column of most recent spectral input. (1:4)
    integer            :: istp     ! Most recent step from swd file in memory
    integer            :: ipt(4,4) ! 'pointer' to which spectral columns correspond to 
                                   ! indices i-1, i, i+1 and i+2 in the interpolation scheme.
    complex(c_float), allocatable :: c_win(:,:,:)  ! Input window of c-spectral components (1dim=-ny:ny, 2dim=0:nx, 3dim=1:4)
    complex(c_float), allocatable :: ct_win(:,:,:) ! Input window of ct-spectral components (1dim=-ny:ny, 2dim=0:nx, 3dim=1:4)
    complex(c_float), allocatable :: h_win(:,:,:)  ! Input window of h-spectral components (1dim=-ny:ny, 2dim=0:nx, 3dim=1:4)
    complex(c_float), allocatable :: ht_win(:,:,:) ! Input window of ht-spectral components (1dim=-ny:ny, 2dim=0:nx, 3dim=1:4)
    complex(wp), allocatable :: c_cur(:,:,:)    ! Spectral c-values at current time (1dim=1:2, 2dim=0:ny, 3dim=0:nx)
    complex(wp), allocatable :: ct_cur(:,:,:)   ! Spectral ct-values at current time (1dim=1:2, 2dim=0:ny, 3dim=0:nx)
    complex(wp), allocatable :: h_cur(:,:,:)    ! Spectral ct-values at current time (1dim=1:2, 2dim=0:ny, 3dim=0:nx)
    complex(wp), allocatable :: ht_cur(:,:,:)   ! Spectral ht-values at current time (1dim=1:2, 2dim=0:ny, 3dim=0:nx)
    type(spectral_interpolation)  :: tpol   ! Temporal interpolation scheme
contains
    procedure :: close              ! Destructor
    procedure :: update_time        ! Obtain spectral data for current time
    procedure :: phi                ! Calculate potential at location for current time
    procedure :: stream             ! Stream function at location for current time
    procedure :: phi_t              ! Calculate d(potential)/dt (Euler) at location for current time
    procedure :: grad_phi           ! Calculate particle velocity at location for current time
    procedure :: grad_phi_2nd       ! Calculate second order spatial gradients of potential
    procedure :: acc_euler          ! Calculate Euler acceleration (grad(phi_t)) at location for current time
    procedure :: acc_particle       ! Calculate particle acceleration at location for current time
    procedure :: elev               ! Calculate surface elevation at location for current time
    procedure :: elev_t             ! Calculate d(surface elevation)/dt (Euler) at location for current time
    procedure :: grad_elev          ! Calculate gradient of surface elevation at location for current time
    procedure :: grad_elev_2nd      ! Calculate second order spatial gradients of elevation
    procedure :: pressure           ! Fully nonlinear Bernoulli pressure
    procedure :: bathymetry         ! Return local depth at application position (x, y)
    procedure :: bathymetry_nvec    ! Unit normal vector of sea floor into the ocean at (x,y)
    procedure :: convergence        ! For a specific (t,x,y,z) return a csv-file on how velocity, elevation
                                    ! and pressure converge as a function of number of spectral components
    procedure :: strip              ! Create a new SWD file containing only the time steps within a window.
    procedure :: get_int            ! Extract a specified int parameter
    procedure :: get_logical        ! Extract a specified logical parameter
    procedure :: get_real           ! Extract a specified real parameter
    procedure :: get_chr            ! Extract a specified char parameter
    procedure :: elev_fft           ! Surface elevation on a regular grid using FFT 
end type spectral_wave_data_shape_4_impl_1

interface spectral_wave_data_shape_4_impl_1
    module procedure constructor
end interface

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

subroutine close(self)
class(spectral_wave_data_shape_4_impl_1) :: self  ! Object to destruct
!
logical opened
!
inquire(unit=self % unit, opened=opened)
if (opened) close(self % unit)
if (allocated(self % cid)) deallocate(self % cid)
if (allocated(self % c_win)) deallocate(self % c_win)
if (allocated(self % ct_win)) deallocate(self % ct_win)
if (allocated(self % h_win)) deallocate(self % h_win)
if (allocated(self % ht_win)) deallocate(self % ht_win)
if (allocated(self % c_cur)) deallocate(self % c_cur)
if (allocated(self % ct_cur)) deallocate(self % ct_cur)
if (allocated(self % h_cur)) deallocate(self % h_cur)
if (allocated(self % ht_cur)) deallocate(self % ht_cur)
!
self % file = '0'
self % unit = 0
self % nx = 0
self % ny = 0
self % nsumx = 0
self % nsumy = 0
!
end subroutine  close
    
!==============================================================================

function constructor(file, x0, y0, t0, beta, rho, nsumx, nsumy, ipol, norder, &
                     dc_bias) result(self)
character(len=*),    intent(in)  :: file  ! File containing HOSM data
real(knd),           intent(in)  :: x0,y0 ! Spatial seed
real(knd),           intent(in)  :: t0    ! Temporal seed >=0
real(knd),           intent(in)  :: beta  ! Wave heading (deg)
real(knd), optional, intent(in)  :: rho   ! Density of water (applied for pressure calculations)
integer, optional,   intent(in)  :: nsumx ! If present and nsumx>-1: apply nsumx number of spectral
                                          ! components in x-direction, else apply all spectral
                                          ! components from SWD file in summation.
integer, optional,   intent(in)  :: nsumy ! If present and nsumy>-1: apply nsumy number of spectral
                                          ! components in y-direction, else apply all spectral
                                          ! components from SWD file in summation.
integer, optional,   intent(in)  :: ipol  ! Index to request actual temporal interpolation scheme
                                          !   0 = Default (C^2 continous scheme)
                                          !   1 = C^1 continous
                                          !   2 = C^3 continous
integer, optional,   intent(in) :: norder ! Expansion order to apply in kinematics for z>0
                                          ! <0: apply exp(kj z)
                                          ! 0:  apply expansion order specified on swd file (default)
                                          ! >0: apply expansion order = norder 
logical, optional,   intent(in):: dc_bias ! True: apply zero frequency amplitudes from SWD file. 
                                          ! False: Suppress contribution from zero frequency amplitudes (Default)
type(spectral_wave_data_shape_4_impl_1) :: self  ! Object to construct
!
integer :: i, ios, ix, iy, err_id
integer(int64) :: ipos1, ipos2
integer(c_int) :: fmt, shp, amp, nx, ny, order, nid, nsteps, nstrip
real(c_float) :: dkx, dky, dt, grav, lscale, magic
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::constructor'
character(len=250) :: err_msg(5)
complex(wp) :: fval, dfval
real(wp) :: dt_tpol
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
!
call self % error % clear()
!
if (present(rho)) then
    self % rho = rho
else
    self % rho = 1025.0_wp
end if

if ( t0 < 0.0_wp ) then
    err_msg(1) = 'The temporal seed t0 should be zero or positive.'
    write(err_msg(2),'(a,f0.8)') 't0 = ', t0
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:2))
    return
end if
self % t0 = t0
self % x0 = x0
self % y0 = y0
self % file = file

call swd_validate_binary_convention(self % file, err_id, err_msg(2))
if (err_msg(2) /= '') then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    call self % error % set_id_msg(err_proc, err_id, err_msg(1:2))
    return
end if

call open_swd_file(newunit=self % unit, file=self % file, status='old', &
                   as_little_endian=.true., iostat=ios)
if ( ios/=0 ) then
    err_msg(1) = 'Not able to open SWD file:'
    err_msg(2) = self % file
    call self % error % set_id_msg(err_proc, 1001, err_msg(1:2))
    return
end if

!------------------------------------------------------------------------------------
! The swd file format assumes 4 bytes for each integer and real and 1 byte for each
! character. As a consequence we apply temporary C variables for reading variables
! from the swd file. (Size of Fortran variables depends on compilers)
!------------------------------------------------------------------------------------

read(self % unit, end=98, err=99) magic
! magic is already validated in 'swd_validate_binary_convention'

read(self % unit, end=98, err=99) fmt
if ( fmt /= 100 ) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    write(err_msg(2),'(a,i0,a)') 'fmt=', fmt, ' is an unknown swd format parameter.'
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
    return
end if
self % fmt = fmt

read(self % unit, end=98, err=99) shp
self % shp = shp

read(self % unit, end=98, err=99) amp
self % amp = amp

read(self % unit, end=98, err=99) cprog
self % prog = trim(cprog)

read(self % unit, end=98, err=99) cdate 
self % date = trim(cdate)

read(self % unit, end=98, err=99) nid
allocate(character(len=int(nid)) :: cid)
allocate(character(len=int(nid)) :: self % cid)

read(self % unit, end=98, err=99) cid
self % cid = cid(:nid)

read(self % unit, end=98, err=99) grav
self % grav = grav

read(self % unit, end=98, err=99) lscale
self % lscale = lscale

read(self % unit, end=98, err=99) nstrip
self % nstrip = nstrip

read(self % unit, end=98, err=99) nsteps
self % nsteps = nsteps

read(self % unit, end=98, err=99) dt
self % dt = dt

read(self % unit, end=98, err=99) order
self % order = order

read(self % unit, end=98, err=99) nx
self % nx = nx

read(self % unit, end=98, err=99) ny
self % ny = ny

read(self % unit, end=98, err=99) dkx
self % dkx = dkx

read(self % unit, end=98, err=99) dky
self % dky = dky

self % norder = self % order
if (present(norder)) then
    if (norder /= 0) then
        self % norder = norder
    end if
end if
        
if (present(dc_bias)) then
    self % dc_bias = dc_bias
else
    self % dc_bias = .false.
end if

if (present(nsumx)) then
    if (nsumx < 0) then
        self % nsumx = nx
    else
        self % nsumx = nsumx
        if (nsumx > nx) then
            write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
            write(err_msg(2),'(a,i0)') 'Number of spectral components in SWD file (x-dir) = ', nx
            write(err_msg(3),'(a,i0)') 'Number of requested spectral components (x-dir) = ', nsumx
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
            return
        end if
    end if
else
    self % nsumx = nx
end if

if (present(nsumy)) then
    if (nsumy < 0) then
        self % nsumy = ny
    else
        self % nsumy = nsumy
        if (nsumy > ny) then
            write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
            write(err_msg(2),'(a,i0)') 'Number of spectral components in SWD file (y-dir) = ', ny
            write(err_msg(3),'(a,i0)') 'Number of requested spectral components (y-dir) = ', nsumy
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
            return
        end if
    end if
else
    self % nsumy = ny
end if

! make object for FFT-based evaluations
self % fft = swd_fft(self % nsumx, self % nsumy, self % dkx, self % dky)

if (self % nsteps == 1) then
    dt_tpol = 1.0_wp
else
    dt_tpol = self % dt
end if

if (present(ipol)) then
    call self % tpol % construct(ischeme=ipol, delta_t=dt_tpol, ierr=i)
else
    call self % tpol % construct(ischeme=0, delta_t=dt_tpol, ierr=i)
end if
self % ipol = self % tpol % ischeme
if (i /= 0) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a)') "ipol is out of bounds."
    write(err_msg(3),'(a,i0)') 'ipol = ', self % ipol
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if

self % sbeta = sin(beta*pi/180.0_wp)
self % cbeta = cos(beta*pi/180.0_wp)
self % tmax = self % dt * (self % nsteps - 1) - self % t0
if (self % tmax < self % dt) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a)') "Constructor parameter t0 is too large."
    write(err_msg(3),'(a,f0.4)') 't0 = ', self % t0
    write(err_msg(4),'(a,f0.4)') 'Max application time = ', self % tmax
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4))
    return
end if

allocate( self % c_win(-self % ny:self % ny, 0:self % nx, 4),  &
          self % ct_win(-self % ny:self % ny, 0:self % nx, 4), &
          self % h_win(-self % ny:self % ny, 0:self % nx, 4),  &
          self % ht_win(-self % ny:self % ny, 0:self % nx, 4), &
          self % c_cur(2, 0:self % ny, 0:self % nx),           &
          self % ct_cur(2, 0:self % ny, 0:self % nx),          &
          self % h_cur(2, 0:self % ny, 0:self % nx),           &
          self % ht_cur(2, 0:self % ny, 0:self % nx), stat=i)
if (i /= 0) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a)') 'Not able to allocate space for storing window of spectral components.'
    write(err_msg(3),'(a,i0,a,i0,a)') 'Number of spectral components (nx,ny) is (', &
                                      self % nx, ',', self % ny, ')'
    call self % error % set_id_msg(err_proc, 1005, err_msg(1:3))
    return
end if

! The first timestep is put into memory.
associate(c => self % c_win, ct => self % ct_win, h => self % h_win, ht => self % ht_win)
    ! request file position where the temporal functions start
    inquire(self % unit, pos=self % ipos0) 

    ! set to zero intially
    h = czero_c
    ht = czero_c
    c = czero_c
    ct = czero_c

    read(self % unit, end=98, err=99) h(:,:,2)
    read(self % unit, end=98, err=99) ht(:,:,2)
    if (self % amp < 3) then
        read(self % unit, end=98, err=99) c(:,:,2)
        read(self % unit, end=98, err=99) ct(:,:,2)
    end if
    ipos1 = self % ipos0
    inquire(self % unit, pos=ipos2)
    ! Storage fortran units per complex (c_float based)
    if (self % amp == 3) then
        self % size_complex = (ipos2 - ipos1) / (2 * (self%nx + 1) * ( 2 * self%ny + 1))
    else
        self % size_complex = (ipos2 - ipos1) / (4 * (self%nx + 1) * ( 2 * self%ny + 1))
    end if
    self % size_step = (ipos2 - ipos1)

end associate
self % istp = 1  ! The most recent physical step in memory

self % icur = 3  ! The column to store next data. Cycles with repetitons from 1 to 4
! self % ipt(1:4, icur) represent i-1, i, i+1 and i+2 in the interpolation scheme
self % ipt(:,1) = [1,2,3,4]
self % ipt(:,2) = [2,3,4,1]
self % ipt(:,3) = [3,4,1,2]
self % ipt(:,4) = [4,1,2,3]

return
!
98 continue
err_msg(1) = 'End of file when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
return
!
99 continue
err_msg(1) = 'Error when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
!
end function constructor

!==============================================================================

subroutine update_time(self, time)
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self  ! Update data in memory (if needed)
real(knd), intent(in) :: time  ! Current time in simulation program
!
integer(int64) :: ipos
integer :: istp_max, i, ix, iy, i1, i2, i3, i4, istp_min, ios, imove
real(wp) :: delta, teps
complex(wp) :: fval, dfval
complex(c_float) :: cdum
complex(wp), parameter :: czero = cmplx(0.0_wp, 0.0_wp, wp)
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::update_time'
character(len=250) :: err_msg(6)
!
teps = spacing(time) * 10.0_wp
if (time > self % tmax + teps) then
    ! time accounted for round-off is definetly too large...
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'Requested time is too large!'
    write(err_msg(3),'(a,f0.5)') 'User time = ', time
    write(err_msg(4),'(a,f0.5)') 'Max user time = ', self % tmax
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4))
    return
end if
! In case time is _very_ close to self % tmax we make it close to
! infinitesimal less than self % tmax due to the interpolation scheme...
self % tswd = self % t0 + min(time, self % tmax - teps)
if (self % tswd < - teps) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'time corresponds to negative swd time!'
    write(err_msg(3),'(a,f0.5)') 'time user = ', time
    write(err_msg(4),'(a,f0.5)') 'time swd = ', self % tswd
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4))
    return
else if (self % tswd < 0.0_wp) then
    ! Is very close to negative zero. Interpreted as zero in the interpolation scheme...
    self % tswd = 0.0_wp
end if
! From now on, we know that self % tswd is within the acceptable interpolation range.
! If any I/O error will occur it must be due less data than expected
! in the SWD file... (The wave_generator ended prematurely)

! We need to store the 4 time steps: istp_min, istp_min+1, ..., istp_max in memory
! tswd=0.0 corresponds to time step 1. The last step in file is nsteps.
! Minimum time step in memory: =0 indicates need of padding below tswd = 0.0
if (self % nsteps == 1) then
    istp_min = 1
    delta = 0.0_wp
    istp_max = 1
    self % icur = 1
else
    istp_min = int((self % tswd - teps) / self % dt)
    delta = self % tswd / self % dt - istp_min  ! delta in [0.0, 1.0] 
    ! Maximum time step in memory: =nsteps+1 indicates padding beyond tswd_max
    istp_max = istp_min + 3
end if  

associate(c => self % c_win, ct => self % ct_win, h => self % h_win, &
          ht => self % ht_win, ic => self % icur, ip => self % ipt)

    ! Calculate how many time steps we need to shift in the SWD file...
    ! imove < 0 is negative time stepping from user. 
    ! imove > 4 is fast forward and we will skip some time steps
    imove = istp_max - self % istp

    if (imove < 0 .or. imove > 4) then
        ! We need to fill the 4 time step buffer from scratch.
        ! We will reposition the file pointer, but we need to make sure that
        ! 3 time steps prior to istp_max will be stored in memory...
        if (istp_min == 0) then
            ! Wind back to tswd < dt_swd. 
            ! We have to pad data later. Consequently we skip the first column...
            self % istp = 0
            ic = 2 ! ic is a position counter 1, 2, 3 or 4.
        else
            self % istp = istp_min - 1
            ic = 1 ! ic is a position counter 1, 2, 3 or 4.
        end if            
        ipos = self % ipos0 + int(self % istp, int64) * self % size_step
        ! Reposition the file pointer by reading a dummy complex number 
        ! just before the goodies...
        read(self % unit, pos=ipos - self % size_complex, iostat=ios) cdum 
        if (ios/=0) then
            err_msg(1) = 'User time beyond what is available from SWD file:'
            err_msg(2) = self % file
            err_msg(3) = 'The file has less content than expected. Not able to recover!!!'
            write(err_msg(4),'(a,f0.5)') 'Requested swd-time = ', self % tswd
            write(err_msg(5),'(a,f0.5)') 'Requested user-time = ', time
            call self % error % set_id_msg(err_proc, 1003, err_msg(1:5))
            return
        end if
    end if

    do i = 1, istp_max - self % istp
        if (self % istp == self % nsteps) then
            ! We apply padding for data at the very end of the file
            ic = ic + 1
            if (ic > 4) ic = 1
            i1 = ip(1,ic)
            i2 = ip(2,ic)
            i3 = ip(3,ic)
            i4 = ip(4,ic)
            do ix = 0, self % nsumx
                do iy = -self % nsumy, self % nsumy
                    ! Potential and d/dt of potential
                    call self % tpol % pad_right(              &
                                cmplx(c(iy,ix,i1),  kind=wp),  &
                                cmplx(c(iy,ix,i2),  kind=wp),  &
                                cmplx(c(iy,ix,i3),  kind=wp),  &
                                cmplx(ct(iy,ix,i1), kind=wp),  &
                                cmplx(ct(iy,ix,i2), kind=wp),  &
                                cmplx(ct(iy,ix,i3), kind=wp),  &
                                fval, dfval)
                    c(iy,ix,i4) = fval
                    ct(iy,ix,i4) = dfval
                    ! Wave height and d/dt of wave height
                    call self % tpol % pad_right(              &
                                cmplx(h(iy,ix,i1),  kind=wp),  &
                                cmplx(h(iy,ix,i2),  kind=wp),  &
                                cmplx(h(iy,ix,i3),  kind=wp),  &
                                cmplx(ht(iy,ix,i1), kind=wp),  &
                                cmplx(ht(iy,ix,i2), kind=wp),  &
                                cmplx(ht(iy,ix,i3), kind=wp),  &
                                fval, dfval)
                    h(iy,ix,i4) = fval
                    ht(iy,ix,i4) = dfval
                end do
            end do
            self % istp = istp_max
        else
            read(self % unit, end=98, err=99) h(:,:,ic)
            read(self % unit, end=98, err=99) ht(:,:,ic)
            if (self % amp < 3) then
                read(self % unit, end=98, err=99) c(:,:,ic)
                read(self % unit, end=98, err=99) ct(:,:,ic)
            else
                c(:,:,ic) = czero_c
                ct(:,:,ic) = czero_c
            end if
            self % istp = self % istp + 1
            ic = ic + 1
            if (ic > 4) ic = 1
        end if
    end do

    if (istp_min == 0) then
        ! Padding in first column because tswd < dt_swd. 
        do ix = 0, self % nsumx
            do iy = -self % nsumy, self % nsumy
                ! Potential and d/dt of potential
                call self % tpol % pad_left(              &
                            cmplx(c(iy,ix,2), kind=wp),   &
                            cmplx(c(iy,ix,3), kind=wp),   &
                            cmplx(c(iy,ix,4), kind=wp),   &
                            cmplx(ct(iy,ix,2), kind=wp),  &
                            cmplx(ct(iy,ix,3), kind=wp),  &
                            cmplx(ct(iy,ix,4), kind=wp),  &
                            fval, dfval)
                c(iy,ix,1) = fval
                ct(iy,ix,1) = dfval
                ! Wave height and d/dt of wave height
                call self % tpol % pad_left(              &
                            cmplx(h(iy,ix,2), kind=wp),   &
                            cmplx(h(iy,ix,3), kind=wp),   &
                            cmplx(h(iy,ix,4), kind=wp),   &
                            cmplx(ht(iy,ix,2), kind=wp),  &
                            cmplx(ht(iy,ix,3), kind=wp),  &
                            cmplx(ht(iy,ix,4), kind=wp),  &
                            fval, dfval)
                h(iy,ix,1) = fval
                ht(iy,ix,1) = dfval
            end do
        end do
    end if
    ! Finally we do the actual interpolation...
    i1 = ip(1,ic)
    i2 = ip(2,ic)
    i3 = ip(3,ic)
    i4 = ip(4,ic)
    do ix = 0, self % nsumx
        do iy = -self % nsumy, -1
            ! Potential and d/dt of potential
            call self % tpol % scheme(delta,       &
                        cmplx(c(iy,ix,i1), kind=wp),   &
                        cmplx(c(iy,ix,i2), kind=wp),   &
                        cmplx(c(iy,ix,i3), kind=wp),   &
                        cmplx(c(iy,ix,i4), kind=wp),   &
                        cmplx(ct(iy,ix,i1), kind=wp),  &
                        cmplx(ct(iy,ix,i2), kind=wp),  &
                        cmplx(ct(iy,ix,i3), kind=wp),  &
                        cmplx(ct(iy,ix,i4), kind=wp),  &
                        self % c_cur(2,-iy,ix), self % ct_cur(2,-iy,ix))
            ! Wave height and d/dt of wave height
            call self % tpol % scheme(delta,       &
                        cmplx(h(iy,ix,i1), kind=wp),   &
                        cmplx(h(iy,ix,i2), kind=wp),   &
                        cmplx(h(iy,ix,i3), kind=wp),   &
                        cmplx(h(iy,ix,i4), kind=wp),   &
                        cmplx(ht(iy,ix,i1), kind=wp),  &
                        cmplx(ht(iy,ix,i2), kind=wp),  &
                        cmplx(ht(iy,ix,i3), kind=wp),  &
                        cmplx(ht(iy,ix,i4), kind=wp),  &
                        self % h_cur(2,-iy,ix), self % ht_cur(2,-iy,ix))
        end do
        iy = 0
        ! Potential and d/dt of potential
        call self % tpol % scheme(delta,       &
                    cmplx(c(iy,ix,i1), kind=wp),   &
                    cmplx(c(iy,ix,i2), kind=wp),   &
                    cmplx(c(iy,ix,i3), kind=wp),   &
                    cmplx(c(iy,ix,i4), kind=wp),   &
                    cmplx(ct(iy,ix,i1), kind=wp),  &
                    cmplx(ct(iy,ix,i2), kind=wp),  &
                    cmplx(ct(iy,ix,i3), kind=wp),  &
                    cmplx(ct(iy,ix,i4), kind=wp),  &
                    self % c_cur(1,iy,ix), self % ct_cur(1,iy,ix))
        self % c_cur(2,iy,ix) = czero
        self % ct_cur(2,iy,ix) = czero
        ! Wave height and d/dt of wave height
        call self % tpol % scheme(delta,       &
                    cmplx(h(iy,ix,i1), kind=wp),   &
                    cmplx(h(iy,ix,i2), kind=wp),   &
                    cmplx(h(iy,ix,i3), kind=wp),   &
                    cmplx(h(iy,ix,i4), kind=wp),   &
                    cmplx(ht(iy,ix,i1), kind=wp),  &
                    cmplx(ht(iy,ix,i2), kind=wp),  &
                    cmplx(ht(iy,ix,i3), kind=wp),  &
                    cmplx(ht(iy,ix,i4), kind=wp),  &
                    self % h_cur(1,iy,ix), self % ht_cur(1,iy,ix))
        self % h_cur(2,iy,ix) = czero
        self % ht_cur(2,iy,ix) = czero
        do iy = 1, self % nsumy
            ! Potential and d/dt of potential
            call self % tpol % scheme(delta,       &
                        cmplx(c(iy,ix,i1), kind=wp),   &
                        cmplx(c(iy,ix,i2), kind=wp),   &
                        cmplx(c(iy,ix,i3), kind=wp),   &
                        cmplx(c(iy,ix,i4), kind=wp),   &
                        cmplx(ct(iy,ix,i1), kind=wp),  &
                        cmplx(ct(iy,ix,i2), kind=wp),  &
                        cmplx(ct(iy,ix,i3), kind=wp),  &
                        cmplx(ct(iy,ix,i4), kind=wp),  &
                        self % c_cur(1,iy,ix), self % ct_cur(1,iy,ix))
            ! Wave height and d/dt of wave height
            call self % tpol % scheme(delta,       &
                        cmplx(h(iy,ix,i1), kind=wp),   &
                        cmplx(h(iy,ix,i2), kind=wp),   &
                        cmplx(h(iy,ix,i3), kind=wp),   &
                        cmplx(h(iy,ix,i4), kind=wp),   &
                        cmplx(ht(iy,ix,i1), kind=wp),  &
                        cmplx(ht(iy,ix,i2), kind=wp),  &
                        cmplx(ht(iy,ix,i3), kind=wp),  &
                        cmplx(ht(iy,ix,i4), kind=wp),  &
                        self % h_cur(1,iy,ix), self % ht_cur(1,iy,ix))
        end do
    end do
end associate

if (.not. self % dc_bias) then
    self % c_cur(:, 0, 0) = czero
    self % ct_cur(:, 0, 0) = czero
    self % h_cur(:, 0, 0) = czero
    self % ht_cur(:, 0, 0) = czero
end if

return
!
98 continue
err_msg(1) = 'End of file when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
return
!
99 continue
err_msg(1) = 'Error when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
!
end subroutine update_time

!==============================================================================

function SfunTaylor(z, kjxjy, order) result(res) ! Value of Sfun based on Taylor expansion
real(wp), intent(in) :: z   ! z-position (>0)
real(wp), intent(in) :: kjxjy ! Actual k
integer,  intent(in) :: order ! expansion order
real(wp)             :: res
!
integer :: p
real(wp) :: ap1, apj
!
ap1 = kjxjy * z
apj = 1.0_wp
res = 1.0_wp
do p = 1, order - 1
    apj = apj * ap1 / p
    res = res + apj
end do
!
end function SfunTaylor

!==============================================================================

function phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self  ! Actual class
real(knd), intent(in) :: x,y,z ! Position application program
real(knd)             :: res   ! Potential at (x,y,z)
!
integer :: jx, jy
real(wp) :: xswd, yswd, Zfun, kjxjy, zwp, reswp, kjx2
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Cfun1
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
zwp = z
reswp = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    kjx2 = (jx * self % dkx)**2
    do jy = 0, self % nsumy
        Yfun = kappa2 * Yfun
        Cfun1 = self % c_cur(1,jy,jx) * Yfun + &
                self % c_cur(2,jy,jx) * conjg(Yfun)
        kjxjy = sqrt(kjx2 + (jy * self % dky)**2)
        if (zwp > 0.0_wp .and. self % norder > 0) then
            Zfun = SfunTaylor(zwp, kjxjy, self % norder)
        else
            Zfun = exp(zwp * kjxjy)
        endif
        reswp = reswp + real(Cfun1 * Xfun) * Zfun
    end do
end do
res = reswp
!
end function phi

!==============================================================================

function stream(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self  ! Actual class
real(knd), intent(in) :: x,y,z ! Position application program
real(knd)             :: res   ! Stream function at (x,y,z)
!
res = 0.0_knd   ! Stream functions make no sense in short crested seas
!
end function stream

!==============================================================================

function phi_t(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self  ! Actual class
real(knd),       intent(in) :: x,y,z ! Position application program
real(knd)                   :: res   ! Euler time derivative of potential at (x,y,z)
!
integer :: jx, jy
real(wp) :: xswd, yswd, Zfun, kjxjy, zwp, reswp, kjx2
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Ctfun1
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
zwp = z
reswp = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    kjx2 = (jx * self % dkx)**2
    do jy = 0, self % nsumy
        Yfun = kappa2 * Yfun
        Ctfun1 = self % ct_cur(1,jy,jx) * Yfun +   &
                 self % ct_cur(2,jy,jx) * conjg(Yfun)
        kjxjy = sqrt(kjx2 +(jy * self % dky)**2)
        if (zwp > 0.0_wp .and. self % norder > 0) then
            Zfun = SfunTaylor(zwp, kjxjy, self % norder)
        else
            Zfun = exp(zwp * kjxjy)
        endif
        reswp = reswp + real(Ctfun1 * Xfun) * Zfun
    end do
end do
res = reswp
!
end function phi_t

!==============================================================================

function grad_phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res(3) ! Particle velocity at (x,y,z)
!
integer :: jx, jy
real(wp) :: xswd, yswd, Zfun, kjxjy
real(wp) :: phi_xswd, phi_yswd, phi_zswd, kjx, kjy, kjx2, zwp
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Cfun1, Cfun2, XZfun, CfunA, CfunB
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
zwp = z
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
!
phi_xswd = 0.0_wp
phi_yswd = 0.0_wp
phi_zswd = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    kjx = jx * self % dkx
    kjx2 = kjx * kjx
    do jy = 0, self % nsumy
        kjy = jy * self % dky
        kjxjy = sqrt(kjx2 + kjy*kjy)
        if (zwp > 0.0_wp .and. self % norder > 0) then
            Zfun = SfunTaylor(zwp, kjxjy, self % norder)
        else
            Zfun = exp(zwp * kjxjy)
        endif
        Yfun = kappa2 * Yfun
        CfunA = self % c_cur(1,jy,jx) * Yfun
        CfunB = self % c_cur(2,jy,jx) * conjg(Yfun)
        XZfun = Xfun * Zfun
        Cfun1 = (CfunA + CfunB) * XZfun
        Cfun2 = (CfunA - CfunB) * XZfun  ! d(conjg(Yfun))/dy = - d(Yfun)/dy
        phi_xswd = phi_xswd + kjx * Cfun1 % im
        phi_yswd = phi_yswd + kjy * Cfun2 % im
        phi_zswd = phi_zswd + kjxjy * Cfun1 % re
    end do
end do
res(1) = phi_xswd * self % cbeta - phi_yswd * self % sbeta
res(2) = phi_xswd * self % sbeta + phi_yswd * self % cbeta
res(3) = phi_zswd
!
end function grad_phi

!==============================================================================

function grad_phi_2nd(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd),       intent(in) :: x,y,z  ! Position application program
real(knd)                   :: res(6) ! Second order gradients of potential at (x,y,z)
                                      ! res(1) = d^2(potential) / dx^2
                                      ! res(2) = d^2(potential) / dx dy
                                      ! res(3) = d^2(potential) / dx dz
                                      ! res(4) = d^2(potential) / dy dy
                                      ! res(5) = d^2(potential) / dy dz
                                      ! res(6) = d^2(potential) / dz dz
!
integer :: jx, jy
real(wp) :: xswd, yswd, Zfun, Zfun_z, kjxjy
real(wp) :: phi_xx_swd, phi_xy_swd, phi_xz_swd, phi_yy_swd, phi_yz_swd, phi_zz_swd
real(wp) :: kjx, kjy, kjx2, zwp, cc, cs, ss
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Cfun1, Cfun2, CfunA, CfunB
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
zwp = z
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
!
phi_xx_swd = 0.0_wp
phi_xy_swd = 0.0_wp
phi_xz_swd = 0.0_wp
phi_yy_swd = 0.0_wp
phi_yz_swd = 0.0_wp
phi_zz_swd = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    kjx = jx * self % dkx
    kjx2 = kjx * kjx
    do jy = 0, self % nsumy
        kjy = jy * self % dky
        kjxjy = sqrt(kjx2 + kjy*kjy)
        if (zwp > 0.0_wp .and. self % norder > 0) then
            Zfun = SfunTaylor(zwp, kjxjy, self % norder)
        else
            Zfun = exp(zwp * kjxjy)
        endif
        Zfun_z = kjxjy * Zfun
        Yfun = kappa2 * Yfun
        CfunA = self % c_cur(1,jy,jx) * Yfun
        CfunB = self % c_cur(2,jy,jx) * conjg(Yfun)
        Cfun1 = (CfunA + CfunB) * Xfun
        Cfun2 = (CfunA - CfunB) * Xfun  ! d(conjg(Yfun))/dy = - d(Yfun)/dy
        phi_xx_swd = phi_xx_swd - kjx*kjx * Cfun1 % re * Zfun
        phi_xy_swd = phi_xy_swd - kjx*kjy * Cfun2 % re * Zfun
        phi_xz_swd = phi_xz_swd + kjx * Cfun1 % im * Zfun_z
        phi_yy_swd = phi_yy_swd - kjy*kjy * Cfun1 % re * Zfun
        phi_yz_swd = phi_yz_swd + kjy * Cfun2 % im * Zfun_z
        phi_zz_swd = phi_zz_swd + kjxjy*kjxjy * Cfun1 % re * Zfun
    end do
end do
cc = self % cbeta * self % cbeta
cs = self % cbeta * self % sbeta
ss = self % sbeta * self % sbeta
res(1) = phi_xx_swd * cc - phi_xy_swd * cs * 2.0_wp + phi_yy_swd * ss
res(2) = phi_xy_swd * (cc - ss) + (phi_xx_swd - phi_yy_swd) * cs
res(3) = phi_xz_swd * self % cbeta - phi_yz_swd * self % sbeta
res(4) = phi_yy_swd * cc + phi_xy_swd * cs * 2.0_wp + phi_xx_swd * ss
res(5) = phi_yz_swd * self % cbeta + phi_xz_swd * self % sbeta
res(6) = phi_zz_swd
!
end function grad_phi_2nd

!==============================================================================

function acc_euler(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res(3) ! Euler acceleration at (x,y,z)
!
integer :: jx, jy
real(wp) :: xswd, yswd, Zfun, kjxjy
real(wp) :: phi_xtswd, phi_ytswd, phi_ztswd, kjx, kjy, kjx2, zwp
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Cfun1, Cfun2, XZfun, CfunA, CfunB
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
zwp = z
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
!
phi_xtswd = 0.0_wp
phi_ytswd = 0.0_wp
phi_ztswd = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    kjx = jx * self % dkx
    kjx2 = kjx * kjx
    do jy = 0, self % nsumy
        kjy = jy * self % dky
        kjxjy = sqrt(kjx2 + kjy*kjy)
        if (zwp > 0.0_wp .and. self % norder > 0) then
            Zfun = SfunTaylor(zwp, kjxjy, self % norder)
        else
            Zfun = exp(zwp * kjxjy)
        endif
        Yfun = kappa2 * Yfun
        CfunA = self % ct_cur(1,jy,jx) * Yfun
        CfunB = self % ct_cur(2,jy,jx) * conjg(Yfun)
        XZfun = Xfun * Zfun
        Cfun1 = (CfunA + CfunB) * XZfun
        Cfun2 = (CfunA - CfunB) * XZfun  ! d(conjg(Yfun))/dy = - d(Yfun)/dy
        phi_xtswd = phi_xtswd + kjx * Cfun1 % im 
        phi_ytswd = phi_ytswd + kjy * Cfun2 % im
        phi_ztswd = phi_ztswd + kjxjy * Cfun1 % re
    end do
end do
res(1) = phi_xtswd * self % cbeta - phi_ytswd * self % sbeta
res(2) = phi_xtswd * self % sbeta + phi_ytswd * self % cbeta
res(3) = phi_ztswd
!
end function acc_euler

!==============================================================================

function acc_particle(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res(3) ! Particle acceleration at (x,y,z)
!
real(wp) :: a_euler(3), vel(3), g2nd(6)
!
! Could be about 3 times faster by merging these 3 functions in one local loop...
a_euler = self % acc_euler(x, y, z)
vel = self % grad_phi(x, y, z)
g2nd = self % grad_phi_2nd(x, y, z)
!
res(1) = a_euler(1) + vel(1) * g2nd(1) + vel(2) * g2nd(2) + vel(3) * g2nd(3)
res(2) = a_euler(2) + vel(1) * g2nd(2) + vel(2) * g2nd(4) + vel(3) * g2nd(5)
res(3) = a_euler(3) + vel(1) * g2nd(3) + vel(2) * g2nd(5) + vel(3) * g2nd(6)
!
end function acc_particle

!==============================================================================

function pressure(self, x, y, z) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res    ! Fully nonlinear pressure
!
integer :: jx, jy
real(wp) :: xswd, yswd, Zfun, kjxjy
real(wp) :: phi_t, phi_xswd, phi_yswd, phi_zswd, kjx, kjy, kjx2, zwp
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Cfun1, Cfun2, XZfun, CfunA, CfunB, Ctfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
zwp = z
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
phi_xswd = 0.0_wp
phi_yswd = 0.0_wp
phi_zswd = 0.0_wp
phi_t = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    kjx = jx * self % dkx
    kjx2 = kjx * kjx
    do jy = 0, self % nsumy
        kjy = jy * self % dky
        kjxjy = sqrt(kjx2 + kjy*kjy)
        if (zwp > 0.0_wp .and. self % norder > 0) then
            Zfun = SfunTaylor(zwp, kjxjy, self % norder)
        else
            Zfun = exp(zwp * kjxjy)
        endif
        Yfun = kappa2 * Yfun
        CfunA = self % c_cur(1,jy,jx) * Yfun
        CfunB = self % c_cur(2,jy,jx) * conjg(Yfun)
        Ctfun = self % ct_cur(1,jy,jx) * Yfun +   &
                self % ct_cur(2,jy,jx) * conjg(Yfun)
        XZfun = Xfun * Zfun
        Cfun1 = (CfunA + CfunB) * XZfun
        Cfun2 = (CfunA - CfunB) * XZfun  ! d(conjg(Yfun))/dy = - d(Yfun)/dy
        phi_xswd = phi_xswd + kjx * Cfun1 % im 
        phi_yswd = phi_yswd + kjy * Cfun2 % im
        phi_zswd = phi_zswd + kjxjy * Cfun1 % re ! Zfun_z = kjxjy * Zfun
        phi_t = phi_t + real(Ctfun * Xfun) * Zfun
    end do
end do
res = (-phi_t - (phi_xswd**2 + phi_yswd**2 + phi_zswd**2) * 0.5_wp - &
       z * self % grav) * self % rho
!
end function pressure

!==============================================================================

subroutine convergence(self, x, y, z, csv)
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self   ! Actual class
real(knd),        intent(in) :: x,y,z  ! Position application program
character(len=*), intent(in) :: csv    ! New output file
!
integer :: j, jx, jy, lucsv
real(wp) :: xswd, yswd, Zfun, kjxjy, elev, phi_x, phi_y, phi_z, prs
real(wp) :: phi_t, phi_xswd, phi_yswd, phi_zswd, kjx, kjy, kjx2, zwp
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Cfun1, Cfun2, XZfun
complex(wp) :: CfunA, CfunB, Ctfun, Hfun1
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::convergence'
character(len=250) :: err_msg(3)
!
open( newunit=lucsv, file=csv, status='replace', iostat=j)
if ( j /= 0 ) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'Not able create csv file:'
    err_msg(3) = csv
    call self % error % set_id_msg(err_proc, 1001, err_msg(1:3))
    return
end if
write(lucsv,'(a)') 'For each jx: accumulated values including all jy is presented'
write(lucsv,'(a,f0.9)') 'x = ', x
write(lucsv,'(a,f0.9)') 'y = ', y
write(lucsv,'(a,f0.9)') 'z = ', z
write(lucsv,'(a,f0.9)') 't = ', self % tswd - self % t0
write(lucsv,'(100(a,:","))') 'jx', 'velx', 'vely', 'velz', 'elev', 'prs'

xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
zwp = z
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
phi_xswd = 0.0_wp
phi_yswd = 0.0_wp
phi_zswd = 0.0_wp
phi_t = 0.0_wp
elev = 0.0_wp
Zfun = 1.0_wp  ! For infinite water Zfun=Sfun (Ufun=1, Vfun=0)
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    ! For each jx we output the total contribution from all y-directions...
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    kjx = jx * self % dkx
    kjx2 = kjx * kjx
    do jy = 0, self % nsumy
        kjy = jy * self % dky
        kjxjy = sqrt(kjx2 + kjy*kjy)
        if (zwp > 0.0_wp .and. self % norder > 0) then
            Zfun = SfunTaylor(zwp, kjxjy, self % norder)
        else
            Zfun = exp(zwp * kjxjy)
        endif
        Yfun = kappa2 * Yfun
        CfunA = self % c_cur(1,jy,jx) * Yfun
        CfunB = self % c_cur(2,jy,jx) * conjg(Yfun)
        Ctfun = self % ct_cur(1,jy,jx) * Yfun +   &
                self % ct_cur(2,jy,jx) * conjg(Yfun)
        Hfun1 = self % h_cur(1,jy,jx) * Yfun +   &
                self % h_cur(2,jy,jx) * conjg(Yfun)
        XZfun = Xfun * Zfun
        Cfun1 = (CfunA + CfunB) * XZfun
        Cfun2 = (CfunA - CfunB) * XZfun  ! d(conjg(Yfun))/dy = - d(Yfun)/dy
        phi_xswd = phi_xswd + kjx * Cfun1 % im 
        phi_yswd = phi_yswd + kjy * Cfun2 % im
        phi_zswd = phi_zswd + kjxjy * Cfun1 % re ! Zfun_z = kjxjy * Zfun
        phi_t = phi_t + real(Ctfun * Xfun) * Zfun
        elev = elev + real(Hfun1 * Xfun)
    end do
    phi_x = phi_xswd * self % cbeta - phi_yswd * self % sbeta
    phi_y = phi_xswd * self % sbeta + phi_yswd * self % cbeta
    phi_z = phi_zswd
    prs = (-phi_t - (phi_x**2 + phi_y**2 + phi_z**2) * 0.5_wp - z * self % grav) * self % rho
    write(lucsv, '(100(f0.13,:","))') real(j,wp), phi_x, phi_y, phi_z, elev, prs
end do
close(lucsv)
!
end subroutine convergence

!==============================================================================

subroutine strip(self, tmin, tmax, file_swd)
! Store part of SWD file into a new SWD file
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self   ! Actual class
real(knd),        intent(in) :: tmin, tmax ! Time window to store
character(len=*), intent(in) :: file_swd   ! Name of new swd file
!
integer(int64) :: ipos, ipos_save
integer :: i, ios, lures, istep_first, istep_last, n_camp_per_time_step
real(wp) :: tmin_swd, tmax_swd, tmax_allow
integer(c_int) fmt, shp, amp, nid, nstrip, nsteps, order, nx, ny
real(c_float) dkx, dky, grav, lscale, dt, magic
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
complex(c_float), allocatable :: camp(:)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::strip'
character(len=250) :: err_msg(5)
!
if (tmin >= tmax) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'tmin should be less than tmax'
    write(err_msg(3),'(a,f0.5)') 'tmin = ', tmin
    write(err_msg(4),'(a,f0.5)') 'tmax = ', tmax
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4))
    return
end if
tmax_allow = self % get_real('tmax')
if (tmax >= tmax_allow) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    write(err_msg(2),'(a,f0.5)') 'tmax should be less than: ', tmax_allow
    write(err_msg(3),'(a,f0.5)') 'tmax = ', tmax
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
    return
end if
tmin_swd = self % t0 + tmin
tmax_swd = self % t0 + tmax
istep_first = floor(tmin_swd / self % dt) + 1 ! t=0.0 is time step 1
istep_last = ceiling(tmax_swd / self % dt) + 1
!
call open_swd_file(newunit=lures, file=file_swd, status='replace', &
                   as_little_endian=.true., iostat=ios)
if ( ios/=0 ) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'Not able to open new SWD file for writing.'
    err_msg(3) = file_swd
    call self % error % set_id_msg(err_proc, 1001, err_msg(1:3))
    return
end if
!
! We first store the current file pointer to be restored when finished
inquire(self % unit, pos=ipos_save)
read(self % unit, pos=1, end=98, err=99) magic
write(lures) magic

read(self % unit, end=98, err=99) fmt
write(lures) fmt

read(self % unit, end=98, err=99) shp
write(lures) shp

read(self % unit, end=98, err=99) amp
write(lures) amp

read(self % unit, end=98, err=99) cprog
write(lures) cprog

read(self % unit, end=98, err=99) cdate
write(lures) cdate

read(self % unit, end=98, err=99) nid
write(lures) nid

allocate(character(len=int(nid)) :: cid)
read(self % unit, end=98, err=99) cid
write(lures) cid

read(self % unit, end=98, err=99) grav
write(lures) grav

read(self % unit, end=98, err=99) lscale
write(lures) lscale

read(self % unit, end=98, err=99) nstrip
nstrip = nstrip + istep_first - 1
write(lures) nstrip

read(self % unit, end=98, err=99) nsteps
nsteps = istep_last - istep_first + 1
write(lures) nsteps

read(self % unit, end=98, err=99) dt
write(lures) dt

read(self % unit, end=98, err=99) order
write(lures) order

read(self % unit, end=98, err=99) nx
write(lures) nx

read(self % unit, end=98, err=99) ny
write(lures) ny

read(self % unit, end=98, err=99) dkx
write(lures) dkx

read(self % unit, end=98, err=99) dky
write(lures) dky

! Only complex numbers in remaining part of swd file...
n_camp_per_time_step = self % size_step / self % size_complex

allocate(camp(n_camp_per_time_step), stat=i)
if (i /= 0) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'Not able to allocate workspace for camp(:).'
    write(err_msg(3),'(a,i0)') 'Requested size is ', n_camp_per_time_step
    call self % error % set_id_msg(err_proc, 1001, err_msg(1:3))
    return
end if

! We don't need the first 'istep_first-1' time steps in the SWD file.
! We will reposition the file pointer
inquire(self % unit, pos=ipos)
! requested new file position
ipos = ipos + int(istep_first - 1, int64) * self % size_step
! Reposition the file pointer by reading a dummy complex number 
! just before the goodies...
read(self % unit, pos=ipos - self % size_complex, end=98, err=99) camp(1)
do i = 1, nsteps
    read(self % unit, end=98, err=99) camp(:)
    write(lures) camp(:)
end do
close(lures)
! Reposition and check of file pointer
read(self % unit, pos=ipos_save - self % size_complex, end=98, err=99) camp(1)
inquire(self % unit, pos=ipos)
if (ipos /= ipos_save) then
    err_msg(1) = 'Unexpected position of file pointer.'
    write(err_msg(2),'(a,i0)') 'ipos = ', ipos
    write(err_msg(3),'(a,i0)') 'ipos_save = ', ipos_save
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:3))
end if
return
!
98 continue
err_msg(1) = 'End of file when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
return

99 continue
err_msg(1) = 'Error when reading data from file:'
err_msg(2) = self % file
call self % error % set_id_msg(err_proc, 1003, err_msg(1:2))
!
end subroutine strip

!==============================================================================

function elev(self, x, y) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self ! Actual class
real(knd), intent(in) :: x,y  ! Position application program
real(knd)             :: res  ! Surface elevation at (x,y)
!
integer :: jx, jy
real(wp) :: xswd, yswd, reswp
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Hfun1
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
reswp = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    do jy = 0, self % nsumy
        Yfun = kappa2 * Yfun
        Hfun1 = self % h_cur(1,jy,jx) * Yfun +   &
                self % h_cur(2,jy,jx) * conjg(Yfun)
        reswp = reswp + real(Hfun1 * Xfun)
    end do
end do
res = reswp
!
end function elev

!==============================================================================

function elev_t(self, x, y) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self ! Actual class
real(knd), intent(in) :: x,y  ! Position application program
real(knd)             :: res  ! d/dt of surface elevation at (x,y)
!
integer :: jx, jy
real(wp) :: xswd, yswd, reswp
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Htfun1
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))

reswp = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    do jy = 0, self % nsumy
        Yfun = kappa2 * Yfun
        Htfun1 = self % ht_cur(1,jy,jx) * Yfun + &
                 self % ht_cur(2,jy,jx) * conjg(Yfun)
        reswp = reswp + real(Htfun1 * Xfun)
    end do
end do
res = reswp
!
end function elev_t

!==============================================================================

function grad_elev(self, x, y) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y    ! Position application program
real(knd)             :: res(3) ! x, y and z gradients of surface elevation at (x,y)
!
integer :: jx, jy
real(wp) :: xswd, yswd, elev_xswd, elev_yswd, kjx, kjy
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Hfun1, Hfun2, HfunA, HfunB
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))

elev_xswd = 0.0_wp
elev_yswd = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    kjx = jx * self % dkx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    do jy = 0, self % nsumy
        kjy = jy * self % dky
        Yfun = kappa2 * Yfun
        HfunA = self % h_cur(1,jy,jx) * Yfun
        HfunB = self % h_cur(2,jy,jx) * conjg(Yfun)
        Hfun1 = (HfunA + HfunB) * Xfun
        Hfun2 = (HfunA - HfunB) * Xfun  ! d(conjg(Yfun))/dy = - d(Yfun)/dy
        elev_xswd = elev_xswd + kjx * Hfun1 % im
        elev_yswd = elev_yswd + kjy * Hfun2 % im
    end do
end do
res(1) = elev_xswd * self % cbeta - elev_yswd * self % sbeta
res(2) = elev_xswd * self % sbeta + elev_yswd * self % cbeta
res(3) = 0.0_knd
!
end function grad_elev

!==============================================================================

function grad_elev_2nd(self, x, y) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y    ! Position application program
real(knd)             :: res(3) ! Second order gradients of surface elevation
                                ! res(1) = d^2(elevation) / dx^2
                                ! res(2) = d^2(elevation) / dx dy
                                ! res(3) = d^2(elevation) / dy dy
!
integer :: jx, jy
real(wp) :: xswd, yswd, elev_xx_swd, elev_xy_swd, elev_yy_swd, kjx, kjy
real(wp) :: cc, cs, ss
complex(wp) :: kappa1, kappa2, Xfun, Yfun, Hfun1, Hfun2, HfunA, HfunB
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
kappa1 = exp(cmplx(0.0_wp, -self % dkx * xswd, kind=wp))
kappa2 = exp(cmplx(0.0_wp, -self % dky * yswd, kind=wp))
elev_xx_swd = 0.0_wp
elev_xy_swd = 0.0_wp
elev_yy_swd = 0.0_wp
Xfun = conjg(kappa1)
do jx = 0, self % nsumx
    kjx = jx * self % dkx
    Xfun = kappa1 * Xfun
    Yfun = conjg(kappa2)
    do jy = 0, self % nsumy
        kjy = jy * self % dky
        Yfun = kappa2 * Yfun
        HfunA = self % h_cur(1,jy,jx) * Yfun
        HfunB = self % h_cur(2,jy,jx) * conjg(Yfun)
        Hfun1 = (HfunA + HfunB) * Xfun
        Hfun2 = (HfunA - HfunB) * Xfun  ! d(conjg(Yfun))/dy = - d(Yfun)/dy
        elev_xx_swd = elev_xx_swd - kjx*kjx * Hfun1 % re
        elev_xy_swd = elev_xy_swd - kjx*kjy * Hfun2 % re
        elev_yy_swd = elev_yy_swd - kjy*kjy * Hfun1 % re
    end do
end do
cc = self % cbeta * self % cbeta
cs = self % cbeta * self % sbeta
ss = self % sbeta * self % sbeta
res(1) = elev_xx_swd * cc - elev_xy_swd * cs * 2.0_wp + elev_yy_swd * ss
res(2) = elev_xy_swd * (cc - ss) + (elev_xx_swd - elev_yy_swd) * cs
res(3) = elev_yy_swd * cc + elev_xy_swd * cs * 2.0_wp + elev_xx_swd * ss
!
end function grad_elev_2nd

!==============================================================================

function bathymetry(self, x, y) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self ! Actual class
real(knd), intent(in) :: x,y  ! Position application program
real(knd)             :: res  ! Local water depth
!
res = -1.0_knd
!
end function bathymetry

!==============================================================================

function bathymetry_nvec(self, x, y) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y    ! Position application program
real(knd)             :: res(3) ! Unit normal vector into ocean at (x,y)
!
res(1) = 0.0_knd
res(2) = 0.0_knd
res(3) = 1.0_knd
!
end function bathymetry_nvec

!==============================================================================

function get_int(self, name) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of int parameter
integer                      :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::get_int'
character(len=250) :: err_msg(1)
!
select case(name)
case('fmt')
    res = self % fmt
case('shp')
    res = self % shp
case('amp')
    res = self % amp
case('nid')
    res = len(self % cid)
case('nstrip')
    res = self % nstrip
case('nsteps')
    res = self % nsteps
case('order')
    res = self % order
case('norder')
    res = self % norder
case('nx')
    res = self % nx
case('ny')
    res = self % ny
case('ipol')
    res = self % ipol
case('impl')
    res = 1
case('nsumx')
    res = self % nsumx
case('nsumy')
    res = self % nsumy
case default
    res = huge(res)
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_int

!==============================================================================

function get_logical(self, name) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of logical parameter
logical                      :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::get_logical'
character(len=250) :: err_msg(1)
!
select case(name)
case('dc_bias')
    res = self % dc_bias
case default
    res = .false.
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_logical

!==============================================================================

function get_real(self, name) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of float parameter
real(knd)                    :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::get_real'
character(len=250) :: err_msg(1)
!
select case(name)
! Constructor parameters...
case('t0')
    res = self % t0
case('x0')
    res = self % x0
case('y0')
    res = self % y0
case('beta')
    res = atan2(self % sbeta, self % cbeta) * 180.0_wp / pi
    if (res < 0.0_wp) res = res + 360.0_wp  ! 0.0 <= beta < 360.0
case('rho')
    res = self % rho
! All actual float scalars from the swd file...
case('magic')
    res = swd_magic_number
case('grav')
    res = self % grav
case('lscale')
    res = self % lscale
case('dt')
    res = self % dt
case('dkx')
    res = self % dkx
case('dky')
    res = self % dky
! Some special requests....
case('tmax')
    res = self % tmax
case('d')
    res = -1.0_knd
case('lmin')
    res = 2 * pi / sqrt((self % dkx * self % nsumx)**2 + (self % dky * self % nsumy)**2)
case('lmax')
    res = 2 * pi / min(self % dkx, self % dky)
case('sizex')
    res = 2 * pi / self % dkx
case('sizey')
    res = 2 * pi / self % dky
case default
    res = huge(res)
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_real

!==============================================================================

function get_chr(self, name) result(res)
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self   ! Actual class
character(len=*),  intent(in) :: name ! Name of char parameter
character(len=:), allocatable :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::get_chr'
character(len=250) :: err_msg(1)
!
select case(name)
case('file', 'file_swd')
    res = self % file
case('version')
    res = version
case('class')
    res = 'spectral_wave_data_shape_4_impl_1'
case('prog')
    res = self % prog
case('date')
    res = self % date
case('cid')
    res = self % cid
case default
    res = 'unknown name specified'
    write(err_msg(1),'(a,i0,a,a,a)') 'For this shape class (shp=', &
            self % shp, ') the key "', trim(name), '" is unknown.'
    call self % error % set_id_msg(err_proc, 1004, err_msg)
end select
!
end function get_chr

!==============================================================================

function elev_fft(self, nx_fft_in, ny_fft_in) result(elev)
class(spectral_wave_data_shape_4_impl_1), intent(inout) :: self ! Actual class
integer, optional, intent(in) :: nx_fft_in, ny_fft_in
real(knd), allocatable :: elev(:, :)
complex(wp) :: c_fft(self % nsumx + 1, 2*self % nsumy + 1)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_4_impl_1::elev_fft'
character(len=:), allocatable :: err_msg(:)

c_fft = self % fft % swd_to_fft_coeffs_2D(self % h_cur(:, 0:self % nsumy, 0:self % nsumx))
elev = self % fft % fft_field_2D(c_fft, nx_fft_in, ny_fft_in)

if (self % fft % error % raised()) then
    err_msg = [self % fft % error % get_msg()]
    call self % error % set_id_msg(err_proc, &
                                   self % fft % error % get_id(), &
                                   err_msg)
end if

end function elev_fft

!==============================================================================

end module spectral_wave_data_shape_4_impl_1_def
