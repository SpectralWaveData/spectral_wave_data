module spectral_wave_data_shape_3_impl_1_def

use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float

use kind_values, only: knd => kind_swd_interface, wp => kind_swd_internal

use open_swd_file_def, only: open_swd_file, swd_validate_binary_convention, &
                             swd_magic_number
use spectral_wave_data_def, only: spectral_wave_data
use spectral_interpolation_def, only: spectral_interpolation
use swd_version, only: version

implicit none
private

! This module provides an implemention of the shape 3 class of the 
! spectral-wave-data API.
!
! Written by Jens Bloch Helmers, August, 2. 2019
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: spectral_wave_data_shape_3_impl_1
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type, extends(spectral_wave_data) :: spectral_wave_data_shape_3_impl_1
    integer            :: n    ! Number of spectral components exp(kz) in x-direction
    integer            :: nh   ! Number of spectral components exp(-kz) in x-direction
    integer            :: nsumx ! Number of spectral components exp(kz) in x-summation (<=n)
    integer            :: nsumxh ! Number of spectral components exp(-kz) in x-summation (<=nh)
    integer            :: isf  ! Flag to indicate geometric description of the sea floor
    integer            :: nsf  ! Number of points defining the sea floor
    real(wp)           :: d    ! Average water depth (<0 if infinite depth)
    real(wp)           :: dk   ! Constant spacing of wave numbers in x-direction
    integer            :: icur     ! 'pointer' to column of most recent spectral input. (1:4)
    integer            :: istp     ! Most recent step from swd file in memory
    integer            :: ipt(4,4) ! 'pointer' to which spectral columns correspond to 
                                   ! indices i-1, i, i+1 and i+2 in the interpolation scheme.
    complex(c_float), allocatable :: c_win(:,:)  ! Input window of c-spectral components 1dim=1:n, 2dim=1:4
    complex(c_float), allocatable :: ct_win(:,:) ! Input window of ct-spectral components 1dim=1:n, 2dim=1:4 (if fmt=2 or 3)
    complex(c_float), allocatable :: ch_win(:,:)  ! Input window of ch-spectral components 1dim=1:n, 2dim=1:4
    complex(c_float), allocatable :: cht_win(:,:) ! Input window of cht-spectral components 1dim=1:n, 2dim=1:4 (if fmt=2 or 3)
    complex(c_float), allocatable :: h_win(:,:)  ! Input window of h-spectral components 1dim=1:n, 2dim=1:4
    complex(c_float), allocatable :: ht_win(:,:) ! Input window of ht-spectral components 1dim=1:n, 2dim=1:4 (if fmt=3)
    complex(wp), allocatable :: c_cur(:)    ! Spectral c-values at current time (1dim=1:n)
    complex(wp), allocatable :: ct_cur(:)   ! Spectral ct-values at current time (1dim=1:n)
    complex(wp), allocatable :: h_cur(:)    ! Spectral ct-values at current time (1dim=1:n)
    complex(wp), allocatable :: ht_cur(:)   ! Spectral ht-values at current time (1dim=1:n)
    complex(wp), allocatable :: ch_cur(:)   ! Spectral ch-values at current time (1dim=1:nh)
    complex(wp), allocatable :: cht_cur(:)  ! Spectral cht-values at current time (1dim=1:nh)
    real(wp), allocatable :: xsf(:)   ! x-locations defining the sea floor
    real(wp), allocatable :: zsf(:)   ! z-locations defining the sea floor
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
end type spectral_wave_data_shape_3_impl_1

interface spectral_wave_data_shape_3_impl_1
    module procedure constructor
end interface

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

subroutine close(self)
class(spectral_wave_data_shape_3_impl_1) :: self  ! Object to destruct
!
logical opened
!
inquire(unit=self % unit, opened=opened)
if (opened) close(self % unit)
if (allocated(self % cid)) deallocate(self % cid)
if (allocated(self % c_win)) deallocate(self % c_win)
if (allocated(self % ct_win)) deallocate(self % ct_win)
if (allocated(self % ch_win)) deallocate(self % ch_win)
if (allocated(self % cht_win)) deallocate(self % cht_win)
if (allocated(self % h_win)) deallocate(self % h_win)
if (allocated(self % ht_win)) deallocate(self % ht_win)
if (allocated(self % h_cur)) deallocate(self % h_cur)
if (allocated(self % ht_cur)) deallocate(self % ht_cur)
if (allocated(self % c_cur)) deallocate(self % c_cur)
if (allocated(self % ct_cur)) deallocate(self % ct_cur)
if (allocated(self % ch_cur)) deallocate(self % ch_cur)
if (allocated(self % cht_cur)) deallocate(self % cht_cur)
!
self % file = '0'
self % unit = 0
self % n = 0
self % nsumx = 0
!
end subroutine  close
    
!==============================================================================

function constructor(file, x0, y0, t0, beta, rho, nsumx, ipol, norder, &
                     dc_bias) result(self)
character(len=*),    intent(in)  :: file  ! File containing HOSM data
real(knd),           intent(in)  :: x0,y0 ! Spatial seed
real(knd),           intent(in)  :: t0    ! Temporal seed >=0
real(knd),           intent(in)  :: beta  ! Wave heading (deg)
real(knd), optional, intent(in)  :: rho   ! Density of water (applied for pressure calculations)
integer, optional,   intent(in)  :: nsumx ! If present and nsumx>-1: apply nsumx number of spectral
                                          ! components in x-directions, else apply all spectral
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
type(spectral_wave_data_shape_3_impl_1) :: self  ! Object to construct
!
integer :: i, ios

integer(int64) :: ipos1, ipos2
integer(c_int) :: fmt, shp, amp, n, nh, order, nid, nsteps, nstrip, isf, nsf
real(c_float) :: dk, dt, grav, lscale, xsf, zsf, magic
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
real(wp) :: sum
integer :: n_cmplx_per_step
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::constructor'
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

call swd_validate_binary_convention(self % file, err_msg(2))
if (err_msg(2) /= '') then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    call self % error % set_id_msg(err_proc, 1002, err_msg(1:2))
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

read(self % unit, end=98, err=99) n
self % n = n

read(self % unit, end=98, err=99) nh
self % nh = nh

read(self % unit, end=98, err=99) dk
self % dk = dk

read(self % unit, end=98, err=99) isf
self % isf = isf
if (self % isf /= 0) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a)') 'In the current implementation only isf=0 is supported.'
    write(err_msg(3),'(a,i0)') 'isf = ', self % isf
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:3))
    return
end if

read(self % unit, end=98, err=99) nsf
self % nsf = nsf

if (nsf > 0) then
    allocate(self % xsf(nsf), self % zsf(nsf))
    do i = 1, self % nsf
        read(self % unit, end=98, err=99) xsf
        self % xsf(i) = xsf
    end do
    do i = 1, self % nsf
        read(self % unit, end=98, err=99) zsf
        self % zsf(i) = zsf
    end do
end if
if (nsf == 0) then
    self % d = -1.0_wp
     self % nsumxh = -1
else if (nsf == 1) then
    if (self % zsf(1) >= 0.0_wp) then
        write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
        write(err_msg(2),'(a)') 'zsf(1) should be negative'
        write(err_msg(3),'(a,f0.8)') 'zsf(1) = ', self % zsf(1)
        call self % error % set_id_msg(err_proc, 1003, err_msg(1:3))
        return
    end if
    self % d = - self % zsf(1)
    self % nsumxh = int(atanh(1.0_wp - 100*epsilon(self % d)) / (self % dk * self % d))
    self % nsumxh = min(self % nh, self % nsumxh)
else if (nsf > 1) then
     self % nsumxh = self % nh
    ! Check if xsf(:) and zsf(:) are sound
    if (abs(self % xsf(1)) > epsilon(self % xsf)) then
        write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
        write(err_msg(2),'(a)') 'xsf(1) should be 0.0'
        write(err_msg(3),'(a,f0.8)') 'xsf(1) = ', self % xsf(1)
        call self % error % set_id_msg(err_proc, 1003, err_msg(1:3))
        return
    end if
    if (abs(self % xsf(nsf) - 2*pi/dk) > 1.0e-4_wp * 2*pi/dk) then
        write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
        write(err_msg(2),'(a)') 'xsf(nsf) should be 2*pi/dk'
        write(err_msg(3),'(a,f0.8)') 'xsf(nsf) = ', self % xsf(nsf)
        call self % error % set_id_msg(err_proc, 1003, err_msg(1:3))
        return
    end if
    do i = 2, self % nsf
        if (self % xsf(i-1) >= self % xsf(i)) then
            write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
            write(err_msg(2),'(a)') 'xsf(1:nsf) should be monotonic increasing'
            write(err_msg(3),'(a,f0.8)') 'xsf(i-1) = ', self % xsf(i-1)
            write(err_msg(4),'(a,f0.8)') 'xsf(i) = ', self % xsf(i)
            write(err_msg(5),'(a,i0)') 'i = ', i
            call self % error % set_id_msg(err_proc, 1003, err_msg(1:5))
            return
        end if
    end do
    do i = 1, self % nsf
        if (self % zsf(i) >= 0.0_wp) then
            write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
            write(err_msg(2),'(a)') 'zsf(1:nsf) should all be negative'
            write(err_msg(3),'(a,f0.8)') 'zsf(i) = ', self % zsf(i)
            write(err_msg(4),'(a,i0)') 'i = ', i
            call self % error % set_id_msg(err_proc, 1003, err_msg(1:4))
            return
        end if
    end do
    ! Calculate average depth
    sum = 0.0_wp
    do i = 2, self % nsf
        sum = sum + 0.5_wp * (self % zsf(i-1) + self % zsf(i)) * &
                             (self % xsf(i-1) + self % xsf(i))
    end do
    self % d = - sum / self % xsf(nsf)
else
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a)') 'nsf should not be negative'
    write(err_msg(3),'(a,i0)') 'nsf = ', self % nsf
    call self % error % set_id_msg(err_proc, 1003, err_msg(1:3))
    return
end if
    
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
        self % nsumx = n
    else
        self % nsumx = nsumx
        if (nsumx > n) then
            write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
            write(err_msg(2),'(a,i0)') 'Number of spectral components in SWD file = ', n
            write(err_msg(3),'(a,i0)') 'Number of requested spectral components   = ', nsumx
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
            return
        end if
    end if
else
    self % nsumx = n
end if

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

allocate( self % c_win(0:self % n, 4),   &
          self % ct_win(0:self % n, 4),  &
          self % ch_win(0:self % nh, 4), &
          self % cht_win(0:self % nh, 4),&
          self % h_win(0:self % n, 4),   &
          self % ht_win(0:self % n, 4),  &
          self % c_cur(0:self % n),      &
          self % ct_cur(0:self % n),     &
          self % ch_cur(0:self % nh),    &
          self % cht_cur(0:self % nh),   &
          self % h_cur(0:self % n),      &
          self % ht_cur(0:self % n),    stat=i)
if (i /= 0) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a)') 'Not able to allocate space for storing window of spectral components.'
    write(err_msg(3),'(a,i0)') 'Number of spectral components = ', self % n
    call self % error % set_id_msg(err_proc, 1005, err_msg(1:3))
    return
end if

! The first timestep is put into memory.
associate(c => self % c_win, ct => self % ct_win, ch => self % ch_win, &
          cht => self % cht_win, h => self % h_win, ht => self % ht_win)
    ! request file position where the temporal functions start
    inquire(self % unit, pos=self % ipos0)

    ! set to zero intially
    h = czero_c
    ht = czero_c
    c = czero_c
    ct = czero_c
    ch = czero_c
    cht = czero_c

    n_cmplx_per_step = 2 * (self%n + 1)
    read(self % unit, end=98, err=99) h(:,2)
    read(self % unit, end=98, err=99) ht(:,2)
    if (self % amp < 3) then
        n_cmplx_per_step = n_cmplx_per_step + 2 * (self%n + 1)
        read(self % unit, end=98, err=99) c(:,2)
        read(self % unit, end=98, err=99) ct(:,2)
        if (self % nsf > 1) then
            n_cmplx_per_step = n_cmplx_per_step + 2 * (self % nh + 1)
            read(self % unit, end=98, err=99) ch(:,2)
            read(self % unit, end=98, err=99) cht(:,2)
        end if
    end if
    ipos1 = self % ipos0
    inquire(self % unit, pos=ipos2)
    ! Storage fortran units per complex (c_float based)
    self % size_complex = (ipos2 - ipos1) / (n_cmplx_per_step)
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
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self  ! Update data in memory (if needed)
real(knd),          intent(in)    :: time  ! Current time in simulation program
!
integer(int64) :: ipos
integer :: istp_max, i, j, i1, i2, i3, i4, istp_min, ios, imove
real(wp) :: delta, teps, gamma_h
complex(wp) :: fval, dfval
complex(c_float) :: cdum
complex(wp), parameter :: czero = cmplx(0.0_wp, 0.0_wp, wp)
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::update_time'
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

associate(c => self % c_win, ct => self % ct_win, ch => self % ch_win,    &
          cht => self % cht_win, h => self % h_win, ht => self % ht_win,  &
          ic => self % icur, ip => self % ipt)

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

    do j  = 1, istp_max - self % istp
        if (self % istp == self % nsteps) then
            ! We apply padding for data at the very end of the file
            ic = ic + 1
            if (ic > 4) ic = 1
            i1 = ip(1,ic)
            i2 = ip(2,ic)
            i3 = ip(3,ic)
            i4 = ip(4,ic)
            do concurrent (i = 0 : self % nsumx)
                ! Potential and d/dt of potential
                call self % tpol % pad_right(          &
                            cmplx(c(i,i1), kind=wp),   &
                            cmplx(c(i,i2), kind=wp),   &
                            cmplx(c(i,i3), kind=wp),   &
                            cmplx(ct(i,i1), kind=wp),  &
                            cmplx(ct(i,i2), kind=wp),  &
                            cmplx(ct(i,i3), kind=wp),  &
                            fval, dfval)
                c(i,i4) = fval
                ct(i,i4) = dfval
                ! Wave height and d/dt of wave height
                call self % tpol % pad_right(          &
                            cmplx(h(i,i1), kind=wp),   &
                            cmplx(h(i,i2), kind=wp),   &
                            cmplx(h(i,i3), kind=wp),   &
                            cmplx(ht(i,i1), kind=wp),  &
                            cmplx(ht(i,i2), kind=wp),  &
                            cmplx(ht(i,i3), kind=wp),  &
                            fval, dfval)
                h(i,i4) = fval
                ht(i,i4) = dfval
            end do
            if (self % nsf > 1) then
                do concurrent (i = 0 : self % nsumxh)
                    call self % tpol % pad_right(           &
                                cmplx(ch(i,i1), kind=wp),   &
                                cmplx(ch(i,i2), kind=wp),   &
                                cmplx(ch(i,i3), kind=wp),   &
                                cmplx(cht(i,i1), kind=wp),  &
                                cmplx(cht(i,i2), kind=wp),  &
                                cmplx(cht(i,i3), kind=wp),  &
                                fval, dfval)
                    ch(i,i4) = fval
                    cht(i,i4) = dfval
                end do
            end if
            self % istp = istp_max
        else
            read(self % unit, end=98, err=99) h(:,ic)
            read(self % unit, end=98, err=99) ht(:,ic)
            if (self % amp < 3) then
                read(self % unit, end=98, err=99) c(:,ic)
                read(self % unit, end=98, err=99) ct(:,ic)
                if (self % nsf > 1) then
                    read(self % unit, end=98, err=99) ch(:,ic)
                    read(self % unit, end=98, err=99) cht(:,ic)
                end if
            else
                c(:,ic) = czero_c
                ct(:,ic) = czero_c
            end if
            self % istp = self % istp + 1
            ic = ic + 1
            if (ic > 4) ic = 1
        end if
    end do

    if (istp_min == 0) then
        ! Padding in first column because tswd < dt_swd. 
        do concurrent (i = 0 : self % nsumx)
            ! Potential and d/dt of potential
            call self % tpol % pad_left(          &
                        cmplx(c(i,2), kind=wp),   &
                        cmplx(c(i,3), kind=wp),   &
                        cmplx(c(i,4), kind=wp),   &
                        cmplx(ct(i,2), kind=wp),  &
                        cmplx(ct(i,3), kind=wp),  &
                        cmplx(ct(i,4), kind=wp),  &
                        fval, dfval)
            c(i,1) = fval
            ct(i,1) = dfval
            ! Wave height and d/dt of wave height
            call self % tpol % pad_left(          &
                        cmplx(h(i,2), kind=wp),   &
                        cmplx(h(i,3), kind=wp),   &
                        cmplx(h(i,4), kind=wp),   &
                        cmplx(ht(i,2), kind=wp),  &
                        cmplx(ht(i,3), kind=wp),  &
                        cmplx(ht(i,4), kind=wp),  &
                        fval, dfval)
            h(i,1) = fval
            ht(i,1) = dfval
        end do
        if (self % nsf > 1) then
            do concurrent (i = 0 : self % nsumxh)
                call self % tpol % pad_left(           &
                            cmplx(ch(i,2), kind=wp),   &
                            cmplx(ch(i,3), kind=wp),   &
                            cmplx(ch(i,4), kind=wp),   &
                            cmplx(cht(i,2), kind=wp),  &
                            cmplx(cht(i,3), kind=wp),  &
                            cmplx(cht(i,4), kind=wp),  &
                            fval, dfval)
                ch(i,1) = fval
                cht(i,1) = dfval
            end do
        end if
    end if

    i1 = ip(1,ic)
    i2 = ip(2,ic)
    i3 = ip(3,ic)
    i4 = ip(4,ic)
    do concurrent (i = 0 : self % nsumx)
        ! Potential and d/dt of potential
        call self % tpol % scheme(delta,       &
                    cmplx(c(i,i1), kind=wp),   &
                    cmplx(c(i,i2), kind=wp),   &
                    cmplx(c(i,i3), kind=wp),   &
                    cmplx(c(i,i4), kind=wp),   &
                    cmplx(ct(i,i1), kind=wp),  &
                    cmplx(ct(i,i2), kind=wp),  &
                    cmplx(ct(i,i3), kind=wp),  &
                    cmplx(ct(i,i4), kind=wp),  &
                    self % c_cur(i), self % ct_cur(i))
        ! Wave height and d/dt of wave height
        call self % tpol % scheme(delta,       &
                    cmplx(h(i,i1), kind=wp),   &
                    cmplx(h(i,i2), kind=wp),   &
                    cmplx(h(i,i3), kind=wp),   &
                    cmplx(h(i,i4), kind=wp),   &
                    cmplx(ht(i,i1), kind=wp),  &
                    cmplx(ht(i,i2), kind=wp),  &
                    cmplx(ht(i,i3), kind=wp),  &
                    cmplx(ht(i,i4), kind=wp),  &
                    self % h_cur(i), self % ht_cur(i))
    end do
    do concurrent (i = 0 : self % nsumxh)
        if (self % nsf > 1) then
            ! Potential and d/dt of potential
            call self % tpol % scheme(delta,        &
                        cmplx(ch(i,i1), kind=wp),   &
                        cmplx(ch(i,i2), kind=wp),   &
                        cmplx(ch(i,i3), kind=wp),   &
                        cmplx(ch(i,i4), kind=wp),   &
                        cmplx(cht(i,i1), kind=wp),  &
                        cmplx(cht(i,i2), kind=wp),  &
                        cmplx(cht(i,i3), kind=wp),  &
                        cmplx(cht(i,i4), kind=wp),  &
                        self % ch_cur(i), self % cht_cur(i))
        else if (self % nsf == 1) then
            gamma_h = exp(- 2 * i * self % dk * self % d)
            self % ch_cur(i) = gamma_h * self % c_cur(i) 
            self % cht_cur(i) = gamma_h * self % ct_cur(i) 
        end if
    end do          
end associate

if (.not. self % dc_bias) then
    self % c_cur(0) = czero
    self % ct_cur(0) = czero
    self % h_cur(0) = czero
    self % ht_cur(0) = czero
    if (self % nsumxh > -1) then
        self % ch_cur(0) = czero
        self % cht_cur(0) = czero
    end if
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

subroutine update_time2(self, time)
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self  ! Update data in memory (if needed)
real(knd),          intent(in)    :: time  ! Current time in simulation program
!
integer(int64) :: ipos
integer :: istp_req, i, i1, i2, i3, i4, itmp, ios, imove
real(wp) :: delta, tmp, gamma_h
complex(c_float) :: cdum
complex(wp), parameter :: czero = cmplx(0.0_wp, 0.0_wp, wp)
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::update_time'
character(len=250) :: err_msg(6)

if (time > self % tmax) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'Requested time is too large!'
    write(err_msg(3),'(a,f0.5)') 'User time = ', time
    write(err_msg(4),'(a,f0.5)') 'Max user time = ', self % tmax
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4))
    return
end if
self % tswd = self % t0 + time
if (self % tswd < 0.0) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(self % file)
    err_msg(2) = 'time corresponds to negative swd time!'
    write(err_msg(3),'(a,f0.5)') 'time user = ', time
    write(err_msg(4),'(a,f0.5)') 'time swd = ', self % tswd
    call self % error % set_id_msg(err_proc, 1004, err_msg(1:4))
    return
end if

tmp = self % tswd / self % dt
itmp = int(tmp - 0.0000001_wp * self % dt)
delta = tmp - itmp     ! delta in [0, 1)
istp_req = itmp + 3    ! The most recent required swd step to be in memory

associate(c => self % c_win, ct => self % ct_win, ch => self % ch_win,    &
          cht => self % cht_win, h => self % h_win, ht => self % ht_win,  &
          ic => self % icur, ip => self % ipt)

    imove = istp_req - self % istp
    if (imove < 0 ) self % eof = .false.
    if (imove < 0 .or. imove > 4) then
        ! imove < 0 is negative time stepping from user. 
        ! imove > 4 is fast forward and we will skip some time steps
        ! We need to fill the 4 time step buffer from scratch.
        ! We will reposition the file pointer, but we need to make sure that
        ! 3 time steps prior to istp_req will be stored in memory...
        ! request old file position:
        inquire(self % unit, pos=ipos)
        if (istp_req == 3) then
            ! Wind back to tswd < dt_swd. 
            ! We have to pad data later. Consequently we skip the first column...
            ipos = ipos + int(imove - 3, int64) * self % size_step
            self % istp = istp_req - 3
            ic = 2 ! ic is a position counter 1,2,3 or 4.
        else
            ipos = ipos + int(imove - 4, int64) * self % size_step
            self % istp = istp_req - 4
            ic = 1 ! ic is a position counter 1,2,3 or 4.
        end if            
        ! Reposition the file pointer by reading a dummy complex number 
        ! just before the goodies...
        read(self % unit, pos=ipos - self % size_complex, end=98, iostat=ios) cdum 
        if (ios/=0) then
            err_msg(1) = 'User time beyond what is available from SWD file:'
            err_msg(2) = self % file
            err_msg(3) = 'The file has less content than expected. Not able to recover!!!'
            write(err_msg(4),'(a,f0.5)') 'Requested swd-time = ', self % tswd
            write(err_msg(5),'(a,f0.5)') 'Requested user-time = ', time
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:5))
            return
        end if
    end if

    if (.not. self % eof) then
        do i  = 1, istp_req - self % istp
            read(self % unit, end=98, iostat=ios) h(:,ic)
            if (ios/=0) then
                self % eof = .true.
                ! We apply padding for data at the end of the file
                ic = ic + 1
                if (ic > 4) ic = 1
                i1 = ip(1,ic)
                i2 = ip(2,ic)
                i3 = ip(3,ic)
                i4 = ip(4,ic)
                ! NOTE: Should be updated to new padding scheme (Not important)
                c(:,i4) = 3*(c(:,i3) - c(:,i2)) + c(:,i1)
                ct(:,i4) = 3*(ct(:,i3) - ct(:,i2)) + ct(:,i1)
                if (self % nsf > 1) then
                    ch(:,i4) = 3*(ch(:,i3) - ch(:,i2)) + ch(:,i1)
                    cht(:,i4) = 3*(cht(:,i3) - cht(:,i2)) + cht(:,i1)
                end if
                h(:,i4) = 3*(h(:,i3) - h(:,i2)) + h(:,i1)
                ht(:,i4) = 3*(ht(:,i3) - ht(:,i2)) + ht(:,i1)
                exit
            end if
            read(self % unit, end=98, err=99) ht(:,ic)
            if (self % amp < 3) then
                read(self % unit, end=98, err=99) c(:,ic)
                read(self % unit, end=98, err=99) ct(:,ic)
                if (self % nsf > 1) then
                    read(self % unit, end=98, err=99) ch(:,ic)
                    read(self % unit, end=98, err=99) cht(:,ic)
                end if
            else
                c(:,ic) = czero_c
                ct(:,ic) = czero_c
                if (self % nsf > 1) then
                    ch(:,ic) = czero_c
                    cht(:,ic) = czero_c
                end if
            end if
            self % istp = self % istp + 1
            ic = ic + 1
            if (ic > 4) ic = 1
        end do
    end if

    if (imove < 0 .and. istp_req == 3) then
        ! Padding in first column because tswd < dt_swd. 
        ! NOTE: Should be updated to the new padding scheme (not important)
        h(:,1) = 3*(h(:,2) - h(:,3)) + h(:,4)
        ht(:,1) = 3*(ht(:,2) - ht(:,3)) + ht(:,4)
        c(:,1) = 3*(c(:,2) - c(:,3)) + c(:,4)
        ct(:,1) = 3*(ct(:,2) - ct(:,3)) + ct(:,4)
    end if

    if (self % eof .and. self % tswd > (self % istp - 1) * self % dt) then
        err_msg(1) = 'User time beyond what is available from SWD file:'
        err_msg(2) = self % file
        err_msg(3) = 'Extrapolation of spectral data is not performed.'
        write(err_msg(4),'(a,f0.5)') 'Max stored swd-time = ', &
                          (self % istp - 1) * self % dt
        write(err_msg(5),'(a,f0.5)') 'Requested swd-time = ', self % tswd
        write(err_msg(6),'(a,f0.5)') 'Requested user-time = ', time
        call self % error % set_id_msg(err_proc, 1004, err_msg(1:6))
        return
    end if
    
    i1 = ip(1,ic)
    i2 = ip(2,ic)
    i3 = ip(3,ic)
    i4 = ip(4,ic)
    do concurrent (i = 0 : self % nsumx)
        ! Potential and d/dt of potential
        call self % tpol % scheme(delta,       &
                    cmplx(c(i,i1), kind=wp),   &
                    cmplx(c(i,i2), kind=wp),   &
                    cmplx(c(i,i3), kind=wp),   &
                    cmplx(c(i,i4), kind=wp),   &
                    cmplx(ct(i,i1), kind=wp),  &
                    cmplx(ct(i,i2), kind=wp),  &
                    cmplx(ct(i,i3), kind=wp),  &
                    cmplx(ct(i,i4), kind=wp),  &
                    self % c_cur(i), self % ct_cur(i))
        ! Wave height and d/dt of wave height
        call self % tpol % scheme(delta,       &
                    cmplx(h(i,i1), kind=wp),   &
                    cmplx(h(i,i2), kind=wp),   &
                    cmplx(h(i,i3), kind=wp),   &
                    cmplx(h(i,i4), kind=wp),   &
                    cmplx(ht(i,i1), kind=wp),  &
                    cmplx(ht(i,i2), kind=wp),  &
                    cmplx(ht(i,i3), kind=wp),  &
                    cmplx(ht(i,i4), kind=wp),  &
                    self % h_cur(i), self % ht_cur(i))
    end do
    do concurrent (i = 0 : self % nsumxh)
        if (self % nsf > 1) then
            ! Potential and d/dt of potential
            call self % tpol % scheme(delta,        &
                        cmplx(ch(i,i1), kind=wp),   &
                        cmplx(ch(i,i2), kind=wp),   &
                        cmplx(ch(i,i3), kind=wp),   &
                        cmplx(ch(i,i4), kind=wp),   &
                        cmplx(cht(i,i1), kind=wp),  &
                        cmplx(cht(i,i2), kind=wp),  &
                        cmplx(cht(i,i3), kind=wp),  &
                        cmplx(cht(i,i4), kind=wp),  &
                        self % ch_cur(i), self % cht_cur(i))
        else if (self % nsf == 1) then
            gamma_h = exp(- 2 * i * self % dk * self % d)
            self % ch_cur(i) = gamma_h * self % c_cur(i) 
            self % cht_cur(i) = gamma_h * self % ct_cur(i) 
        end if
    end do          
end associate

if (.not. self % dc_bias) then
    self % c_cur(0) = czero
    self % ct_cur(0) = czero
    self % h_cur(0) = czero
    self % ht_cur(0) = czero
    if (self % nsumxh > -1) then
        self % ch_cur(0) = czero
        self % cht_cur(0) = czero
    end if
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
end subroutine update_time2

!==============================================================================

function ZfunTaylor(z, j, dk, order) result(res) ! Value of Zfun(j) based on Taylor expansion
real(knd),       intent(in) :: z   ! z-position (>0)
integer,         intent(in) :: j   ! Index of Sfun
real(wp),        intent(in) :: dk  ! self % dk
integer,         intent(in) :: order ! self % order
real(wp)                    :: res
!
integer :: p
real(wp) :: ap1, apj
!
ap1 = dk * z * j
apj = 1.0_wp
res = 1.0_wp
do p = 1, order - 1
    apj = apj * ap1 / p
    res = res + apj
end do
!
end function ZfunTaylor

!==============================================================================

function ZhfunTaylor(z, j, dk, order) result(res) ! Value of Zhfun(j) based on Taylor expansion
real(knd),       intent(in) :: z   ! z-position (>0)
integer,         intent(in) :: j   ! Index of Sfun
real(wp),        intent(in) :: dk  ! self % dk
integer,         intent(in) :: order ! self % order
real(wp)                    :: res
!
integer :: p
real(wp) :: ap1, apj
!
ap1 = dk * z * j
apj = 1.0_wp
res = 1.0_wp
do p = 1, order - 1
    apj = apj * ap1 / p
    res = res + apj
end do
!
end function ZhfunTaylor

!==============================================================================

function phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self  ! Actual class
real(knd), intent(in) :: x,y,z ! Position application program
real(knd)             :: res   ! Potential at (x,y,z)
!
integer :: j
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
! Includes contributions from j=0
if (self % nsumxh < 0) then
    res = self % c_cur(0) % re
else
    res = self % c_cur(0) % re + self % ch_cur(0) % re
end if
if (z > 0 .and. self % norder > 0) then
    Xfun = 1.0_wp
    do j = 1, self % nsumx
        Xfun = kappa1 * Xfun
        Zfun = ZfunTaylor(z, j, self %dk, self % norder) 
        res = res + real(self % c_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    do j = 1, self % nsumxh
        Xfun = kappa1 * Xfun
        Zhfun = ZhfunTaylor(z, j, self %dk, self % norder) 
        res = res + real(self % ch_cur(j) * Xfun) * Zhfun
    end do
else
    kappa2 = exp(self % dk * z)
    kappa3 = 1.0_wp / kappa2
    Xfun = 1.0_wp
    Zfun = 1.0_wp
    do j = 1, self % nsumx
        Xfun = kappa1 * Xfun
        Zfun = kappa2 * Zfun
        res = res + real(self % c_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    Zhfun = 1.0_wp
    do j = 1, self % nsumxh
        Xfun = kappa1 * Xfun
        Zhfun = kappa3 * Zhfun
        res = res + real(self % ch_cur(j) * Xfun) * Zhfun
    end do
end if
!
end function phi

!==============================================================================

function stream(self, x, y, z) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self  ! Actual class
real(knd), intent(in) :: x,y,z ! Position application program
real(knd)             :: res   ! Stream function at (x,y,z)
!
integer :: j
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
! Includes contributions from j=0
if (self % nsumxh < 0) then
    res = self % c_cur(0) % im
else
    res = self % c_cur(0) % im - self % ch_cur(0) % im
end if
if (z > 0 .and. self % norder > 0) then
    Xfun = 1.0_wp
    do j = 1, self % nsumx
        Xfun = kappa1 * Xfun
        Zfun = ZfunTaylor(z, j, self %dk, self % norder) 
        res = res + aimag(self % c_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    do j = 1, self % nsumxh
        Xfun = kappa1 * Xfun
        Zhfun = ZhfunTaylor(z, j, self %dk, self % norder) 
        res = res - aimag(self % ch_cur(j) * Xfun) * Zhfun
    end do
else
    Xfun = 1.0_wp
    kappa2 = exp(self % dk * z)
    Zfun = 1.0_wp
    do j = 1, self % nsumx
        Xfun = kappa1 * Xfun
        Zfun = kappa2 * Zfun
        res = res + aimag(self % c_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    kappa3 = 1.0_wp / kappa2
    Zhfun = 1.0_wp
    do j = 1, self % nsumxh
        Xfun = kappa1 * Xfun
        Zhfun = kappa3 * Zhfun
        res = res - aimag(self % ch_cur(j) * Xfun) * Zhfun
    end do
end if
!
end function stream

!==============================================================================

function phi_t(self, x, y, z) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self  ! Actual class
real(knd),       intent(in) :: x,y,z ! Position application program
real(knd)                   :: res   ! Euler time derivative of potential at (x,y,z)
!
integer :: j
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
! Includes contributions from j=0
if (self % nsumxh < 0) then
    res = self % ct_cur(0) % re
else
    res = self % ct_cur(0) % re + self % cht_cur(0) % re
end if
if (z > 0 .and. self % norder > 0) then
    Xfun = 1.0_wp
    do j = 1, self % nsumx
        Xfun = kappa1 * Xfun
        Zfun = ZfunTaylor(z, j, self %dk, self % norder) 
        res = res + real(self % ct_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    do j = 1, self % nsumxh
        Xfun = kappa1 * Xfun
        Zhfun = ZhfunTaylor(z, j, self %dk, self % norder) 
        res = res + real(self % cht_cur(j) * Xfun) * Zhfun
    end do
else
    kappa2 = exp(self % dk * z)
    kappa3 = 1.0_wp / kappa2
    Xfun = 1.0_wp
    Zfun = 1.0_wp
    do j = 1, self % nsumx
        Xfun = kappa1 * Xfun
        Zfun = kappa2 * Zfun
        res = res + real(self % ct_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    Zhfun = 1.0_wp
    do j = 1, self % nsumxh
        Xfun = kappa1 * Xfun
        Zhfun = kappa3 * Zhfun
        res = res + real(self % cht_cur(j) * Xfun) * Zhfun
    end do
end if
!
end function phi_t

!==============================================================================

function grad_phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res(3) ! Particle velocity at (x,y,z)
!
integer :: j
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
real(wp) :: phi_xswd, phi_z, kval
complex(wp) :: kappa1, Xfun, cxz
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
phi_xswd = 0.0_wp  ! Includes contribution from j=0
phi_z = 0.0_wp     ! Includes contribution from j=0
if (z > 0 .and. self % norder > 0) then
    Xfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = ZfunTaylor(z, j, self %dk, self % norder) 
        cxz = (self % c_cur(j) * Xfun) * (kval * Zfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z + cxz % re
    end do
    Xfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = ZhfunTaylor(z, j, self %dk, self % norder) 
        cxz = (self % ch_cur(j) * Xfun) * (kval * Zhfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z - cxz % re
    end do
else
    kappa2 = exp(self % dk * z)
    Xfun = 1.0_wp
    kval = 0.0_wp
    Zfun = 1.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = kappa2 * Zfun
        cxz = (self % c_cur(j) * Xfun) * (kval * Zfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z + cxz % re
    end do
    kappa3 = 1.0_wp / kappa2
    Xfun = 1.0_wp
    kval = 0.0_wp
    Zhfun = 1.0_wp
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = kappa3 * Zhfun
        cxz = (self % ch_cur(j) * Xfun) * (kval * Zhfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z - cxz % re
    end do
end if
!
res(1) = phi_xswd * self % cbeta
res(2) = phi_xswd * self % sbeta
res(3) = phi_z
!
end function grad_phi

!==============================================================================

function grad_phi_2nd(self, x, y, z) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
real(knd),       intent(in) :: x,y,z  ! Position application program
real(knd)                   :: res(6) ! Second order gradients of potential at (x,y,z)
                                      ! res(1) = d^2(potential) / dx^2
                                      ! res(2) = d^2(potential) / dx dy
                                      ! res(3) = d^2(potential) / dx dz
                                      ! res(4) = d^2(potential) / dy dy
                                      ! res(5) = d^2(potential) / dy dz
                                      ! res(6) = d^2(potential) / dz dz
!
integer :: j
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
real(wp) :: phi_xx_swd, phi_xz_swd, kval
complex(wp) :: kappa1, Xfun, cof
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
phi_xx_swd = 0.0_wp  ! Includes contribution from j=0
phi_xz_swd = 0.0_wp  ! Includes contribution from j=0
if (z > 0 .and. self % norder > 0) then
    Xfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = ZfunTaylor(z, j, self %dk, self % norder) 
        cof = (self % c_cur(j) * Xfun) * (kval * kval * Zfun)
        phi_xx_swd = phi_xx_swd - cof % re
        phi_xz_swd = phi_xz_swd + cof % im
    end do
    Xfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = ZhfunTaylor(z, j, self %dk, self % norder) 
        cof = (self % ch_cur(j) * Xfun) * (kval * kval * Zhfun)
        phi_xx_swd = phi_xx_swd - cof % re
        phi_xz_swd = phi_xz_swd - cof % im
    end do
else
    kappa2 = exp(self % dk * z)
    Xfun = 1.0_wp
    Zfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = kappa2 * Zfun
        cof = (self % c_cur(j) * Xfun) * (kval * kval * Zfun)
        phi_xx_swd = phi_xx_swd - cof % re
        phi_xz_swd = phi_xz_swd + cof % im
    end do
    kappa3 = 1.0_wp / kappa2
    Xfun = 1.0_wp
    Zhfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = kappa3 * Zhfun
        cof = (self % ch_cur(j) * Xfun) * (kval * kval * Zhfun)
        phi_xx_swd = phi_xx_swd - cof % re
        phi_xz_swd = phi_xz_swd - cof % im
    end do
end if
!
res(1) = phi_xx_swd * self % cbeta * self % cbeta
res(2) = phi_xx_swd * self % sbeta * self % cbeta
res(3) = phi_xz_swd * self % cbeta
res(4) = phi_xx_swd * self % sbeta * self % sbeta
res(5) = phi_xz_swd * self % sbeta
res(6) = - phi_xx_swd
!
end function grad_phi_2nd

!==============================================================================

function acc_euler(self, x, y, z) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res(3) ! Euler acceleration at (x,y,z)
!
integer :: j
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
real(wp) :: phit_xswd, phit_z, kval
complex(wp) :: kappa1, Xfun, cof
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
phit_xswd = 0.0_wp  ! Includes contribution from j=0
phit_z = 0.0_wp     ! Includes contribution from j=0
kval = 0.0_wp
if (z > 0 .and. self % norder > 0) then
    Xfun = 1.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = ZfunTaylor(z, j, self % dk, self % norder)
        cof = (self % ct_cur(j) * Xfun) * (kval * Zfun)
        phit_xswd = phit_xswd + cof % im
        phit_z = phit_z + cof % re
    end do
    Xfun = 1.0_wp
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = ZhfunTaylor(z, j, self % dk, self % norder)
        cof = (self % cht_cur(j) * Xfun) * (kval * Zhfun)
        phit_xswd = phit_xswd + cof % im
        phit_z = phit_z - cof % re
    end do
else
    kappa2 = exp(self % dk * z)
    Xfun = 1.0_wp
    Zfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = kappa2 * Zfun
        cof = (self % ct_cur(j) * Xfun) * (kval * Zfun)
        phit_xswd = phit_xswd + cof % im
        phit_z = phit_z + cof % re
    end do
    kappa3 = 1.0_wp / kappa2
    Xfun = 1.0_wp
    Zhfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = kappa3 * Zhfun
        cof = (self % cht_cur(j) * Xfun) * (kval * Zhfun)
        phit_xswd = phit_xswd + cof % im
        phit_z = phit_z - cof % re
    end do
end if
res(1) = phit_xswd * self % cbeta
res(2) = phit_xswd * self % sbeta
res(3) = phit_z
!
end function acc_euler

!==============================================================================

function acc_particle(self, x, y, z) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
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
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res    ! Fully nonlinear pressure
!
integer :: j
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
real(wp) :: phi_xswd, phi_z, phi_t, kval
complex(wp) :: kappa1, Xfun, cxz
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
kappa2 = exp(self % dk * z)
phi_xswd = 0.0_wp   ! Includes contribution from j=0
phi_z = 0.0_wp      ! Includes contribution from j=0
! Includes contribution from j=0
if (self % nsumxh < 0) then
    phi_t = self % ct_cur(0) % re
else
    phi_t = self % ct_cur(0) % re + self % cht_cur(0) % re
end if
if (z > 0 .and. self % norder > 0) then
    Xfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = ZfunTaylor(z, j, self % dk, self % norder) 
        cxz = (self % c_cur(j) * Xfun) * (kval * Zfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z + cxz % re
        phi_t = phi_t + real(self % ct_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    kval = 0.0_wp
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = ZhfunTaylor(z, j, self % dk, self % norder) 
        cxz = (self % ch_cur(j) * Xfun) * (kval * Zhfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z - cxz % re
        phi_t = phi_t + real(self % cht_cur(j) * Xfun) * Zhfun
    end do
else
    Xfun = 1.0_wp
    kval = 0.0_wp
    Zfun = 1.0_wp
    do j = 1, self % nsumx
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zfun = kappa2 * Zfun
        cxz = (self % c_cur(j) * Xfun) * (kval * Zfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z + cxz % re
        phi_t = phi_t + real(self % ct_cur(j) * Xfun) * Zfun
    end do
    Xfun = 1.0_wp
    kval = 0.0_wp
    Zhfun = 1.0_wp
    kappa3 = 1.0_wp / kappa2
    do j = 1, self % nsumxh
        kval = kval + self % dk
        Xfun = kappa1 * Xfun
        Zhfun = kappa3 * Zhfun
        cxz = (self % ch_cur(j) * Xfun) * (kval * Zhfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z - cxz % re
        phi_t = phi_t + real(self % cht_cur(j) * Xfun) * Zhfun
    end do
end if
res = (-phi_t - (phi_xswd**2 + phi_z**2) * 0.5_wp - z * self % grav) * self % rho
!
end function pressure

!==============================================================================

subroutine convergence(self, x, y, z, csv)
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self   ! Actual class
real(knd),        intent(in) :: x,y,z  ! Position application program
character(len=*), intent(in) :: csv    ! New output file
!
integer :: j, lucsv
real(wp) :: xswd, kappa2, kappa3, Zfun, Zhfun
real(wp) :: phi_xswd, phi_x, phi_y, phi_z, phi_t, kval, elev, prs
complex(wp) :: kappa1, Xfun, cxz
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::convergence'
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
write(lucsv,'(a,f0.9)') 'x = ', x
write(lucsv,'(a,f0.9)') 'y = ', y
write(lucsv,'(a,f0.9)') 'z = ', z
write(lucsv,'(a,f0.9)') 't = ', self % tswd - self % t0
write(lucsv,'(100(a,:","))') 'jx', 'velx', 'vely', 'velz', 'elev', 'prs'

xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
kappa2 = exp(self % dk * z)
kappa3 = 1.0_wp / kappa2
Xfun = 1.0_wp
phi_xswd = 0.0_wp  ! Includes contribution from j=0
phi_z = 0.0_wp     ! Includes contribution from j=0
! Includes contribution from j=0
if (self % nsumxh < 0) then
    phi_t = self % ct_cur(0) % re
else
    phi_t = self % ct_cur(0) % re + self % cht_cur(0) % re
end if
elev = self % h_cur(0) % re    ! Includes contribution from j=0
phi_x = phi_xswd * self % cbeta
phi_y = phi_xswd * self % sbeta
prs = (-phi_t - (phi_xswd**2 + phi_z**2) * 0.5_wp - z * self % grav) * self % rho
write(lucsv, '(100(f0.13,:","))') real(0,wp), phi_x, phi_y, phi_z, elev, prs
kval = 0.0_wp
Zfun = 1.0_wp
Zhfun = 1.0_wp
do j = 1, self % nsumx
    kval = kval + self % dk
    Xfun = kappa1 * Xfun
    if (z > 0 .and. self % norder > 0) then
        Zfun = ZfunTaylor(z, j, self % dk, self % norder) 
    else
        Zfun = kappa2 * Zfun
    end if
    cxz = (self % c_cur(j) * Xfun) * (kval * Zfun)
    phi_xswd = phi_xswd + cxz % im
    phi_z = phi_z + cxz % re
    phi_t = phi_t + real(self % ct_cur(j) * Xfun) * Zfun
    elev = elev + real(self % h_cur(j) * Xfun)
    if (j <= self % nsumxh) then
        if (z > 0 .and. self % norder > 0) then
            Zhfun = ZhfunTaylor(z, j, self % dk, self % norder) 
        else
            Zhfun = kappa3 * Zhfun
        end if
        cxz = (self % ch_cur(j) * Xfun) * (kval * Zhfun)
        phi_xswd = phi_xswd + cxz % im
        phi_z = phi_z - cxz % re
        phi_t = phi_t + real(self % cht_cur(j) * Xfun) * Zhfun
    end if
    ! It is assumed that self % nsumxh <= self % nsumx
    phi_x = phi_xswd * self % cbeta
    phi_y = phi_xswd * self % sbeta
    prs = (-phi_t - (phi_xswd**2 + phi_z**2) * 0.5_wp - z * self % grav) * self % rho
    write(lucsv, '(100(f0.13,:","))') real(j,wp), phi_x, phi_y, phi_z, elev, prs
end do
close(lucsv)
!
end subroutine convergence

!==============================================================================

subroutine strip(self, tmin, tmax, file_swd)
! Store part of SWD file into a new SWD file
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self   ! Actual class
real(knd),        intent(in) :: tmin, tmax ! Time window to store
character(len=*), intent(in) :: file_swd   ! Name of new swd file
!
integer(int64) :: ipos, ipos_save
integer :: i, ios, lures, istep_first, istep_last, n_camp_per_time_step
real(wp) :: tmin_swd, tmax_swd, tmax_allow
integer(c_int) fmt, shp, amp, nid, nstrip, nsteps, order, n, nh, isf, nsf
real(c_float) dk, grav, lscale, dt, xzsf, magic
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
complex(c_float), allocatable :: camp(:)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::strip'
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

read(self % unit, end=98, err=99) n
write(lures) n

read(self % unit, end=98, err=99) nh
write(lures) nh

read(self % unit, end=98, err=99) dk
write(lures) dk

read(self % unit, end=98, err=99) isf
write(lures) isf

read(self % unit, end=98, err=99) nsf
write(lures) nsf

do i = 1, 2 * nsf !Read/write both xsf and zsf
    read(self % unit, end=98, err=99) xzsf
    write(lures) xzsf
end do

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
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self ! Actual class
real(knd),       intent(in) :: x,y  ! Position application program
real(knd)                   :: res  ! Surface elevation at (x,y)
!
integer :: j
real(wp) :: xswd
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
res = self % h_cur(0) % re ! Includes contribution from j=0
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    res = res + real(self % h_cur(j) * Xfun)
end do
!
end function elev

!==============================================================================

function elev_t(self, x, y) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self ! Actual class
real(knd),       intent(in) :: x,y  ! Position application program
real(knd)                   :: res  ! d/dt of surface elevation at (x,y)
!
integer :: j
real(wp) :: xswd
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
res = self % ht_cur(0) % re    ! Includes contribution from j=0
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    res = res + real(self % ht_cur(j) * Xfun, wp)
end do
!
end function elev_t

!==============================================================================

function grad_elev(self, x, y) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
real(knd),          intent(in) :: x,y    ! Position application program
real(knd)                      :: res(3) ! x, y and z gradients of surface elevation at (x,y)
!
integer :: j
real(wp) :: xswd, elev_x_swd, kval
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
elev_x_swd = 0.0_wp  ! Includes contribution from j=0
kval = 0.0_wp
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    kval = kval + self % dk
    elev_x_swd = elev_x_swd + kval * aimag(self % h_cur(j) * Xfun)
end do
res(1) = elev_x_swd * self % cbeta
res(2) = elev_x_swd * self % sbeta
res(3) = 0.0_knd
!
end function grad_elev

!==============================================================================

function grad_elev_2nd(self, x, y) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
real(knd),          intent(in) :: x,y    ! Position application program
real(knd)                      :: res(3) ! Second order gradients of surface elevation
                                         ! res(1) = d^2(elevation) / dx^2
                                         ! res(2) = d^2(elevation) / dx dy
                                         ! res(3) = d^2(elevation) / dy dy
!
integer :: j
real(wp) :: xswd, elev_xx_swd, kval
complex(wp) :: kappa1, Xfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
kappa1 = exp(cmplx(0.0_wp, -self % dk * xswd, kind=wp))
Xfun = 1.0_wp
elev_xx_swd = 0.0_wp  ! Includes contribution from j=0
kval = 0.0_wp
do j = 1, self % nsumx
    Xfun = kappa1 * Xfun
    kval = kval + self % dk
    elev_xx_swd = elev_xx_swd - kval * kval * real(self % h_cur(j) * Xfun)
end do
res(1) = elev_xx_swd * self % cbeta * self % cbeta
res(2) = elev_xx_swd * self % sbeta * self % cbeta
res(3) = elev_xx_swd * self % sbeta * self % sbeta
!
end function grad_elev_2nd

!==============================================================================

function bathymetry(self, x, y) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self ! Actual class
real(knd), intent(in) :: x,y  ! Position application program
real(knd)             :: res  ! Local water depth
!
integer :: i
real(wp) :: xswd, xdom
!
if (self % isf == 0) then
    if (self % nsf == 0) then
        res = -1.0_wp
        return
    else if (self % nsf == 1) then
        res = self % d
    else
        ! The bathymetry is periodic
        xswd = self % x0 + x * self % cbeta + y * self % sbeta
        xdom = 2 * pi / self % dk
        xswd = modulo(xswd, xdom)   ! xswd in [0, xdom)
        do i = 2, self % nsf
            if (self % xsf(i) > xswd) then
                ! Linear interpolation for local depth (sign shift)
                res = -self % zsf(i-1) - (xswd - self % xsf(i-1)) *  &
                       (self % zsf(i) - self % zsf(i-1)) /           &
                       (self % xsf(i) - self % xsf(i-1))
                return
            end if
        end do
    end if
end if
!
end function bathymetry

!==============================================================================

function bathymetry_nvec(self, x, y) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y    ! Position application program
real(knd)             :: res(3) ! Unit normal vector into ocean at (x,y)
!
integer :: i
real(wp) :: xswd, xdom, dzdx_swd, nvec_swd(3)
!
if (self % isf == 0) then
    if (self % nsf < 2) then
        res = [0.0_wp, 0.0_wp, 1.0_wp]
        return
    else
        ! The bathymetry is periodic
        xswd = self % x0 + x * self % cbeta + y * self % sbeta
        xdom = 2 * pi / self % dk
        xswd = modulo(xswd, xdom)   ! xswd in [0, xdom)
        do i = 2, self % nsf
            if (self % xsf(i) > xswd) then
                ! Exact for linear interpolation
                dzdx_swd = (self % zsf(i) - self % zsf(i-1)) / &
                           (self % xsf(i) - self % xsf(i-1))
                nvec_swd = [-dzdx_swd, 0.0_wp, 1.0_wp] / &
                           sqrt(1.0_wp + dzdx_swd**2)
                res = [nvec_swd(1) * self % cbeta, &
                       nvec_swd(1) * self % sbeta, &
                       nvec_swd(3)]  ! nvec_swd(2)=0
                return
            end if
        end do
    end if
end if
!
end function bathymetry_nvec

!==============================================================================

function get_int(self, name) result(res)
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of int parameter
integer                      :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::get_int'
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
case('n')
    res = self % n
case('nh')
    res = self % nh
case('nsumxh')
    res = self % nsumxh
case('isf')
    res = self % isf
case('nsf')
    res = self % nsf
case('ipol')
    res = self % ipol
case('impl')
    res = 1
case('nsumx')
    res = self % nsumx
case('nsumy')
    res = -1
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
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of logical parameter
logical                      :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::get_logical'
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
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of float parameter
real(knd)                    :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::get_real'
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
case('dk')
    res = self % dk
! Some special requests....
case('tmax')
    res = self % tmax
case('d')
    res = self % d
case('lmin')
    res = 2 * pi / (self % dk * self % nsumx)
case('lmax')
    res = 2 * pi / self % dk
case('sizex')
    res = 2 * pi / self % dk
case('sizey')
    res = 0.0_knd
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
class(spectral_wave_data_shape_3_impl_1), intent(inout) :: self   ! Actual class
character(len=*),  intent(in) :: name ! Name of char parameter
character(len=:), allocatable :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_3_impl_1::get_chr'
character(len=250) :: err_msg(1)
!
select case(name)
case('file', 'file_swd')
    res = self % file
case('version')
    res = version
case('class')
    res = 'spectral_wave_data_shape_3_impl_1'
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

end module spectral_wave_data_shape_3_impl_1_def
