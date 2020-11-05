module spectral_wave_data_shape_6_impl_1_def

use, intrinsic :: iso_fortran_env, only: int64
use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float

use kind_values, only: knd => kind_swd_interface, wp => kind_swd_internal

use open_swd_file_def, only: open_swd_file, swd_validate_binary_convention, &
                             swd_magic_number
use spectral_wave_data_def, only: spectral_wave_data
use swd_version, only: version

implicit none
private

! This module provides an implemention of the shape 6 class of the 
! spectral-wave-data API.
!
! Written by Jens Bloch Helmers, November, 16. 2019
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: spectral_wave_data_shape_6_impl_1
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type wave
    real(wp)  :: omg ! Wave frequencies (rad/s)
    real(wp)  :: kw  ! Wave numbers
    real(wp)  :: kwx ! Wave numbers x-direction
    real(wp)  :: kwy ! Wave numbers y-direction
    real(wp)  :: amp ! Single amplitudes
    real(wp)  :: phs ! Phase (radians)
    real(wp)  :: bfac ! exp(-2*kw*d)
    real(wp)  :: fracbfac ! 1/(1 + bfac)
    real(wp)  :: tanhkd   ! tanh(kw*d)
    complex(wp) :: c_cur  ! Input window of c-spectral components (1:n)
    complex(wp) :: ct_cur ! Input window of ct-spectral components (1:n)
    complex(wp) :: h_cur  ! Input window of h-spectral components (1:n)
    complex(wp) :: ht_cur ! Input window of ht-spectral components (1:n)
end type wave

type, extends(spectral_wave_data) :: spectral_wave_data_shape_6_impl_1
    integer   :: n      ! Number of spectral components
    logical   :: lcrest ! True if long-crested wave field
    real(wp)  :: d      ! Constant water depth (<0 is infinite depth)
    type(wave), allocatable :: wd(:) ! Wave component data
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
end type spectral_wave_data_shape_6_impl_1

interface spectral_wave_data_shape_6_impl_1
    module procedure constructor
end interface

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

subroutine close(self)
class(spectral_wave_data_shape_6_impl_1) :: self  ! Object to destruct
!
logical opened
!
inquire(unit=self % unit, opened=opened)
if (opened) close(self % unit)
if (allocated(self % wd)) deallocate(self % wd)
!
self % file = '0'
self % unit = 0
self % n = 0
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
type(spectral_wave_data_shape_6_impl_1) :: self  ! Object to construct
!
integer :: i, j, ios, err_id

integer(c_int) :: fmt, shp, amp, n, order, nid, nsteps, nstrip
real(c_float) :: kw, gam, gam_prev, wamp, phs, dt, grav, lscale, d, magic
real(wp) :: gamwp
character(kind=c_char, len=:), allocatable :: cid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
complex(c_float), parameter :: czero_c = cmplx(0.0_c_float, 0.0_c_float, c_float)
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_6_impl_1::constructor'
character(len=250) :: err_msg(5)
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
        
if (present(dc_bias)) then
    self % dc_bias = dc_bias
else
    self % dc_bias = .false.
end if

self % sbeta = sin(beta*pi/180.0_wp)
self % cbeta = cos(beta*pi/180.0_wp)
self % tmax = huge(self % tmax) / 100

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
self % nstrip = 0  ! Not relevant

read(self % unit, end=98, err=99) nsteps
self % nsteps = 0  ! Not relevant

read(self % unit, end=98, err=99) dt
self % dt = huge(self % dt) * 0.1_wp ! Not relevant

read(self % unit, end=98, err=99) order ! Not relevant
self % order = 0  ! By default linear waves should be consistent between wave generators
self % norder = self % order
if (present(norder)) then
    if (norder < 3) then
        self % norder = norder
    else
        write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
        write(err_msg(2),'(a)') 'For shp=6, norder should be less than 3.'
        write(err_msg(3),'(a,i0)') 'norder = ', norder
        call self % error % set_id_msg(err_proc, 1004, err_msg(1:3))
        return
    end if
end if

read(self % unit, end=98, err=99) n
self % n = n

read(self % unit, end=98, err=99) d
self % d = d

allocate( self % wd(n), stat=i)
if (i /= 0) then
    write(err_msg(1),'(a,a)') 'Input file: ', trim(self % file)
    write(err_msg(2),'(a)') 'Not able to allocate space'
    write(err_msg(3),'(a,i0)') 'Number of spectral components (n) is ', self % n
    call self % error % set_id_msg(err_proc, 1005, err_msg(1:3))
    return
end if

self % lcrest = .true.
associate(wd => self % wd)
    do j = 1, self % n
        read(self % unit, end=98, err=99) wamp, kw, gam, phs
        if (j > 1) then
            if (abs(gam - gam_prev) > 0.00001_wp) then
                self % lcrest = .false.
            end if
        end if
        wd(j) % amp = wamp
        wd(j) % phs = phs
        gamwp = gam
        wd(j) % kw = kw
        wd(j) % kwx = wd(j) % kw * cos(gamwp)
        wd(j) % kwy = wd(j) % kw * sin(gamwp)
        if (self % d < 0.0_wp) then
            wd(j) % omg = sqrt(wd(j) % kw * self % grav)
            wd(j) % bfac = 0.0_wp
            wd(j) % fracbfac = 1.0_wp
            wd(j) % tanhkd = 1.0_wp
        else
            wd(j) % omg = sqrt(wd(j) % kw * self % grav * tanh(wd(j) % kw * self % d))
            wd(j) % bfac = exp(-2.0_wp * wd(j) % kw * self % d)
            wd(j) % fracbfac = 1.0_wp / (1.0_wp + wd(j) % bfac)
            wd(j) % tanhkd = tanh(wd(j) % kw * self % d)
        end if
        gam_prev = gam
    end do
end associate

close(self % unit)

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
end function constructor

!==============================================================================

subroutine update_time(self, time)
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self
real(knd), intent(in) :: time  ! Current time in simulation program
!
integer :: j
complex(wp), parameter :: imag = cmplx(0.0_wp, 1.0_wp, wp)
!
self % tswd = self % t0 + time
associate(wd => self % wd)
    do concurrent (j = 1 : self % n)
        wd(j) % h_cur  = wd(j) % amp * exp(cmplx(0.0_wp,   &
                         wd(j) % omg * self % tswd + wd(j) % phs, wp))
        wd(j) % ht_cur = imag * wd(j) % omg * wd(j) % h_cur
        wd(j) % c_cur  = imag * self % grav * wd(j) % h_cur / wd(j) % omg
        wd(j) % ct_cur = - self % grav * wd(j) % h_cur
    end do
end associate
!
end subroutine update_time

!==============================================================================

function phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self  ! Actual class
real(knd), intent(in) :: x,y,z ! Position application program
real(knd)             :: res   ! Potential at (x,y,z)
!
integer :: j
real(wp) :: xswd, yswd, Zfun, zwp, reswp, a, c, elev
complex(wp) :: Efun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
if (self % norder == 0) then
    zwp = min(z, 0.0_wp)
else if (self % norder < 2) then
    zwp = z
else
    ! Wheeler stretching
    elev = self % elev(x, y)
    zwp = z - elev
    if (self % d > 0.0_wp) then
        zwp = zwp / (1.0_wp + elev / self % d)
    end if
end if
reswp = 0.0_wp
if (self % d < 0.0_wp) then
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &
                      c_cur => self % wd(j) % c_cur)
                Zfun = exp(kw * zwp)
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(c_cur * Efun) * Zfun
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      c_cur => self % wd(j) % c_cur)
                Zfun = 1.0_wp + kw * zwp
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(c_cur * Efun) * Zfun
            end associate
        end do
    end if
else
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,              &
                      kwx => self % wd(j) % kwx,            &
                      kwy => self % wd(j) % kwy,            &  
                      bfac => self % wd(j) % bfac,          &
                      fracbfac => self % wd(j) % fracbfac,  &
                      c_cur => self % wd(j) % c_cur)
                a = exp(kw * zwp)
                c = bfac / a**2
                Zfun = a * (1.0_wp + c) * fracbfac
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(c_cur * Efun) * Zfun
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,       &
                      kwx => self % wd(j) % kwx,     &
                      kwy => self % wd(j) % kwy,     &  
                      c_cur => self % wd(j) % c_cur, &
                      tanhkd => self % wd(j) % tanhkd)
                Zfun = 1.0_wp + tanhkd * kw * zwp
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(c_cur * Efun) * Zfun
            end associate
        end do
    end if
end if
res = reswp
!
end function phi

!==============================================================================

function stream(self, x, y, z) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self  ! Actual class
real(knd), intent(in) :: x,y,z ! Position application program
real(knd)             :: res   ! Stream function at (x,y,z)
integer :: j
real(wp) :: xswd, yswd, Zhfun, zwp, reswp, a, c, elev
complex(wp) :: Efun
!
if (.not. self % lcrest) then
    res = 0.0_knd   ! Stream functions make no sense in short crested seas
    return
end if
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
if (self % norder == 0) then
    zwp = min(z, 0.0_wp)
else if (self % norder < 2) then
    zwp = z
else
    ! Wheeler stretching
    elev = self % elev(x, y)
    zwp = z - elev
    if (self % d > 0.0_wp) then
        zwp = zwp / (1.0_wp + elev / self % d)
    end if
end if
reswp = 0.0_wp
if (self % d < 0.0_wp) then
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      c_cur => self % wd(j) % c_cur)
                Zhfun = exp(kw * zwp)
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + aimag(c_cur * Efun) * Zhfun
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      c_cur => self % wd(j) % c_cur)
                Zhfun = 1.0_wp + kw * zwp
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + aimag(c_cur * Efun) * Zhfun
            end associate
        end do
    end if
else
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,              &
                      kwx => self % wd(j) % kwx,            &
                      kwy => self % wd(j) % kwy,            &  
                      bfac => self % wd(j) % bfac,          &
                      fracbfac => self % wd(j) % fracbfac,  &
                      c_cur => self % wd(j) % c_cur)
                a = exp(kw * zwp)
                c = bfac / a**2
                Zhfun = a * (1.0_wp - c) * fracbfac
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + aimag(c_cur * Efun) * Zhfun
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,       &
                      kwx => self % wd(j) % kwx,     &
                      kwy => self % wd(j) % kwy,     &  
                      c_cur => self % wd(j) % c_cur, &
                      tanhkd => self % wd(j) % tanhkd)
                Zhfun = tanhkd + kw * zwp
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + aimag(c_cur * Efun) * Zhfun
            end associate
        end do
    end if
end if
res = reswp
!
end function stream

!==============================================================================

function phi_t(self, x, y, z) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self  ! Actual class
real(knd),       intent(in) :: x,y,z ! Position application program
real(knd)                   :: res   ! Euler time derivative of potential at (x,y,z)
!
integer :: j
real(wp) :: xswd, yswd, Zfun, zwp, reswp, a, c, elev
complex(wp) :: Efun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
if (self % norder == 0) then
    zwp = min(z, 0.0_wp)
else if (self % norder < 2) then
    zwp = z
else
    ! Wheeler stretching
    elev = self % elev(x, y)
    zwp = z - elev
    if (self % d > 0.0_wp) then
        zwp = zwp / (1.0_wp + elev / self % d)
    end if
end if
reswp = 0.0_wp
if (self % d < 0.0_wp) then
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      ct_cur => self % wd(j) % ct_cur)
                Zfun = exp(kw * zwp)
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(ct_cur * Efun) * Zfun
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      ct_cur => self % wd(j) % ct_cur)
                Zfun = 1.0_wp + kw * zwp
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(ct_cur * Efun) * Zfun
            end associate
        end do
    end if
else
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,              &
                      kwx => self % wd(j) % kwx,            &
                      kwy => self % wd(j) % kwy,            &  
                      bfac => self % wd(j) % bfac,          &
                      fracbfac => self % wd(j) % fracbfac,  &
                      ct_cur => self % wd(j) % ct_cur)
                a = exp(kw * zwp)
                c = bfac / a**2
                Zfun = a * (1.0_wp + c) * fracbfac
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(ct_cur * Efun) * Zfun
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,         &
                      kwx => self % wd(j) % kwx,       &
                      kwy => self % wd(j) % kwy,       &  
                      ct_cur => self % wd(j) % ct_cur, &
                      tanhkd => self % wd(j) % tanhkd)
                Zfun = 1.0_wp + tanhkd * kw * zwp
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                reswp = reswp + real(ct_cur * Efun) * Zfun
            end associate
        end do
    end if
end if
res = reswp
!
end function phi_t

!==============================================================================

function grad_phi(self, x, y, z) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res(3) ! Particle velocity at (x,y,z)
!
integer :: j
real(wp) :: xswd, yswd, Zfun, Zhfun, Zfun_z, zwp, reswp, a, c, elev
real(wp) :: phi_xswd, phi_yswd, phi_zswd
complex(wp) :: Efun, Cfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
if (self % norder == 0) then
    zwp = min(z, 0.0_wp)
else if (self % norder < 2) then
    zwp = z
else
    ! Wheeler stretching
    elev = self % elev(x, y)
    zwp = z - elev
    if (self % d > 0.0_wp) then
        zwp = zwp / (1.0_wp + elev / self % d)
    end if
end if
phi_xswd = 0.0_wp
phi_yswd = 0.0_wp
phi_zswd = 0.0_wp
if (self % d < 0.0_wp) then
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &
                      c_cur => self % wd(j) % c_cur)
                Zfun = exp(kw * zwp)
                Zfun_z = kw * Zfun  ! Zhfun = Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      c_cur => self % wd(j) % c_cur)
                Zfun = 1.0_wp + kw * zwp
                Zfun_z = kw * Zfun  ! Zhfun = Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
        end do
    end if
else
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,              &
                      kwx => self % wd(j) % kwx,            &
                      kwy => self % wd(j) % kwy,            &  
                      bfac => self % wd(j) % bfac,          &
                      fracbfac => self % wd(j) % fracbfac,  &
                      c_cur => self % wd(j) % c_cur)
                a = exp(kw * zwp)
                c = bfac / a**2
                Zfun = a * (1.0_wp + c) * fracbfac
                Zfun_z = kw * a * (1.0_wp - c) * fracbfac
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,       &
                      kwx => self % wd(j) % kwx,     &
                      kwy => self % wd(j) % kwy,     &  
                      c_cur => self % wd(j) % c_cur, &
                      tanhkd => self % wd(j) % tanhkd)
                Zfun = 1.0_wp + tanhkd * kw * zwp
                Zhfun = tanhkd + kw * zwp
                Zfun_z = kw * Zhfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
        end do
    end if
end if
res(1) = phi_xswd * self % cbeta - phi_yswd * self % sbeta
res(2) = phi_xswd * self % sbeta + phi_yswd * self % cbeta
res(3) = phi_zswd
!
end function grad_phi

!==============================================================================

function grad_phi_2nd(self, x, y, z) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
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
real(wp) :: xswd, yswd, zwp, Zfun, Zhfun, Zfun_z, Zfun_zz, elev, a, c, cc, cs, ss
real(wp) :: phi_xx_swd, phi_xy_swd, phi_xz_swd, phi_yy_swd, phi_yz_swd, phi_zz_swd
complex(wp) :: Efun, Cfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
if (self % norder == 0) then
    zwp = min(z, 0.0_wp)
else if (self % norder < 2) then
    zwp = z
else
    ! Wheeler stretching
    elev = self % elev(x, y)
    zwp = z - elev
    if (self % d > 0.0_wp) then
        zwp = zwp / (1.0_wp + elev / self % d)
    end if
end if
phi_xx_swd = 0.0_wp
phi_xy_swd = 0.0_wp
phi_xz_swd = 0.0_wp
phi_yy_swd = 0.0_wp
phi_yz_swd = 0.0_wp
phi_zz_swd = 0.0_wp
if (self % d < 0.0_wp) then
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &
                      c_cur => self % wd(j) % c_cur)
                Zfun = exp(kw * zwp)
                Zfun_z = kw * Zfun  ! Zhfun = Zfun
                Zfun_zz = kw * Zfun_z
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xx_swd = phi_xx_swd - kwx*kwx * Cfun % re * Zfun
                phi_xy_swd = phi_xy_swd - kwx*kwy * Cfun % re * Zfun
                phi_xz_swd = phi_xz_swd + kwx * Cfun % im * Zfun_z
                phi_yy_swd = phi_yy_swd - kwy*kwy * Cfun % re * Zfun
                phi_yz_swd = phi_yz_swd + kwy * Cfun % im * Zfun_z
                phi_zz_swd = phi_zz_swd + Cfun % re * Zfun_zz
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      c_cur => self % wd(j) % c_cur)
                Zfun = 1.0_wp + kw * zwp
                Zfun_z = kw * Zfun  ! Zhfun = Zfun
                Zfun_zz = kw * Zfun_z
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xx_swd = phi_xx_swd - kwx*kwx * Cfun % re * Zfun
                phi_xy_swd = phi_xy_swd - kwx*kwy * Cfun % re * Zfun
                phi_xz_swd = phi_xz_swd + kwx * Cfun % im * Zfun_z
                phi_yy_swd = phi_yy_swd - kwy*kwy * Cfun % re * Zfun
                phi_yz_swd = phi_yz_swd + kwy * Cfun % im * Zfun_z
                phi_zz_swd = phi_zz_swd + Cfun % re * Zfun_zz
            end associate
        end do
    end if
else
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,              &
                      kwx => self % wd(j) % kwx,            &
                      kwy => self % wd(j) % kwy,            &  
                      bfac => self % wd(j) % bfac,          &
                      fracbfac => self % wd(j) % fracbfac,  &
                      c_cur => self % wd(j) % c_cur)
                a = exp(kw * zwp)
                c = bfac / a**2
                Zfun = a * (1.0_wp + c) * fracbfac
                Zfun_z = kw * a * (1.0_wp - c) * fracbfac
                Zfun_zz = kw**2 * Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xx_swd = phi_xx_swd - kwx*kwx * Cfun % re * Zfun
                phi_xy_swd = phi_xy_swd - kwx*kwy * Cfun % re * Zfun
                phi_xz_swd = phi_xz_swd + kwx * Cfun % im * Zfun_z
                phi_yy_swd = phi_yy_swd - kwy*kwy * Cfun % re * Zfun
                phi_yz_swd = phi_yz_swd + kwy * Cfun % im * Zfun_z
                phi_zz_swd = phi_zz_swd + Cfun % re * Zfun_zz
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,       &
                      kwx => self % wd(j) % kwx,     &
                      kwy => self % wd(j) % kwy,     &  
                      c_cur => self % wd(j) % c_cur, &
                      tanhkd => self % wd(j) % tanhkd)
                Zfun = 1.0_wp + tanhkd * kw * zwp
                Zhfun = tanhkd + kw * zwp
                Zfun_z = kw * Zhfun
                Zfun_zz = kw * kw * Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_xx_swd = phi_xx_swd - kwx*kwx * Cfun % re * Zfun
                phi_xy_swd = phi_xy_swd - kwx*kwy * Cfun % re * Zfun
                phi_xz_swd = phi_xz_swd + kwx * Cfun % im * Zfun_z
                phi_yy_swd = phi_yy_swd - kwy*kwy * Cfun % re * Zfun
                phi_yz_swd = phi_yz_swd + kwy * Cfun % im * Zfun_z
                phi_zz_swd = phi_zz_swd + Cfun % re * Zfun_zz
            end associate
        end do
    end if
end if
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
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res(3) ! Euler acceleration at (x,y,z)
!
integer :: j
real(wp) :: xswd, yswd, Zfun, Zhfun, Zfun_z, zwp, reswp, a, c, elev
real(wp) :: phi_xtswd, phi_ytswd, phi_ztswd
complex(wp) :: Efun, Ctfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
if (self % norder == 0) then
    zwp = min(z, 0.0_wp)
else if (self % norder < 2) then
    zwp = z
else
    ! Wheeler stretching
    elev = self % elev(x, y)
    zwp = z - elev
    if (self % d > 0.0_wp) then
        zwp = zwp / (1.0_wp + elev / self % d)
    end if
end if
phi_xtswd = 0.0_wp
phi_ytswd = 0.0_wp
phi_ztswd = 0.0_wp
if (self % d < 0.0_wp) then
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &
                      ct_cur => self % wd(j) % ct_cur)
                Zfun = exp(kw * zwp)
                Zfun_z = kw * Zfun  ! Zhfun = Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Ctfun = ct_cur * Efun
                phi_xtswd = phi_xtswd + kwx * Ctfun % im * Zfun
                phi_ytswd = phi_ytswd + kwy * Ctfun % im * Zfun
                phi_ztswd = phi_ztswd + Ctfun % re * Zfun_z
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,     &
                      kwx => self % wd(j) % kwx,   &
                      kwy => self % wd(j) % kwy,   &  
                      ct_cur => self % wd(j) % ct_cur)
                Zfun = 1.0_wp + kw * zwp
                Zfun_z = kw * Zfun  ! Zhfun = Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Ctfun = ct_cur * Efun
                phi_xtswd = phi_xtswd + kwx * Ctfun % im * Zfun
                phi_ytswd = phi_ytswd + kwy * Ctfun % im * Zfun
                phi_ztswd = phi_ztswd + Ctfun % re * Zfun_z
            end associate
        end do
    end if
else
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,              &
                      kwx => self % wd(j) % kwx,            &
                      kwy => self % wd(j) % kwy,            &  
                      bfac => self % wd(j) % bfac,          &
                      fracbfac => self % wd(j) % fracbfac,  &
                      ct_cur => self % wd(j) % ct_cur)
                a = exp(kw * zwp)
                c = bfac / a**2
                Zfun = a * (1.0_wp + c) * fracbfac
                Zfun_z = kw * a * (1.0_wp - c) * fracbfac
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Ctfun = ct_cur * Efun
                phi_xtswd = phi_xtswd + kwx * Ctfun % im * Zfun
                phi_ytswd = phi_ytswd + kwy * Ctfun % im * Zfun
                phi_ztswd = phi_ztswd + Ctfun % re * Zfun_z
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,       &
                      kwx => self % wd(j) % kwx,     &
                      kwy => self % wd(j) % kwy,     &  
                      ct_cur => self % wd(j) % ct_cur, &
                      tanhkd => self % wd(j) % tanhkd)
                Zfun = 1.0_wp + tanhkd * kw * zwp
                Zhfun = tanhkd + kw * zwp
                Zfun_z = kw * Zhfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Ctfun = ct_cur * Efun
                phi_xtswd = phi_xtswd + kwx * Ctfun % im * Zfun
                phi_ytswd = phi_ytswd + kwy * Ctfun % im * Zfun
                phi_ztswd = phi_ztswd + Ctfun % re * Zfun_z
            end associate
        end do
    end if
end if
res(1) = phi_xtswd * self % cbeta - phi_ytswd * self % sbeta
res(2) = phi_xtswd * self % sbeta + phi_ytswd * self % cbeta
res(3) = phi_ztswd
!
end function acc_euler

!==============================================================================

function acc_particle(self, x, y, z) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
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
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y,z  ! Position application program
real(knd)             :: res    ! Fully nonlinear pressure
!
integer :: j
real(wp) :: xswd, yswd, Zfun, Zhfun, Zfun_z, elev, a, c
real(wp) :: phi_t, phi_xswd, phi_yswd, phi_zswd, zwp
complex(wp) :: Efun, Cfun, Ctfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
if (self % norder == 0) then
    zwp = min(z, 0.0_wp)
else if (self % norder < 2) then
    zwp = z
else
    ! Wheeler stretching
    elev = self % elev(x, y)
    zwp = z - elev
    if (self % d > 0.0_wp) then
        zwp = zwp / (1.0_wp + elev / self % d)
    end if
end if
phi_xswd = 0.0_wp
phi_yswd = 0.0_wp
phi_zswd = 0.0_wp
phi_t = 0.0_wp
if (self % d < 0.0_wp) then
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,       &
                      kwx => self % wd(j) % kwx,     &
                      kwy => self % wd(j) % kwy,     &
                      c_cur => self % wd(j) % c_cur, &
                      ct_cur => self % wd(j) % ct_cur)
                Zfun = exp(kw * zwp)
                Zfun_z = kw * Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_t = phi_t + real(ct_cur * Efun) * Zfun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
       end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,       &
                      kwx => self % wd(j) % kwx,     &
                      kwy => self % wd(j) % kwy,     &
                      c_cur => self % wd(j) % c_cur, &
                      ct_cur => self % wd(j) % ct_cur)
                Zfun = 1.0_wp + kw * zwp
                Zfun_z = kw * Zfun  ! Zhfun = Zfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_t = phi_t + real(ct_cur * Efun) * Zfun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
        end do
    end if
else
    if (zwp <= 0.0_wp .or. self % norder /= 1) then
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,              &
                      kwx => self % wd(j) % kwx,            &
                      kwy => self % wd(j) % kwy,            &  
                      bfac => self % wd(j) % bfac,          &
                      fracbfac => self % wd(j) % fracbfac,  &
                      c_cur => self % wd(j) % c_cur,        &
                      ct_cur => self % wd(j) % ct_cur)
                a = exp(kw * zwp)
                c = bfac / a**2
                Zfun = a * (1.0_wp + c) * fracbfac
                Zfun_z = kw * a * (1.0_wp - c) * fracbfac
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_t = phi_t + real(ct_cur * Efun) * Zfun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
        end do
    else
        do j = 1, self % n
            associate(kw => self % wd(j) % kw,         &
                      kwx => self % wd(j) % kwx,       &
                      kwy => self % wd(j) % kwy,       &
                      c_cur => self % wd(j) % c_cur,   &
                      ct_cur => self % wd(j) % ct_cur, &
                      tanhkd => self % wd(j) % tanhkd)
                Zfun = 1.0_wp + tanhkd * kw * zwp
                Zhfun = tanhkd + kw * zwp
                Zfun_z = kw * Zhfun
                Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
                Cfun = c_cur * Efun
                phi_t = phi_t + real(ct_cur * Efun) * Zfun
                phi_xswd = phi_xswd + kwx * Cfun % im * Zfun
                phi_yswd = phi_yswd + kwy * Cfun % im * Zfun
                phi_zswd = phi_zswd + Cfun % re * Zfun_z
            end associate
        end do
    end if
end if
! Note that hydrostatic component is always evaluated at physical location:
res = (-phi_t - (phi_xswd**2 + phi_yswd**2 + phi_zswd**2) * 0.5_wp - &
       z * self % grav) * self % rho
!
end function pressure

!==============================================================================

subroutine convergence(self, x, y, z, csv)
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self   ! Actual class
real(knd),        intent(in) :: x,y,z  ! Position application program
character(len=*), intent(in) :: csv    ! New output file
!
print*, 'spectral_wave_data::CONVERGENCE method is not relevant for shp=6.'
!
end subroutine convergence

!==============================================================================

subroutine strip(self, tmin, tmax, file_swd)
! Store part of SWD file into a new SWD file
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self   ! Actual class
real(knd),        intent(in) :: tmin, tmax ! Time window to store
character(len=*), intent(in) :: file_swd   ! Name of new swd file
!
print*, 'spectral_wave_data::STRIP method is not relevant for shp=6.'
!
end subroutine strip

!==============================================================================

function elev(self, x, y) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self ! Actual class
real(knd), intent(in) :: x,y  ! Position application program
real(knd)             :: res  ! Surface elevation at (x,y)
!
integer :: j
real(wp) :: xswd, yswd, reswp
complex(wp) :: Efun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
reswp = 0.0_wp
do j = 1, self % n
    associate(kwx => self % wd(j) % kwx,  &
              kwy => self % wd(j) % kwy,  &
              h_cur => self % wd(j) % h_cur)
        Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
        reswp = reswp + real(h_cur * Efun)
    end associate
end do
res = reswp
!
end function elev

!==============================================================================

function elev_t(self, x, y) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self ! Actual class
real(knd), intent(in) :: x,y  ! Position application program
real(knd)             :: res  ! d/dt of surface elevation at (x,y)
!
integer :: j
real(wp) :: xswd, yswd, reswp
complex(wp) :: Efun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
reswp = 0.0_wp
do j = 1, self % n
    associate(kwx => self % wd(j) % kwx,  &
              kwy => self % wd(j) % kwy,  &
              ht_cur => self % wd(j) % ht_cur)
        Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
        reswp = reswp + real(ht_cur * Efun)
    end associate
end do
res = reswp
!
end function elev_t

!==============================================================================

function grad_elev(self, x, y) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y    ! Position application program
real(knd)             :: res(3) ! x, y and z gradients of surface elevation at (x,y)
!
integer :: j
real(wp) :: xswd, yswd, elev_xswd, elev_yswd, HfunIm
complex(wp) :: Efun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
elev_xswd = 0.0_wp
elev_yswd = 0.0_wp
do j = 1, self % n
    associate(kwx => self % wd(j) % kwx,  &
              kwy => self % wd(j) % kwy,  &
              h_cur => self % wd(j) % h_cur)
        Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
        HfunIm = aimag(h_cur * Efun)
        elev_xswd = elev_xswd + kwx * HfunIm
        elev_yswd = elev_yswd + kwy * HfunIm
    end associate
end do
res(1) = elev_xswd * self % cbeta - elev_yswd * self % sbeta
res(2) = elev_xswd * self % sbeta + elev_yswd * self % cbeta
res(3) = 0.0_knd
!
end function grad_elev

!==============================================================================

function grad_elev_2nd(self, x, y) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
real(knd), intent(in) :: x,y    ! Position application program
real(knd)             :: res(3) ! Second order gradients of surface elevation
                                ! res(1) = d^2(elevation) / dx^2
                                ! res(2) = d^2(elevation) / dx dy
                                ! res(3) = d^2(elevation) / dy dy
!
integer :: j
real(wp) :: xswd, yswd, elev_xx_swd, elev_xy_swd, elev_yy_swd
real(wp) :: HfunRe, cc, cs, ss
complex(wp) :: Efun, Hfun
!
xswd = self % x0 + x * self % cbeta + y * self % sbeta
yswd = self % y0 - x * self % sbeta + y * self % cbeta
elev_xx_swd = 0.0_wp
elev_xy_swd = 0.0_wp
elev_yy_swd = 0.0_wp
do j = 1, self % n
    associate(kwx => self % wd(j) % kwx,  &
              kwy => self % wd(j) % kwy,  &
              h_cur => self % wd(j) % h_cur)
        Efun = exp(cmplx(0.0_wp, -kwx * xswd - kwy * yswd, kind=wp))
        HfunRe = real(h_cur * Efun)
        elev_xx_swd = elev_xx_swd - kwx * kwx * HfunRe
        elev_xy_swd = elev_xy_swd - kwx * kwy * HfunRe
        elev_yy_swd = elev_yy_swd - kwy * kwy * HfunRe
    end associate
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
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self ! Actual class
real(knd), intent(in) :: x,y  ! Position application program
real(knd)             :: res  ! Local water depth
!
res = self % d
!
end function bathymetry

!==============================================================================

function bathymetry_nvec(self, x, y) result(res)
class(spectral_wave_data_shape_6_impl_1), intent(in) :: self   ! Actual class
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
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of int parameter
integer                      :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_6_impl_1::get_int'
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
    res = count(abs(self % wd(:) % kwy) < 0.000001)
case('ny')
    res = count(abs(self % wd(:) % kwx) < 0.000001)
case('n')
    res = self % n
case('ipol')
    res = -1
case('impl')
    res = 1
case('nsumx')
    res = -1
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
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of logical parameter
logical                      :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_6_impl_1::get_logical'
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
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self   ! Actual class
character(len=*), intent(in) :: name ! Name of float parameter
real(knd)                    :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_6_impl_1::get_real'
character(len=250) :: err_msg(1)
integer :: i
real(wp) :: kmin, kmax
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
    res = -1.0_wp
case('dkx')
    res = -1.0_wp
case('dky')
    res = -1.0_wp
case('d')
    res = self % d
! Some special requests....
case('tmax')
    res = self % tmax
case('lmin')
    kmax = self % wd(1) % kw
    do i = 2, self % n
        kmax = max(kmax, self % wd(i) % kw)
    end do
    res = 2 * pi / kmax
case('lmax')
    kmin = self % wd(1) % kw
    do i = 2, self % n
        kmin = min(kmin, self % wd(i) % kw)
    end do
    res = 2 * pi / kmin
case('sizex')
    res = -1.0_wp
case('sizey')
    res = -1.0_wp
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
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self   ! Actual class
character(len=*),  intent(in) :: name ! Name of char parameter
character(len=:), allocatable :: res  ! Value of name parameter
!
character(len=*), parameter :: err_proc = 'spectral_wave_data_shape_6_impl_1::get_chr'
character(len=250) :: err_msg(1)
!
select case(name)
case('file', 'file_swd')
    res = self % file
case('version')
    res = version
case('class')
    res = 'spectral_wave_data_shape_6_impl_1'
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
class(spectral_wave_data_shape_6_impl_1), intent(inout) :: self ! Actual class
integer, optional, intent(in) :: nx_fft_in, ny_fft_in
real(knd), allocatable :: elev(:, :)

allocate(elev(nx_fft_in, ny_fft_in))
elev = 0.0_knd

end function elev_fft

!==============================================================================

end module spectral_wave_data_shape_6_impl_1_def
