module spectral_wave_data_def

use, intrinsic :: iso_fortran_env, only: int64

use kind_values, only: knd => kind_swd_interface, wp => kind_swd_internal

use spectral_wave_data_error, only: swd_error
use swd_fft_def, only: swd_fft

implicit none
private

! This module provides the abstract base class for spectral_wave_data_X objects.
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
public :: spectral_wave_data
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type, abstract :: spectral_wave_data
    ! Common attributes for all shape classes
    character(len=30)  :: prog ! Name of the program who created this swd file including version.
    character(len=20)  :: date ! Date and time this swd file was created
    character(len=200) :: file ! Name of swd file
    integer            :: unit ! Unit number associated with swd file
    integer            :: fmt  ! Code to identify format of swd file.
    integer            :: shp  ! Index of actual spectral shape class
    integer            :: amp  ! Index of which spectral amplitudes are available
    character(len=:), allocatable :: cid  ! Identification text in swd file
    integer            :: nstrip ! Number of initial time steps removed from original simulation
    integer            :: nsteps ! Total number of time steps in swd file.
    integer            :: order  ! Order of perturbation (<0 if fully nonlinear) applied in wave generator
    integer            :: norder  ! Expansion order to apply in kinematics for z>0
                                  ! <0: apply exp(kj z)
                                  ! 0:  apply expansion order specified on swd file
                                  ! >0: apply expansion order = norder 
    integer            :: ipol ! Index defining the temporal interpolation scheme
    real(wp)           :: dt   ! Constant time step in swd file
    real(wp)           :: t0   ! Input seed for time (t0>=0)
    real(wp)           :: x0,y0 ! Input seed for spatial location
    real(wp)           :: tswd ! Current swd time
    real(wp)           :: grav ! Acceleration of gravity
    real(wp)           :: lscale ! Number of length units in wave generator per meter.
    real(wp)           :: rho  ! Density of water
    real(wp)           :: cbeta ! cos(beta), beta=angle between swd and application x-axis
    real(wp)           :: sbeta ! sin(beta), beta=angle between swd and application x-axis
    real(wp)           :: tmax  ! Maximum allowed simulation time (user system)
    integer            :: size_complex ! On most systems size_complex=8 for c_float based numbers
    integer            :: size_step ! Fortran storage size per time step
    integer(int64)     :: ipos0 ! File postion where temporal functions starts
    logical            :: eof  ! End-of-file detected for SWD file
    logical            :: dc_bias ! True: apply zero frequency amplitudes from SWD file. 
                                  ! False: Suppress contribution from zero frequency amplitudes (Default)
    type(swd_error)    :: error   ! Abort free error handler
    type(swd_fft)      :: fft     ! FFT-based evaluations
contains
    procedure(update_time), deferred :: update_time       ! Obtain spectral data for current time
    procedure(phi),         deferred :: phi               ! Calculate potential at location for current time
    procedure(stream),      deferred :: stream            ! Calculate stream function
    procedure(phi_t),       deferred :: phi_t             ! Calculate d(potential)/dt (Euler) at location for current time
    procedure(grad_phi),    deferred :: grad_phi          ! Calculate particle velocity at location for current time
    procedure(grad_phi_2nd),deferred :: grad_phi_2nd      ! Calculate second order spatial gradients of potential
    procedure(acc_euler),   deferred :: acc_euler         ! Calculate Euler acceleration (grad(phi_t)) at location for current time
    procedure(acc_particle),deferred :: acc_particle      ! Calculate particle acceleration at location for current time
    procedure(elev),        deferred :: elev              ! Calculate surface elevation at location for current time
    procedure(elev_t),      deferred :: elev_t            ! Calculate d(surface elevation)/dt (Euler) at location for current time
    procedure(grad_elev),   deferred :: grad_elev         ! Calculate gradient of surface elevation at location for current time
    procedure(grad_elev_2nd),deferred:: grad_elev_2nd     ! Calculate second order spatial gradients of elevation
    procedure(pressure),    deferred :: pressure          ! Fully nonlinear Bernoulli pressure
    procedure(bathymetry),  deferred :: bathymetry        ! Return local depth at application position (x, y)
    procedure(bathymetry_nvec),deferred :: bathymetry_nvec ! Unit normal vector of sea floor into the ocean at (x,y)
    procedure(convergence), deferred :: convergence       ! For a specific (t,x,y,z) return a csv-file on how particle velocity, elevation
                                                          ! and pressure converge as a function of number of spectral components
    procedure(strip),       deferred :: strip             ! Create a new SWD file based on a time window of current SWD file.
    procedure(get_int),     deferred :: get_int           ! Extract a specified int parameter
    procedure(get_logical), deferred :: get_logical       ! Extract a specified logical parameter
    procedure(get_real),    deferred :: get_real          ! Extract a specified float parameter
    procedure(get_chr),     deferred :: get_chr           ! Extract a specified char parameter
    procedure(close),       deferred :: close             ! Manual destructor
    procedure(elev_fft),    deferred :: elev_fft          ! Surface elevation on a regular grid using FFT 
    procedure(grad_phi_fft),deferred :: grad_phi_fft      ! Grad phi on a regular grid using FFT 
    procedure :: x_fft                                    ! x-vector corresponding to output from FFT routines
    procedure :: y_fft                                    ! y-vector corresponding to output from FFT routines
    procedure :: error_raised                             ! Return .true. if error has been signaled
    procedure :: error_id                                 ! Return error id
    procedure :: error_msg                                ! Return error message
    procedure :: error_clear                              ! Clear error signal (id=0)
end type spectral_wave_data

abstract interface
    !
    subroutine update_time(self, time)
        import
        class(spectral_wave_data), intent(inout) :: self  ! Update data in memory (if needed)
        real(knd),                 intent(in)    :: time  ! Current time in simulation program
    end subroutine update_time

    function phi(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self  ! Actual class
        real(knd),                 intent(in) :: x,y,z ! Position application program
        real(knd)                             :: res   ! Potential at (x,y,z)
    end function phi
    
    function stream(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self  ! Actual class
        real(knd),                 intent(in) :: x,y,z ! Position application program
        real(knd)                             :: res   ! Stream function value at (x,y,z)
    end function stream

    function phi_t(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self  ! Actual class
        real(knd),                 intent(in) :: x,y,z ! Position application program
        real(knd)                             :: res   ! Euler time derivative of potential at (x,y,z)
    end function phi_t

    function grad_phi(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y,z  ! Position application program
        real(knd)                             :: res(3) ! Particle velocity at (x,y,z)
    end function grad_phi

    function grad_phi_2nd(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y,z  ! Position application program
        real(knd)                             :: res(6) ! Second order gradients of potential at (x,y,z)
                                                        ! res(1) = d^2(potential) / dx^2
                                                        ! res(2) = d^2(potential) / dx dy
                                                        ! res(3) = d^2(potential) / dx dz
                                                        ! res(4) = d^2(potential) / dy^2
                                                        ! res(5) = d^2(potential) / dy dz
                                                        ! res(6) = d^2(potential) / dz^2
    end function grad_phi_2nd

    function acc_euler(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y,z  ! Position application program
        real(knd)                             :: res(3) ! Euler acceleration at (x,y,z)
    end function acc_euler

    function acc_particle(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y,z  ! Position application program
        real(knd)                             :: res(3) ! Particle acceleration at (x,y,z)
    end function acc_particle

    function elev(self, x, y) result(res)
        import
        class(spectral_wave_data), intent(in) :: self ! Actual class
        real(knd),                 intent(in) :: x,y  ! Position application program
        real(knd)                             :: res  ! Surface elevation at (x,y)
    end function elev
    
    function elev_t(self, x, y) result(res)
        import
        class(spectral_wave_data), intent(in) :: self ! Actual class
        real(knd),                 intent(in) :: x,y  ! Position application program
        real(knd)                             :: res  ! d/dt of surface elevation at (x,y)
    end function elev_t

    function grad_elev(self, x, y) result(res)
        import
        class(spectral_wave_data), intent(in) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y    ! Position application program
        real(knd)                             :: res(3) ! x, y and z gradients of surface elevation at (x,y)
    end function grad_elev

    function grad_elev_2nd(self, x, y) result(res)
        import
        class(spectral_wave_data), intent(in) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y    ! Position application program
        real(knd)                             :: res(3) ! Second order gradients of surface elevation
                                                        ! res(1) = d^2(elevation) / dx^2
                                                        ! res(2) = d^2(elevation) / dx dy
                                                        ! res(3) = d^2(elevation) / dy^2
    end function grad_elev_2nd

    function pressure(self, x, y, z) result(res)
        import
        class(spectral_wave_data), intent(in) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y,z  ! Position application program
        real(knd)                             :: res    ! Nonlinear pressure
    end function pressure

    function bathymetry(self, x, y) result(res)
        import
        class(spectral_wave_data), intent(in) :: self ! Actual class
        real(knd),                 intent(in) :: x,y  ! Position application program
        real(knd)                             :: res  ! Local depth at (x,y)
    end function bathymetry
    
    function bathymetry_nvec(self, x, y) result(res)
        import
        class(spectral_wave_data), intent(in) :: self ! Actual class
        real(knd),                 intent(in) :: x,y  ! Position application program
        real(knd)                             :: res(3) ! Unit normal vector of sea floor into the ocean at (x,y)
    end function bathymetry_nvec

    subroutine convergence(self, x, y, z, csv)
        import
        class(spectral_wave_data), intent(inout) :: self   ! Actual class
        real(knd),                 intent(in) :: x,y,z  ! Position application program
        character(len=*),          intent(in) :: csv    ! Name of output csv-file
    end subroutine convergence

    subroutine strip(self, tmin, tmax, file_swd)
        ! Create a new swd file containing the spectral information limited
        ! to the application time window:  tmin <= t <= tmax.
        import
        class(spectral_wave_data), intent(inout) :: self       ! Actual class
        real(knd),                 intent(in) :: tmin, tmax ! Time window application program
        character(len=*),          intent(in) :: file_swd   ! Name of new swd file
    end subroutine strip

    function get_int(self, name) result(res)
        import
        class(spectral_wave_data), intent(inout) :: self ! Actual class
        character(len=*),          intent(in) :: name ! Name of int parameter
        integer                               :: res  ! Value of name parameter
    end function get_int

    function get_logical(self, name) result(res)
        import
        class(spectral_wave_data), intent(inout) :: self ! Actual class
        character(len=*),          intent(in) :: name ! Name of logical parameter
        logical                               :: res  ! Value of name parameter
    end function get_logical

    function get_real(self, name) result(res)
        import
        class(spectral_wave_data), intent(inout) :: self ! Actual class
        character(len=*), intent(in) :: name ! Name of real parameter
        real(knd)                    :: res  ! Value of name parameter
    end function get_real

    function get_chr(self, name) result(res)
        import
        class(spectral_wave_data), intent(inout) :: self ! Actual class
        character(len=*),          intent(in) :: name ! Name of char parameter
        character(len=:), allocatable         :: res  ! Value of name parameter
    end function get_chr

    subroutine close(self)
        import
        class(spectral_wave_data) :: self  ! Object to destruct
    end subroutine close

    function elev_fft(self, nx_fft_in, ny_fft_in) result(elev)
        import
        class(spectral_wave_data), intent(inout) :: self ! Actual class
        integer, optional, intent(in) :: nx_fft_in, ny_fft_in
        real(knd), allocatable :: elev(:, :)
    end function elev_fft

    function grad_phi_fft(self, z, nx_fft_in, ny_fft_in) result(grad_phi)
        import
        class(spectral_wave_data), intent(inout) :: self ! Actual class
        real(wp), intent(in) :: z
        integer, optional, intent(in) :: nx_fft_in, ny_fft_in
        real(knd), allocatable :: grad_phi(:, :, :)
    end function grad_phi_fft

end interface

contains

!==============================================================================

function error_raised(self) result(res)
class(spectral_wave_data), intent(in) :: self ! Error handler
logical                               :: res  ! .true. if error has been signaled
!
res = self % error % raised()
!
end function error_raised

!==============================================================================

function error_id(self) result(res)
class(spectral_wave_data), intent(in) :: self ! Error handler
integer                               :: res  ! Return error code
!
res = self % error % id
!
end function error_id

!==============================================================================

function error_msg(self) result(res)
class(spectral_wave_data), intent(in)   :: self ! Error handler
character(len=len_trim(self%error%msg)) :: res  ! Return error code
!
res = self % error % msg
!
end function error_msg

!==============================================================================

subroutine error_clear(self)
class(spectral_wave_data), intent(inout) :: self ! Error handler
!
call self % error % clear()
!
end subroutine error_clear

!==============================================================================

function x_fft(self, nx_fft_in)
class(spectral_wave_data), intent(inout) :: self ! Actual class
integer, optional, intent(in) :: nx_fft_in
real(knd), allocatable :: x_fft(:)
character(len=*), parameter :: err_proc = 'spectral_wave_data::x_fft'
character(len=:), allocatable :: err_msg(:)

x_fft = self % fft % x_fft(nx_fft_in)

if (self % fft % error % raised()) then
    err_msg = [self % fft % error % get_msg()]
    call self % error % set_id_msg(err_proc, &
                                   self % fft % error % get_id(), &
                                   err_msg)
end if

end function x_fft

!==============================================================================

function y_fft(self, ny_fft_in)
class(spectral_wave_data), intent(inout) :: self ! Actual class
integer, optional, intent(in) :: ny_fft_in
real(knd), allocatable :: y_fft(:)
character(len=*), parameter :: err_proc = 'spectral_wave_data::y_fft'
character(len=:), allocatable :: err_msg(:)

y_fft = self % fft % y_fft(ny_fft_in)

if (self % fft % error % raised()) then
    err_msg = [self % fft % error % get_msg()]
    call self % error % set_id_msg(err_proc, &
                                   self % fft % error % get_id(), &
                                   err_msg)
end if

end function y_fft

!==============================================================================
  
end module spectral_wave_data_def
