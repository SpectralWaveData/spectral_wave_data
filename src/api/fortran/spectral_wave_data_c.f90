module spectral_wave_data_c_def

use, intrinsic :: iso_fortran_env, only: output_unit
use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_double, c_ptr, c_loc, &
                                         c_f_pointer, c_null_char, c_bool

use kind_values, only: wp => kind_swd_interface, c_wp => kind_swd_c

use spectral_wave_data_def,           only: spectral_wave_data
use spectral_wave_data_allocate_def,  only: spectral_wave_data_allocate

implicit none
private

! This module provides a standard C-compatible interface to the
! fortran spectral_wave_data class.
!
! Written by Jens Bloch Helmers, July, 20. 2019
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: spectral_wave_data_c  ! Target for C-pointer. Fortran SWD object is in target.
!
! C-compatible structs
!
public :: vector        ! Vector [x,y,z]
public :: vector_phi_2nd ! Vector [xx,xy,xz,yy,yz,zz]
public :: vector_elev_2nd ! Vector [xx,xy,yy]
!
! Public fortran functions with pure C-interface  
! All C-compatible names are defined by bind(c, name=...) in the function definitions
!
public :: constructor   ! Construct a C-pointer to a new specific SWD object
public :: update_time   ! Set current user time
public :: phi           ! Wave potential
public :: stream        ! Stream function
public :: phi_t         ! d(potential)/dt (Euler derivative)
public :: grad_phi      ! Particle velocity
public :: grad_phi_2nd  ! 2nd order gradients of potential
public :: acc_euler     ! Euler acceleration (grad(phi_t))
public :: acc_particle  ! Particle acceleration
public :: elev          ! Surface elevation
public :: elev_t        ! d(surface elevation)/dt (Euler derivative)
public :: grad_elev     ! Gradient of surface elevation
public :: grad_elev_2nd ! 2nd order gradients of elevation
public :: pressure      ! Fully nonlinear Bernoulli pressure
public :: bathymetry    ! Vertical distance from z=0 to sea floor (<0 if infinite)
public :: bathymetry_nvec ! Unit normal vector of sea floor into the ocean
public :: convergence   ! For a specific (t,x,y,z) return a csv-file on how velocity, elevation
                        ! and pressure converge as a function of number of spectral components
public :: strip         ! Create an new SWD file based on a time window of current SWD file.
public :: get_int       ! Extract a specified int parameter from the swd file.
public :: get_bool      ! Extract a specified bool parameter from the swd file.
public :: get_real      ! Extract a specified real parameter from the swd file.
public :: get_chr       ! Extract a specified char parameter from the swd file.
public :: error_raised  ! Return true if an error has been raised, otherwise false
public :: error_get_id  ! Return error id
public :: error_get_msg ! Return error message
public :: error_clear   ! Clear the error flag
public :: close         ! Destructor
public :: elev_fft      ! Surface elevation FFT
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type :: spectral_wave_data_c
    class(spectral_wave_data), allocatable :: obj
end type spectral_wave_data_c

type, bind(c) :: vector
    real(c_wp) :: x, y, z
end type vector

type, bind(c) :: vector_phi_2nd
    real(c_wp) :: xx, xy, xz, yy, yz, zz
end type vector_phi_2nd

type, bind(c) :: vector_elev_2nd
    real(c_wp) :: xx, xy, yy
end type vector_elev_2nd

contains

!==============================================================================

function constructor(file_swd, x0, y0, t0, beta, rho, nsumx, nsumy, impl, &
                     ipol, norder, dc_bias) bind(c, name='swd_api_allocate')
!------------------------------------------------------------------------------
! Provides a C pointer to an object containing the actual Fortran SWD object.
!------------------------------------------------------------------------------
type(c_ptr) :: constructor  ! Pointer to construct

! Name of swd file
character(c_char), intent(in) :: file_swd(*)

! Relation between SWD and application coordinates (beta is in degree)
real(c_wp), value, intent(in) :: x0, y0, t0, beta

! Density of water (applied for pressure calculations)
real(c_wp), value, intent(in) :: rho      

! Number of applied spectral comp. in x and y-dir. 
! If negative: apply all spectral components from the SWD file 
integer(c_int), value, intent(in) :: nsumx
integer(c_int), value, intent(in) :: nsumy

! Index to determine actual derived class
!    0 = Default based on content of SWD file 
!   <0 = In-house and experimental implementations
!   >0 = Validated implementations available open software
integer(c_int), value, intent(in) :: impl     

! Index to request actual temporal interpolation scheme
!   0 = C^2 continous scheme (default)
!   1 = C^1 continous
!   2 = C^3 continous
integer(c_int), value, intent(in) :: ipol

! Expansion order to apply in kinematics for z>0
!   0  = Apply expansion order specified on swd file (default)
!   <0 = Apply exp(kj z)
!   >0 = Apply expansion order = norder 
integer(c_int), value, intent(in) :: norder

! Control application of zero-frequency bias present in SWD file
!  False = Suppress contribution from zero frequency amplitudes (default)
!  True  = Apply zero frequency amplitudes from SWD file. 
logical(c_bool), value, intent(in):: dc_bias
!------------------------------------------------------------------------------

integer :: impl_f, nsumx_f, nsumy_f, ipol_f, norder_f
logical :: dc_bias_f
type(spectral_wave_data_c), pointer :: swdc
real(wp) :: x0_f, y0_f, t0_f, beta_f, rho_f
character(len=:), allocatable :: file_swd_f

! Convert to fortran types
file_swd_f = str_c2f(file_swd)
x0_f = x0
y0_f = y0
t0_f = t0
beta_f = beta
rho_f = rho
nsumx_f = nsumx
nsumy_f = nsumy
impl_f = impl
ipol_f = ipol
norder_f = norder
dc_bias_f = dc_bias
!
allocate(swdc)
call spectral_wave_data_allocate(swdc % obj, file_swd_f, x0_f, y0_f,    &
                                 t0_f, beta_f, rho_f, nsumx_f, nsumy_f, &
                                 impl_f, ipol_f, norder_f, dc_bias_f)

constructor = c_loc(swdc)
!
end function constructor

!==============================================================================

subroutine update_time(this, time) bind(c, name='swd_api_update_time')
type(c_ptr), value             :: this  ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: time  ! Current application time
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: time_f
time_f = time
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
call swd % obj % update_time(time_f)
!
end subroutine update_time

!==============================================================================

function phi(this, x, y, z) bind(c, name='swd_api_phi')
type(c_ptr), value             :: this  ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y,z ! Position application program
real(c_wp)                     :: phi   ! Potential at (x,y,z)
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
phi = swd % obj % phi(x_f, y_f, z_f)
!
end function phi

!==============================================================================

function stream(this, x, y, z) bind(c, name='swd_api_stream')
type(c_ptr), value             :: this   ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y,z  ! Position application program
real(c_wp)                     :: stream ! Stream function at at (x,y,z)
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
stream = swd % obj % stream(x_f, y_f, z_f)
!
end function stream

!==============================================================================

function phi_t(this, x, y, z) bind(c, name='swd_api_phi_t')
type(c_ptr), value             :: this  ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y,z ! Position application program
real(c_wp)                     :: phi_t ! d(potential)/dt (Euler) at (x,y,z)
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
phi_t = swd % obj % phi_t(x_f, y_f, z_f)
!
end function phi_t

!==============================================================================

function grad_phi(this, x, y, z) bind(c, name='swd_api_grad_phi')
type(c_ptr), value             :: this     ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y,z    ! Position application program
type(vector)                   :: grad_phi ! Particle velocity at (x,y,z)
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f, vec(3)
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
vec = swd % obj % grad_phi(x_f, y_f, z_f)
grad_phi = vector(vec(1), vec(2), vec(3))
!
end function grad_phi

!==============================================================================

function grad_phi_2nd(this, x, y, z) bind(c, name='swd_api_grad_phi_2nd')
type(c_ptr), value             :: this     ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y,z    ! Position application program
type(vector_phi_2nd)           :: grad_phi_2nd ! [xx,xy,xz,yy,yz,zz] at (x,y,z)
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f, vec(6)
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
vec = swd % obj % grad_phi_2nd(x_f, y_f, z_f)
grad_phi_2nd = vector_phi_2nd(vec(1), vec(2), vec(3), vec(4), vec(5), vec(6))
!
end function grad_phi_2nd

!==============================================================================

function acc_euler(this, x, y, z) bind(c, name='swd_api_acc_euler')
!--------------------------------------------------------------------------    
! Actual spectral_wave_data_c object
type(c_ptr), value :: this
    
! Position application program
real(c_wp), value, intent(in) :: x,y,z
    
! Euler acceleration (grad(phi_t)) at (x,y,z) for current time)
type(vector) :: acc_euler
!--------------------------------------------------------------------------    
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f, acc(3)
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
acc = swd % obj % acc_euler(x_f, y_f, z_f)
acc_euler = vector(acc(1), acc(2), acc(3))
!
end function acc_euler

!==============================================================================

function acc_particle(this, x, y, z) bind(c, name='swd_api_acc_particle')
!--------------------------------------------------------------------------    
! Actual spectral_wave_data_c object
type(c_ptr), value :: this
    
! Position application program
real(c_wp), value, intent(in) :: x,y,z
    
! Euler acceleration (grad(phi_t)) at (x,y,z) for current time)
type(vector) :: acc_particle
!--------------------------------------------------------------------------    
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f, acc(3)
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
acc = swd % obj % acc_particle(x_f, y_f, z_f)
acc_particle = vector(acc(1), acc(2), acc(3))
!
end function acc_particle

!==============================================================================

function elev(this, x, y) bind(c, name='swd_api_elev')
type(c_ptr), value             :: this  ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y   ! Position application program
real(c_wp)                     :: elev  ! Wave elevation at (x,y)
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f
x_f = x
y_f = y
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
elev = swd % obj % elev(x_f, y_f)
!
end function elev

!==============================================================================

function elev_t(this, x, y) bind(c, name='swd_api_elev_t')
type(c_ptr), value             :: this   ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y    ! Position application program
real(c_wp)                     :: elev_t ! d(surface elevation)/dt (Euler derivative)
                                         ! at (x,y) for current time
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f
x_f = x
y_f = y
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
elev_t = swd % obj % elev_t(x_f, y_f)
!
end function elev_t

!==============================================================================

function grad_elev(this, x, y) bind(c, name='swd_api_grad_elev')
type(c_ptr), value             :: this      ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y       ! Position application program
type(vector)                   :: grad_elev ! Gradient of surface elevation at (x,y) for current time
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, grad(3)
x_f = x
y_f = y
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
grad = swd % obj % grad_elev(x_f, y_f)
grad_elev = vector(grad(1), grad(2), grad(3))
!
end function grad_elev

!==============================================================================

function grad_elev_2nd(this, x, y) bind(c, name='swd_api_grad_elev_2nd')
type(c_ptr), value             :: this      ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y       ! Position application program
type(vector_elev_2nd)          :: grad_elev_2nd ! Gradient of surface elevation at (x,y) for current time
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, grad(3)
x_f = x
y_f = y
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
grad = swd % obj % grad_elev_2nd(x_f, y_f)
grad_elev_2nd = vector_elev_2nd(grad(1), grad(2), grad(3))
!
end function grad_elev_2nd

!==============================================================================

function pressure(this, x, y, z) bind(c, name='swd_api_pressure')
type(c_ptr), value             :: this     ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y,z    ! Position application program
real(c_wp)                     :: pressure ! Fully nonlinear Bernoulli pressure
                                            ! at (x,y,z) for current time
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, z_f
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
pressure = swd % obj % pressure(x_f, y_f, z_f)
!
end function pressure

!==============================================================================

function bathymetry(this, x, y) bind(c, name='swd_api_bathymetry')
type(c_ptr), value             :: this       ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y        ! Position application program
real(c_wp)                     :: bathymetry ! Local water depth at
                                                ! (x,y) for current time
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f
x_f = x
y_f = y
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
bathymetry = swd % obj % bathymetry(x_f, y_f)
!
end function bathymetry

!==============================================================================

function bathymetry_nvec(this, x, y) bind(c, name='swd_api_bathymetry_nvec')
type(c_ptr), value :: this  ! Actual spectral_wave_data_c object
real(c_wp), value, intent(in) :: x,y   ! Position application program
type(vector) :: bathymetry_nvec ! Unit normal vector of sea floor into the ocean at(x, y)
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: x_f, y_f, res(3)
x_f = x
y_f = y
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
res = swd % obj % bathymetry_nvec(x_f, y_f)
bathymetry_nvec = vector(res(1), res(2), res(3))
!
end function bathymetry_nvec

!==============================================================================

subroutine convergence(this, x, y, z, csv) bind(c, name='swd_api_convergence')
!--------------------------------------------------------------------------
! For a specific(x, y, z) at this t, return a CSV file on how particle velocity,
! elevation and pressure converge as a function of number of spectral components
!--------------------------------------------------------------------------
type(c_ptr), value  :: this ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: x,y,z    ! Position application program
character(c_char),  intent(in) :: csv(*)   ! Name of output CSV file
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
character(len=:), allocatable :: csv_f
real(wp) :: x_f, y_f, z_f
x_f = x
y_f = y
z_f = z
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
csv_f = str_c2f(csv)  ! Convert to fortran string
!
call swd % obj % convergence(x_f, y_f, z_f, csv_f)
!
end subroutine convergence

!==============================================================================

subroutine strip(this, tmin, tmax, swd_file) bind(c, name='swd_api_strip')
!--------------------------------------------------------------------------
! Create a new SWD file containing only the time steps within 
! the range [tmin, tmax]
!--------------------------------------------------------------------------
type(c_ptr), value :: this ! Actual spectral_wave_data_c object
real(c_wp),  value, intent(in) :: tmin, tmax ! Time window
character(c_char),  intent(in) :: swd_file(*)   ! Name of new SWD file
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
real(wp) :: tmin_f, tmax_f
character(len=:), allocatable :: swd_file_f
!
swd_file_f = str_c2f(swd_file)  ! Convert to fortran string
tmin_f = tmin
tmax_f = tmax
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
call swd % obj % strip(tmin_f, tmax_f, swd_file_f)
!
end subroutine strip

!==============================================================================

function get_int(this, name) bind(c, name='swd_api_get_int')
!--------------------------------------------------------------------------
! Extract 'name' from the SWD file assuming it is an int
! Error signal raised if 'name' is not sound (error_get_id()==1004)
!--------------------------------------------------------------------------
type(c_ptr), value :: this ! Actual spectral_wave_data_c object
integer(c_int) :: get_int  ! The value of the requested parameter
character(c_char), intent(in) :: name(*)  ! Name of character variable
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
character(len=:), allocatable :: name_f
!
name_f = str_c2f(name)  ! Convert to fortran string
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
get_int = swd % obj % get_int(name_f)
!
end function get_int

!==============================================================================

function get_bool(this, name) bind(c, name='swd_api_get_bool')
!--------------------------------------------------------------------------
! Extract 'name' from the SWD file assuming it is a bool
! Error signal raised if 'name' is not sound (error_get_id()==1004)
!--------------------------------------------------------------------------
type(c_ptr), value :: this ! Actual spectral_wave_data_c object
logical(c_bool) :: get_bool ! The value of the requested parameter
character(c_char), intent(in) :: name(*)  ! Name of character variable
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
character(len=:), allocatable :: name_f
!
name_f = str_c2f(name)  ! Convert to fortran string
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
get_bool = swd % obj % get_logical(name_f)
!
end function get_bool

!==============================================================================

function get_real(this, name) bind(c, name='swd_api_get_real')
!--------------------------------------------------------------------------
! Extract 'name' from the SWD file assuming it is a float
! Error signal raised if 'name' is not sound (error_get_id()==1004)
!--------------------------------------------------------------------------
type(c_ptr), value :: this ! Actual spectral_wave_data_c object
real(c_wp) :: get_real  ! The value of the requested parameter
character(c_char), intent(in) :: name(*)  ! Name of character variable
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
character(len=:), allocatable :: name_f
!
name_f = str_c2f(name)  ! Convert to fortran string
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
get_real = swd % obj % get_real(name_f)
!
end function get_real

!==============================================================================

function get_chr(this, name) bind(c, name='swd_api_get_chr')
!--------------------------------------------------------------------------
! Extract 'name' from the SWD file assuming it is a character array
! Error signal raised if 'name' is not sound (error_get_id()==1004)
!--------------------------------------------------------------------------
type(c_ptr), value  :: this     ! Actual spectral_wave_data_c object
type(c_ptr)         :: get_chr  ! The value of the requested parameter
character(c_char),  intent(in) :: name(*)  ! Name of character variable
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
character(kind=c_char, len=:), allocatable, save, target :: value_c
character(len=:), allocatable :: name_f, value_f
!
name_f = str_c2f(name)  ! Convert to fortran string
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
value_f = swd % obj % get_chr(name_f)
value_c = trim(value_f) // c_null_char
get_chr = c_loc(value_c)
!
end function get_chr

!==============================================================================

function error_raised(this) bind(c, name='swd_api_error_raised')
!--------------------------------------------------------------------------
! Check if an error has been raised in the SWD object
!--------------------------------------------------------------------------
type(c_ptr), value :: this ! Actual spectral_wave_data_c object
logical(c_bool) :: error_raised  ! True if an error is raised (false if not)
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
!
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
error_raised = swd % obj % error % raised()
!
end function error_raised

!==============================================================================

function error_get_id(this) bind(c, name='swd_api_error_get_id')
!--------------------------------------------------------------------------
! Extract error id from SWD object
!--------------------------------------------------------------------------
type(c_ptr), value :: this ! Actual spectral_wave_data_c object
integer(c_int) :: error_get_id  ! The error id
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
!
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
error_get_id = swd % obj % error % get_id()
!
end function error_get_id

!==============================================================================

function error_get_msg(this) bind(c, name='swd_api_error_get_msg')
!--------------------------------------------------------------------------
! Extract error message from SWD object
!--------------------------------------------------------------------------
type(c_ptr), value  :: this     ! Actual spectral_wave_data_c object
type(c_ptr)         :: error_get_msg  ! The error message
!--------------------------------------------------------------------------
!
type(spectral_wave_data_c), pointer :: swd
character(kind=c_char, len=:), allocatable, save, target :: value_c
character(len=:), allocatable :: value_f
!
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
value_f = swd % obj % error % get_msg()
value_c = trim(value_f) // c_null_char
error_get_msg = c_loc(value_c)
!
end function error_get_msg

!==============================================================================

subroutine error_clear(this) bind(c, name='swd_api_error_clear')
type(c_ptr), value :: this  ! Actual spectral_wave_data_c object
!
type(spectral_wave_data_c), pointer :: swd
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
call swd % obj % error % clear()
!
end subroutine error_clear

!==============================================================================

subroutine close(this) bind(c, name='swd_api_close')
type(c_ptr), value :: this  ! Actual spectral_wave_data_c object to destroy
!
type(spectral_wave_data_c), pointer :: swd
!
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)
!
call swd % obj % close()
deallocate(swd % obj)
!
end subroutine close

!==============================================================================

function str_c2f(c_string) result (f_string)    ! Local helper function
character(c_char), intent(in) :: c_string(*)    ! Classical array of C-characters
character(len=:), allocatable :: f_string       ! Corresponding Fortran string
!
integer :: i, n
i = 0
do
    i = i + 1
    if (c_string(i)==c_null_char) exit
    if (i>1000000) then
        ! Seems that c_string is not null terminated!"...
        i = 20 ! Make it short and sweet to indicate bug...
    end if
end do
n = i - 1
allocate(character(len=n) :: f_string)
do i = 1, n
    f_string(i:i) = c_string(i)
end do
!
end function str_c2f

!==============================================================================

subroutine elev_fft(this, nx_fft, ny_fft, elev_arr) bind(c, name='swd_api_elev_fft_')
type(c_ptr), value              :: this        ! Actual spectral_wave_data_c object
integer(c_int), value           :: nx_fft, ny_fft       ! Dimensions of output-grid
real(c_wp), allocatable, intent(out) :: elev_arr(:, :)
type(spectral_wave_data_c), pointer :: swd

integer :: nx_fft_f, ny_fft_f
        
nx_fft_f = nx_fft
ny_fft_f = ny_fft
call c_f_pointer(this, swd) ! Find corresponding Fortran object (swd)

allocate(elev_arr(nx_fft, ny_fft))
elev_arr = swd % obj % elev_fft(nx_fft_f, ny_fft_f)

end subroutine elev_fft

!==============================================================================

end module spectral_wave_data_c_def
