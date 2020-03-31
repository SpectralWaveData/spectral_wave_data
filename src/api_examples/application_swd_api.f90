program application_swd_api

use, intrinsic :: iso_fortran_env, only: output_unit

use kind_values, only: wp => kind_swd_interface

use spectral_wave_data_def, only: spectral_wave_data
use spectral_wave_data_allocate_def, only: spectral_wave_data_allocate

implicit none
    
integer :: norder, impl, ipol, nsumx, nsumy, nx, ny, shp
real(wp) ::  x0, y0, t0, t, dt, tend, beta, x, y, z, rho
logical :: dc_bias
character(len=200) :: file_swd
character(len=:), allocatable :: cid

! Polymorphic object (type unknown at compile time) It can deal with any shape classes
class(spectral_wave_data), allocatable :: swd

! First we define some constructor parameters...

! Name of actual swd file
file_swd = "stokesTest_shallowWater.swd"

! Relate the application and SWD coordinate systems
x0 = 0.0_wp;  y0 = 0.0_wp;  t0 = 0.0_wp;  beta = 0.0_wp

! Optional: impl=0 means apply recommended implementation based on header of swd file.
impl = 0

! Optional: Density of water needed for pressure calculations
rho = 1025.0_wp

! Optional: number of spectral components to apply in calculations.
nsumx = -1;  nsumy = -1   ! Negative: apply all components from swd-file.

! Optional: Select temporal interpolation scheme
ipol = 0   ! C^2 continous

! Optional: For z>0 apply same series expansion as in wave generator
norder = 0

! Optional: Control zero-frequency signals from the wave generator
dc_bias = .false.   ! Do not include zero-frequency signals

!---------------------------------------------------------------------------
! Allocate and construct actual type of swd object...
!---------------------------------------------------------------------------
call spectral_wave_data_allocate(swd, file_swd, x0, y0, t0, beta, rho=rho,  &
                                 nsumx=nsumx, nsumy=nsumy, impl=impl,       &
                                 ipol=ipol, norder=norder, dc_bias=dc_bias)
if (swd % error % raised()) then
    print*, swd % error % get_msg()    ! Stop if any errors occured.
    stop
end if
!---------------------------------------------------------------------------
! Time domain simulation....
!---------------------------------------------------------------------------
tend = 5.0_wp;  dt = 1.0_wp;  t = 0.0_wp
do
    if (t > tend) exit
    
    ! Tell the swd object current application time...
    call swd % update_time(t)
    if (swd % error % raised()) then
        print*, swd % error % get_msg()    ! Stop if t is out of bounds
        stop
    end if
   
    ! Application coordinate where we need kinematics...
    x=4.3_wp;  y=5.4_wp;  z=-1.7_wp

    ! Complete set of kinematic functions....
    
    !  Velocity potential
    print*, 'phi = ', swd % phi(x,y,z)

    !  Stream function
    print*, 'varphi = ', swd % stream(x,y,z)

    !  partial(phi) / partial t
    print*, 'phi_t = ', swd % phi_t(x,y,z)
    
    ! particle velocity
    print*, 'grad_phi = ', swd % grad_phi(x,y,z)
    
    ! 2nd order gradients of phi
    print*, 'phi_xx, phi_xy, phi_xz, phi_yy, phi_yz, phi_zz = ', swd % grad_phi_2nd(x,y,z)
    
    ! Local (Euler) acceleration
    print*, 'acc_euler = ', swd % acc_euler(x,y,z)

    ! Particle acceleration
    print*, 'acc_particle = ', swd % acc_particle(x,y,z)

    ! Surface elevation
    print*, 'elev = ', swd % elev(x,y)

    ! local time derivative of wave elevation
    print*, 'elev_t = ', swd % elev_t(x,y)
    
    ! spatial gradient of wave elevation
    print*, 'grad_elev = ', swd % grad_elev(x,y)

    ! 2nd order gradients of wave elevation
    print*, 'elev_xx, elev_xy, elev_yy = ', swd % grad_elev_2nd(x,y)
 
    ! Total pressure
    print*, 'pressure = ', swd % pressure(x,y,z)

    ! Vertical distance from z=0 to sea floor (<0 if infinite)
    print*, 'local depth = ', swd % bathymetry(x,y)

    ! Unit normal vector of sea floor into the ocean at (x,y)
    print*, 'nx, ny, nz = ', swd % bathymetry_nvec(x,y)

    ! For a specific (x,y,z) at this t, return a CSV file on how particle velocity,
    ! elevation and pressure converge as a function of number of spectral components
    call swd % convergence(x, y, z, 'dump.csv')
    
    t = t + dt
end do

! To save storage for an interesting event you may create a new SWD file
! containing only the time steps within a specified time window.
! In this case we only care about the time interval [100.0, 200.0]
call swd % strip(tmin=100.0_wp, tmax=200.0_wp, file_swd='stripped.swd')

!---------------------------------------------------------------------------    
! There are 3 methods returning a value from the swd file...
!   swd.get_chr(swd, name)
!   swd.get_int(swd, name)
!   swd.get_real(swd, name)
! where name is the requested parameter in the swd file.
! Three examples are given below...
!---------------------------------------------------------------------------

! Extract the cid string from the SWD file
! (contains typical the content of the input file applied in the wave generator)
cid = swd % get_chr("cid")

! The shp parameter from the swd file
shp = swd % get_int("shp")

! Time step in SWD file
dt = swd % get_real("dt")

! Name of actual implementation class
print*, "cls_name = ", swd % get_chr('class')

! Close SWD file and release memory
call swd % close
!
end program application_swd_api
