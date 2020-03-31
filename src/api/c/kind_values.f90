module kind_values

use, intrinsic :: iso_fortran_env, only: real32, real64, real128
use, intrinsic :: iso_c_binding,   only: c_float, c_double

implicit none

! Written by Jens Bloch Helmers, 19. march 1995.

! This module defines some kind values.
!----------------------------------------------------------------------------

integer, private, parameter :: single = kind(0.0e0)
integer, private, parameter :: double = kind(0.0d0)

integer, public, parameter :: kind_swd_c = c_double        ! C interface
integer, public, parameter :: kind_swd_interface = real64  ! Fortran interface
integer, public, parameter :: kind_swd_internal = real64   ! Internal computations

end module kind_values
