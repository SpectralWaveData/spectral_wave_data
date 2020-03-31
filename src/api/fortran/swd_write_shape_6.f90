module swd_write_shape_6_def

use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float, c_null_char

implicit none
private

! Please adjust wp to the working precision applied in the wave generator
integer, parameter :: wp = kind(1.0d0) ! (double precision)


! This module provides a class for writing SWD files for shape 6.
!
! Written by Jens Bloch Helmers, November, 17. 2019
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: swd_write_shape_6
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type :: swd_write_shape_6
    integer :: amp    ! Type of spectral amplitudes to be provided
    integer :: lu     ! Unit number associated with swd file
    integer :: nsteps ! Expected number of time steps
    integer :: it     ! Current time steps
    integer :: n      ! Number of spectral components in x-dir
    integer :: nh     ! Number of auxilary shapes (hat(n) in theory manual)
    integer :: nsf    ! Number of sea floor points
contains
    procedure :: build   ! Build the file in one go....
end type swd_write_shape_6

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

subroutine build(self, file, prog, cid, grav, lscale, depth, nwaves, awaves, kwaves, &
                 dwaves, pwaves)
class(swd_write_shape_6),  intent(out) :: self  ! Object to create
character(len=*),  intent(in)  :: file   ! Requested name for new SWD file
character(len=*),  intent(in)  :: prog   ! Name of actual wave generator including version number
character(len=*),  intent(in)  :: cid    ! Text describing actual input parameters to the
                                         ! wave generator. Hint: write(cid, nml=wave_input_namelist)
real(wp),          intent(in)  :: grav   ! Acceleration of gravity applied in wave generator
real(wp),          intent(in)  :: lscale ! Number of length units per meter applied in wave generator
                                         ! E.g. if US feet was applied: Lscale = 1/0.3048 = 3.2808  
real(wp),          intent(in)  :: depth  ! Constant water depth (<0 indicates infinite depth)
integer,           intent(in)  :: nwaves ! n as defined in the theory documentation
real(wp),          intent(in)  :: awaves(:) ! Single wave amplitude for each component
real(wp),          intent(in)  :: kwaves(:) ! Wave number for each component
real(wp),          intent(in)  :: dwaves(:) ! Direction gamma as defined in doc. for each wave component (rad)
real(wp),          intent(in)  :: pwaves(:) ! Phase delta as defined in doc. (rad)
!
integer :: nid, i
integer, parameter :: fmt = 100
integer, parameter :: shp = 6
integer, parameter :: amp = 1
integer, parameter :: order = 0
integer, parameter :: nstrip = 0
integer, parameter :: nsteps = 0
real(wp), parameter :: dt = -1.0_wp
! We need some C compatible characters
character(kind=c_char, len=:), allocatable :: ccid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
!
open( newunit=self % lu, file=file, access='stream', form='unformatted', &
      status='replace', action='write' )
!
!
nid = len_trim(cid) + 1  ! one extra for C null character
allocate(character(len=nid) :: ccid)
ccid = trim(cid) // c_null_char
cdate = timestamp() // c_null_char
cprog = trim(prog) // c_null_char 
!
write(self % lu) 37.0221_c_float
write(self % lu) int(fmt, c_int)
write(self % lu) int(shp, c_int)
write(self % lu) int(amp, c_int)
write(self % lu) cprog
write(self % lu) cdate
write(self % lu) int(nid, c_int)
write(self % lu) ccid(:nid)
write(self % lu) real(grav, c_float)
write(self % lu) real(lscale, c_float)
write(self % lu) int(nstrip, c_int)
write(self % lu) int(nsteps, c_int)
write(self % lu) real(dt, c_float)
write(self % lu) int(order, c_int)
write(self % lu) int(nwaves, c_int)
write(self % lu) real(depth, c_float)
do i = 1, nwaves
    write(self % lu) real(awaves(i), c_float)
    write(self % lu) real(kwaves(i), c_float)
    write(self % lu) real(dwaves(i), c_float)
    write(self % lu) real(pwaves(i), c_float)
end do
close(self % lu)
!  
contains

  function timestamp() result(timestr)
    ! Return a string with the current time in the form YYYY-MM-DD hh-mm-ss
    character (len = 19) :: timestr

    integer :: y, m, d, h, n, s
    integer :: values(8)
    
    call date_and_time (values = values)
    
    y = values(1)
    m = values(2)
    d = values(3)
    h = values(5)
    n = values(6)
    s = values(7)
    
    write (timestr, '(i4.4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2)' ) &
         y,'-',m,'-',d,' ',h,':',n,':',s
    
  end function timestamp
!
end subroutine build

!==============================================================================

end module swd_write_shape_6_def
