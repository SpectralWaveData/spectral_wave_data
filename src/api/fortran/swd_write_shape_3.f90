module swd_write_shape_3_def

use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float, c_null_char

implicit none
private

! Please adjust wp to the working precision applied in the wave generator
integer, parameter :: wp = kind(1.0d0) ! (double precision)


! This module provides a class for writing SWD files for shape 3.
!
! Written by Jens Bloch Helmers, August, 1. 2019
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: swd_write_shape_3
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type :: swd_write_shape_3
    integer :: amp    ! Type of spectral amplitudes to be provided
    integer :: lu     ! Unit number associated with swd file
    integer :: nsteps ! Expected number of time steps
    integer :: it     ! Current time steps
    integer :: n      ! Number of spectral components in x-dir
    integer :: nh     ! Number of auxilary shapes (hat(n) in theory manual)
    integer :: nsf    ! Number of sea floor points
contains
    procedure :: init    ! Open the swd file and write the header section
    procedure :: update  ! Add spectral data at current time step
    procedure :: close   ! close the file
end type swd_write_shape_3

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

subroutine init(self, file, prog, cid, grav, Lscale, amp, &
                n, nh, order, dk, dt, nsteps, isf, nsf, xsf, zsf)
class(swd_write_shape_3),  intent(out) :: self  ! Object to create
character(len=*),  intent(in)  :: file  ! Requested name for new SWD file
character(len=*),  intent(in)  :: prog  ! Name of actual wave generator including version number
character(len=*),  intent(in)  :: cid   ! Text describing actual input parameters to the
                                        ! wave generator. Hint: write(cid, nml=wave_input_namelist)
real(wp),          intent(in)  :: grav  ! Acceleration of gravity applied in wave generator
real(wp),          intent(in)  :: lscale! Number of length units per meter applied in wave generator
                                        ! E.g. if US feet was applied: Lscale = 1/0.3048 = 3.2808  
integer,           intent(in)  :: amp   ! Flag to indicate type of provided spectral amplitudes
integer,           intent(in)  :: n     ! n as defined in the theory documentation
integer,           intent(in)  :: nh    ! hat(n) as defined in the theory documentation
integer,           intent(in)  :: order ! Expansion order applied in wave generator. (<0 for fully non-linear)
real(wp),          intent(in)  :: dk    ! Spacing of k-numbers in x-direction
real(wp),          intent(in)  :: dt    ! Spacing of time steps in the swd file
integer,           intent(in)  :: nsteps ! Total number of time steps to be stored in the swd file.
integer,           intent(in)  :: isf    ! Flag to indicate geometric description of the sea floor
                                         ! 0: Piecewise linear sea floor
integer,           intent(in)  :: nsf    ! Number of offset points describing the sea floor
real(wp),          intent(in)  :: xsf(:) ! X-positions of sea floor points
real(wp),          intent(in)  :: zsf(:) ! Z-positions of sea floor points
!
integer :: nid, i
integer, parameter :: fmt = 100
integer, parameter :: shp = 3
integer, parameter :: nstrip = 0
! We need some C compatible characters
character(kind=c_char, len=:), allocatable :: ccid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate

call input_sanity_checks

open( newunit=self % lu, file=file, access='stream', form='unformatted', &
      status='replace', action='write', convert='little_endian')
! The CONVERT option is an Intel, gfortran, HP and IBM extension. 
! Modify this parameter for other compilers if your compiler complains.

!
self % amp = amp
self % nsteps = nsteps
self % it = 0
self % n = n
self % nh = nh
self % nsf = nsf
!
ccid = trim(cid) // c_null_char
nid = len(ccid)
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
write(self % lu) ccid
write(self % lu) real(grav, c_float)
write(self % lu) real(lscale, c_float)
write(self % lu) int(nstrip, c_int)
write(self % lu) int(nsteps, c_int)
write(self % lu) real(dt, c_float)
write(self % lu) int(order, c_int)
write(self % lu) int(n, c_int)
write(self % lu) int(nh, c_int)
write(self % lu) real(dk, c_float)
write(self % lu) int(isf, c_int)
write(self % lu) int(nsf, c_int)
do i = 1, nsf
    write(self % lu) real(xsf(i), c_float)
end do
do i = 1, nsf
    write(self % lu) real(zsf(i), c_float)
end do

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
  subroutine input_sanity_checks
    if (nsf > 0) then
        ! Sea floor is completely submerged...
        if (any(zsf(1:nsf) >= 0.0_wp)) then
            print*, "ERROR: from swd_write_shape_3 :: init"
            print*, "All zsf(:) should be negative"
            print*, "zsf(1:nsf) = ", zsf(1:nsf)
            stop
        end if
        ! Sea floor x-range is in domain [0, 2pi/dk]
        if (abs(xsf(1)) > epsilon(xsf)) then
            print*, 'ERROR from swd_write_shape_3 :: init'
            print*, 'xsf(1) should be 0.0'
            print*, 'xsf(1) = ', xsf(1)
            stop
        end if
        if (abs(xsf(nsf) - 2*pi/dk) > 1.0e-4_wp * 2*pi/dk) then
            print*, 'ERROR from swd_write_shape_3 :: init'
            print*, 'xsf(nsf) should be 2*pi/dk'
            print*, 'xsf(nsf) = ', xsf(nsf)
            stop
        end if
        ! Sea floor x-values should be monotonic increasing
        do i = 2, nsf
            if (xsf(i-1) >= xsf(i)) then
                print*, 'ERROR from swd_write_shape_3 :: init'
                print*, 'xsf(1:nsf) should be monotonic increasing'
                print*, 'xsf(1:nsf) = ', xsf(1:nsf)
                stop
            end if
        end do
    end if
    if (nsf==0 .and. nh/=0) then  
        print*, 'ERROR from swd_write_shape_3 :: init'
        print*, 'nh should be 0 if nsf is 0'
        print*, 'nsf, nh = ', nsf, nh
        stop
    end if
    if (nsf>0 .and. nh<=0) then  
        print*, 'ERROR from swd_write_shape_3 :: init'
        print*, 'nh should be positive if nsf is positive'
        print*, 'nsf, nh = ', nsf, nh
        stop
    end if
  
  end subroutine input_sanity_checks
  
end subroutine init

!==============================================================================

subroutine update(self, h, ht, c, ct, ch, cht)
! Output of temporal functions as defined in the shape 3 class
class(swd_write_shape_3), intent(inout) :: self   ! Object to update
complex(wp),           intent(in) :: h(0:)   ! h(0:n) temporal amp
complex(wp),           intent(in) :: ht(0:)  ! dh/dt(0:n) temporal amp
complex(wp), optional, intent(in) :: c(0:)   ! c(0:n) temporal amp
complex(wp), optional, intent(in) :: ct(0:)  ! dc/dt(0:n) temporal amp
complex(wp), optional, intent(in) :: ch(0:)  ! hat(c)(0:nh) temporal amp
complex(wp), optional, intent(in) :: cht(0:) ! d hat(c)/dt(0:nh) temporal amp
!
self % it = self % it + 1
call dump(h, self % n)
call dump(ht, self % n)
if (self % amp < 3) then
    call dump(c, self % n)
    call dump(ct, self % n)
    if (self % nsf > 1) then
        call dump(ch, self % nh)
        call dump(cht, self % nh)
    end if
end if

contains

    subroutine dump(array, n)
        complex(wp), intent(in) :: array(0:)
        integer,     intent(in) :: n
        integer :: j
        do j = 0, n
            write(self % lu) cmplx(array(j), kind=c_float)
        end do
    end subroutine dump
!
end subroutine update

!==============================================================================

subroutine close(self)
class(swd_write_shape_3), intent(inout) :: self   ! Object to update
!
if (self % it /= self % nsteps) then
    print*, "WARNING: from swd_write_shape_3 :: close"
    print*, "Specified number of time steps = ", self % nsteps
    print*, "Number of provided time steps = ", self % it
end if
!
close(self % lu)
!
end subroutine close

!==============================================================================

end module swd_write_shape_3_def
