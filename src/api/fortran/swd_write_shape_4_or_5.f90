module swd_write_shape_4_or_5_def

use, intrinsic :: iso_c_binding,   only: c_char, c_int, c_float, c_null_char

implicit none
private

! Please adjust wp to the working precision applied in the wave generator
integer, parameter :: wp = kind(1.0d0) ! (double precision)


! This module provides a class for writing SWD files for shape 4 and 5.
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
public :: swd_write_shape_4_or_5
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type swd_write_shape_4_or_5
    integer :: amp    ! Type of spectral amplitudes to be provided
    integer :: lu     ! Unit number associated with swd file
    integer :: nsteps ! Expected number of time steps
    integer :: it     ! Current time steps
    integer :: nx     ! Number of points in physical space (x-dir)
    integer :: ny     ! Number of points in physical space (y-dir)
    integer :: nx_swd ! nx in swd theory description (=int(nx/2))
    integer :: ny_swd ! ny in swd theory description (=int(ny/2))
contains
    procedure :: init    ! Open the swd file and write the header section
    procedure :: update  ! Add spectral data at current time step
    procedure :: update_fft  ! Alternative output method based on FFTW r2c data
    procedure :: close   ! close the file
end type swd_write_shape_4_or_5
  
contains

!==============================================================================

subroutine init(self, file, prog, cid, grav, Lscale, amp, &
                nx, ny, order, dkx, dky, dt, nsteps, d)
class(swd_write_shape_4_or_5),  intent(out) :: self  ! Object to create
character(len=*),  intent(in)  :: file  ! Requested name for new SWD file
character(len=*),  intent(in)  :: prog  ! Name of actual wave generator including version number
character(len=*),  intent(in)  :: cid   ! Text describing actual input parameters to the
                                        ! wave generator. Hint: write(cid, nml=wave_input_namelist)
real(wp),          intent(in)  :: grav  ! Acceleration of gravity applied in wave generator
real(wp),          intent(in)  :: lscale! Number of length units per meter applied in wave generator
                                        ! E.g. if US feet was applied: Lscale = 1/0.3048 = 3.2808  
integer,           intent(in)  :: amp   ! Flag to indicate type of provided spectral amplitudes
integer,           intent(in)  :: nx    ! nx as defined in the theory documentation
integer,           intent(in)  :: ny    ! ny as defined in the theory documentation
integer,           intent(in)  :: order ! Expansion order applied in wave generator. (<0 for fully non-linear)
real(wp),          intent(in)  :: dkx   ! Spacing of k-numbers in x-direction
real(wp),          intent(in)  :: dky   ! Spacing of k-numbers in y-direction
real(wp),          intent(in)  :: dt    ! Spacing of time steps in the swd file
integer,           intent(in)  :: nsteps ! Total number of time steps to be stored in the swd file.
real(wp),          intent(in), optional  :: d ! Water depth (if not present, deep water (shape 4) is assumed)
!
integer :: nid
integer, parameter :: fmt = 100
integer :: shp
integer, parameter :: nstrip = 0
! We need some C compatible characters
character(kind=c_char, len=:), allocatable :: ccid
character(kind=c_char, len=30) :: cprog
character(kind=c_char, len=20) :: cdate
!
open( newunit=self % lu, file=file, access='stream', form='unformatted', &
      status='replace', action='write', convert='little_endian')
! The CONVERT option is an Intel, GFortran, HP and IBM extension. 
! Modify this parameter for other compilers if your compiler complains.

!
self % amp = amp
self % nsteps = nsteps
self % it = 0
self % nx = nx
self % ny = ny
self % nx_swd = nx/2
self % ny_swd = ny/2
!
shp = 4
if (present(d)) then
   if (d > 0.0_wp) then
      shp = 5
   end if
end if
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
write(self % lu) int(self % nx_swd, c_int)
write(self % lu) int(self % ny_swd, c_int)
write(self % lu) real(dkx, c_float)
write(self % lu) real(dky, c_float)
if (shp == 5) then
   write(self % lu) real(d, c_float)
end if

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
end subroutine init

!==============================================================================

subroutine update(self, h, ht, c, ct)
! Output of temporal functions as defined in the shape 4 and 5 classes
class(swd_write_shape_4_or_5), intent(inout) :: self   ! Object to update
complex(wp),           intent(in) :: h(-self % ny_swd:, 0:)   ! h(-ny_swd:ny_swd, 0:nx_swd) temporal amp
complex(wp),           intent(in) :: ht(-self % ny_swd:, 0:)  ! dh/dt(-ny_swd:ny_swd, 0:nx_swd) temporal amp
complex(wp), optional, intent(in) :: c(-self % ny_swd:, 0:)   ! c(-ny_swd:ny_swd, 0:nx_swd) temporal amp
complex(wp), optional, intent(in) :: ct(-self % ny_swd:, 0:)  ! dc/dt(-ny_swd:ny_swd, 0:nx_swd) temporal amp
!
self % it = self % it + 1
call dump(h)
call dump(ht)
if (self % amp < 3) then
    call dump(c)
    call dump(ct)
end if

contains

    subroutine dump(array)
        complex(wp), intent(in) :: array(-self % ny_swd:, 0:)
        integer :: jx, jy
        !---------------------------------------------------------------------------
        ! Explicit loops to demonstrate element ordering for all languages.
        !---------------------------------------------------------------------------
        do jx = 0, self % nx_swd
            do jy = -self % ny, self % ny_swd
                write(self % lu) cmplx(array(jy, jx), kind=c_float)
            end do
        end do
    end subroutine dump
!
end subroutine update

!==============================================================================

subroutine update_fft(self, h_fft, ht_fft, c_fft, ct_fft)
! Alternative to using the 'update' routine.
! Output of temporal functions as defined in the shape 4 and 5 classes
! Input arrays are in the form returned by FFTW r2c 2D transforms
class(swd_write_shape_4_or_5), intent(inout) :: self   ! Object to update
complex(wp),           intent(in) :: h_fft(:, :)   ! h_fft(1:nx/2+1, 1:ny) temporal amp
complex(wp),           intent(in) :: ht_fft(:, :)  ! dh/dt(1:nx/2+1, 1:ny) temporal amp
complex(wp), optional, intent(in) :: c_fft(:, :)   ! c(1:nx/2+1, 1:ny) temporal amp
complex(wp), optional, intent(in) :: ct_fft(:, :)  ! dc/dt(1:nx/2+1, 1:ny) temporal amp
!
self % it = self % it + 1
call dump_fft(h_fft)
call dump_fft(ht_fft)
if (self % amp < 3) then
    call dump_fft(c_fft)
    call dump_fft(ct_fft)
end if

contains
  
  subroutine dump_fft(array)
    complex(wp), intent(in) :: array(:,:)
    integer :: jx, jy
    real(wp) :: sc
    !---------------------------------------------------------------------------
    ! Explicit loops to demonstrate element ordering for all languages.
    !---------------------------------------------------------------------------
    do jx = 1, self % nx/2 + 1
       sc = 2.0_wp ! scaling factor between fft and swd format
       if (jx == 1) sc = 1.0_wp ! special case for zero-frequency
       if (mod(self % nx, 2) == 0 .and. jx == self % nx/2 + 1) sc = 1.0_wp ! special case for nyquist-frequency for even nx
       ! the first negative freq (-ny_swd)
       if (mod(self % ny, 2) == 0) then
          write(self % lu) cmplx(0.5_wp*sc*conjg(array(jx,self % ny - self % ny/2 + 1)), kind=c_float) 
       else
          write(self % lu) cmplx(sc*conjg(array(jx,self % ny - self % ny/2 + 1)), kind=c_float)
       end if
       ! all other negative freqs up to -1
       do jy = self % ny - self % ny/2 + 2, self % ny 
          write(self % lu) cmplx(sc*conjg(array(jx,jy)), kind=c_float)
       end do
       ! 0-freq and all positive freqs
       do jy = 1, self % ny - self % ny/2 
          write(self % lu) cmplx(sc*conjg(array(jx,jy)), kind=c_float)
       end do
       ! if even ny, nyquist freq must be on positive side as well
       if (mod(self % ny, 2) == 0) then
          write(self % lu) cmplx(0.5_wp*sc*conjg(array(jx,self % ny - self % ny/2 + 1)), kind=c_float)
       end if
    end do
  end subroutine dump_fft
!
end subroutine update_fft

!==============================================================================

subroutine close(self)
class(swd_write_shape_4_or_5), intent(inout) :: self   ! Object to update
!
if (self % it /= self % nsteps) then
    print*, "WARNING: from swd_write_shape_4_or_5 :: close"
    print*, "Specified number of time steps = ", self % nsteps
    print*, "Number of provided time steps = ", self % it
end if
!
close(self % lu)
!
end subroutine close

!==============================================================================

end module swd_write_shape_4_or_5_def
