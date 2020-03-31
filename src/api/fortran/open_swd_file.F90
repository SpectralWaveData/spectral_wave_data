module open_swd_file_def

use, intrinsic :: iso_c_binding, only: c_int, c_float, c_double

implicit none
private

! This module provides functions for open and testing endianess of a swd file. 
!
! Written by Jens Bloch Helmers, January, 10. 2020
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: swd_magic_number       ! The magic number
public :: swd_validate_binary_convention  ! Check binary convention of SWD file
public :: open_swd_file          ! Open swd file using requested endian
public :: ok_swd_magic_number    ! Return .true. if magic number is correct
                                 ! Used for checking byte soundness of SWD file 
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

!                     King Harald V was borned 37.0221
real(c_float), parameter :: swd_magic_number = 37.0221_c_float

contains

!==============================================================================

subroutine open_swd_file(newunit, file, status, as_little_endian, iostat)
integer,           intent(out) :: newunit ! Associated file unit number (output)
character(len=*),  intent(in)  :: file    ! Name of SWD file to open with correct endian
character(len=*),  intent(in)  :: status  ! 'new', 'replace' or 'old'
logical,           intent(in)  :: as_little_endian ! True if little_endian is applied
integer,           intent(out) :: iostat  ! =0 if succefull open statement
!
integer :: lu
character(len=5) :: action
character(len=20) :: endian
!
if (status == 'old' .or. status == 'OLD') then
    action = 'read'
else
    action = 'write'
end if

#if (__GFORTRAN__ || __INTEL_COMPILER)

! The convert argument in the open statement is a GFortran and Intel language extension
if (as_little_endian) then
    endian = 'little_endian'
else
    endian = 'big_endian'
end if
open(newunit=newunit, file=file, access='stream', form='unformatted', &
     status=status, action=action, convert=endian, iostat=iostat)

#else

! Some other compilers (like IBM) also support the convert argument.
! Consequently, we try the same recipe as above. If this crash
! please add a similar preprocessor section for your specific compiler.
if (as_little_endian) then
    endian = 'little_endian'
else
    endian = 'big_endian'
end if
open(newunit=newunit, file=file, access='stream', form='unformatted', &
    status=status, action=action, convert=endian, iostat=iostat)

#endif
!
end subroutine open_swd_file

!==============================================================================

subroutine swd_validate_binary_convention(file_swd, msg)
character(len=*), intent(in)  :: file_swd ! Name of actual SWD file.
character(len=*), intent(out) :: msg      ! Error message. Empty if no errors.
!
integer :: lu, ios
real(c_float) :: value
real(c_double) :: dvalue
!
! SWD should be little_endian...
call open_swd_file(newunit=lu, file=file_swd, status='old', &
                   as_little_endian=.true., iostat=ios)
if (ios /= 0) then
    msg = "Not able to open SWD file."
    return
end if

read(lu) value
close(lu)

if (ok_swd_magic_number(value)) then
    msg = ""
else
    ! Try big endian...
    call open_swd_file(newunit=lu, file=file_swd, status='old', &
                       as_little_endian=.false., iostat=ios)
    read(lu) value
    close(lu)
    if (ok_swd_magic_number(value)) then
        msg = "This SWD file is written in 'big_endian' byte convention. "// &
              "It should be in 'little_endian'."
    else
        ! Could it be double precision in SWD file?
        call open_swd_file(newunit=lu, file=file_swd, status='old', &
                           as_little_endian=.true., iostat=ios)
        read(lu) dvalue
        close(lu)
        if (ok_swd_magic_number(real(dvalue, c_float))) then
            msg = "This SWD file seems to written using double precision. "// &
                  "It should be in single precision."
        else
            ! Try big endian using double...
            call open_swd_file(newunit=lu, file=file_swd, status='old', &
                               as_little_endian=.false., iostat=ios)
            read(lu) dvalue
            close(lu)
            if (ok_swd_magic_number(real(dvalue, c_float))) then
                msg = "This SWD file seems to written using double precision and big " // &
                      "endian. It should be in single precision and little endian."
            else
                msg = "This SWD file has the wrong initial magic number."
            end if
        end if
    end if
end if
!
end subroutine swd_validate_binary_convention

!==============================================================================
  
function ok_swd_magic_number(val) result(res)
real(c_float), intent(in) :: val
logical :: res  ! True if val and swd_magic_number are very close
!
res = abs(val - swd_magic_number) / swd_magic_number < 1.e-6_c_float
!
end function ok_swd_magic_number

!==============================================================================

end module open_swd_file_def
