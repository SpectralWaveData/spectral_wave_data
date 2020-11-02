module spectral_wave_data_error

implicit none
private

! This module a simple error handler for spectral_wave_data to ensure that
! spectral_wave_data will never abort. Consequently, it is the end user who
! is responsible for checking if any error has occured when applying this
! library.
!
! Written by Jens Bloch Helmers, Januar, 12. 2020
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: swd_error  ! Class for handling errors
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type swd_error
    integer  :: id=0  ! Code to flag kind of error
                      ! 1001: SWD file not able to open
                      ! 1002: SWD file has wrong binary convention
                      ! 1003: SWD file data read error
                      ! 1004: invalid input parameters
                      ! 1005: Not able to allocate memory
    character(len=:), allocatable :: msg ! Characters to describe the error
contains
    procedure :: raised     ! Return .true. if error has been signaled
    procedure :: set_id_msg ! Set error code and message characters
    procedure :: get_id     ! Return error id
    procedure :: get_msg    ! Return error message
    procedure :: clear      ! Clear error signal (id=0)
end type swd_error

contains

!==============================================================================

function raised(self) result(res)
class(swd_error), intent(in) :: self ! Error handler
logical                      :: res  ! .true. if error has been signaled
!
res = self % id /= 0
!
end function raised

!==============================================================================

subroutine set_id_msg(self, proc, id, msg, fmt_in)
class(swd_error), intent(out) :: self   ! Error handler
integer,          intent(in)  :: id     ! Error code
character(len=*), intent(in)  :: proc   ! Name of calling procedure
character(len=*), intent(in)  :: msg(:) ! Error messages
integer, optional, intent(in)  :: fmt_in ! Set to /=0 for format more suitable for 
                                         ! "forwarding" errors from other modules/routines
!
integer :: i, n, fmt
character(len=*), parameter :: header = "Error detected in Spectral_Wave_Data:"
character(len=24) :: cid
!
call self % clear()
!
fmt = 0
if (present(fmt_in)) then
    fmt = fmt_in
end if

self % id = id
write(cid,'(a,i0)') 'SWD error code = #', id

! Count size of error string
if (fmt == 0) then
    n = len_trim(header) + 1 + len_trim(proc) + 1 + len_trim(cid)
else
    n = len_trim(proc)
end if

do i = 1, size(msg)
    if (len_trim(msg(i)) == 0) cycle
    n = n + 1 ! Insert end-of-line character
    n = n + len_trim(msg(i))
end do
allocate(character(len=n) :: self % msg)

! Assign to self %msg
if (fmt == 0) then
    self % msg = trim(header) // new_line('fortran') // trim(proc) // new_line('fortran') // trim(cid)
else
    self % msg = trim(proc)
end if
do i = 1, size(msg)
    if (len_trim(msg(i)) == 0) cycle
    self % msg = trim(self % msg) // new_line('fortran') // trim(msg(i))
end do
!
end subroutine set_id_msg

!==============================================================================

function get_id(self) result(res)
class(swd_error), intent(in) :: self ! Error handler
integer                      :: res  ! Return error code
!
res = self % id
!
end function get_id

!==============================================================================

function get_msg(self) result(res)
class(swd_error), intent(in)  :: self ! Error handler
character(len=:), allocatable :: res  ! Return error code
!
res = trim(self % msg)
!
end function get_msg

!==============================================================================

subroutine clear(self)
class(swd_error), intent(inout) :: self   ! Error handler
!
self % id = 0
if (allocated(self % msg)) deallocate(self % msg)
!
end subroutine clear

!==============================================================================

end module spectral_wave_data_error
