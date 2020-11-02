module swd_fft_def

use kind_values, only: wp => kind_swd_internal
use swd_fft_fftw3, only : irfft2  ! alternatively, use/write another FFT-module that implements a irfft2-routine
use spectral_wave_data_error, only: swd_error

implicit none
private

! This module provides various routines used for the FFT evaluation.
!
! Written by Odin Gramstad, November, 02. 2020
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: swd_fft  ! Class for handling FFT-based evaluation
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type :: swd_fft
    private
    integer :: sizex, sizey
    integer :: nsumx, nsumy
    type(swd_error), public :: error
contains
    private
    procedure, public :: fft_field_1D
    procedure, public :: x_fft
    procedure, public :: y_fft
    procedure :: nx_fft
    procedure :: ny_fft
    procedure :: swd_to_fft_coeffs_1D
end type swd_fft

interface swd_fft
    module procedure constructor
end interface swd_fft

contains

!==============================================================================

function constructor(nsumx, nsumy, sizex, sizey) result(self)
integer, intent(in) :: nsumx, nsumy
real(wp), intent(in) :: sizex, sizey

type(swd_fft) :: self

self % nsumx = nsumx
self % sizex = sizex

! it is assumed nsumy <= 1 means long-crested waves
if (nsumy <= 1) then
    self % nsumy = 1
    self % sizey = 0.0_wp
else
    self % nsumy = nsumy
    self % sizey = sizey
end if

call self % error % clear()

end function constructor

!==============================================================================

function x_fft(self, nx_fft_in)
class(swd_fft), intent(inout) :: self
integer, optional, intent(in) :: nx_fft_in
real(wp), allocatable :: x_fft(:)
integer :: nx_fft
integer :: i

nx_fft = self % nx_fft(nx_fft_in)
x_fft = [(real(i, wp)*self % sizex/real(nx_fft, wp), i=0, nx_fft - 1)]

end function x_fft

!==============================================================================

function y_fft(self, ny_fft_in)
class(swd_fft), intent(inout) :: self
integer, optional, intent(in) :: ny_fft_in
real(wp), allocatable :: y_fft(:)
integer :: ny_fft
integer :: i

ny_fft = self % ny_fft(ny_fft_in)
y_fft = [(real(i, wp)*self % sizey/real(ny_fft, wp), i=0, ny_fft - 1)] 

end function y_fft

!==============================================================================

function fft_field_1D(self, swd_coeffs, nx_fft_in) result(f)
class(swd_fft), intent(inout) :: self
complex(wp), intent(in) :: swd_coeffs(0:)
integer, optional, intent(in) :: nx_fft_in
real(wp), allocatable :: f(:, :)
integer :: nx_fft

nx_fft = self % nx_fft(nx_fft_in)

if (self % error % raised()) then
    allocate(f(1, 1))
    f = huge(f)
    return
end if

f = irfft2(self % swd_to_fft_coeffs_1D(swd_coeffs), 2*self % nsumx, 1, nx_fft, 1)

end function fft_field_1D

!==============================================================================

function swd_to_fft_coeffs_1D(self, swd_coeffs) result(fft_coeffs)
class(swd_fft), intent(inout) :: self
complex(wp), intent(in) :: swd_coeffs(0:self % nsumx)
complex(wp), dimension(self % nsumx + 1, 1) :: fft_coeffs

fft_coeffs = cmplx(0.0_wp, 0.0_wp, kind=wp)
fft_coeffs(1:self % nsumx + 1, 1) = 0.5_wp*conjg(swd_coeffs(0:self % nsumx))

end function swd_to_fft_coeffs_1D

!==============================================================================

function nx_fft(self, nx_fft_in)
class(swd_fft), intent(inout) :: self
integer, optional, intent(in) :: nx_fft_in
integer :: nx_fft

character(len=*), parameter :: err_proc = 'swd_fft::nx_fft'
character(len=2500) :: err_msg(5)

if (present(nx_fft_in)) then
    if (nx_fft_in < 0) then
        nx_fft = -2*nx_fft_in*self % nsumx
    elseif (nx_fft_in >= 2*self % nsumx) then
        nx_fft = nx_fft_in
    else
        write(err_msg(1),'(a)') "Invalid grid size nx_fft."
        write(err_msg(2),'(a)') "nx_fft must either be a negative integer, or a positive"
        write(err_msg(3),'(a)') "integer larger than or equal to the size of the" 
        write(err_msg(4),'(a)') "smallest grid resolving all coefficients in the swd-file."
        write(err_msg(5),'(a, I0)') 'Smallest possible nx_fft = ', 2*self % nsumx
        call self % error % set_id_msg(err_proc, 1004, err_msg)
        nx_fft = 1
        return
    end if
else
    nx_fft = 2*self % nsumx
end if

end function nx_fft

!==============================================================================

function ny_fft(self, ny_fft_in)
class(swd_fft), intent(inout) :: self
integer, optional, intent(in) :: ny_fft_in
integer :: ny_fft

character(len=*), parameter :: err_proc = 'swd_fft::ny_fft'
character(len=2500) :: err_msg(5)

if (self % nsumy == 1) then  ! 1D
    if (present(ny_fft_in)) then
        if (abs(ny_fft_in) /= 1) then
            write(err_msg(1),'(a)') "Invalid grid size ny_fft."
            write(err_msg(2),'(a)') "ny_fft must be 1 for unidirectional waves."
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:2))
            ny_fft = 1
            return
        end if
    end if
    ny_fft = 1
else if (present(ny_fft_in)) then
    if (ny_fft_in < 0) then
        ny_fft = -2*ny_fft_in*self % nsumy
    elseif (ny_fft_in >= 2*self % nsumy) then
        ny_fft = ny_fft_in
    else
        write(err_msg(1),'(a)') "Invalid grid size ny_fft."
        write(err_msg(2),'(a)') "ny_fft must either be a negative integer, or a positive"
        write(err_msg(3),'(a)') "integer larger than or equal to the size of the" 
        write(err_msg(4),'(a)') "smallest grid resolving all coefficients in the swd-file."
        write(err_msg(5),'(a, I0)') 'Smallest possible ny_fft = ', 2*self % nsumy
        call self % error % set_id_msg(err_proc, 1004, err_msg)
        ny_fft = 1
        return
    end if
else
    ny_fft = 2*self % nsumy
end if

end function ny_fft

!==============================================================================

end module swd_fft_def
