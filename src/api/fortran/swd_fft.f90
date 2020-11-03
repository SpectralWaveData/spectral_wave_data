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
    real(wp) :: dkx, dky
    real(wp) :: sizex, sizey
    integer :: nsumx, nsumy
    type(swd_error), public :: error
contains
    private
    procedure, public :: fft_field_1D
    procedure, public :: fft_field_2D
    procedure, public :: x_fft
    procedure, public :: y_fft
    procedure :: nx_fft
    procedure :: ny_fft
    procedure :: swd_to_fft_coeffs_1D
    procedure :: swd_to_fft_coeffs_2D
end type swd_fft

interface swd_fft
    module procedure constructor
end interface swd_fft

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

function constructor(nsumx, nsumy, dkx, dky) result(self)
integer, intent(in) :: nsumx, nsumy
real(wp), intent(in) :: dkx, dky

type(swd_fft) :: self

self % nsumx = nsumx
self % dkx = dkx
self % sizex = 2.0_wp*pi/dkx

! it is assumed nsumy <= 1 means long-crested waves
if (nsumy <= 1) then
    self % nsumy = 1
    self % dky = 0.0_wp
    self % sizey = 0.0_wp
else
    self % nsumy = nsumy
    self % dky = dky
    self % sizey = 2.0_wp*pi/dky
end if

call self % error % clear()

end function constructor

!==============================================================================

function x_fft(self, nx_fft_in)
class(swd_fft), intent(inout) :: self
integer, optional, intent(in) :: nx_fft_in
real(wp), allocatable :: x_fft(:)
real(wp) :: dx
integer :: nx_fft
integer :: i

nx_fft = self % nx_fft(nx_fft_in)
dx = self % sizex/real(nx_fft, wp)
x_fft = [(i*dx, i=0, nx_fft - 1)]

end function x_fft

!==============================================================================

function y_fft(self, ny_fft_in)
class(swd_fft), intent(inout) :: self
integer, optional, intent(in) :: ny_fft_in
real(wp), allocatable :: y_fft(:)
real(wp) :: dy
integer :: ny_fft
integer :: i

ny_fft = self % ny_fft(ny_fft_in)
dy = self % sizey/real(ny_fft, wp)
y_fft = [(i*dy, i=0, ny_fft - 1)] 

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

function fft_field_2D(self, swd_coeffs, nx_fft_in, ny_fft_in) result(f)
class(swd_fft), intent(inout) :: self
complex(wp), intent(in) :: swd_coeffs(:, 0:, 0:)
integer, optional, intent(in) :: nx_fft_in, ny_fft_in
real(wp), allocatable :: f(:, :)
integer :: nx_fft, ny_fft

nx_fft = self % nx_fft(nx_fft_in)
ny_fft = self % ny_fft(ny_fft_in)

if (self % error % raised()) then
    allocate(f(1, 1))
    f = huge(f)
    return
end if

f = irfft2(self % swd_to_fft_coeffs_2D(swd_coeffs), 2*self % nsumx, 2*self % nsumy, nx_fft, ny_fft)

end function fft_field_2D

!==============================================================================

function swd_to_fft_coeffs_1D(self, swd_coeffs) result(fft_coeffs)
class(swd_fft), intent(inout) :: self
complex(wp), intent(in) :: swd_coeffs(0:self % nsumx)
complex(wp), dimension(self % nsumx + 1, 1) :: fft_coeffs

fft_coeffs = cmplx(0.0_wp, 0.0_wp, kind=wp)
fft_coeffs(1:self % nsumx + 1, 1) = 0.5_wp*conjg(swd_coeffs(0:self % nsumx))
fft_coeffs(self % nsumx + 1, 1) = 2.0_wp*fft_coeffs(self % nsumx + 1, 1) ! nyquist wavenumber

end function swd_to_fft_coeffs_1D

!==============================================================================

function swd_to_fft_coeffs_2D(self, swd_coeffs) result(fft_coeffs)
class(swd_fft), intent(in) :: self  ! Actual class
complex(wp), intent(in) :: swd_coeffs(2, 0:self % nsumy, 0:self % nsumx)
complex(wp), dimension(self % nsumx + 1, 2*self % nsumy) :: fft_coeffs
integer :: ix, iy
real(wp) :: sc

fft_coeffs = cmplx(0.0_wp, 0.0_wp, kind=wp)

do ix = 0, self % nsumx
    sc = 0.5_wp ! scaling factor between fft and swd format
    if (ix == 0) sc = 1.0_wp ! special case for zero-wavenumber
    if (ix == self % nsumx) sc = 1.0_wp ! special case for nyquist wavenumber

    ! 0-freq and all positive wavenumbers
    do iy = 1, self % nsumy
        fft_coeffs(ix + 1, iy) = sc*conjg(swd_coeffs(1, iy - 1, ix))
    end do

    ! the largest negative wavenumber (nyquist)
    fft_coeffs(ix + 1, self % nsumy + 1) = sc*(conjg(swd_coeffs(1, self % nsumy, ix) + swd_coeffs(2, self % nsumy, ix)))

    ! all other negative wavenumbers up to -1
    do iy = self % nsumy + 2, 2*self % nsumy
        fft_coeffs(ix + 1, iy) = sc*conjg(swd_coeffs(2, 2*self % nsumy + 1 - iy, ix))
    end do   
end do

end function swd_to_fft_coeffs_2D

!==============================================================================

function nx_fft(self, nx_fft_in)
class(swd_fft), intent(inout) :: self
integer, optional, intent(in) :: nx_fft_in
integer :: nx_fft

character(len=*), parameter :: err_proc = 'swd_fft::nx_fft'
character(len=250) :: err_msg(5)

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
        call self % error % set_id_msg(err_proc, 1004, err_msg, 1)
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
character(len=250) :: err_msg(5)

if (self % nsumy == 1) then  ! 1D
    if (present(ny_fft_in)) then
        if (abs(ny_fft_in) /= 1) then
            write(err_msg(1),'(a)') "Invalid grid size ny_fft."
            write(err_msg(2),'(a)') "ny_fft must be 1 for unidirectional waves."
            call self % error % set_id_msg(err_proc, 1004, err_msg(1:2), 1)
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
        call self % error % set_id_msg(err_proc, 1004, err_msg, 1)
        ny_fft = 1
        return
    end if
else
    ny_fft = 2*self % nsumy
end if

end function ny_fft

!==============================================================================

end module swd_fft_def
