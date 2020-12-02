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
    real(wp) :: d
    integer :: nsumx, nsumy
    real(wp), allocatable, public :: kx(:,:), ky(:,:), k(:,:), tanhkh(:,:)
    type(swd_error), public :: error
contains
    private
    procedure, public :: fft_field_1D
    procedure, public :: fft_field_2D
    procedure, public :: impl2_to_impl1
    procedure, public :: x_fft
    procedure, public :: y_fft
    procedure, public :: swd_to_fft_coeffs_1D
    procedure, public :: swd_to_fft_coeffs_2D
    procedure, public :: close
    procedure :: nx_fft
    procedure :: ny_fft

end type swd_fft

interface swd_fft
    module procedure constructor
end interface swd_fft

real(wp), parameter :: pi = 3.14159265358979323846264338327950288419716939937510582097494_wp

contains

!==============================================================================

function constructor(nsumx, nsumy, dkx, dky, d) result(self)
integer, intent(in) :: nsumx, nsumy
real(wp), intent(in) :: dkx, dky
real(wp), intent(in) :: d ! water depth

type(swd_fft) :: self

integer :: i

self % nsumx = nsumx
self % dkx = dkx
self % sizex = 2.0_wp*pi/dkx
self % d = d

! it is assumed nsumy <= 1 means long-crested waves
if (nsumy <= 1) then
    self % nsumy = 0
    self % dky = 0.0_wp
    self % sizey = 0.0_wp
else
    self % nsumy = nsumy
    self % dky = dky
    self % sizey = 2.0_wp*pi/dky
end if

call self % error % clear()

! build matrices for derivatives, etc
allocate(self % kx(self % nsumx + 1, 2*self % nsumy + 1))
allocate(self % ky(self % nsumx + 1, 2*self % nsumy + 1))
allocate(self % k(self % nsumx + 1, 2*self % nsumy + 1))
allocate(self % tanhkh(self % nsumx + 1, 2*self % nsumy + 1))

! kx-matrix
do i = 1, self % nsumx + 1
    self % kx(i, :) = self%dkx*(i - 1)
end do

! ky-matrix and k-matrix
if (self % nsumy > 1) then
    do i = 1, self % nsumy + 1
        self % ky(:, i) = self%dky*(i - 1)
    end do
    do i = self % nsumy + 2, 2*self % nsumy + 1
        self % ky(:, i) = self%dky*(i - 2*self % nsumy - 2)
    end do
    self % k = sqrt(self % kx**2 + self % ky**2)
else
    self % ky = 0.0_wp
    self % k = self % kx
end if

! tanh(kh)-matrix
if (self % d > 0.0_wp) then
    self % tanhkh = tanh(self % k*self % d)
else
    self % tanhkh = 1.0_wp
end if 

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

function fft_field_1D(self, fft_coeffs, nx_fft_in) result(f)
class(swd_fft), intent(inout) :: self
complex(wp), intent(in) :: fft_coeffs(:, :)
integer, optional, intent(in) :: nx_fft_in
real(wp), allocatable :: f(:, :)
integer :: nx_fft

nx_fft = self % nx_fft(nx_fft_in)

if (self % error % raised()) then
    allocate(f(1, 1))
    f = huge(f)
    return
end if

f = irfft2(fft_coeffs, 2*self % nsumx, 1, nx_fft, 1)

end function fft_field_1D

!==============================================================================

function fft_field_2D(self, fft_coeffs, nx_fft_in, ny_fft_in) result(f)
class(swd_fft), intent(inout) :: self
complex(wp), intent(in) :: fft_coeffs(:, :)
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

f = irfft2(fft_coeffs, 2*self % nsumx, 2*self % nsumy + 1, nx_fft, ny_fft)

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
complex(wp), dimension(self % nsumx + 1, 2*self % nsumy + 1) :: fft_coeffs
integer :: ix, iy
real(wp) :: sc

fft_coeffs = cmplx(0.0_wp, 0.0_wp, kind=wp)

do ix = 0, self % nsumx
    sc = 0.5_wp ! scaling factor between fft and swd format
    if (ix == 0) sc = 1.0_wp ! special case for zero-wavenumber
    if (ix == self % nsumx) sc = 1.0_wp ! special case for nyquist wavenumber

    ! 0-freq and all positive wavenumbers
    do iy = 1, self % nsumy + 1
        fft_coeffs(ix + 1, iy) = sc*conjg(swd_coeffs(1, iy - 1, ix))
    end do

    ! all other negative wavenumbers up to -1
    do iy = self % nsumy + 2, 2*self % nsumy + 1
        fft_coeffs(ix + 1, iy) = sc*conjg(swd_coeffs(2, 2*self % nsumy + 2 - iy, ix))
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

if (self % nsumy == 0) then  ! 1D
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
    if (ny_fft_in <= -1) then
        ny_fft = -2*ny_fft_in*self % nsumy
    elseif (ny_fft_in >= 2*self % nsumy) then ! arbitrary grid size
        ny_fft = ny_fft_in
    else ! error for "invalid" downsampling (use nsumx/nsumy instead)
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

function impl2_to_impl1(self, impl2_coeffs) result(impl1_coeffs)
class(swd_fft), intent(in) :: self  ! Actual class
complex(wp), intent(in) :: impl2_coeffs(4, (self%nsumx + 1)*(self%nsumx + 2)/2)
complex(wp) :: impl1_coeffs(2, 0:self % nsumx, 0:self % nsumx)

integer ii, ix, iy

ii = 0
do ix = 0, self % nsumx
    iy = 0
    ii = ii + 1
    if (ix == 0) then
        impl1_coeffs(1,iy,ix) = impl2_coeffs(1,ii)
        cycle
    end if
    impl1_coeffs(1,iy,ix) = impl2_coeffs(1,ii)
    impl1_coeffs(1,ix,iy) = impl2_coeffs(3,ii)
    impl1_coeffs(2,ix,iy) = impl2_coeffs(4,ii)
    do iy = 1, ix - 1
        ii = ii + 1
        impl1_coeffs(1,iy,ix) = impl2_coeffs(1,ii)
        impl1_coeffs(2,iy,ix) = impl2_coeffs(2,ii)
        impl1_coeffs(1,ix,iy) = impl2_coeffs(3,ii)
        impl1_coeffs(2,ix,iy) = impl2_coeffs(4,ii)
    end do
    iy = ix
    ii = ii + 1
    impl1_coeffs(1,iy,ix) = impl2_coeffs(1,ii)
    impl1_coeffs(2,iy,ix) = impl2_coeffs(2,ii)
end do

end function impl2_to_impl1

!=============================================================================

subroutine close(self)
class(swd_fft) :: self  ! Object to destruct
!
if (allocated(self % kx)) deallocate(self % kx)
if (allocated(self % ky)) deallocate(self % ky)
if (allocated(self % k)) deallocate(self % k)
if (allocated(self % tanhkh)) deallocate(self % tanhkh)
!
end subroutine  close

!=============================================================================

end module swd_fft_def
