module swd_fft_fftw3

use, intrinsic :: iso_c_binding

use kind_values, only: wp => kind_swd_interface

implicit none

include 'fftw3.f03'

! This module provides fft routines using the FFTW3 library.
!
! Written by Odin Gramstad, October, 26. 2020
!
!------------------------------------------------------------------------------

private
public :: irfft2

contains

!==============================================================================

function irfft2(fh, nx, ny, nxE, nyE) result(fE)
    ! 2D IFFT of fh, possibly with padding to outputsize nxE x nyE 
    complex(wp), dimension(nx/2 + 1, ny), intent(in) :: fh
    integer, intent(in) :: nx, ny, nxE, nyE
    real(wp), dimension(nxE, nyE) :: fE

    complex(wp), dimension(nxE/2 + 1, nyE) :: fhE
    type(C_PTR) :: plan_IFFT
    integer, parameter :: flag = FFTW_ESTIMATE


    if (nx == nxE .and. ny == nyE) then
        fhE = fh
    else if (nx <= nxE .and. ny <= nyE) then
        fhE = zeropad(fh, nx, ny, nxE, nyE)
    else
        error stop "Can only resample on a finer grid." ! This can never happen due to checking in elev_fft
    end if

    call dfftw_plan_dft_c2r_2d(plan_IFFT, nxE, nyE, fhE, fE, flag)
    call dfftw_execute_dft_c2r(plan_IFFT, fhE, fE)
    call dfftw_destroy_plan(plan_IFFT)
    call fftw_cleanup()

end function irfft2

!==============================================================================

function zeropad(fh, nx, ny, nxE, nyE) result(fhE)
    ! Zero padding of fh, corresponding to a FFT interpolation
    ! of f = ifft(fh) from a grid nx x ny to a new grid nxE x nyE.
    ! Tested to agree with scipy.resample

    complex(wp), dimension(nx/2 + 1, ny), intent(in) :: fh
    integer, intent(in) :: nx, ny, nxE, nyE
    complex(wp), dimension(nxE/2 + 1, nyE) :: fhE
    integer nxh, nyh, nxhE, nyhE

    nxh = nx/2 + 1
    nxhE = nxE/2 + 1
    nyh = ny/2 + 1
    nyhE = nyE/2 + 1

    fhE(1:nxhE, 1:nyE) = cmplx(0.0_wp, 0.0_wp, kind=wp)

    if (ny == 1 .and. nyE == 1) then
        fhE(1:nxh, 1) = fh(1:nxh, 1)
        if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
            fhE(nxh, 1) = 0.5_wp*fhE(nxh, 1)
        end if
    else if (ny > 1 .and. nyE > 1) then
        fhE(1:nxh, 1:nyh) = fh(1:nxh, 1:nyh)
        if (mod(ny, 2) == 0) then
            fhE(1:nxh, nyE - nyh + 3:nyE) = fh(1:nxh, nyh + 1:ny)
        else
            fhE(1:nxh, nyE - nyh + 2:nyE) = fh(1:nxh, nyh + 1:ny)
        end if
        ! handeling of nyquist frequency. This is tested carefully against
        ! scipy.signal.resample and is therefore very likely correct, although
        ! a bit difficult to understand
        if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
            fhE(nxh, :) = 0.5_wp*fhE(nxh, :)
        end if
        if (mod(ny, 2) == 0 .and. ny /= nyE) then
            fhE(:, nyh) = 0.5_wp*fhE(:, nyh)
            fhE(:, nyE - nyh + 2) = fhE(:, nyh)
        end if
    else
        error stop "Invalid input."  ! This can never happen due to checking in elev_fft
    end if

end function zeropad

!==============================================================================

end module swd_fft_fftw3
