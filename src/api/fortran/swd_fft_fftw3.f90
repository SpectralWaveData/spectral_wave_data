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

function irfft2(fh, nx_in, ny_in, nx_out, ny_out) result(fE)
    ! 2D real IFFT of fh, possibly with padding to give outputsize nx_out x ny_out 
    complex(wp), dimension(nx_in/2 + 1, ny_in), intent(in) :: fh
    integer, intent(in) :: nx_in, ny_in, nx_out, ny_out
    real(wp), dimension(nx_out, ny_out) :: fE

    integer :: nx_int, ny_int
    complex(wp), dimension(nx_out/2 + 1, ny_out) :: fhE
    type(C_PTR) :: plan_IFFT
    integer, parameter :: flag = FFTW_ESTIMATE

    if (nx_in == nx_out .and. ny_in == ny_out) then ! same size
        fhE = fh
    else if (nx_in <= nx_out .and. ny_in <= ny_out) then ! upsample in both dimensions
        fhE = zeropad(fh, nx_in, ny_in, nx_out, ny_out)
    else if (nx_in >= nx_out .and. ny_in >= ny_out) then ! downsample in both dimensions
        fhE = truncate(fh, nx_in, ny_in, nx_out, ny_out)
    else ! valid in all possible cases
        nx_int = max(nx_out, nx_in)
        ny_int = max(ny_out, ny_in)
        fhE = truncate(zeropad(fh, nx_in, ny_in, nx_int, ny_int), nx_int, ny_int, nx_out, ny_out)
    end if

    call dfftw_plan_dft_c2r_2d(plan_IFFT, nx_out, ny_out, fhE, fE, flag)
    call dfftw_execute_dft_c2r(plan_IFFT, fhE, fE)
    call dfftw_destroy_plan(plan_IFFT)
    call fftw_cleanup()

end function irfft2

!==============================================================================

function zeropad(fh, nx, ny, nxE, nyE) result(fhE)
    ! Zero padding of fh, corresponding to an upsampling
    ! of f = ifft(fh) from a grid (nx x ny) to a new grid (nxE x nyE),
    ! where nxE >= nx and nyE >= ny.
    ! Tested to agree with scipy.signal.resample

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
        ! This can never happen due to checking in nx_fft and ny_fft
        error stop "swd_fft_fftw::zeropad: Invalid input."
    end if

end function zeropad

!==============================================================================

function truncate(fhE, nxE, nyE, nx, ny) result(fh)
    ! Truncation of fhE, corresponding to a downsampling
    ! of fE = ifft(fhE) from a grid (nxE x nyE) to a new grid (nx x ny),
    ! where nxE >= nx and nyE >= ny.
    ! Tested to agree with scipy.signal.resample

    complex(wp), dimension(nxE/2 + 1, nyE), intent(in) :: fhE
    integer, intent(in) :: nxE, nyE, nx, ny
    complex(wp), dimension(nx/2 + 1, ny) :: fh
    integer nxh, nyh, nxhE, nyhE

    nxh = nx/2 + 1
    nxhE = nxE/2 + 1
    nyh = ny/2 + 1
    nyhE = nyE/2 + 1
  
    fh = cmplx(0.0_wp, 0.0_wp, kind=wp)

    if (ny == 1 .and. nyE == 1) then
        fh(1:nxh, 1) = fhE(1:nxh, 1)
        if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
            fh(nxh, 1) = 2.0_wp*fh(nxh, 1)%re
        end if
    else if (ny > 1 .and. nyE > 1) then
        fh(1:nxh, 1:nyh) = fhE(1:nxh, 1:nyh)
        if (mod(ny, 2) == 0) then
            fh(1:nxh, nyh + 1:ny) = fhE(1:nxh, nyE - nyh + 3:nyE)
        else
            fh(1:nxh, nyh + 1:ny) = fhE(1:nxh, nyE - nyh + 2:nyE)
        end if
        ! handeling of nyquist frequency. This is tested carefully against
        ! scipy.signal.resample and is therefore very likely correct, although
        ! a bit difficult to understand
        if (mod(ny, 2) == 0 .and. ny /= nyE) then
            fh(1, nyh) = 2.0_wp*fh(1, nyh)%re
            fh(2:nxh, nyh) = fh(2:nxh, nyh) + fhE(2:nxh, nyE - nyh + 2)
        end if
        if (mod(nx, 2) == 0 .and. nxh /= nxhE) then
            fh(nxh, 1) = 2.0_wp*fh(nxh, 1)%re
            fh(nxh, 2:ny) = fh(nxh, 2:ny) + conjg(fh(nxh, ny:2:-1))
        end if
    else
        ! This can never happen due to checking in nx_fft and ny_fft
        error stop "swd_fft_fftw::truncate: Invalid input."  
    end if

end function truncate

!==============================================================================

end module swd_fft_fftw3
