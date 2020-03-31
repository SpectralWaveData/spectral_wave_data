module spectral_wave_data_allocate_def

use, intrinsic :: iso_c_binding,   only: c_int, c_float, c_char

! In kind_values.f90 adjust kind_swd_interface to match applied 
! precision in the application program. (real, double precision etc.)
use kind_values, only: knd => kind_swd_interface

! Open file with correct endianess
use open_swd_file_def, only: open_swd_file, swd_validate_binary_convention

! The abstract base class
use spectral_wave_data_def, only: spectral_wave_data

! Include support for all spectral_wave_data implementations
use spectral_wave_data_shape_1_impl_1_def, only: spectral_wave_data_shape_1_impl_1
use spectral_wave_data_shape_2_impl_1_def, only: spectral_wave_data_shape_2_impl_1
use spectral_wave_data_shape_3_impl_1_def, only: spectral_wave_data_shape_3_impl_1
use spectral_wave_data_shape_4_impl_1_def, only: spectral_wave_data_shape_4_impl_1
use spectral_wave_data_shape_4_impl_2_def, only: spectral_wave_data_shape_4_impl_2
use spectral_wave_data_shape_5_impl_1_def, only: spectral_wave_data_shape_5_impl_1
use spectral_wave_data_shape_6_impl_1_def, only: spectral_wave_data_shape_6_impl_1

implicit none
private

! This module provides a procedure for selecting and allocating
! the proper SWD class based on the header of the actual SWD-file
! and optional requested implementation.
!
! Developers may add their own implementations which can
! be selected using this procedure.
!
! Written by Jens Bloch Helmers, August, 2. 2019
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: spectral_wave_data_allocate
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

contains

!==============================================================================

subroutine spectral_wave_data_allocate(swd, file_swd, x0, y0, t0, beta, rho, &
                                       nsumx, nsumy, impl, ipol, norder,     &
                                       dc_bias)
!------------------------------------------------------------------------------
! This procedure selects and allocates the proper SWD class based on 
! the header of the actual SWD-file and an optional requested implementation.
!------------------------------------------------------------------------------

! The actual class to be allocated
class(spectral_wave_data), allocatable, intent(out) :: swd

! Name of actual swd file
character(len=*), intent(in) :: file_swd

! Relation between SWD and application coordinates
real(knd), intent(in) :: x0, y0, t0, beta 

! The remaining parameters are optional...

! Density of water (applied for pressure calculations)
real(knd), optional, intent(in) :: rho

! If present and nsumx>-1: apply nsumx number of spectral components
!  in x-direction, else apply all spectral components from the SWD file.
integer, optional, intent(in) :: nsumx 


! If present and nsumy>-1: apply nsumy number of spectral components 
! in y-direction, else apply all spectral components from the SWD file.
integer, optional, intent(in) :: nsumy

! Index to determine actual derived class: 
!   0  = Default
!   <0 = In-house and experimental implementations 
!   >0 = Validated implementations available open-software
integer, optional, intent(in) :: impl  

! Index to request actual temporal interpolation scheme
!   0 = Default (C^2 continous scheme)
!   1 = C^1 continous
!   2 = C^3 continous
integer, optional, intent(in) :: ipol

! Expansion order to apply in kinematics for z>0
!   0  = Apply expansion order specified on swd file (default)
!   <0 = Apply exp(kj z)
!   >0 = Apply expansion order = norder 
integer, optional, intent(in) :: norder

! Control application of zero-frequency bias present in SWD file
!  True  = Apply zero frequency amplitudes from SWD file. 
!  False = Suppress contribution from zero frequency amplitudes (Default)
logical, optional, intent(in):: dc_bias
!------------------------------------------------------------------------------

integer :: impl_swd, nsumx_swd, nsumy_swd, ios, shp, amp, ipol_swd
integer :: norder_swd, err_id, lu
real(knd) :: rho_swd
logical :: select_ok, dc_bias_swd, equal_nx_ny
character(len=*), parameter :: err_proc = 'spectral_wave_data_allocate_def::spectral_wave_data_allocate'
character(len=250) :: err_msg(5)

! Assign default values if not present...
if (present(rho)) then
    rho_swd = rho
else
    rho_swd = 1025.0
end if
if (present(nsumx)) then
    nsumx_swd = nsumx
else
    nsumx_swd = -1
end if
if (present(nsumy)) then
    nsumy_swd = nsumy
else
    nsumy_swd = -1
end if
if (present(impl)) then
    impl_swd = impl
else
    impl_swd = 0
end if

if (present(ipol)) then
    ipol_swd = ipol
else
    ipol_swd = 0
end if

if (present(norder)) then
    norder_swd = norder
else
    norder_swd = 0
end if

if (present(dc_bias)) then
    dc_bias_swd = dc_bias
else
    dc_bias_swd = .false.
end if

! Check if it is possible to open input SWD file...
call open_swd_file(newunit=lu, file=file_swd, status='old', &
                   as_little_endian=.true., iostat=ios)
if (ios /= 0) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(file_swd)
    write(err_msg(2),'(a)') "Not able to open SWD file."
    allocate(spectral_wave_data_shape_1_impl_1 :: swd)  ! Empty object to store errors
    call swd % error % set_id_msg(err_proc, 1001, err_msg(1:2))
    return
end if
close(lu)

! Before we decode the SWD file we check if
! it applies the correct binary convention...
call swd_validate_binary_convention(file_swd, err_msg(2))
if (err_msg(2) /= '') then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(file_swd)
    allocate(spectral_wave_data_shape_1_impl_1 :: swd)  ! Empty object to store errors
    call swd % error % set_id_msg(err_proc, 1002, err_msg(1:2))
    return
end if

! Obtain shp and amp from actual SWD file
call swd_sniff_shp_amp()
if (err_msg(2) /= '') then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(file_swd)
    allocate(spectral_wave_data_shape_1_impl_1 :: swd)  ! Empty object to store errors
    call swd % error % set_id_msg(err_proc, err_id, err_msg(1:2))
    return
end if

select_ok = .false.
select case(shp)
case(1)
    if (impl_swd == 0 .or. impl_swd == 1) then
        if (amp/=2) then
            select_ok = .true.
            allocate(swd,                                      &
                source=spectral_wave_data_shape_1_impl_1(      &
                file_swd, x0, y0, t0, beta, rho=rho_swd,       &
                nsumx=nsumx_swd, ipol=ipol_swd, norder=norder_swd, &
                dc_bias=dc_bias_swd), stat=ios)
        end if
    end if
case(2)
    if (impl_swd == 0 .or. impl_swd == 1) then
        if (amp/=2) then
            select_ok = .true.
            allocate(swd,                                      &
                source=spectral_wave_data_shape_2_impl_1(      &
                file_swd, x0, y0, t0, beta, rho=rho_swd,       &
                nsumx=nsumx_swd, ipol=ipol_swd, norder=norder_swd, &
                dc_bias=dc_bias_swd), stat=ios)
        end if
    end if
case(3)
    if (impl_swd == 0 .or. impl_swd == 1) then
        if (amp/=2) then
            select_ok = .true.
            allocate(swd,                                      &
                source=spectral_wave_data_shape_3_impl_1(      &
                file_swd, x0, y0, t0, beta, rho=rho_swd,       &
                nsumx=nsumx_swd, ipol=ipol_swd, norder=norder_swd, &
                dc_bias=dc_bias_swd), stat=ios)
        end if
    end if
case(4)
    if (impl_swd == 0) then
        equal_nx_ny = swd_sniff_equal_nx_ny()
        if (err_msg(2) /= '') then
            write(err_msg(1),'(a,a)') 'SWD file: ', trim(file_swd)
            allocate(spectral_wave_data_shape_4_impl_1 :: swd)  ! Empty object to store errors
            call swd % error % set_id_msg(err_proc, err_id, err_msg(1:2))
            return
        end if
        if (equal_nx_ny) then
            impl_swd = 2
        else
            impl_swd = 1
        end if
    end if
    if (impl_swd == 1) then
        if (amp/=2) then
            select_ok = .true.
            allocate(swd,                                            &
                source=spectral_wave_data_shape_4_impl_1(            &
                file_swd, x0, y0, t0, beta, rho=rho_swd,             &
                nsumx=nsumx_swd, nsumy=nsumy_swd, norder=norder_swd, &
                dc_bias=dc_bias_swd, ipol=ipol_swd), stat=ios)
        end if
    else if (impl_swd == 2) then
        if (amp/=2) then
            select_ok = .true.
            allocate(swd,                                            &
                source=spectral_wave_data_shape_4_impl_2(            &
                file_swd, x0, y0, t0, beta, rho=rho_swd,             &
                nsumx=nsumx_swd, nsumy=nsumy_swd, norder=norder_swd, &
                dc_bias=dc_bias_swd, ipol=ipol_swd), stat=ios)
        end if
    end if
case(5)
    if (impl_swd == 0 .or. impl_swd == 1) then
        if (amp/=2) then
            select_ok = .true.
            allocate(swd,                                            &
                source=spectral_wave_data_shape_5_impl_1(            &
                file_swd, x0, y0, t0, beta, rho=rho_swd,             &
                nsumx=nsumx_swd, nsumy=nsumy_swd, norder=norder_swd, &
                dc_bias=dc_bias_swd, ipol=ipol_swd), stat=ios)
       end if
    end if
case(6)
    if (impl_swd == 0 .or. impl_swd == 1) then
        select_ok = .true.
        allocate(swd,                                      &
            source=spectral_wave_data_shape_6_impl_1(      &
            file_swd, x0, y0, t0, beta, rho=rho_swd,       &
            nsumx=nsumx_swd, ipol=ipol_swd, norder=norder_swd, &
            dc_bias=dc_bias_swd), stat=ios)
    end if
end select

if (select_ok) then
    if (swd % error % raised()) then
        return
    end if
else
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(file_swd)
    write(err_msg(2),'(a)')  "Unsupported combination of shp/amp/impl"
    write(err_msg(3),'(a,i0)') "shp = ", shp
    write(err_msg(4),'(a,i0)') "amp = ", amp
    write(err_msg(5),'(a,i0)') "impl = ", impl
    call swd % error % set_id_msg(err_proc, 1003, err_msg(1:5))
    return
end if

if (ios /= 0) then
    write(err_msg(1),'(a,a)') 'SWD file: ', trim(file_swd)
    write(err_msg(2),'(a)')  'Not possible to allocate actual class!'
    write(err_msg(3),'(a,i0)') "shp = ", shp
    write(err_msg(4),'(a,i0)') "amp = ", amp
    write(err_msg(5),'(a,i0)') "impl = ", impl
    allocate(spectral_wave_data_shape_1_impl_1 :: swd)  ! Empty object to store errors
    call swd % error % set_id_msg(err_proc, 1005, err_msg(1:5))
    return
end if

! SWD now has a data type and content...

contains

    subroutine swd_sniff_shp_amp()
    ! Extract shp and amp from swd file
    integer :: luswd
    integer(c_int) :: c_fmt, c_shp, c_amp
    real(c_float) :: c_magic
    err_msg(2) = ''
    call open_swd_file(newunit=luswd, file=file_swd, status='old', &
                       as_little_endian=.true., iostat=ios)
    if ( ios/=0 ) then
        err_id = 1001
        err_msg(2) = 'Not able to open SWD file.'
        return
    end if
    read(luswd) c_magic
    read(luswd) c_fmt
    if (c_fmt /= int(100, c_int)) then
        err_id = 1003
        write(err_msg(2),'(a,i0)') 'Unexpected version number of SWD-file. fmt = ', c_fmt
        return
    end if
    read(luswd) c_shp
    read(luswd) c_amp
    shp = c_shp
    amp = c_amp
    close(luswd)        
    end subroutine swd_sniff_shp_amp
    
    logical function swd_sniff_equal_nx_ny() result(res)
    ! Return true if dkx==dky and nx==ny
    integer :: luswd
    integer(c_int) :: c_fmt, c_shp, c_amp
    integer(c_int) :: c_nx, c_ny, c_order, c_nid, c_nsteps, c_nstrip
    real(c_float) :: c_dkx, c_dky, c_dt, c_grav, c_lscale, c_magic
    character(kind=c_char, len=:), allocatable :: cid
    character(kind=c_char, len=30) :: c_prog
    character(kind=c_char, len=20) :: c_date
    err_msg(2) = ''
    res = .false.
    call open_swd_file(newunit=luswd, file=file_swd, status='old', &
                       as_little_endian=.true., iostat=ios)
    if ( ios/=0 ) then
        err_id = 1001
        err_msg(2) = 'Not able to open SWD file.'
        return
    end if
    read(luswd) c_magic
    read(luswd) c_fmt
    if (c_fmt /= int(100, c_int)) then
        err_id = 1003
        write(err_msg(2),'(a,i0)') 'Unexpected version number of SWD-file. fmt = ', c_fmt
        return
    end if
    read(luswd) c_shp
    if (c_shp == 4 .or. c_shp == 5) then
        read(luswd) c_amp
        read(luswd) c_prog
        read(luswd) c_date 
        read(luswd) c_nid
        allocate(character(len=int(c_nid)) :: cid)
        read(luswd) cid
        read(luswd) c_grav
        read(luswd) c_lscale
        read(luswd) c_nstrip
        read(luswd) c_nsteps
        read(luswd) c_dt
        read(luswd) c_order
        read(luswd) c_nx
        read(luswd) c_ny
        read(luswd) c_dkx
        read(luswd) c_dky
        res = (c_nx == c_ny) .and. (abs(c_dkx - c_dky) < epsilon(c_dkx))
    end if
    close(luswd)        
    end function swd_sniff_equal_nx_ny

end subroutine spectral_wave_data_allocate

!==============================================================================

end module spectral_wave_data_allocate_def

