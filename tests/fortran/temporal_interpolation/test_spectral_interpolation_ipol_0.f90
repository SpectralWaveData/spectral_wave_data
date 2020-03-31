program test_spectral_interpolation_ipol_0
    
use kind_values, only:  wp => kind_swd_internal
use spectral_interpolation_def, only: spectral_interpolation

implicit none

! This program tests the ipol=0 temporal scheme

! Dummy polynomials: f(t) = a0 + a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5
complex(wp), parameter :: a0 = cmplx(0.156_wp, 0.17_wp, wp)
complex(wp), parameter :: a1 = cmplx(1.234_wp, -0.54_wp, wp)
complex(wp), parameter :: a2 = cmplx(-2.34_wp, 0.154_wp, wp)
complex(wp), parameter :: a3 = cmplx(-0.18_wp, -0.25_wp, wp)
complex(wp), parameter :: a4 = cmplx(0.58_wp, -0.74_wp, wp)
complex(wp), parameter :: a5 = cmplx(0.23_wp, 0.09_wp, wp)

real(wp), parameter :: delta_t = 0.34_wp
real(wp), parameter :: t_im1 = 0.15_wp
real(wp), parameter :: t_i = t_im1 + delta_t
real(wp), parameter :: t_ip1 = t_i + delta_t
real(wp), parameter :: t_ip2 = t_ip1 + delta_t

type(spectral_interpolation) :: si
real(wp) :: t, delta
complex(wp) :: c1, dc1, ddc1, c2, dc2, ddc2
complex(wp) :: f_im1, f_i, f_ip1, f_ip2
complex(wp) :: df_im1, df_i, df_ip1, df_ip2
complex(wp) :: ddf_i, ddf_ip1

integer :: ipol, i, ierr, n_error = 0
real(wp), parameter :: deltas(5) = [0.0_wp, 0.29_wp, 0.64_wp, 0.84_wp, 1.0_wp]

print*, " Test of ipol=0..."
ipol = 0
call si % construct(ipol, delta_t, ierr)
if (ierr /= 0) then
    print *, "Unknown interpolation scheme: ", ipol
    stop
end if

f_im1 = f5(t_im1)
f_i = f5(t_i)
f_ip1 = f5(t_ip1)
f_ip2 = f5(t_ip2)
df_im1 = df5(t_im1)
df_i = df5(t_i)
df_ip1 = df5(t_ip1)
df_ip2 = df5(t_ip2)
! Adjust f_im1 make ddf_i based on finite-difference match polynomial at t_i
ddf_i = ddf5(t_i)
f_im1 = 0.5_wp * ddf_i * delta_t**2 + 2 * f_i - f_ip1 + (df_ip1 - df_im1) * delta_t / 4
! Adjust f_ip2 to make ddf_ip1 based on finite-difference match polynomial at t_ip1
ddf_ip1 = ddf5(t_ip1)
f_ip2 = 0.5_wp * ddf_ip1 * delta_t**2 + 2 * f_ip1 - f_i + (df_ip2 - df_i) * delta_t / 4

do i = 1, size(deltas)
    delta = deltas(i)
    t = t_i + delta * delta_t
    c1 = f5(t)
    dc1 = df5(t)
    ddc1 = ddf5(t)
    call si % scheme(delta, &
                     f_im1, f_i, f_ip1, f_ip2, &
                     df_im1, df_i, df_ip1, df_ip2, &
                     c2, dc2)
    print*, "Delta = ", delta
    print*, "Compare function values: ", c1, c2, compare(c1, c2)
    print*, "Compare function derivatives: ", dc1, dc2, compare(dc1, dc2)
    if (i==1) then
        print*, "f_i, df_i: ", f_i, df_i
        print*, "compare(c2, f_i)", compare(c2, f_i)
        print*, "compare(dc2, df_i)", compare(dc2, df_i)
        print*, "compare(ddc1, ddf_i)", compare(ddc1, ddf_i)
    else if (i==size(deltas)) then
        print*, "f_ip1, df_ip1: ", f_ip1, df_ip1
        print*, "compare(c2, f_ip1)", compare(c2, f_ip1)
        print*, "compare(dc2, df_ip1)", compare(dc2, df_ip1)
        print*, "compare(ddc1, ddf_ip1)", compare(ddc1, ddf_ip1)
    end if
end do

print*, "Number of detected ERRORs = ", n_error
  
contains

function f5(t) result(res)
real(wp), intent(in) :: t
complex(wp) res
res = a0 + a1*t + a2*t**2 + a3*t**3 + a4*t**4 + a5*t**5
end function f5

function df5(t) result(res)
real(wp), intent(in) :: t
complex(wp) res
res = a1 + 2.0_wp * a2*t + 3.0_wp * a3*t**2 + 4.0_wp * a4*t**3 + 5.0_wp * a5*t**4
end function df5

function ddf5(t) result(res)
real(wp), intent(in) :: t
complex(wp) res
res = 2.0_wp * a2 + 6.0_wp * a3*t + 12.0_wp * a4*t**2 + 20.0_wp * a5*t**3
end function ddf5

function compare(c1, c2) result(res)
complex(wp), intent(in) :: c1, c2
character(len=10) res
real(wp), parameter :: eps = 100*epsilon(eps)
if (abs(c1 % re - c2 % re) > eps * max(abs(c1 % re), abs(c2 % re)) .or. &
    abs(c1 % im - c2 % im) > eps * max(abs(c1 % im), abs(c2 % im))) then
    res = "ERROR!"
    n_error = n_error + 1
else
    res = "OK"
end if
end function compare

end program test_spectral_interpolation_ipol_0
