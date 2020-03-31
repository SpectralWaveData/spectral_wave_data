module spectral_interpolation_def

use kind_values, only: wp => kind_swd_internal

implicit none
private

! This module provides the actual spectral interpolation scheme.
!
! Notation follows the document:
! "A Method for Handling Nonlinear Ocean Waves in Time Domain Seakeeping Programs"
! Jens B. Helmers, December 2015
!
! Written by Jens Bloch Helmers, December, 7. 2015
!
!------------------------------------------------------------------------------

!##############################################################################
!
!              B E G I N    P U B L I C    Q U A N T I T I E S
!
!------------------------------------------------------------------------------
!
public :: spectral_interpolation  ! Class for handling spectral interpolation
!
!------------------------------------------------------------------------------
!
!                E N D    P U B L I C    Q U A N T I T I E S
!
!##############################################################################

type spectral_interpolation
    integer  :: ischeme  ! 0 = 5th order interpolation
                         ! 1 = 3rd order interpolation
                         ! 2 = 7th order interpolation
    real(wp) :: delta_t  ! Spacing time samples
    procedure(scheme_interface), pointer :: scheme  ! Applied scheme
    procedure(pad_left_interface), pointer :: pad_left ! Applied scheme
    procedure(pad_right_interface), pointer :: pad_right ! Applied scheme
  contains
    procedure :: construct          ! Select actual interpolation scheme
    procedure :: scheme_1           ! 5th order spline interpolation
    procedure :: scheme_2           ! 3rd order spline interpolation
    procedure :: scheme_3           ! 7th order spline interpolation
    procedure :: pad_left_1         ! Left padding for scheme_1
    procedure :: pad_left_2         ! Left padding for scheme_2
    procedure :: pad_left_3         ! Left padding for scheme_3
    procedure :: pad_right_1        ! Right padding for scheme_1
    procedure :: pad_right_2        ! Right padding for scheme_2
    procedure :: pad_right_3        ! Right padding for scheme_3
end type spectral_interpolation

abstract interface
    pure subroutine scheme_interface(self, delta, &
                                f_im1, f_i, f_ip1, f_ip2, &
                                df_im1, df_i, df_ip1, df_ip2, &
                                f, df)
       import spectral_interpolation, wp
       class(spectral_interpolation), intent(in) :: self
       real(wp),    intent(in)  :: delta  ! (t - t_i)/delta_t
       complex(wp), intent(in)  :: f_im1  ! f(t_(i-1))
       complex(wp), intent(in)  :: f_i    ! f(t_i)
       complex(wp), intent(in)  :: f_ip1  ! f(t_(i+1))
       complex(wp), intent(in)  :: f_ip2  ! f(t_(i+2))
       complex(wp), intent(in)  :: df_im1 ! df/dt(t_(i-1))
       complex(wp), intent(in)  :: df_i   ! df/dt(t_i)
       complex(wp), intent(in)  :: df_ip1 ! df/dt(t_(i+1))
       complex(wp), intent(in)  :: df_ip2 ! df/dt(t_(i+2))
       complex(wp), intent(out) :: f     ! f(t)
       complex(wp), intent(out) :: df    ! df/dt(t)
    end subroutine scheme_interface

    pure subroutine pad_left_interface(self, f_0, f_1, f_2, &
                                       df_0, df_1, df_2, f_m1, df_m1)
       import spectral_interpolation, wp
       class(spectral_interpolation), intent(in) :: self
       complex(wp), intent(in)  :: f_0    ! f(t=0)
       complex(wp), intent(in)  :: f_1    ! f(t=dt)
       complex(wp), intent(in)  :: f_2    ! f(t=2dt)
       complex(wp), intent(in)  :: df_0   ! df(t=0)/dt
       complex(wp), intent(in)  :: df_1   ! df(t=dt)/dt
       complex(wp), intent(in)  :: df_2   ! df(t=2dt)/dt
       complex(wp), intent(out) :: f_m1   ! f(t=-dt)
       complex(wp), intent(out) :: df_m1  ! df(t=-dt)/dt
    end subroutine pad_left_interface

    pure subroutine pad_right_interface(self, f_nm2, f_nm1, f_n, &
                                      df_nm2, df_nm1, df_n, f_np1, df_np1)
       import spectral_interpolation, wp
       class(spectral_interpolation), intent(in) :: self
       complex(wp), intent(in)  :: f_nm2  ! f(t=tmax-2dt)
       complex(wp), intent(in)  :: f_nm1  ! f(t=tmax-dt)
       complex(wp), intent(in)  :: f_n    ! f(t=tmax)
       complex(wp), intent(in)  :: df_nm2 ! df(t=tmax-2dt)/dt
       complex(wp), intent(in)  :: df_nm1 ! df(t=tmax-dt)/dt
       complex(wp), intent(in)  :: df_n   ! df(t=tmax)/dt
       complex(wp), intent(out) :: f_np1  ! f(t=tmax+dt)
       complex(wp), intent(out) :: df_np1 ! df(t=tmax+dt)/dt
    end subroutine pad_right_interface

end interface

contains

!==============================================================================

subroutine construct(self, ischeme, delta_t, ierr)
class(spectral_interpolation), intent(out) :: self    ! Update data in memory (if needed)
integer, optional,             intent(in)  :: ischeme ! 0 = C^2 continous scheme
                                                      ! 1 = C^1 continous
                                                      ! 2 = C^3 continous
real(wp),                      intent(in)  :: delta_t ! Time step (constant spacing)
integer,                       intent(out) :: ierr    ! Error if ierr/=0 
self % delta_t = delta_t
if (present(ischeme)) then
    self % ischeme = ischeme
else
    self % ischeme = 0
end if
!
ierr = 0
select case(self % ischeme)
case(0)
    self % scheme => scheme_1
    self % pad_left => pad_left_1
    self % pad_right => pad_right_1
case(1)
    self % scheme => scheme_2
    self % pad_left => pad_left_2
    self % pad_right => pad_right_2
case(2)
    self % scheme => scheme_3
    self % pad_left => pad_left_3
    self % pad_right => pad_right_3
case default
    ! ischeme is out of bounds...
    ierr = 1
end select
!
end subroutine construct
    
!==============================================================================

pure subroutine scheme_1(self, delta, &
                    f_im1, f_i, f_ip1, f_ip2, &
                    df_im1, df_i, df_ip1, df_ip2, &
                    f, df)
class(spectral_interpolation), intent(in) :: self
real(wp),    intent(in)  :: delta  ! (t - t_i)/delta_t
complex(wp), intent(in)  :: f_im1  ! f(t_(i-1))
complex(wp), intent(in)  :: f_i    ! f(t_i)
complex(wp), intent(in)  :: f_ip1  ! f(t_(i+1))
complex(wp), intent(in)  :: f_ip2  ! f(t_(i+2))
complex(wp), intent(in)  :: df_im1 ! df/dt(t_(i-1))
complex(wp), intent(in)  :: df_i   ! df/dt(t_i)
complex(wp), intent(in)  :: df_ip1 ! df/dt(t_(i+1))
complex(wp), intent(in)  :: df_ip2 ! df/dt(t_(i+2))
complex(wp), intent(out) :: f     ! f(t)
complex(wp), intent(out) :: df    ! df/dt(t)
!
real(wp) :: dt4
complex(wp) :: q1, q2, q3, q4, q5
!
dt4 = 0.25_wp * self % delta_t
q1 = self % delta_t * df_i
q2 = f_im1 - 2.0_wp * f_i + f_ip1 + (df_im1 - df_ip1) * dt4
q3 = -3.0_wp * f_im1 - 3.0_wp * f_i + 5 * f_ip1 + f_ip2 - &
     (3.0_wp * df_im1 + 23.0_wp * df_i + 13.0_wp * df_ip1 + df_ip2) * dt4
q4 = 3.0_wp * f_im1 + 7.0_wp * f_i - 8.0_wp * f_ip1 - 2.0_wp * f_ip2 + &
     (3.0_wp * df_im1 + 30.0_wp * df_i + 25.0_wp * df_ip1 + 2.0_wp * df_ip2) * dt4
q5 = -f_im1 - 3.0_wp * f_i + 3.0_wp * f_ip1 + f_ip2 - &
    (df_im1 + 11.0_wp * df_i + 11.0_wp * df_ip1 + df_ip2) * dt4
!
f = f_i + (q1 + (q2  + (q3 + (q4 + q5 * delta) * delta) * delta) * delta) * delta
df = df_i + (2.0_wp * q2 + (3.0_wp * q3 + (4.0_wp * q4 + 5.0_wp * q5 * delta) &
    * delta) * delta) * delta / self % delta_t
!
end subroutine scheme_1
                    
!==============================================================================

pure subroutine scheme_2(self, delta, &
                    f_im1, f_i, f_ip1, f_ip2, &
                    df_im1, df_i, df_ip1, df_ip2, &
                    f, df)
class(spectral_interpolation), intent(in) :: self
real(wp),    intent(in)  :: delta  ! (t - t_i)/delta_t
complex(wp), intent(in)  :: f_im1  ! f(t_(i-1))
complex(wp), intent(in)  :: f_i    ! f(t_i)
complex(wp), intent(in)  :: f_ip1  ! f(t_(i+1))
complex(wp), intent(in)  :: f_ip2  ! f(t_(i+2))
complex(wp), intent(in)  :: df_im1 ! df/dt(t_(i-1))
complex(wp), intent(in)  :: df_i   ! df/dt(t_i)
complex(wp), intent(in)  :: df_ip1 ! df/dt(t_(i+1))
complex(wp), intent(in)  :: df_ip2 ! df/dt(t_(i+2))
complex(wp), intent(out) :: f     ! f(t)
complex(wp), intent(out) :: df    ! df/dt(t)
!
real(wp) :: c
complex(wp) :: a, b, d
!
a = self % delta_t * df_i - (f_ip1 - f_i)
b = - self % delta_t * df_ip1 + (f_ip1 - f_i)
c = 1.0_wp - delta
d = a * c + b * delta
!
f = c * f_i + delta * f_ip1 + delta * c * d
df = ((f_ip1 - f_i) + (1.0_wp - 2.0_wp * delta) * d +  &
     delta * c * (b - a)) / self % delta_t
!
end subroutine scheme_2
                    
!==============================================================================

pure subroutine scheme_3(self, delta, &
                    f_im1, f_i, f_ip1, f_ip2, &
                    df_im1, df_i, df_ip1, df_ip2, &
                    f, df)
class(spectral_interpolation), intent(in) :: self
real(wp),    intent(in)  :: delta  ! (t - t_i)/delta_t
complex(wp), intent(in)  :: f_im1  ! f(t_(i-1))
complex(wp), intent(in)  :: f_i    ! f(t_i)
complex(wp), intent(in)  :: f_ip1  ! f(t_(i+1))
complex(wp), intent(in)  :: f_ip2  ! f(t_(i+2))
complex(wp), intent(in)  :: df_im1 ! df/dt(t_(i-1))
complex(wp), intent(in)  :: df_i   ! df/dt(t_i)
complex(wp), intent(in)  :: df_ip1 ! df/dt(t_(i+1))
complex(wp), intent(in)  :: df_ip2 ! df/dt(t_(i+2))
complex(wp), intent(out) :: f     ! f(t)
complex(wp), intent(out) :: df    ! df/dt(t)
!
complex(wp) :: q1, q2, q3, q4, q5, q6, q7
!
q1 = self % delta_t * df_i
q2 = (df_ip1 - df_im1) * (0.25_wp * self % delta_t) 
q3 = (df_im1 - 2.0_wp * df_i + df_ip1) * (self % delta_t / 6.0_wp) 
q4 = 35.0_wp * (f_ip1 - f_i) + &
    (22.0_wp * df_im1 - 241.0_wp * df_i - 214.0_wp * df_ip1 + 13.0_wp * df_ip2) * (self % delta_t / 12.0_wp) 
q5 = 84.0_wp * (f_i - f_ip1) + &
    (-4.0_wp * df_im1 + 47.0_wp * df_i + 44.0_wp * df_ip1 - 3.0_wp * df_ip2) * self % delta_t
q6 = 70.0_wp * (f_ip1 - f_i) + &
    (37.0_wp * df_im1 - 461.0_wp * df_i - 449.0_wp * df_ip1 + 33.0_wp * df_ip2) * (self % delta_t / 12.0_wp) 
q7 = 20.0_wp * (f_i - f_ip1) + &
    (-5.0_wp * df_im1 + 65.0_wp * df_i + 65.0_wp * df_ip1 - 5.0_wp * df_ip2) * (self % delta_t / 6.0_wp)
!
f = f_i + (q1 + (q2  + (q3 + (q4 + (q5 + (q6 + q7 * delta) * &
    delta) * delta) * delta) * delta) * delta) * delta
df = df_i + (2.0_wp * q2 + (3.0_wp * q3 + (4.0_wp * q4 + (5.0_wp * q5 + (6.0_wp * q6 + 7.0_wp * q7 *  &
    delta) * delta) * delta) * delta) * delta) * delta / self % delta_t
!
end subroutine scheme_3
                    
!==============================================================================

pure subroutine pad_left_1(self, f_0, f_1, f_2, df_0, df_1, df_2, f_m1, df_m1)
class(spectral_interpolation), intent(in) :: self
complex(wp), intent(in)  :: f_0    ! f(t=0)
complex(wp), intent(in)  :: f_1    ! f(t=dt)
complex(wp), intent(in)  :: f_2    ! f(t=2dt)
complex(wp), intent(in)  :: df_0   ! df(t=0)/dt
complex(wp), intent(in)  :: df_1   ! df(t=dt)/dt
complex(wp), intent(in)  :: df_2   ! df(t=2dt)/dt
complex(wp), intent(out) :: f_m1   ! f(t=-dt)
complex(wp), intent(out) :: df_m1  ! d(f(t=-dt))/dt
!
f_m1 = f_0 + (df_1 - 3.0_wp * df_0) * self % delta_t * 0.5_wp
df_m1 = 2.0_wp * df_0 - df_1
!
end subroutine pad_left_1

!==============================================================================

pure subroutine pad_left_2(self, f_0, f_1, f_2, df_0, df_1, df_2, f_m1, df_m1)
class(spectral_interpolation), intent(in) :: self
complex(wp), intent(in)  :: f_0    ! f(t=0)
complex(wp), intent(in)  :: f_1    ! f(t=dt)
complex(wp), intent(in)  :: f_2    ! f(t=2dt)
complex(wp), intent(in)  :: df_0   ! df(t=0)/dt
complex(wp), intent(in)  :: df_1   ! df(t=dt)/dt
complex(wp), intent(in)  :: df_2   ! df(t=2dt)/dt
complex(wp), intent(out) :: f_m1   ! f(t=-dt)
complex(wp), intent(out) :: df_m1  ! d(f(t=-dt))/dt
!
f_m1 = f_0 + (df_1 - 3.0_wp * df_0) * self % delta_t * 0.5_wp
df_m1 = 2.0_wp * df_0 - df_1

! Legacy scheme: quadratic extrapolation of first 3 time steps
!f_m1 = 3.0_wp * (f_0 - f_1) + f_2
!df_m1 = 3.0_wp * (df_0 - df_1) + df_2
!
end subroutine pad_left_2

!==============================================================================

pure subroutine pad_left_3(self, f_0, f_1, f_2, df_0, df_1, df_2, f_m1, df_m1)
class(spectral_interpolation), intent(in) :: self
complex(wp), intent(in)  :: f_0    ! f(t=0)
complex(wp), intent(in)  :: f_1    ! f(t=dt)
complex(wp), intent(in)  :: f_2    ! f(t=2dt)
complex(wp), intent(in)  :: df_0   ! df(t=0)/dt
complex(wp), intent(in)  :: df_1   ! df(t=dt)/dt
complex(wp), intent(in)  :: df_2   ! df(t=2dt)/dt
complex(wp), intent(out) :: f_m1   ! f(t=-dt)
complex(wp), intent(out) :: df_m1  ! d(f(t=-dt))/dt
!
f_m1 = f_0 + (df_1 - 3.0_wp * df_0) * self % delta_t * 0.5_wp
df_m1 = 2.0_wp * df_0 - df_1

! Legacy scheme: quadratic extrapolation of first 3 time steps
!f_m1 = 3.0_wp * (f_0 - f_1) + f_2
!df_m1 = 3.0_wp * (df_0 - df_1) + df_2
!
end subroutine pad_left_3

!==============================================================================

pure subroutine pad_right_1(self, f_nm2, f_nm1, f_n, &
                            df_nm2, df_nm1, df_n, f_np1, df_np1)
class(spectral_interpolation), intent(in) :: self
complex(wp), intent(in)  :: f_nm2  ! f(t=tmax-2dt)
complex(wp), intent(in)  :: f_nm1  ! f(t=tmax-dt)
complex(wp), intent(in)  :: f_n    ! f(t=tmax)
complex(wp), intent(in)  :: df_nm2 ! df(t=tmax-2dt)/dt
complex(wp), intent(in)  :: df_nm1 ! df(t=tmax-dt)/dt
complex(wp), intent(in)  :: df_n   ! df(t=tmax)/dt
complex(wp), intent(out) :: f_np1  ! f(t=tmax+dt)
complex(wp), intent(out) :: df_np1 ! df(t=tmax+dt)/dt
!
f_np1 = f_n - (df_nm1 - 3.0_wp * df_n) * self % delta_t * 0.5_wp
df_np1 = 2.0_wp * df_n - df_nm1
!
end subroutine pad_right_1

!==============================================================================

pure subroutine pad_right_2(self, f_nm2, f_nm1, f_n, &
                            df_nm2, df_nm1, df_n, f_np1, df_np1)
class(spectral_interpolation), intent(in) :: self
complex(wp), intent(in)  :: f_nm2  ! f(t=tmax-2dt)
complex(wp), intent(in)  :: f_nm1  ! f(t=tmax-dt)
complex(wp), intent(in)  :: f_n    ! f(t=tmax)
complex(wp), intent(in)  :: df_nm2 ! df(t=tmax-2dt)/dt
complex(wp), intent(in)  :: df_nm1 ! df(t=tmax-dt)/dt
complex(wp), intent(in)  :: df_n   ! df(t=tmax)/dt
complex(wp), intent(out) :: f_np1  ! f(t=tmax+dt)
complex(wp), intent(out) :: df_np1 ! df(t=tmax+dt)/dt
!
f_np1 = f_n - (df_nm1 - 3.0_wp * df_n) * self % delta_t * 0.5_wp
df_np1 = 2.0_wp * df_n - df_nm1

! Legacy scheme: quadratic extrapolation of last 3 time steps
!f_np1 = 3.0_wp * (f_n - f_nm1) + f_nm2
!df_np1 = 3.0_wp * (df_n - df_nm1) + df_nm2
!
end subroutine pad_right_2

!==============================================================================

pure subroutine pad_right_3(self, f_nm2, f_nm1, f_n, &
                            df_nm2, df_nm1, df_n, f_np1, df_np1)
class(spectral_interpolation), intent(in) :: self
complex(wp), intent(in)  :: f_nm2  ! f(t=tmax-2dt)
complex(wp), intent(in)  :: f_nm1  ! f(t=tmax-dt)
complex(wp), intent(in)  :: f_n    ! f(t=tmax)
complex(wp), intent(in)  :: df_nm2 ! df(t=tmax-2dt)/dt
complex(wp), intent(in)  :: df_nm1 ! df(t=tmax-dt)/dt
complex(wp), intent(in)  :: df_n   ! df(t=tmax)/dt
complex(wp), intent(out) :: f_np1  ! f(t=tmax+dt)
complex(wp), intent(out) :: df_np1 ! df(t=tmax+dt)/dt
!
f_np1 = f_n - (df_nm1 - 3.0_wp * df_n) * self % delta_t * 0.5_wp
df_np1 = 2.0_wp * df_n - df_nm1

! Legacy scheme: quadratic extrapolation of last 3 time steps
!f_np1 = 3.0_wp * (f_n - f_nm1) + f_nm2
!df_np1 = 3.0_wp * (df_n - df_nm1) + df_nm2
!
end subroutine pad_right_3

!==============================================================================

end module spectral_interpolation_def
