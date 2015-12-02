!Subroutine to create the shape and strength distribution of the vorticity in a VAWT wake
!To be used with VAWT Wake Model.py


!Creating the fit of the vorticity distribution
subroutine skew(y,e,w,a,i,gam_skew)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: y,e,w,a,i

    ! out
    real(dp), intent(out) :: gam_skew

    !local
    ! real(dp) :: gp
    intrinsic exp
    intrinsic erf
    intrinsic sqrt

    ! Exponentially Modified Gaussian Distribution
    gam_skew = i*a/2.0_dp*exp(a/2.0_dp*(2.0_dp*e+a*w**2.0_dp-2.0_dp*y))*(1.0_dp-erf((e + a*w**2.0_dp -y)/(sqrt(2.0_dp)*w)))

end subroutine skew


!Calculating vorticity strength in the s and t directions
subroutine gamma(x,y,dia,e1,e2,e3,w1,w2,a1,a2,i1,i2,i3,gam_lat)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x,y,dia
    real(dp), intent(in) :: e1,e2,e3,w1,w2,a1,a2
    real(dp), intent(in) :: i1,i2,i3

    ! out
    real(dp), intent(out) :: gam_lat

    ! local
    real(dp) :: xd,yd,e,w,a,i,g1,g2,gam
    intrinsic exp

    xd = x/dia ! normalizing s by the diameter
    yd = y/dia ! normalizing t by the diameter

    e = e1*xd**2 + e2*xd + e3
    w = w1*xd + w2
    a = a1*xd + a2
    i = i1/(1.0_dp + exp(i2*(xd - i3)))

    ! Limiting the fits to maximum values the Weibull distribution can handle
    if (e < 0.0_dp) then
        e = 0.0_dp
    end if
    if (a > 0.0_dp) then
       a = 0.0_dp
    end if


    call skew(yd,e,w,a,i,g1)
    call skew(yd,-e,-w,-a,-i,g2)

    gam = (g1 - g2)

    gam_lat = gam ! readjusting the vorticity strength to remove normalization

end subroutine gamma

!Creating the integral for VAWT Wake Model.py
subroutine integrand(y,x,x0,y0,dia,e1,e2,e3,w1,w2,a1,a2,i1,i2,i3,inte)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: y,x,x0,y0,dia
    real(dp), intent(in) :: e1,e2,e3,w1,w2,a1,a2
    real(dp), intent(in) :: i1,i2,i3

    ! out
    real(dp), intent(out) :: inte

    ! local
    real(dp) :: gammav


    !Specifying the strength of the vorticity
    call gamma(x,y,dia,e1,e2,e3,w1,w2,a1,a2,i1,i2,i3,gammav)

    ! inte = gammav*((td - yd)/((xd - xd)**2 + (td - yd)**2))
    inte = gammav*((y - y0)/((x - x0)**2 + (y - y0)**2))

end subroutine integrand


! To build for Python interface: f2py -c  --opt=-O2 -m _vortrun vorticity.f90
