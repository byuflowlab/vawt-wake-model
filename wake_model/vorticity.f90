! Subroutine to create the shape and strength distribution of the vorticity in a VAWT wake
! To be used with VAWT_Wake_Model.py


! Creating the fit of the vorticity distribution
subroutine dist(y,loc,spr,skw,scl,gam_skew)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: y,loc,spr,skw,scl

    ! out
    real(dp), intent(out) :: gam_skew

    !local
    ! real(dp) :: gp
    intrinsic exp
    intrinsic erf
    intrinsic sqrt

    ! Exponentially Modified Gaussian Distribution
    gam_skew = scl*skw/2.0_dp*exp(skw/2.0_dp*(2.0_dp*loc+skw*spr**2.0_dp-2.0_dp*y))&
    *(1.0_dp-erf((loc + skw*spr**2.0_dp -y)/(sqrt(2.0_dp)*spr)))

end subroutine dist


! Calculating vorticity strength in the x and y directions
subroutine gamma(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gam_lat)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x,y,dia
    real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), intent(out) :: gam_lat

    ! local
    real(dp) :: xd,yd,loc,spr,skw,scl,g1,g2
    intrinsic exp

    xd = x/dia ! normalizing x by the diameter
    yd = y/dia ! normalizing y by the diameter

    loc = loc1*xd**2 + loc2*xd + loc3
    spr = spr1*xd + spr2
    skw = skw1*xd + skw2
    scl = scl1/(1.0_dp + exp(scl2*(xd - scl3)))

    ! Limiting the fits to maximum values the Weibull distribution can handle
    if (loc < 0.0_dp) then
        loc = 0.0_dp
    end if
    if (skw > 0.0_dp) then
       skw = 0.0_dp
    end if

    call dist(yd,loc,spr,skw,scl,g1)
    call dist(yd,-loc,-spr,-skw,-scl,g2)

    gam_lat = (g1 - g2)

end subroutine gamma

! Creating the integral for VAWT_Wake_Model.py
subroutine integrand(y,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,inte)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: y,x,x0,y0,dia
    real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), intent(out) :: inte

    ! local
    real(dp) :: gammav

    ! Specifying the strength of the vorticity
    call gamma(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gammav)

    inte = gammav*((y - y0)/((x - x0)**2 + (y - y0)**2))

end subroutine integrand


! To build for Python interface: f2py -c  --opt=-O2 -m _vortmodel vorticity.f90
