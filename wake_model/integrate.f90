! Subroutine to create the shape and strength distribution of the vorticity in a VAWT wake
! To be used with VAWT_Wake_Model.py

include '_quadpack.so'

! Creating the fit of the vorticity distribution
subroutine quad(func, a, b, args, full_output=0, epsabs=1.49e-8, epsrel=1.49e-8,limit=50,)
    implicit none
    integer, parameter :: dp = kind(0.d0)
    ! in
    real(dp), intent(in) :: func, a, b, args
    ! out
    real(dp), intent(out) :: gam_skew
    !local
    
    ! real(dp) :: gp
    intrinsic exp
    intrinsic erf
    intrinsic sqrt

    retval = _quadpack._qagse(func,a,b,args,full_output,epsabs,epsrel,limit)

end subroutine quad


! Calculating vorticity strength in the x and y directions
subroutine infunc()
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

    ! Limiting the fits to maximum values the Modified Gaussian Distribution can handle
    if (loc < 0.2_dp) then
        loc = 0.2_dp
    end if
    if (skw > 0.0_dp) then
       skw = 0.0_dp
    end if

    call dist(yd,loc,spr,skw,scl,g1)
    call dist(yd,-loc,-spr,-skw,-scl,g2)

    gam_lat = (g1 - g2)

end subroutine infunc

! Creating the integral for VAWT_Wake_Model.py
subroutine dblquad(y,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,inte)
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

end subroutine dblquad


! To build for Python interface: f2py -c  --opt=-O2 -m _vortmodel vorticity.f90
