! Parameterized VAWT Wake Model Fortran Routines
! Developed by Eric Tingey at Brigham Young University
! This code models the wake behind a vertical-axis wind turbine based on
! tip-speed ratio, solidity and wind speed by converting the vorticity of
! the wake into velocity information. The model uses CFD data obtained
! from STAR-CCM+ turbine simulations serve as the basis of the initial
! wake model.
! Only valid for tip-speed ratios between 1.5 and 7.0 and solidities between
! 0.15 and 1.0. Reynolds numbers should also be around the range of 600,000 to
! 6,000,000.
! In this code, the x and y coordinates (looking down on the turbine) are
! made according to:
! --------------->--------------------------------------------------------
! --------------->--------------------------------------------------------
! --------------->---------=====--------#################-----------Y-----
! --------------->------//       \\#############################----|-----
! -FREE-STREAM--->-----|| TURBINE ||########## WAKE ###############-|___X-
! ----WIND------->-----||         ||###############################-------
! --------------->------\\       //#############################----------
! --------------->---------=====--------#################-----------------
! --------------->--------------------------------------------------------
! --------------->--------------------------------------------------------
! The imported vorticity data also assumes symmetry in the wake and therefore
! rotation direction for wake velocity calculation is irrelevant.


! trapezoidal integration
subroutine trapz(n,x,y,integral) ! integrate y with respect to x
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: x,y

    ! out
    real(dp), intent(out) :: integral

    ! local
    integer :: i
    real(dp) :: inte

    inte = 0.0_dp
    do i = 1,n-1
      inte = inte + (x(i+1)-x(i))*0.5_dp*(y(i) + y(i+1))
    end do

    integral = inte

end subroutine trapz


! integration for a periodic function where end points don't reach ends (uses trapezoidal method)
subroutine pInt(n,theta,f,integral)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: theta,f

    ! out
    real(dp), intent(out) :: integral

    ! local
    real(dp) :: inte,dtheta

    call trapz(n,theta,f,inte)

    ! add end points
    dtheta = 2.0_dp*theta(1)  ! assumes equally spaced, starts at 0
    integral = inte + dtheta * 0.5_dp*(f(1) + f(n))

end subroutine pInt


! linear interpolation (specifically for extracting airfoil data)
subroutine interpolate(n,x,y,xval,yval)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: x,y
    real(dp), intent(in) :: xval

    ! out
    real(dp), intent(out) :: yval

    ! local
    integer :: i

    ! assuming the values of x are in accending order
    do i = 1,n
      if (xval < x(i)) then
        yval = y(i-1) + (xval - x(i-1))*((y(i)-y(i-1))/(x(i)-x(i-1)))
        exit
      else if (xval == x(i)) then
        yval = y(i)
        exit
      end if

    end do

end subroutine interpolate


! cubic spline interpolation setup (specifically for extracting airfoil data)
subroutine splineint(n,x,y,xval,yval)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n
    real(dp), dimension(n), intent(in) :: x,y
    real(dp), intent(in) :: xval

    ! out
    real(dp), intent(out) :: yval

    ! local
    integer :: i
    real(dp) :: x1,x2,x3,y1,y2,y3

    ! assuming the values of x are in accending order
    do i = 1,n
      if (xval < x(i)) then ! check that given x value is below current x point
        if (i == 2) then ! x value is at the beginning of the data set
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)
          y1 = y(1)
          y2 = y(2)
          y3 = y(3)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else if (i == n) then ! x value is at the end of the data set
          x1 = x(n-2)
          x2 = x(n-1)
          x3 = x(n)
          y1 = y(n-2)
          y2 = y(n-1)
          y3 = y(n)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else
          if (xval <= (x(i)+x(i-1))/2.0_dp) then ! interpolate on beginning half
            x1 = x(i-2)
            x2 = x(i-1)
            x3 = x(i)
            y1 = y(i-2)
            y2 = y(i-1)
            y3 = y(i)
            call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
          else ! interpolate on ending half
            x1 = x(i-1)
            x2 = x(i)
            x3 = x(i+1)
            y1 = y(i-1)
            y2 = y(i)
            y3 = y(i+1)
            call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
          end if
        end if
        exit
      else if (xval == x(i)) then ! no interpolation needed for value in data set
        yval = y(i)
        exit
      end if

    end do

end subroutine splineint


! cubic spline interpolation (specifically for extracting airfoil data)
subroutine cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x1,x2,x3,y1,y2,y3,xval

    ! out
    real(dp), intent(out) :: yval

    ! local
    real(dp) :: a11,a12,a13,a21,a22,a23,a31,a32,a33,b1,b2,b3
    real(dp) :: bot,xtop,ytop,ztop,k1,k2,k3,a,b,t

    ! solving tridiagonal linear equation system
    a11 = 2.0_dp/(x2-x1)
    a12 = 1.0_dp/(x2-x1)
    a13 = 0.0_dp
    a21 = 1.0_dp/(x2-x1)
    a22 = 2.0_dp*((1.0_dp/(x2-x1))+(1.0_dp/(x3-x2)))
    a23 = 1.0_dp/(x3-x2)
    a31 = 0.0_dp
    a32 = 1.0_dp/(x3-x2)
    a33 = 2.0_dp/(x3-x2)
    b1 = 3.0_dp*(y2-y1)/(x2-x1)**2
    b2 = 3.0_dp*(((y2-y1)/(x2-x1)**2)+((y3-y2)/(x3-x2)**2))
    b3 = 3.0_dp*(y3-y2)/(x3-x2)**2

    ! solving using inverse matrix method
    bot = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32
    if (xval < x2) then
      xtop = b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - a13*a22*b3 - a12*b2*a33 - b1*a23*a32
      ytop = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - a13*b2*a31 - b1*a21*a33 - a11*a23*b3

      k1 = xtop/bot
      k2 = ytop/bot

      a = k1*(x2-x1) - (y2-y1)
      b = -k2*(x2-x1) + (y2-y1)
      t = (xval-x1)/(x2-x1)

      yval = (1.0_dp - t)*y1 + t*y2 + t*(1.0_dp - t)*(a*(1.0_dp - t) + b*t)
    else
      ytop = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - a13*b2*a31 - b1*a21*a33 - a11*a23*b3
      ztop = a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - b1*a22*a31 - a12*a21*b3 - a11*b2*a32

      k2 = ytop/bot
      k3 = ztop/bot

      a = k2*(x3-x2) - (y3-y2)
      b = -k3*(x3-x2) + (y3-y2)
      t = (xval-x2)/(x3-x2)

      yval = (1.0_dp - t)*y2 + t*y3 + t*(1.0_dp - t)*(a*(1.0_dp - t) + b*t)
    end if

end subroutine cubspline


! Calculating EMG parameter values based on given polynomial surface
subroutine parameterval(tsr,sol,coef,val)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: tsr,sol
    real(dp), dimension(10), intent(in) :: coef

    ! out
    real(dp), intent(out) :: val

    ! local
    real(dp) :: a,b,c,d,e,f,g,h,i,j

    a = coef(1)
    b = coef(2)
    c = coef(3)
    d = coef(4)
    e = coef(5)
    f = coef(6)
    g = coef(7)
    h = coef(8)
    i = coef(9)
    j = coef(10)

    val = a + b*tsr + c*sol + d*tsr**2 + e*tsr*sol + f*sol**2 + g*tsr**3 + &
    h*tsr**2*sol + i*tsr*sol**2 + j*sol**3

end subroutine parameterval


! Creating the EMG fit of the vorticity distribution
subroutine EMGdist(y,loc,spr,skw,scl,gam_skew)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: y,loc,spr,skw,scl

    ! out
    real(dp), intent(out) :: gam_skew

    ! local
    intrinsic exp
    intrinsic erf
    intrinsic sqrt

    ! Exponentially Modified Gaussian Distribution
    gam_skew = scl*skw/2.0_dp*exp(skw/2.0_dp*(2.0_dp*loc+skw*spr**2.0_dp-2.0_dp*y))&
    *(1.0_dp-erf((loc + skw*spr**2.0_dp - y)/(sqrt(2.0_dp)*spr)))

end subroutine EMGdist


! Calculating vorticity strength in the x and y directions
subroutine vorticitystrength(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gam_lat)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x,y,dia
    real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), intent(out) :: gam_lat

    ! local
    real(dp) :: loc1d,loc2d,loc3d,spr1d,spr2d,skw1d,skw2d,scl1d,scl2d,scl3d
    real(dp) :: xd,yd,loc,spr,skw,scl,g1,g2
    intrinsic exp

    xd = x/dia ! normalizing x by the diameter
    yd = y/dia ! normalizing y by the diameter

    ! Limiting the parameter components to create expected behavior
    if (loc1 > -0.001_dp) then ! ensure concave down
      loc1d = -0.001_dp
    else
      loc1d = loc1
    end if
    if (loc2 < 0.01_dp) then ! ensure slight increase moving downstream
      loc2d = 0.01_dp
    else
      loc2d = loc2
    end if
    if (loc3 < 0.48_dp) then ! ensure wake originating from edge of turbine
      loc3d = 0.48_dp
    else
      loc3d = loc3
    end if

    loc = loc1d*xd*xd + loc2d*xd + loc3d ! EMG Location

    if (spr1 > -0.001_dp) then ! ensure decrease in value (more spread downstream)
      spr1d = -0.001_dp
    else
      spr1d = spr1
    end if
    if (spr2 > 0.0_dp) then ! ensure value does not begin positive
      spr2d = 0.0_dp
    else
      spr2d = spr2
    end if

    spr = spr1d*xd + spr2d ! EMG Spread

    skw1d = skw1 ! no limitations necessary
    if (skw2 > 0.0_dp) then ! ensure value does not begin positive
      skw2d = 0.0_dp
    else
      skw2d = skw2
    end if

    skw = skw1d*xd + skw2d ! EMG Skew

    if (scl1 < 0.0_dp) then ! ensure positive maximum vorticity strength
      scl1d = 0.0_dp
    else
      scl1d = scl1
    end if
    if (scl2 < 0.05_dp) then ! ensure decay moving downstream
      scl2d = 0.05_dp
    else
      scl2d = scl2
    end if
    if (scl3 < 0.0_dp) then ! ensure decay occurs downstream
      scl3d = 0.0_dp
    else
      scl3d = scl3
    end if

    scl = scl1d/(1.0_dp + exp(scl2d*(xd - scl3d))) ! EMG Scale

    ! Limiting the parameters to the maximum values the EMG distribution can handle
    if (loc < 0.2_dp) then
      loc = 0.2_dp
    end if
    if (spr < -0.5_dp) then
      spr = -0.5_dp
    else if (spr > -0.001_dp) then
      spr = -0.001_dp
    end if
    if (skw > 0.0_dp) then
      skw = 0.0_dp
    end if

    call EMGdist(yd,loc,spr,skw,scl,g1)
    call EMGdist(yd,-loc,-spr,-skw,-scl,g2)

    gam_lat = (g1 - g2)

end subroutine vorticitystrength


! Calculating vorticity strength in the x and y directions
subroutine vorticitystrengthx(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gam_lat)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x,y,dia
    real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), intent(out) :: gam_lat

    ! local
    real(dp) :: pi,loc1d,loc2d,loc3d,spr1d,spr2d,skw1d,skw2d,scl1d,scl2d,scl3d
    real(dp) :: xd,yd,loc,spr,skw,scl
    intrinsic exp
    intrinsic erf
    intrinsic sqrt
    pi = 3.1415926535897932_dp

    xd = x/dia ! normalizing x by the diameter
    yd = y/dia ! normalizing y by the diameter

    ! Limiting the parameter components to create expected behavior
    if (loc1 > -0.001_dp) then ! ensure concave down
      loc1d = -0.001_dp
    else
      loc1d = loc1
    end if
    if (loc2 < 0.01_dp) then ! ensure slight increase moving downstream
      loc2d = 0.01_dp
    else
      loc2d = loc2
    end if
    if (loc3 < 0.48_dp) then ! ensure wake originating from edge of turbine
      loc3d = 0.48_dp
    else
      loc3d = loc3
    end if

    loc = loc1d*xd*xd + loc2d*xd + loc3d ! EMG Location

    if (spr1 > -0.001_dp) then ! ensure decrease in value (more spread downstream)
      spr1d = -0.001_dp
    else
      spr1d = spr1
    end if
    if (spr2 > 0.0_dp) then ! ensure value does not begin positive
      spr2d = 0.0_dp
    else
      spr2d = spr2
    end if

    spr = spr1d*xd + spr2d ! EMG Spread

    skw1d = skw1 ! no limitations necessary
    if (skw2 > 0.0_dp) then ! ensure value does not begin positive
      skw2d = 0.0_dp
    else
      skw2d = skw2
    end if

    skw = skw1d*xd + skw2d ! EMG Skew

    if (scl1 < 0.0_dp) then ! ensure positive maximum vorticity strength
      scl1d = 0.0_dp
    else
      scl1d = scl1
    end if
    if (scl2 < 0.05_dp) then ! ensure decay moving downstream
      scl2d = 0.05_dp
    else
      scl2d = scl2
    end if
    if (scl3 < 0.0_dp) then ! ensure decay occurs downstream
      scl3d = 0.0_dp
    else
      scl3d = scl3
    end if

    scl = scl1d/(1.0_dp + exp(scl2d*(xd - scl3d))) ! EMG Scale

    ! Limiting the parameters to the maximum values the EMG distribution can handle
    if (loc < 0.2_dp) then
      loc = 0.2_dp
    end if
    if (spr < -0.5_dp) then
      spr = -0.5_dp
    else if (spr > -0.001_dp) then
      spr = -0.001_dp
    end if
    if (skw > 0.0_dp) then
      skw = 0.0_dp
    end if

    ! Exponentially Modified Gaussian Distribution
    gam_lat = (1.0_dp/(2.0_dp*spr))*scl*skw*(exp(-(loc-yd)**2/(2.0_dp*spr**2))*&
    sqrt(2.0_dp/pi) + exp(-(loc+yd)**2/(2.0_dp*spr**2))*sqrt(2.0_dp/pi) + &
    exp(0.5_dp*skw*(2.0_dp*loc + skw*spr**2 - 2.0_dp*y)*skw*spr*(-1.0_dp + &
    erf((loc + skw*spr**2 - yd)/(sqrt(2.0_dp)*spr)))) + exp(0.5_dp*skw*(2.0_dp*&
    loc + skw*spr**2 + 2.0_dp*y)*skw*spr*(-1.0_dp + erf((loc + skw*spr**2 + &
    yd)/(sqrt(2.0_dp)*spr)))))


end subroutine vorticitystrengthx


! Calculating vorticity strength in the x and y directions
subroutine vorticitystrengthy(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gam_lat)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x,y,dia
    real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), intent(out) :: gam_lat

    ! local
    real(dp) :: loc1d,loc2d,loc3d,spr1d,spr2d,skw1d,skw2d,scl1d,scl2d,scl3d
    real(dp) :: xd,yd,loc,spr,skw,scl,g1,g2
    intrinsic exp
    intrinsic erf
    intrinsic sqrt

    xd = x/dia ! normalizing x by the diameter
    yd = y/dia ! normalizing y by the diameter

    ! Limiting the parameter components to create expected behavior
    if (loc1 > -0.001_dp) then ! ensure concave down
      loc1d = -0.001_dp
    else
      loc1d = loc1
    end if
    if (loc2 < 0.01_dp) then ! ensure slight increase moving downstream
      loc2d = 0.01_dp
    else
      loc2d = loc2
    end if
    if (loc3 < 0.48_dp) then ! ensure wake originating from edge of turbine
      loc3d = 0.48_dp
    else
      loc3d = loc3
    end if

    loc = loc1d*xd*xd + loc2d*xd + loc3d ! EMG Location

    if (spr1 > -0.001_dp) then ! ensure decrease in value (more spread downstream)
      spr1d = -0.001_dp
    else
      spr1d = spr1
    end if
    if (spr2 > 0.0_dp) then ! ensure value does not begin positive
      spr2d = 0.0_dp
    else
      spr2d = spr2
    end if

    spr = spr1d*xd + spr2d ! EMG Spread

    skw1d = skw1 ! no limitations necessary
    if (skw2 > 0.0_dp) then ! ensure value does not begin positive
      skw2d = 0.0_dp
    else
      skw2d = skw2
    end if

    skw = skw1d*xd + skw2d ! EMG Skew

    if (scl1 < 0.0_dp) then ! ensure positive maximum vorticity strength
      scl1d = 0.0_dp
    else
      scl1d = scl1
    end if
    if (scl2 < 0.05_dp) then ! ensure decay moving downstream
      scl2d = 0.05_dp
    else
      scl2d = scl2
    end if
    if (scl3 < 0.0_dp) then ! ensure decay occurs downstream
      scl3d = 0.0_dp
    else
      scl3d = scl3
    end if

    scl = scl1d/(1.0_dp + exp(scl2d*(xd - scl3d))) ! EMG Scale

    ! Limiting the parameters to the maximum values the EMG distribution can handle
    if (loc < 0.2_dp) then
      loc = 0.2_dp
    end if
    if (spr < -0.5_dp) then
      spr = -0.5_dp
    else if (spr > -0.001_dp) then
      spr = -0.001_dp
    end if
    if (skw > 0.0_dp) then
      skw = 0.0_dp
    end if

    call EMGdist(yd,loc,spr,skw,scl,g1)
    call EMGdist(yd,-loc,-spr,-skw,-scl,g2)

    gam_lat = (g1 - g2)

end subroutine vorticitystrengthy


! Creating the integral for induced x-velocity
subroutine integrandx(y,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,inte)
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
    call vorticitystrength(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gammav)

    inte = gammav*((y - y0)/((x - x0)*(x - x0) + (y - y0)*(y - y0)))

end subroutine integrandx


! Creating the integral for induced y-velocity
subroutine integrandy(y,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,inte)
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
    call vorticitystrength(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gammav)

    inte = gammav*((x0 - x)/((x - x0)*(x - x0) + (y - y0)*(y - y0)))

end subroutine integrandy


! Performing integration to convert vorticity into velocity
subroutine vel_field(xt,yt,x0t,y0t,dia,rot,chord,blades,Vinf,loc1d,loc2d,loc3d,spr1d,spr2d,&
  skw1d,skw2d,scl1d,scl2d,scl3d,m_in,n_in,inte,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: xt,yt,x0t,y0t,dia,rot,chord,Vinf
    real(dp), dimension(10), intent(in) :: loc1d,loc2d,loc3d,spr1d,spr2d,skw1d,skw2d,scl1d,scl2d,scl3d
    integer, intent(in) :: blades,m_in,n_in,inte

    ! out
    real(dp), intent(out) :: velx,vely

    ! local
    integer :: i,j,m,n
    real(dp) :: x0,y0,h,k,velxi,velyi,xsum,ysum,intlim,pi2,tsr,sol,a,b,c,d
    real(dp) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
    real(dp), dimension(m_in) :: xdiv
    real(dp), dimension(n_in) :: ydiv
    real(dp) :: xval1,xval2,xval3,xval4,yval1,yval2,yval3,yval4
    real(dp) :: xsum1,xsum2,xsum3,xsum4,xsum5,xsum6,xsum7,xsum8,xsum9,xsum10,xsum11,xsum12
    real(dp) :: ysum1,ysum2,ysum3,ysum4,ysum5,ysum6,ysum7,ysum8,ysum9,ysum10,ysum11,ysum12
    intrinsic sqrt
    intrinsic abs
    intrinsic exp
    pi2 = 6.28318530718_dp

    tsr = (dia/2.0_dp)*abs(rot)/Vinf
    sol = blades*chord/(dia/2.0_dp)

    call parameterval(tsr,sol,loc1d,loc1)
    call parameterval(tsr,sol,loc2d,loc2)
    call parameterval(tsr,sol,loc3d,loc3)
    call parameterval(tsr,sol,spr1d,spr1)
    call parameterval(tsr,sol,spr2d,spr2)
    call parameterval(tsr,sol,skw1d,skw1)
    call parameterval(tsr,sol,skw2d,skw2)
    call parameterval(tsr,sol,scl1d,scl1)
    call parameterval(tsr,sol,scl2d,scl2)
    call parameterval(tsr,sol,scl3d,scl3)

    ! Bounds of integration
    a = 0.0_dp ! starting at turbine
    b = (scl3 + 5.0_dp)*dia ! ending at the inflection point of the vorticity (when it decays)
    c = -1.0_dp*dia ! one diameter laterally
    d = 1.0_dp*dia ! one diameter laterally

    ! Translating the turbine position (placing turbine at 0,0)
    x0 = x0t - xt
    y0 = y0t - yt

    velxi = 0.0_dp
    velyi = 0.0_dp

    if (inte == 1) then
      ! Using 2D Simpson's Rule to integrate****************************************
      h = (b - a)/(m_in)
      k = (d - c)/(n_in)

      do i = 1,m_in
        xdiv(i) = a + i*h
      end do

      do j = 1,n_in
        ydiv(j) = c + j*k
      end do

      m = m_in/2
      n = n_in/2

      call integrandx(c,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval1)
      call integrandx(d,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval2)
      call integrandx(c,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval3)
      call integrandx(d,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval4)

      call integrandy(c,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval1)
      call integrandy(d,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval2)
      call integrandy(c,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval3)
      call integrandy(d,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval4)

      xsum1 = 0.0_dp
      xsum2 = 0.0_dp
      xsum3 = 0.0_dp
      xsum4 = 0.0_dp
      xsum5 = 0.0_dp
      xsum6 = 0.0_dp
      xsum7 = 0.0_dp
      xsum8 = 0.0_dp
      xsum9 = 0.0_dp
      xsum10 = 0.0_dp
      xsum11 = 0.0_dp
      xsum12 = 0.0_dp

      ysum1 = 0.0_dp
      ysum2 = 0.0_dp
      ysum3 = 0.0_dp
      ysum4 = 0.0_dp
      ysum5 = 0.0_dp
      ysum6 = 0.0_dp
      ysum7 = 0.0_dp
      ysum8 = 0.0_dp
      ysum9 = 0.0_dp
      ysum10 = 0.0_dp
      ysum11 = 0.0_dp
      ysum12 = 0.0_dp

      intlim = 0.0_dp

      do j = 1,n
        if ((ydiv(2*j-1) > intlim) .or. (ydiv(2*j-1) < -intlim)) then
          call integrandx(ydiv(2*j-1),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum1 = xsum1 + xsum
          call integrandy(ydiv(2*j-1),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum1 = ysum1 + ysum
        end if
      end do
      do j = 1,n-1
        if ((ydiv(2*j) > intlim) .or. (ydiv(2*j) < -intlim)) then
          call integrandx(ydiv(2*j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum2 = xsum2 + xsum
          call integrandy(ydiv(2*j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum2 = ysum2 + ysum
        end if
      end do
      do j = 1,n
        if ((ydiv(2*j-1) > intlim) .or. (ydiv(2*j-1) < -intlim)) then
          call integrandx(ydiv(2*j-1),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum3 = xsum3 + xsum
          call integrandy(ydiv(2*j-1),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum3 = ysum3 + ysum
        end if
      end do
      do j = 1,n-1
        if ((ydiv(2*j) > intlim) .or. (ydiv(2*j) < -intlim)) then
          call integrandx(ydiv(2*j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum4 = xsum4 + xsum
          call integrandy(ydiv(2*j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum4 = ysum4 + ysum
        end if
      end do
      do i = 1,m
        call integrandx(c,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        xsum5 = xsum5 + xsum
        call integrandy(c,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ysum5 = ysum5 + ysum
      end do
      do i = 1,m
        call integrandx(d,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        xsum6 = xsum6 + xsum
        call integrandy(d,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ysum6 = ysum6 + ysum
      end do
      do i = 1,m-1
        call integrandx(c,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        xsum7 = xsum7 + xsum
        call integrandy(c,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ysum7 = ysum7 + ysum
      end do
      do i = 1,m-1
        call integrandx(d,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        xsum8 = xsum8 + xsum
        call integrandy(d,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ysum8 = ysum8 + ysum
      end do

      do j = 1,n
        do i = 1,m
          if ((ydiv(2*j-1) > intlim) .or. (ydiv(2*j-1) < -intlim)) then
            call integrandx(ydiv(2*j-1),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum9 = xsum9 + xsum
            call integrandy(ydiv(2*j-1),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum9 = ysum9 + ysum
          end if
        end do
      end do
      do j = 1,n-1
        do i = 1,m
          if ((ydiv(2*j) > intlim) .or. (ydiv(2*j) < -intlim)) then
            call integrandx(ydiv(2*j),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum10 = xsum10 + xsum
            call integrandy(ydiv(2*j),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum10 = ysum10 + ysum
          end if
        end do
      end do
      do j = 1,n
        do i = 1,m-1
          if ((ydiv(2*j-1) > intlim) .or. (ydiv(2*j-1) < -intlim)) then
            call integrandx(ydiv(2*j-1),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum11 = xsum11 + xsum
            call integrandy(ydiv(2*j-1),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum11 = ysum11 + ysum
          end if
        end do
      end do
      do j = 1,n-1
        do i = 1,m-1
          if ((ydiv(2*j) > intlim) .or. (ydiv(2*j) < -intlim)) then
            call integrandx(ydiv(2*j),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum12 = xsum12 + xsum
            call integrandy(ydiv(2*j),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum12 = ysum12 + ysum
          end if
        end do
      end do

      velxi = (h*k*(xval1 + xval2 + xval3 + xval4 + 4.0_dp*xsum1 + 2.0_dp*xsum2 +&
      4.0_dp*xsum3 + 2.0_dp*xsum4 + 4.0_dp*xsum5 + 4.0_dp*xsum6 + 2.0_dp*xsum7 + 2.0_dp*xsum8 +&
      16.0_dp*xsum9 + 8.0_dp*xsum10 + 8.0_dp*xsum11 + 4.0_dp*xsum12)/9.0_dp)*(abs(rot)/pi2)

      velyi = (h*k*(yval1 + yval2 + yval3 + yval4 + 4.0_dp*ysum1 + 2.0_dp*ysum2 +&
      4.0_dp*ysum3 + 2.0_dp*ysum4 + 4.0_dp*ysum5 + 4.0_dp*ysum6 + 2.0_dp*ysum7 + 2.0_dp*ysum8 +&
      16.0_dp*ysum9 + 8.0_dp*ysum10 + 8.0_dp*ysum11 + 4.0_dp*ysum12)/9.0_dp)*(abs(rot)/pi2)

    else if (inte == 2) then
      ! Using 2D Trapezoidal Rule to integrate****************************************
      h = (b - a)/m_in
      k = (d - c)/n_in

      do i = 1,m_in
        xdiv(i) = a + i*h
      end do

      do j = 1,n_in
        ydiv(j) = c + j*k
      end do

      call integrandx(c,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval1)
      call integrandx(c,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval2)
      call integrandx(d,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval3)
      call integrandx(d,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval4)

      call integrandy(c,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval1)
      call integrandy(c,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval2)
      call integrandy(d,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval3)
      call integrandy(d,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval4)

      xsum1 = 0.0_dp
      xsum2 = 0.0_dp
      xsum3 = 0.0_dp
      xsum4 = 0.0_dp
      xsum5 = 0.0_dp

      ysum1 = 0.0_dp
      ysum2 = 0.0_dp
      ysum3 = 0.0_dp
      ysum4 = 0.0_dp
      ysum5 = 0.0_dp

      intlim = 0.0_dp

      do i = 1,m_in-1
        call integrandx(c,xdiv(i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        xsum1 = xsum1 + xsum
        call integrandy(c,xdiv(i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ysum1 = ysum1 + ysum
      end do
      do i = 1,m_in-1
        call integrandx(d,xdiv(i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        xsum2 = xsum2 + xsum
        call integrandy(d,xdiv(i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ysum2 = ysum2 + ysum
      end do
      do j = 1,n_in-1
        if ((ydiv(j) > intlim) .or. (ydiv(j) < -intlim)) then
          call integrandx(ydiv(j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum3 = xsum3 + xsum
          call integrandy(ydiv(j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum3 = ysum3 + ysum
        end if
      end do
      do j = 1,n_in-1
        if ((ydiv(j) > intlim) .or. (ydiv(j) < -intlim)) then
          call integrandx(ydiv(j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum4 = xsum4 + xsum
          call integrandy(ydiv(j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum4 = ysum4 + ysum
        end if
      end do
      do j = 1,n_in-1
        do i = 1,m_in-1
          if ((ydiv(j) > intlim) .or. (ydiv(j) < -intlim)) then
            call integrandx(ydiv(j),xdiv(i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum5 = xsum5 + xsum
            call integrandy(ydiv(j),xdiv(i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum5 = ysum5 + ysum
          end if
        end do
      end do

      velxi = (h*k*(xval1 + xval2 + xval3 + xval4 + 2.0_dp*xsum1 + 2.0_dp*xsum2 +&
      2.0_dp*xsum3 + 2.0_dp*xsum4 + 4.0_dp*xsum5)/4.0_dp)*(abs(rot)/pi2)

      velyi = (h*k*(yval1 + yval2 + yval3 + yval4 + 2.0_dp*ysum1 + 2.0_dp*ysum2 +&
      2.0_dp*ysum3 + 2.0_dp*ysum4 + 4.0_dp*ysum5)/4.0_dp)*(abs(rot)/pi2)

    end if

    ! ****************************************************************************

    velx = velxi/Vinf
    vely = velyi/Vinf

end subroutine vel_field


! Calculating vorticity strength for polynomial surface fitting
subroutine sheet_vort(ndata,xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,&
  coef5,coef6,coef7,coef8,coef9,dia,vort)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: ndata
  real(dp), dimension(10), intent(in) :: coef0,coef1,coef2,coef3,coef4,coef5
  real(dp), dimension(10), intent(in) :: coef6,coef7,coef8,coef9
  real(dp), dimension(ndata), intent(in) :: xttr,ystr,posdn,poslt
  real(dp), intent(in) :: dia
  ! out
  real(dp), dimension(ndata), intent(out) :: vort
  ! local
  integer :: i
  real(dp) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

  do i = 1,ndata
    call parameterval(xttr(i),ystr(i),coef0,loc1)
    call parameterval(xttr(i),ystr(i),coef1,loc2)
    call parameterval(xttr(i),ystr(i),coef2,loc3)
    call parameterval(xttr(i),ystr(i),coef3,spr1)
    call parameterval(xttr(i),ystr(i),coef4,spr2)
    call parameterval(xttr(i),ystr(i),coef5,skw1)
    call parameterval(xttr(i),ystr(i),coef6,skw2)
    call parameterval(xttr(i),ystr(i),coef7,scl1)
    call parameterval(xttr(i),ystr(i),coef8,scl2)
    call parameterval(xttr(i),ystr(i),coef9,scl3)

    call vorticitystrength(posdn(i),poslt(i),dia,loc1,loc2,loc3,spr1,spr2,&
    skw1,skw2,scl1,scl2,scl3,vort(i))

  end do

end subroutine sheet_vort


! Calculating vorticity strength for polynomial surface fitting
subroutine sheet_vel(ndata,xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,&
  coef5,coef6,coef7,coef8,coef9,dia,Vinf,m_in,n_in,inte,vel)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: ndata,m_in,n_in,inte
  real(dp), dimension(10), intent(in) :: coef0,coef1,coef2,coef3,coef4,coef5
  real(dp), dimension(10), intent(in) :: coef6,coef7,coef8,coef9
  real(dp), dimension(ndata), intent(in) :: xttr,ystr,posdn,poslt
  real(dp), intent(in) :: dia,Vinf
  ! out
  real(dp), dimension(ndata), intent(out) :: vel
  ! local
  integer :: i
  real(dp) :: rot,chord,velx,vely

  do i = 1,ndata
    rot = xttr(i)*Vinf/(dia/2.0_dp)
    chord = ystr(i)*(dia/2.0_dp)/3

    call vel_field(0.0_dp,0.0_dp,posdn(i),poslt(i),dia,rot,chord,3,Vinf,&
    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m_in,n_in,&
    inte,velx,vely)

    vel(i) = (velx*Vinf + Vinf)/Vinf

  end do

end subroutine sheet_vel


! Calculating loading forces on VAWT blades and power of the turbine
! Aerodynamics of multiple vertical-axis wind turbines using a
! modified Actuator Cylinder approach
! Developed by Andrew Ning at Brigham Young University
! https://github.com/byuflowlab/vawt-ac
! https://doi.org/10.5281/zenodo.165183
subroutine radialforce(n,f,uvec,vvec,thetavec,af_data,cl_data,cd_data,r,chord,&
  twist,delta,B,Omega,Vinf,Vinfx,Vinfy,rho,interp,q,k,Cp,Tp,Vn,Vt)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n,f,B,interp
    real(dp), dimension(n), intent(in) :: uvec,vvec,thetavec,Vinfx,Vinfy
    real(dp), dimension(f), intent(in) :: af_data,cl_data,cd_data
    real(dp), intent(in) :: r,chord,twist,delta,Omega,Vinf,rho

    ! out
    real(dp), intent(out) :: k,Cp
    real(dp), dimension(n), intent(out) :: q,Tp,Vn,Vt

    ! local
    integer :: i
    real(dp) :: pi,rotation,sigma,Ctend,Cto,a,H,As,P,Pend
    real(dp), dimension(n) :: W2,phi,alpha,qdyn,cn,ct,cl,cd,integrand,Rp,Zp,Qp
    intrinsic sin
    intrinsic cos
    intrinsic tan
    intrinsic abs
    intrinsic sqrt
    intrinsic atan2
    pi = 3.1415926535897932_dp

    ! set the rotation direction
    if (Omega >= 0.0_dp) then
      rotation = 1.0_dp
    else
      rotation = -1.0_dp
    end if

    sigma = B*chord/r
    do i = 1,n
      ! velocity components and angles
      Vn(i) = (Vinf*(1.0_dp + uvec(i)) + Vinfx(i))*sin(thetavec(i)) - &
      (Vinf*vvec(i) + Vinfy(i))*cos(thetavec(i))
      Vt(i) = rotation*((Vinf*(1.0_dp + uvec(i)) + Vinfx(i))*cos(thetavec(i)) + &
      (Vinf*vvec(i) + Vinfy(i))*sin(thetavec(i))) + abs(Omega)*r

      ! Original Code (without ability to add wake velocity deficits)
      ! Vn = Vinf*(1.0_dp + uvec(i))*sin(thetavec(i)) - Vinf*vvec(i)*cos(thetavec(i))
      ! Vt = rotation*(Vinf*(1.0_dp + uvec(i))*cos(thetavec(i)) + &
      ! Vinf*vvec(i)*sin(thetavec(i))) + abs(Omega)*r

      W2(i) = Vn(i)**2 + Vt(i)**2
      phi(i) = atan2(Vn(i), Vt(i))
      alpha(i) = phi(i) - twist

      ! airfoil
      if (interp == 1) then
        call interpolate(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
        call interpolate(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))
      else if (interp == 2) then
        call splineint(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
        call splineint(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))
      end if

      ! rotate force coefficients
      cn(i) = cl(i)*cos(phi(i)) + cd(i)*sin(phi(i))
      ct(i) = cl(i)*sin(phi(i)) - cd(i)*cos(phi(i))

      ! radial force
      q(i) = sigma/(4.0_dp*pi)*cn(i)*(W2(i)/Vinf**2)

      ! instantaneous forces
      qdyn(i) = 0.5_dp*rho*W2(i)
      Rp(i) = -cn(i)*qdyn(i)*chord
      Tp(i) = ct(i)*qdyn(i)*chord/cos(delta)
      Zp(i) = -cn(i)*qdyn(i)*chord*tan(delta)

      ! nonlinear correction factor
      integrand(i) = (W2(i)/Vinf**2)*(cn(i)*sin(thetavec(i)) - &
      rotation*ct(i)*cos(thetavec(i))/cos(delta))
    end do

    call pInt(n,thetavec,integrand,Ctend)
    Cto = sigma/(4.0_dp*pi)*Ctend

    if (Cto > 2.0_dp) then ! propeller brake
      a = 0.5_dp*(1.0_dp + sqrt(1.0_dp + Cto))
      k = 1.0_dp/(a-1.0_dp)
    else if (Cto > 0.96) then ! empirical
      a = 1.0_dp/7.0_dp*(1.0_dp + 3.0_dp*sqrt(7.0_dp/2.0_dp*Cto - 3.0_dp))
      k = 18.0_dp*a/(7.0_dp*a**2 - 2.0_dp*a + 4.0_dp)
    else ! momentum
      a = 0.5_dp*(1.0_dp - sqrt(1.0_dp - Cto))
      k = 1.0_dp/(1.0_dp-a)
    end if

    ! power coefficient
    H = 1.0_dp ! per unit height
    As = 2.0_dp*r*H
    Qp = r*Tp

    call pInt(n,thetavec,Qp,Pend)
    P = abs(Omega)*B/(2.0_dp*pi)*Pend
    Cp = P/(0.5_dp*rho*Vinf**3*As)

end subroutine radialforce


! Calculating effective velocities around given turbine due to wake interaction
! For use only with Simpson's method
subroutine overlap(t,p,xt,yt,diat,rott,chord,blades,x0,y0,dia,Vinf,loc1,loc2,loc3,spr1,spr2,&
  skw1,skw2,scl1,scl2,scl3,m,n,inte,pointcalc,velx,vely)

    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t,p,m,n,blades,inte,pointcalc
    real(dp), dimension(t), intent(in) :: xt,yt,diat,rott
    real(dp), intent(in) :: x0,y0,dia,Vinf,chord
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), dimension(p), intent(out) :: velx,vely

    ! local
    integer :: i,j,k,l
    real(dp) :: pi,theta
    real(dp), dimension(p) :: xd,yd,velxi,velyi,velx_int,vely_int,intex,intey
    intrinsic sin
    intrinsic cos
    intrinsic sqrt
    intrinsic abs
    pi = 3.1415926535897932_dp

    intex = 0.0_dp
    intey = 0.0_dp

    ! finding points around the flight path of the blades
    do i = 1,p
      if (pointcalc == 1) then
        theta = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
        xd(i) = x0 - sin(theta)*(dia/2.0_dp)
        yd(i) = y0 + cos(theta)*(dia/2.0_dp)
      else if (pointcalc == 0) then
        if (i == 1) then
          xd(i) = x0
          yd(i) = y0
        else
          xd(i) = 0.0_dp
          yd(i) = 0.0_dp
        end if
      end if
    end do

    if (t == 1) then ! coupled configuration (only two VAWTs)
      if (pointcalc == 1) then
        do j = 1,p
          call vel_field(xt(1),yt(1),xd(j),yd(j),diat(1),rott(1),chord,blades,Vinf,&
          loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
          velxi(j),velyi(j))
          velx(j) = velxi(j)*Vinf
          vely(j) = velyi(j)*Vinf
        end do
      else if (pointcalc == 0) then
        do j = 1,p
          if (j == 1) then
            call vel_field(xt(1),yt(1),xd(j),yd(j),diat(1),rott(1),chord,blades,Vinf,&
            loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
            velxi(j),velyi(j))
            velx(j) = velxi(j)*Vinf
            vely(j) = velyi(j)*Vinf
          else
            velx(j) = 0.0_dp
            vely(j) = 0.0_dp
          end if
        end do
      end if
    else ! multiple turbine wake overlap
      do j = 1,t
        if (pointcalc == 1) then
          do k = 1,p
            call vel_field(xt(j),yt(j),xd(k),yd(k),diat(j),rott(j),chord,blades,&
            Vinf,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
            velx_int(k),vely_int(k))
            velx_int(k) = -velx_int(k)

            ! sum of squares of velocity deficits
            if (velx_int(k) >= 0.0_dp) then
              intex(k) = intex(k) + (velx_int(k))**2
            else
              intex(k) = intex(k) - (velx_int(k))**2
            end if

            if (vely_int(k) >= 0.0_dp) then
              intey(k) = intey(k) + (vely_int(k))**2
            else
              intey(k) = intey(k) - (vely_int(k))**2
            end if
          end do
        else if (pointcalc == 0) then
          do k = 1,p
            if (k == 1) then
              call vel_field(xt(j),yt(j),xd(k),yd(k),diat(j),rott(j),chord,blades,&
              Vinf,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
              velx_int(k),vely_int(k))
              velx_int(k) = -velx_int(k)
            else
              velx_int(k) = 0.0_dp
            end if

            ! sum of squares of velocity deficits
            if (velx_int(k) >= 0.0_dp) then
              intex(k) = intex(k) + (velx_int(k))**2
            else
              intex(k) = intex(k) - (velx_int(k))**2
            end if

            if (vely_int(k) >= 0.0_dp) then
              intey(k) = intey(k) + (vely_int(k))**2
            else
              intey(k) = intey(k) - (vely_int(k))**2
            end if
          end do
        end if
      end do

      ! square root of sum of squares
      do l = 1,p
        if (intex(l) >= 0.0_dp) then
          velx(l) = -Vinf*(sqrt(intex(l)))
        else
          velx(l) = Vinf*(sqrt(abs(intex(l))))
        end if

        if (intey(l) >= 0.0_dp) then
          vely(l) = Vinf*(sqrt(intey(l)))
        else
          vely(l) = -Vinf*(sqrt(abs(intey(l))))
        end if
      end do
    end if

end subroutine overlap


! Calculating power and coefficient of power using velocity vector summation
subroutine powercalc(n,f,thetavec,Vinf,wake_x,wake_y,Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,&
  Omega,r,H,af_data,cl_data,cd_data,twist,rho,interp,P,Cp)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n,f,interp
    real(dp), dimension(n), intent(in) :: thetavec,wake_x,wake_y
    real(dp), dimension(n), intent(in) :: Vnp,Vnn,Vtp,Vtn,Cpp,Cpn
    real(dp), dimension(f), intent(in) :: af_data,cl_data,cd_data
    real(dp), intent(in) :: Vinf,Omega,r,H,twist,rho

    ! out
    real(dp), intent(out) :: P,Cp

    ! local
    integer :: i
    real(dp) :: pi,As
    real(dp), dimension(n) :: W20,phi0,alpha0,cl0,cd0,ct0
    real(dp), dimension(n) :: Vn,Vt,W2,phi,alpha,cl,cd,ct,Cpl
    intrinsic sin
    intrinsic cos
    intrinsic atan2
    intrinsic abs
    pi = 3.1415926535897932_dp

    if (Omega >= 0.0_dp) then
      do i = 1,n
        ! calculate baseline values
        W20(i) = Vnp(i)**2 + Vtp(i)**2
        phi0(i) = atan2(Vnp(i),Vtp(i))
        alpha0(i) = phi0(i) - twist
        if (interp == 1) then
          call interpolate(f,af_data,cl_data,alpha0(i)*180.0_dp/pi,cl0(i))
          call interpolate(f,af_data,cd_data,alpha0(i)*180.0_dp/pi,cd0(i))
        else if (interp == 2) then
          call splineint(f,af_data,cl_data,alpha0(i)*180.0_dp/pi,cl0(i))
          call splineint(f,af_data,cd_data,alpha0(i)*180.0_dp/pi,cd0(i))
        end if
        ct0(i) = cl0(i)*sin(alpha0(i)) - cd0(i)*cos(alpha0(i))

        ! correct normal/tangential velocities with wake velocities
        Vn(i) = Vnp(i) + wake_x(i)*sin(thetavec(i)) - wake_y(i)*cos(thetavec(i))
        Vt(i) = Vtp(i) + wake_x(i)*cos(thetavec(i)) + wake_y(i)*sin(thetavec(i))
        W2(i) = Vn(i)**2 + Vt(i)**2

        ! compute new inflow angle
        phi(i) = atan2(Vn(i),Vt(i))
        alpha(i) = phi(i) - twist

        ! airfoil
        if (interp == 1) then
          call interpolate(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
          call interpolate(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))
        else if (interp == 2) then
          call splineint(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
          call splineint(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))
        end if

        ! compute new tangential force coefficient
        ct(i) = cl(i)*sin(alpha(i)) - cd(i)*cos(alpha(i))

        ! provide relative correction to power coefficient
        Cpl(i) = Cpp(i)*W2(i)/W20(i)*ct(i)/ct0(i)
      end do
    else
      do i = 1,n
        ! calculate baseline values
        W20(i) = Vnn(i)**2 + Vtn(i)**2
        phi0(i) = atan2(Vnn(i),Vtn(i))
        alpha0(i) = phi0(i) - twist
        if (interp == 1) then
          call interpolate(f,af_data,cl_data,alpha0(i)*180.0_dp/pi,cl0(i))
          call interpolate(f,af_data,cd_data,alpha0(i)*180.0_dp/pi,cd0(i))
        else if (interp == 2) then
          call splineint(f,af_data,cl_data,alpha0(i)*180.0_dp/pi,cl0(i))
          call splineint(f,af_data,cd_data,alpha0(i)*180.0_dp/pi,cd0(i))
        end if
        ct0(i) = cl0(i)*sin(alpha0(i)) - cd0(i)*cos(alpha0(i))

        ! correct normal/tangential velocities with wake velocities
        Vn(i) = Vnn(i) + wake_x(i)*sin(thetavec(i)) - wake_y(i)*cos(thetavec(i))
        Vt(i) = Vtn(i) - wake_x(i)*cos(thetavec(i)) - wake_y(i)*sin(thetavec(i))
        W2(i) = Vn(i)**2 + Vt(i)**2

        ! compute new inflow angle
        phi(i) = atan2(Vn(i),Vt(i))
        alpha(i) = phi(i) - twist

        ! airfoil
        if (interp == 1) then
          call interpolate(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
          call interpolate(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))
        else if (interp == 2) then
          call splineint(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
          call splineint(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))
        end if

        ! compute new tangential force coefficient
        ct(i) = cl(i)*sin(phi(i)) - cd(i)*cos(phi(i))

        ! provide relative correction to power coefficient
        Cpl(i) = Cpn(i)*W2(i)/W20(i)*ct(i)/ct0(i)

      end do
    end if

    ! integrating local power coefficients
    call pInt(n,thetavec,Cpl,Cp)

    ! power calculation
    As = 2.0_dp*r*H ! swept area
    P = Cp*0.5_dp*rho*Vinf**3*As

end subroutine powercalc


! To build for Python interface:
! f2py -c  --opt=-O2 -m _vawtwake VAWT_Wake_Model.f90
! python C:\Python27\Scripts\f2py.py -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _vawtwake VAWT_Wake_Model.f90
