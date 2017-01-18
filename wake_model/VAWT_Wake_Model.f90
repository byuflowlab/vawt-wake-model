! To build for Python interface (Mac): f2py -c  --opt=-O2 -m _vawtwake VAWT_Wake_Model.f90


! trapezoidal integration
subroutine trapz(n,x,y,integral) ! integrate y w.r.t. x
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
      if (xval < x(i)) then
        if (i == 2) then
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)
          y1 = y(1)
          y2 = y(2)
          y3 = y(3)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else if (i == n) then
          x1 = x(n-2)
          x2 = x(n-1)
          x3 = x(n)
          y1 = y(n-2)
          y2 = y(n-1)
          y3 = y(n)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else
          if (xval <= (x(i)+x(i-1))/2.0_dp) then
            x1 = x(i-2)
            x2 = x(i-1)
            x3 = x(i)
            y1 = y(i-2)
            y2 = y(i-1)
            y3 = y(i)
            call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
          else
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
      else if (xval == x(i)) then
        yval = y(i)
        exit
      end if

    end do

end subroutine splineint


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


! Calculating EMG parameter values based on given CFD database
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


! Creating the fit of the vorticity distribution
subroutine EMGdist(y,loc,spr,skw,scl,gam_skew)
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
    real(dp) :: xd,yd,loc,spr,skw,scl,g1,g2
    intrinsic exp

    xd = x/dia ! normalizing x by the diameter
    yd = y/dia ! normalizing y by the diameter

    ! Limiting the parameter components to create expected behavior
    if (loc1 > -0.001_dp) then ! ensure concave down
      if (loc2 < 0.01_dp) then ! ensure slight increase moving downstream
        loc = -0.001_dp*xd*xd + 0.01_dp*xd + loc3
      else
        loc = -0.001_dp*xd*xd + loc2*xd + loc3
      end if
    else
      if (loc2 < 0.01_dp) then ! ensure slight increase moving downstream
        loc = loc1*xd*xd + 0.01_dp*xd + loc3
      else
        loc = loc1*xd*xd + loc2*xd + loc3
      end if
    end if

    if (spr1 > -0.001_dp) then ! ensure decrease in value (more spread downstream)
      spr = -0.001_dp*xd + spr2
    else
      spr = spr1*xd + spr2
    end if

    skw = skw1*xd + skw2 ! no limitations necessary

    if (scl2 < 0.05_dp) then ! ensure decay moving downstream
      scl = scl1/(1.0_dp + exp(0.05_dp*(xd - scl3)))
    else
      scl = scl1/(1.0_dp + exp(scl2*(xd - scl3)))
    end if

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


! Parameterized VAWT Wake Model using CFD vorticity data
! Developed by Eric Tingey at Brigham Young University
! This code models the wake behind a vertical-axis wind turbine based on
! parameters like tip-speed ratio, solidity and wind speed by converting the
! vorticity of the wake into velocity information. The model uses CFD data
! obtained from STAR-CCM+ of simulated turbines to make the wake model as
! accurate as possible.
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

! Performing integration to convert vorticity into velocity
subroutine vel_field(xt,yt,x0t,y0t,dia,rot,chord,blades,velf,loc1d,loc2d,loc3d,spr1d,spr2d,&
  skw1d,skw2d,scl1d,scl2d,scl3d,m_in,n_in,inte,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: xt,yt,x0t,y0t,dia,rot,chord,blades,velf
    real(dp), dimension(10), intent(in) :: loc1d,loc2d,loc3d,spr1d,spr2d,skw1d,skw2d,scl1d,scl2d,scl3d
    integer, intent(in) :: m_in,n_in,inte

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

    tsr = (dia/2.0_dp)*abs(rot)/velf
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
    b = (scl3+5.0_dp)*dia ! ending at the inflection point of the vorticity (when it decays)
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

    velx = velxi/velf
    vely = velyi/velf

end subroutine vel_field


! Aerodynamics of multiple vertical axis wind turbines using a modified Actuator Cylinder approach
! Developed by Andrew Ning at Brigham Young University
! https://github.com/byuflowlab/vawt-ac
! https://doi.org/10.5281/zenodo.165183

! Calculating loading forces on VAWT blades and power of the turbine
subroutine radialforce(n,f,uvec,vvec,thetavec,af_data,cl_data,cd_data,r,chord,&
  twist,delta,B,Omega,Vinf,Vinfx,Vinfy,rho,mu,interp,q,ka,CTo,CPo,Rp,Tp,Zp)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n,f,interp
    real(dp), dimension(n), intent(in) :: uvec,vvec,thetavec,Vinfx,Vinfy
    real(dp), dimension(f), intent(in) :: af_data,cl_data,cd_data
    real(dp), intent(in) :: r,chord,twist,delta,B,Omega,Vinf,rho,mu

    ! out
    real(dp), intent(out) :: ka,CTo,CPo
    real(dp), dimension(n), intent(out) :: q,Rp,Tp,Zp

    ! local
    integer :: i
    real(dp) :: pi,rotation,sigma,CTend,a,H,Sref,P,Pend
    real(dp), dimension(n) :: Vn,Vt,W,phi,alpha,cl,cd,cn,ct,qdyn,integrand,Qp
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

      W(i) = sqrt(Vn(i)**2 + Vt(i)**2)
      phi(i) = atan2(Vn(i), Vt(i))
      alpha(i) = phi(i) - twist
      ! Re = rho*W*chord/mu  ! currently no Re dependence

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
      q(i) = sigma/(4.0_dp*pi)*cn(i)*(W(i)/Vinf)**2

      ! instantaneous forces
      qdyn(i) = 0.5_dp*rho*W(i)**2
      Rp(i) = -cn(i)*qdyn(i)*chord
      Tp(i) = ct(i)*qdyn(i)*chord/cos(delta)
      Zp(i) = -cn(i)*qdyn(i)*chord*tan(delta)

      ! nonlinear correction factor
      integrand(i) = (W(i)/Vinf)**2*(cn(i)*sin(thetavec(i)) - &
      rotation*ct(i)*cos(thetavec(i))/cos(delta))
    end do

    call pInt(n,thetavec,integrand,CTend)
    CTo = sigma/(4.0_dp*pi)*CTend
    if (CTo > 2.0_dp) then
      a = 0.5_dp*(1.0_dp + sqrt(1.0_dp + CTo))
      ka = 1.0_dp/(a-1.0_dp)
    else if (CTo > 0.96) then
      a = 1.0_dp/7.0_dp*(1.0_dp + 3.0_dp*sqrt(7.0_dp/2.0_dp*CTo - 3.0_dp))
      ka = 18.0_dp*a/(7.0_dp*a**2 - 2.0_dp*a + 4.0_dp)
    else
      a = 0.5_dp*(1.0_dp - sqrt(1.0_dp - CTo))
      ka = 1.0_dp/(1.0_dp-a)
    end if

    ! power coefficient
    H = 1.0_dp ! per unit height
    Sref = 2.0_dp*r*H
    Qp = r*Tp

    call pInt(n,thetavec,Qp,Pend)
    P = abs(Omega)*B/(2.0_dp*pi)*Pend
    CPo = P/(0.5_dp*rho*Vinf**3*Sref)

end subroutine radialforce


! Calculating effective velocities around given turbine due to wake interaction
subroutine overlap(t,p,xt,yt,diat,rott,chord,blades,x0,y0,dia,velf,loc1,loc2,loc3,spr1,spr2,&
  skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)

    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t,p,m,n,inte
    real(dp), dimension(t), intent(in) :: xt,yt,diat,rott
    real(dp), intent(in) :: x0,y0,dia,velf,chord,blades
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
      theta = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
      xd(i) = x0 - sin(theta)*(dia/2.0_dp)
      yd(i) = y0 + cos(theta)*(dia/2.0_dp)
    end do

    if (t == 1) then ! coupled configuration (only two VAWTs)
      do j = 1,p
        call vel_field(xt(1),yt(1),xd(j),yd(j),diat(1),rott(1),chord,blades,velf,&
        loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
        velxi(j),velyi(j))
        velx(j) = velxi(j)*velf
        vely(j) = velyi(j)*velf
      end do
    else ! multiple turbine wake overlap
      do j = 1,t
        do k = 1,p
          call vel_field(xt(j),yt(j),xd(k),yd(k),diat(j),rott(j),chord,blades,velf,&
          loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
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
      end do

      ! square root of sum of squares
      do l = 1,p
        if (intex(l) >= 0.0_dp) then
          velx(l) = -velf*(sqrt(intex(l)))
        else
          velx(l) = velf*(sqrt(abs(intex(l))))
        end if

        if (intey(l) >= 0.0_dp) then
          vely(l) = velf*(sqrt(intey(l)))
        else
          vely(l) = -velf*(sqrt(abs(intey(l))))
        end if
      end do
    end if

end subroutine overlap


! Calculating effective velocities around given turbine due to wake interaction at (x0,y0)
subroutine overlappoint(t,xt,yt,diat,rot,chord,blades,x0,y0,velf,&
  loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t,m,n,inte
    real(dp), dimension(t), intent(in) :: xt,yt,diat,rot
    real(dp), intent(in) :: x0,y0,velf,chord,blades
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), intent(out) :: velx,vely

    ! local
    integer :: i
    real(dp) :: velx_int,vely_int,intex,intey
    intrinsic sqrt
    intrinsic abs

    intex = 0.0_dp
    intey = 0.0_dp

    ! sum of squares of velocity deficits
    do i = 1,t
      call vel_field(xt(i),yt(i),x0,y0,diat(i),rot(i),chord,blades,velf,&
      loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
      velx_int,vely_int)
      velx_int = -velx_int

      if (velx_int >= 0.0_dp) then
        intex = intex + (velx_int)**2
      else
        intex = intex - (velx_int)**2
      end if

      if (vely_int >= 0.0_dp) then
        intey = intey + (vely_int)**2
      else
        intey = intey - (vely_int)**2
      end if
    end do

    ! square root of sum of squares
    if (intex >= 0.0_dp) then
      velx = -velf*(sqrt(intex))
    else
      velx = velf*(sqrt(abs(intex)))
    end if

    if (intey >= 0.0_dp) then
      vely = velf*(sqrt(intey))
    else
      vely = -velf*(sqrt(abs(intey)))
    end if

end subroutine overlappoint


! Calculating power and coefficient of power of single or multiple turbines
subroutine powercalc(t,f,p,x,y,dia,rot,velf,loc1,loc2,loc3,spr1,spr2,&
  skw1,skw2,scl1,scl2,scl3,af_data,cl_data,cd_data,chord,twist,delta,blades,&
  H,rho,mu,uvec,vvec,interp,power_tot,CPo)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t,f,p,interp
    real(dp), dimension(t), intent(in) :: x,y,dia,rot
    real(dp), dimension(p), intent(in) :: uvec,vvec
    real(dp), dimension(f), intent(in) :: af_data,cl_data,cd_data
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
    real(dp), intent(in) :: velf,chord,twist,delta,blades,H,rho,mu

    ! out
    real(dp), intent(out) :: power_tot,CPo

    ! local
    integer :: m,n,i,j,inte
    real(dp) :: x0,y0,dia0,rot0,pi,ka,CTo!,CPo1,CPo2
    real(dp), dimension(t) :: xt,yt,diat,rott
    real(dp), dimension(p) :: velx,vely,thetavec,q,Rp,Tp,Zp!,velx1,velx2,vely1,vely2
    pi = 3.1415926535897932_dp

    m = 220
    n = 200
    inte = 1 ! 2D Simpson's Rule
    ! inte = 2 ! 2D Trapezoidal Rule

    do i = 1,p
      thetavec(i) = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
    end do

    ! Velocity and power measurements made at first turbine
    if (t >= 2) then ! needing to calculate wake interaction using sum of squares
      x0 = x(1)
      y0 = y(1)
      dia0 = dia(1)
      rot0 = rot(1)
      do j = 1,t-1
        xt(j) = x(j+1)
        yt(j) = y(j+1)
        diat(j) = dia(j+1)
        rott(j) = rot(j+1)
      end do
      call overlap(t-1,p,xt,yt,diat,rott,chord,blades,x0,y0,dia0,velf,loc1,loc2,loc3,&
      spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
    else ! no need to calculate wake interaction
      x0 = x(1)
      y0 = y(1)
      dia0 = dia(1)
      rot0 = rot(1)
      velx = 0.0_dp
      vely = 0.0_dp
    end if

    call radialforce(p,f,uvec,vvec,thetavec,af_data,cl_data,cd_data,(dia0/2.0_dp),chord,&
      twist,delta,blades,rot0,velf,velx,vely,rho,mu,interp,q,ka,CTo,CPo,Rp,Tp,Zp)

    power_tot = (0.5_dp*rho*velf**3)*(dia0*H)*CPo

end subroutine powercalc


! To build for Python interface:
! f2py -c  --opt=-O2 -m _vawtwake VAWT_Wake_Model.f90
! python C:\Python27\Scripts\f2py.py -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _vawtwake VAWT_Wake_Model.f90
