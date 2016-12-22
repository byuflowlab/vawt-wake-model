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
      if (xval .lt. x(i)) then
        yval = y(i-1) + (xval - x(i-1))*((y(i)-y(i-1))/(x(i)-x(i-1)))
        exit
      else if (xval .eq. x(i)) then
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
      if (xval .lt. x(i)) then
        if (i .eq. 2) then
          x1 = x(1)
          x2 = x(2)
          x3 = x(3)
          y1 = y(1)
          y2 = y(2)
          y3 = y(3)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else if (i .eq. n) then
          x1 = x(n-2)
          x2 = x(n-1)
          x3 = x(n)
          y1 = y(n-2)
          y2 = y(n-1)
          y3 = y(n)
          call cubspline(x1,x2,x3,y1,y2,y3,xval,yval)
        else
          if (xval .le. (x(i)+x(i-1))/2.0_dp) then
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
      else if (xval .eq. x(i)) then
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
    if (xval .lt. x2) then
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

    if (inte .eq. 1) then
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
        if ((ydiv(2*j-1) .gt. intlim) .or. (ydiv(2*j-1) .lt. -intlim)) then
          call integrandx(ydiv(2*j-1),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum1 = xsum1 + xsum
          call integrandy(ydiv(2*j-1),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum1 = ysum1 + ysum
        end if
      end do
      do j = 1,n-1
        if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
          call integrandx(ydiv(2*j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum2 = xsum2 + xsum
          call integrandy(ydiv(2*j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum2 = ysum2 + ysum
        end if
      end do
      do j = 1,n
        if ((ydiv(2*j-1) .gt. intlim) .or. (ydiv(2*j-1) .lt. -intlim)) then
          call integrandx(ydiv(2*j-1),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum3 = xsum3 + xsum
          call integrandy(ydiv(2*j-1),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum3 = ysum3 + ysum
        end if
      end do
      do j = 1,n-1
        if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
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
          if ((ydiv(2*j-1) .gt. intlim) .or. (ydiv(2*j-1) .lt. -intlim)) then
            call integrandx(ydiv(2*j-1),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum9 = xsum9 + xsum
            call integrandy(ydiv(2*j-1),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum9 = ysum9 + ysum
          end if
        end do
      end do
      do j = 1,n-1
        do i = 1,m
          if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
            call integrandx(ydiv(2*j),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum10 = xsum10 + xsum
            call integrandy(ydiv(2*j),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum10 = ysum10 + ysum
          end if
        end do
      end do
      do j = 1,n
        do i = 1,m-1
          if ((ydiv(2*j-1) .gt. intlim) .or. (ydiv(2*j-1) .lt. -intlim)) then
            call integrandx(ydiv(2*j-1),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
            xsum11 = xsum11 + xsum
            call integrandy(ydiv(2*j-1),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
            ysum11 = ysum11 + ysum
          end if
        end do
      end do
      do j = 1,n-1
        do i = 1,m-1
          if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
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

    else if (inte .eq. 2) then
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
        if ((ydiv(j) .gt. intlim) .or. (ydiv(j) .lt. -intlim)) then
          call integrandx(ydiv(j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum3 = xsum3 + xsum
          call integrandy(ydiv(j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum3 = ysum3 + ysum
        end if
      end do
      do j = 1,n_in-1
        if ((ydiv(j) .gt. intlim) .or. (ydiv(j) .lt. -intlim)) then
          call integrandx(ydiv(j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          xsum4 = xsum4 + xsum
          call integrandy(ydiv(j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ysum4 = ysum4 + ysum
        end if
      end do
      do j = 1,n_in-1
        do i = 1,m_in-1
          if ((ydiv(j) .gt. intlim) .or. (ydiv(j) .lt. -intlim)) then
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
    if (Omega .ge. 0.0_dp) then
      rotation = 1.0_dp
    else
      rotation = -1.0_dp
    end if

    sigma = B*chord/r
    do i = 1,n
      ! velocity components and angles
      ! Vn(i) = (Vinf*(1.0_dp + uvec(i)) + Vinfx(i))*sin(thetavec(i)) - &
      ! (Vinf*vvec(i) + Vinfy(i))*cos(thetavec(i))
      ! Vt(i) = rotation*((Vinf*(1.0_dp + uvec(i)) + Vinfx(i))*cos(thetavec(i)) + &
      ! (Vinf*vvec(i) + Vinfy(i))*sin(thetavec(i))) + abs(Omega)*r

      ! Original Code (without ability to add wake velocity deficits)
      Vn = Vinf*(1.0_dp + uvec(i))*sin(thetavec(i)) - Vinf*vvec(i)*cos(thetavec(i))
      Vt = rotation*(Vinf*(1.0_dp + uvec(i))*cos(thetavec(i)) + &
      Vinf*vvec(i)*sin(thetavec(i))) + abs(Omega)*r

      W(i) = sqrt(Vn(i)**2 + Vt(i)**2)
      phi(i) = atan2(Vn(i), Vt(i))
      alpha(i) = phi(i) - twist
      ! Re = rho*W*chord/mu  ! currently no Re dependence

      ! airfoil
      if (interp .eq. 1) then
        call interpolate(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
        call interpolate(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))
      else if (interp .eq. 2) then
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

    ! print *, 'Vn', Vn
    ! print *, 'Vt', Vt
    ! print *, 'W', W
    ! print *, 'phi', phi
    ! print *, 'alpha', alpha
    ! print *, 'cn', cn
    ! print *, 'ct', ct
    ! print *, 'q', q
    ! print *, 'qdyn', qdyn
    ! print *, 'Rp', Rp
    ! print *, 'Tp', Tp
    ! print *, 'Zp', Zp


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
subroutine overlap(t,p,xt,yt,diat,rott,chord,blades,x0,y0,dia,rot,velf,loc1,loc2,loc3,spr1,spr2,&
  skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t,p,m,n,inte
    real(dp), dimension(t), intent(in) :: xt,yt,diat,rott
    real(dp), intent(in) :: x0,y0,dia,rot,velf,chord,blades
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), dimension(p), intent(out) :: velx,vely

    ! local
    integer :: i,j,k,l
    real(dp) :: pi,theta
    real(dp), dimension(p) :: xd,yd,velx_int,vely_int,intex,intey
    intrinsic sin
    intrinsic cos
    intrinsic sqrt
    intrinsic abs
    pi = 3.1415926535897932_dp

    intex = 0.0_dp
    intey = 0.0_dp

    ! finding points around the flight path of the blades
    if (rot .ge. 0.0_dp) then
      do i = 1,p
        theta = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
        xd(i) = x0 - sin(theta)*(dia/2.0_dp)
        yd(i) = y0 + cos(theta)*(dia/2.0_dp)
      end do
    else if (rot .lt. 0.0_dp) then
      do i = 1,p
        theta = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
        xd(i) = x0 + sin(theta)*(dia/2.0_dp)
        yd(i) = y0 + cos(theta)*(dia/2.0_dp)
      end do
    end if

    ! sum of squares of velocity deficits
    do j = 1,t
      do k = 1,p
        call vel_field(xt(j),yt(j),xd(k),yd(k),diat(j),rott(j),chord,blades,velf,&
        loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
        velx_int(k),vely_int(k))
        velx_int(k) = -velx_int(k)

        if (velx_int(k) .ge. 0.0_dp) then
          intex(k) = intex(k) + (velx_int(k))**2
        else
          intex(k) = intex(k) - (velx_int(k))**2
        end if

        if (vely_int(k) .ge. 0.0_dp) then
          intey(k) = intey(k) + (vely_int(k))**2
        else
          intey(k) = intey(k) - (vely_int(k))**2
        end if
      end do
    end do

    ! square root of sum of squares
    do l = 1,p
      if (intex(l) .ge. 0.0_dp) then
        velx(l) = -velf*(sqrt(intex(l)))
      else
        velx(l) = velf*(sqrt(abs(intex(l))))
      end if

      if (intey(l) .ge. 0.0_dp) then
        vely(l) = velf*(sqrt(intey(l)))
      else
        vely(l) = -velf*(sqrt(abs(intey(l))))
      end if
    end do

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

      if (velx_int .ge. 0.0_dp) then
        intex = intex + (velx_int)**2
      else
        intex = intex - (velx_int)**2
      end if

      if (vely_int .ge. 0.0_dp) then
        intey = intey + (vely_int)**2
      else
        intey = intey - (vely_int)**2
      end if
    end do

    ! square root of sum of squares
    if (intex .ge. 0.0_dp) then
      velx = -velf*(sqrt(intex))
    else
      velx = velf*(sqrt(abs(intex)))
    end if

    if (intey .ge. 0.0_dp) then
      vely = velf*(sqrt(intey))
    else
      vely = -velf*(sqrt(abs(intey)))
    end if

end subroutine overlappoint


! Calculating velocities around given turbine due to wake interaction of coupled turbine
subroutine coupled(p,xt,yt,diat,rott,chord,blades,x0,y0,dia,rot,velf,loc1,loc2,loc3,&
  spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: p,m,n,inte
    real(dp), intent(in) :: x0,y0,dia,rot,chord,blades,velf
    real(dp), dimension(1), intent(in) :: xt,yt,diat,rott
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3

    ! out
    real(dp), dimension(p), intent(out) :: velx,vely

    ! local
    integer :: i,j
    real(dp) :: pi,theta
    real(dp), dimension(p) :: xd,yd,velxi,velyi
    intrinsic sin
    intrinsic cos
    pi = 3.1415926535897932_dp

    if (rot .ge. 0.0_dp) then
      do i = 1,p
        theta = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
        xd(i) = x0 - sin(theta)*(dia/2.0_dp)
        yd(i) = y0 + cos(theta)*(dia/2.0_dp)
      end do
    else if (rot .lt. 0.0_dp) then
      do i = 1,p
        theta = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
        xd(i) = x0 + sin(theta)*(dia/2.0_dp)
        yd(i) = y0 + cos(theta)*(dia/2.0_dp)
      end do
    end if

    do j = 1,p
      call vel_field(xt(1),yt(1),xd(j),yd(j),diat(1),rott(1),chord,blades,velf,&
      loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,&
      velxi(j),velyi(j))
      velx(j) = velxi(j)*velf
      vely(j) = velyi(j)*velf
    end do

end subroutine coupled


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

    ! print *, 'thetavec', thetavec*180.0_dp/pi

    ! Velocity and power measurements made at first turbine
    ! if (t .ge. 3) then ! needing to calculate wake interaction using sum of squares
    !   x0 = x(1)
    !   y0 = y(1)
    !   dia0 = dia(1)
    !   rot0 = rot(1)
    !   do j = 1,t-1
    !     xt(j) = x(j+1)
    !     yt(j) = y(j+1)
    !     diat(j) = dia(j+1)
    !     rott(j) = rot(j+1)
    !   end do
    !   call overlap(t-1,p,xt,yt,diat,rott,chord,blades,x0,y0,dia0,rot0,velf,loc1,loc2,loc3,&
    !   spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
    if (t .eq. 2) then ! calculating wake effect from one other turbine
      x0 = x(1)
      y0 = y(1)
      dia0 = dia(1)
      rot0 = rot(1)
      ! xt = x(2)
      ! yt = y(2)
      ! diat = dia(2)
      ! rott = rot(2)
      do j = 1,t
        xt(j) = x(j)
        yt(j) = y(j)
        diat(j) = dia(j)
        rott(j) = rot(j)
      end do
      ! print *, 'x', x0
      ! print *, 'y', y0
      ! ! ! print *, 'dia0', dia0
      ! print *, 'rot', rot0
      ! ! print *, 'xt', xt
      ! ! print *, 'yt', yt
      ! ! print *, 'diat', diat
      ! ! print *, 'rott', rott
      ! call coupled(p,xt,yt,diat,rott,chord,blades,x0,y0,dia0,rot0,velf,loc1,loc2,loc3,&
      !   spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
      call overlap(t,p,xt,yt,diat,rott,chord,blades,x0,y0,dia0,rot0,velf,loc1,loc2,loc3,&
        spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
      ! print *, 'velx', velx/velf+uvec
      ! print *, 'vely', vely/velf+vvec
      ! print *, 'velxw', velx/velf
      ! print *, 'velyw', vely/velf
    else ! no need to calculate wake interaction
      x0 = x(1)
      y0 = y(1)
      dia0 = dia(1)
      rot0 = rot(1)
      velx = 0.0_dp
      vely = 0.0_dp
    end if

    ! velx1 = (/0.16130340694612294,   4.7209492911345255e-002, -6.9198342809759783e-002, &
    ! -0.17266795077267935, -0.25705141118121771, -0.30149461209635875, -0.30208997229695350, &
    ! -0.28682010310347461, -0.26489139726461386, -0.17181964813148370, -0.15751629738982814, &
    ! -0.14059071850988003, -0.12358012988728773, -0.15177988980686594, -0.12689982467597716,  &
    ! -9.3035647675653421e-002, -3.5370883694315639e-002,  1.7940974286311369e-002, &
    ! -2.7342717791507022e-002, -0.17853279296819533, -0.32751267863635319, -0.41887725968713396, &
    ! -0.48505545837396574, -0.45262292862125436, -0.49249101938913242, -0.52644672032700035, &
    ! -0.55088653700710388, -0.68760758709454550, -0.71717005047912985, -0.73261256979626688, &
    ! -0.72515016630734186, -0.65344803171320098, -0.52163574228247211, -0.35101586526196871, &
    ! -0.14756568012472623,   8.2117459990049935e-002/)
    ! vely2 = (/0.32644446300958391     ,  0.32814712351894748    ,   0.30892810286263545   , &
    !    0.27546322867321482    ,   0.23382903113186398   ,    0.21074325824887130   ,    &
    !    0.20686367685664581   ,    0.19698227628947818  ,     0.19695495663507523  ,     &
    !    0.16595561301114545   ,     9.7310869916610560E-002  , 2.6598980994991743E-002 , &
    !    -5.7616996415355892E-002, -0.14532401637448478   ,   -0.21393684432495433   ,  &
    !     -0.25230222056395968    ,  -0.25509737481875117    ,  -0.21496656743990661   ,   &
    !     -0.13366156641203258   ,    -4.4785262809848222E-002  , 1.9853355428332289E-002  , &
    !     6.0102522288489862E-002 ,  7.7888141340459099E-002  , 7.9433550425731386E-002 ,  &
    !     7.9622633580537241E-002 ,  8.5814176540163295E-002  , 9.1252775789093146E-002 , &
    !     0.11338264910130991    ,   0.14944549946972727    ,   0.17806893764500994     ,  &
    !     0.20672332718015207   ,    0.22998847489370358    ,   0.24927371464681622   ,    &
    !     0.27203041502484610    ,   0.29490837094618799  ,     0.31365299754937614      /)
    ! velx2 = (/6.8978536799285339e-002,  1.3488551074846453e-002, -4.7405258956748755e-002, &
    ! -8.5718230757169084e-002, -0.11647266929260776,  -9.5757641686378459e-002, &
    ! -0.12192593841657884, -0.14948674303975545, -0.17534317906368060, -0.28000101557799612, &
    ! -0.31257925448747892, -0.33682458809367427, -0.34320167146225455, -0.30376618484689849, &
    ! -0.22267105799942125, -0.12106806361487750,  -5.3589556959241658e-003,  0.10903625752222976,&
    ! 3.1079897477083698e-002, -0.19642511489387499, -0.39664625398086401, -0.56281733620128549,&
    !  -0.68875525222745715, -0.75297265450832429, -0.75127734988957906, -0.72519960482924606,&
    !   -0.68408405616235313, -0.53577691869372646, -0.50068756894301458, -0.45775640359242509, &
    !   -0.41091586925547391, -0.43834068470827919, -0.36887415246039551, -0.27564295783122483, &
    !   -0.12596434436091719,   2.4924431632394599e-002/)
    ! vely1 = (/0.20283421431456705   ,    0.23504593588595254  ,     0.22442885638108578  ,  &
    !    0.17840920304124411  ,     0.10249922842399739  ,      8.2027928921899473E-003 , &
    !    -8.1364953622737560E-002, -0.15553473206646140   ,   -0.22512992737955984   ,   &
    !    -0.25425460326596039   ,   -0.24974236666407099    ,  -0.25300861658555046   ,  &
    !     -0.24896633974339624   ,   -0.26351208699193518   ,   -0.29646327088809094  ,   &
    !      -0.32135727532434050  ,    -0.33220848608006170   ,   -0.32233456760710577   ,  &
    !       -0.30152064442403570   ,   -0.27485693201339212    ,  -0.24415705084198139    , &
    !        -0.21374607336311752   ,   -0.18716368694322957   ,   -0.15730912365699171   ,  &
    !         -0.12330296501726076   ,    -9.1221637319875065E-002 , -5.4208334732890960E-002, &
    !          -3.3953129158204920E-002 , -3.3054086165578984E-002,  -3.3477693851654014E-002 , &
    !          -4.1210468931201311E-002 , -4.8205085480361655E-002,  -3.9102480073602473E-002 , &
    !          -7.4241829666203130E-003  , 4.8846625370967539E-002 , 0.12955167100955539     /)

    ! velx1 = (/0.0679113,-0.0506645,-0.165054,-0.259275,-0.323055,-0.351683,-0.34497,&
    ! -0.307675,-0.246814,-0.20387,-0.179822,-0.161136,-0.161364,-0.155433,-0.125498,&
    ! -0.0739787,-0.0113531,0.0418549,0.00165605,-0.140511,-0.288316,-0.400703,-0.46464,&
    ! -0.482198,-0.486472,-0.517742,-0.554531,-0.616727,-0.701399,-0.749291,-0.752127,&
    ! -0.703099,-0.599818,-0.444353,-0.244792,-0.0117645/)
    ! vely1 = (/0.230144,0.255665,0.239088,0.189377,0.117213,0.0345924,-0.0460268,&
    ! -0.112563,-0.150632,-0.16313,-0.169162,-0.174522,-0.187745,-0.217349,-0.252865,&
    ! -0.27929,-0.286706,-0.277196,-0.259792,-0.234537,-0.202028,-0.16985,-0.141694,&
    ! -0.114553,-0.0835294,-0.0505825,-0.0183377,0.0106717,0.0259373,0.0253131,&
    ! 0.0176296,0.0118827,0.0172345,0.0414225,0.0899208,0.163296/)

    !isolated
    ! velx1 = (/0.130053,0.0172343,-0.100554,-0.206186,-0.285872,-0.331198,&
    ! -0.339091,-0.311855,-0.257649,-0.198857,-0.17342,-0.153704,-0.144239,&
    ! -0.15129,-0.130448,-0.0854709,-0.0249459,0.0319706,-0.0110843,-0.163447,&
    ! -0.310304,-0.414919,-0.467859,-0.468693,-0.493591,-0.531481,-0.574983,&
    ! -0.664414,-0.742709,-0.779144,-0.763469,-0.691039,-0.562558,-0.385705,&
    ! -0.175121,0.0539504/)
    ! vely1 = (/0.307092,0.348335,0.345646,0.304658,0.234351,0.146861,0.0558436,&
    ! -0.0246643,-0.0807229,-0.104103,-0.111845,-0.119991,-0.129019,-0.153226,-0.191164,&
    ! -0.222886,-0.236586,-0.230104,-0.212217,-0.187211,-0.156823,-0.127404,-0.101221,&
    ! -0.0729178,-0.0406995,-0.00845896,0.0243131,0.0494383,0.057264,0.0517809,0.0422214,&
    ! 0.0387825,0.0504856,0.0839942,0.143082,0.227117/)
    !
    ! velx2 = (/0.0418549,-0.0113531,-0.0739787,-0.125498,-0.155433,-0.161364,&
    ! -0.161136,-0.179822,-0.20387,-0.246814,-0.307676,-0.34497,-0.351683,-0.323055,&
    ! -0.259275,-0.165054,-0.0506645,0.0679113,-0.0117645,-0.244792,-0.444353,-0.599818,&
    ! -0.703099,-0.752127,-0.749291,-0.701399,-0.616727,-0.554531,-0.517742,-0.486472,&
    ! -0.482198,-0.46464,-0.400703,-0.288316,-0.140511,0.00165605/)
    ! vely2 = (/0.277196,0.286706,0.27929,0.252865,0.217349,0.187745,0.174522,&
    ! 0.169162,0.16313,0.150632,0.112563,0.0460268,-0.0345924,-0.117213,-0.189377,&
    ! -0.239088,-0.255665,-0.230144,-0.163296,-0.0899208,-0.0414225,-0.0172345,&
    ! -0.0118827,-0.0176296,-0.0253131,-0.0259374,-0.0106718,0.0183378,0.0505825,&
    ! 0.0835294,0.114553,0.141694,0.16985,0.202028,0.234537,0.259792/)

    call radialforce(p,f,uvec+velx/velf,vvec+vely/velf,thetavec,af_data,cl_data,cd_data,(dia0/2.0_dp),chord,&
      twist,delta,blades,rot0,velf,velx,vely,rho,mu,interp,q,ka,CTo,CPo,Rp,Tp,Zp)
    ! print *, 'Cp', CPo

    ! print *, 'Power Calculation Here'
    ! call radialforce(p,f,velx1,vely1,thetavec,af_data,cl_data,cd_data,(dia0/2.0_dp),chord,&
    !   twist,delta,blades,abs(rot),velf,velx*0.0_dp,vely*0.0_dp,rho,mu,interp,q,ka,CTo,CPo1,Rp,Tp,Zp)
    ! ! call radialforce(p,f,velx2,vely2,thetavec,af_data,cl_data,cd_data,(dia0/2.0_dp),chord,&
    ! !   twist,delta,blades,-abs(rot0),velf,velx*0.0_dp,vely*0.0_dp,rho,mu,interp,q,ka,CTo,CPo2,Rp,Tp,Zp)
    ! CPo = CPo1
    ! ! print *, 'Cp1', CPo1
    ! ! print *, 'Cp2', CPo2

    power_tot = (0.5_dp*rho*velf**3)*(dia0*H)*CPo

end subroutine powercalc


! ! Calculating power and coefficient of power of single or multiple turbines
! subroutine velcalc(t,f,x,y,dia,rot,velf,loc1,loc2,loc3,spr1,spr2,&
!   skw1,skw2,scl1,scl2,scl3,af_data,cl_data,cd_data,chord,twist,delta,blades,&
!   H,rho,mu,uvec,vvec,interp,power_tot,CPo)
!     implicit none
!
!     integer, parameter :: dp = kind(0.d0)
!
!     ! in
!     integer, intent(in) :: t,f,p,interp
!     real(dp), dimension(t), intent(in) :: x,y,dia,rot
!     real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
!     real(dp), intent(in) :: velf,chord,twist,delta,blades,H,rho,mu
!
!     ! out
!     real(dp), intent(out) :: veleff
!
!     ! local
!     integer :: m,n,i,j,inte
!     real(dp) :: x0,y0,dia0,rot0,pi,ka,CTo
!     real(dp), dimension(t-1) :: xt,yt,diat,rott
!     real(dp), dimension(p) :: velx,vely,thetavec,q,Rp,Tp,Zp
!     pi = 3.1415926535897932_dp
!
!     m = 220
!     n = 200
!     inte = 1 ! 2D Simpson's Rule
!     ! inte = 2 ! 2D Trapezoidal Rule
!
!     do i = 1,p
!       thetavec(i) = (2.0_dp*pi/p)*i-(2.0_dp*pi/p)/2.0_dp
!     end do
!
!     ! Velocity and power measurements made at first turbine
!     if (t .ge. 3) then ! needing to calculate wake interaction using sum of squares
!       x0 = x(1)
!       y0 = y(1)
!       dia0 = dia(1)
!       rot0 = rot(1)
!       do j = 1,t-1
!         xt(j) = x(j+1)
!         yt(j) = y(j+1)
!         diat(j) = dia(j+1)
!         rott(j) = rot(j+1)
!       end do
!       call overlap(t-1,p,xt,yt,diat,rott,chord,blades,x0,y0,dia0,velf,loc1,loc2,loc3,&
!       spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
!     else if (t .eq. 2) then ! calculating wake effect from one other turbine
!       x0 = x(1)
!       y0 = y(1)
!       dia0 = dia(1)
!       rot0 = rot(1)
!       xt = x(2)
!       yt = y(2)
!       diat = dia(2)
!       rott = rot(2)
!       call coupled(p,xt,yt,diat,rott,chord,blades,x0,y0,dia0,velf,loc1,loc2,loc3,&
!         spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,inte,velx,vely)
!     else ! no need to calculate wake interaction
!       x0 = x(1)
!       y0 = y(1)
!       dia0 = dia(1)
!       rot0 = rot(1)
!       velx = 0.0_dp
!       vely = 0.0_dp
!     end if
!
!     ! call radialforce(p,f,uvec,vvec,thetavec,af_data,cl_data,cd_data,(dia0/2.0_dp),chord,&
!     !   twist,delta,blades,rot0,velf,velx,vely,rho,mu,interp,q,ka,CTo,CPo,Rp,Tp,Zp)
!
!     power_tot = (0.5_dp*rho*velf**3)*(dia0*H)*CPo
!
! end subroutine velcalc


! To build for Python interface (Mac): f2py -c  --opt=-O2 -m _vawtwake VAWT_Wake_Model.f90
