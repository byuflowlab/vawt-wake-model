! To build for Python interface: f2py -c  --opt=-O2 -m _vawtwake VAWT_Wake_Model.f90



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


! Performing integration to convert vorticity into velocity
subroutine vel_field(xt,yt,x0t,y0t,dia,rot,chord,blades,loc1d,loc2d,loc3d,spr1d,spr2d,&
  skw1d,skw2d,scl1d,scl2d,scl3d,velf,m_in,n_in,inte,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: xt,yt,x0t,y0t,dia,rot,chord,blades
    real(dp), dimension(10), intent(in) :: loc1d,loc2d,loc3d,spr1d,spr2d,skw1d,skw2d,scl1d,scl2d,scl3d
    real(dp), intent(in) :: velf
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
    pi2 = 6.28318530718_dp

    tsr = (dia/2.0_dp)*rot/velf
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
      16.0_dp*xsum9 + 8.0_dp*xsum10 + 8.0_dp*xsum11 + 4.0_dp*xsum12)/9.0_dp)*(rot/pi2)

      velyi = (h*k*(yval1 + yval2 + yval3 + yval4 + 4.0_dp*ysum1 + 2.0_dp*ysum2 +&
      4.0_dp*ysum3 + 2.0_dp*ysum4 + 4.0_dp*ysum5 + 4.0_dp*ysum6 + 2.0_dp*ysum7 + 2.0_dp*ysum8 +&
      16.0_dp*ysum9 + 8.0_dp*ysum10 + 8.0_dp*ysum11 + 4.0_dp*ysum12)/9.0_dp)*(rot/pi2)

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
      2.0_dp*xsum3 + 2.0_dp*xsum4 + 4.0_dp*xsum5)/4.0_dp)*(rot/pi2)

      velyi = (h*k*(yval1 + yval2 + yval3 + yval4 + 2.0_dp*ysum1 + 2.0_dp*ysum2 +&
      2.0_dp*ysum3 + 2.0_dp*ysum4 + 4.0_dp*ysum5)/4.0_dp)*(rot/pi2)

    end if

    ! ****************************************************************************

    velx = velxi/velf
    vely = velyi/velf

end subroutine vel_field


! Calculating loading forces on VAWT blades and power of the turbine
subroutine radialforce(n,f,uvec,vvec,thetavec,af_data,cl_data,cd_data,r,chord,&
  twist,delta,B,Omega,Vinf,rho,mu,q,ka,CTo,CPo,Rp,Tp,Zp)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: n,f
    real(dp), dimension(n), intent(in) :: uvec,vvec,thetavec
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
      Vn(i) = Vinf*(1.0_dp + uvec(i))*sin(thetavec(i)) - Vinf*vvec(i)*cos(thetavec(i))
      Vt(i) = rotation*(Vinf*(1.0_dp + uvec(i))*cos(thetavec(i)) + Vinf*vvec(i)*&
      sin(thetavec(i))) + abs(Omega)*r
      W(i) = sqrt(Vn(i)**2 + Vt(i)**2)
      phi(i) = atan2(Vn(i), Vt(i))
      alpha(i) = phi(i) - twist
      ! Re = rho*W*chord/mu  ! currently no Re dependence

      ! airfoil
      call interpolate(f,af_data,cl_data,alpha(i)*180.0_dp/pi,cl(i))
      call interpolate(f,af_data,cd_data,alpha(i)*180.0_dp/pi,cd(i))

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
      integrand(i) = (W(i)/Vinf)**2 * (cn(i)*sin(thetavec(i)) - rotation*ct(i)*&
      cos(thetavec(i))/cos(delta))
    end do
    ! print *,'uvec',uvec
    ! print *,'vvec',vvec
    ! print *,'thetavec',thetavec
    ! print *,'Vn',Vn
    ! print *,'Vt',Vt
    ! print *,'W',W
    ! print *,'phi',phi
    ! print *,'alpha',alpha
    ! print *,'cl',cl
    ! print *,'cd',cd
    ! print *,'cn',cn
    ! print *,'ct',ct
    ! print *,'q',q
    ! print *,'qdyn',qdyn
    ! print *,'Rp',Rp
    ! print *,'Tp',Tp
    ! print *,'Zp',Zp
    ! print *,'integrand',integrand


    call pInt(n,thetavec,integrand,CTend)
    CTo = sigma/(4.0_dp*pi) * CTend
    if (CTo .gt. 2.0_dp) then
      a = 0.5_dp*(1.0_dp + sqrt(1.0_dp + CTo))
      ka = 1.0_dp / (a-1.0_dp)
    elseif (CTo .gt. 0.96_dp) then
      a = 1.0_dp/7.0_dp*(1.0_dp + 3.0_dp*sqrt(7.0_dp/2.0_dp*CTo - 3.0_dp))
      ka = 18.0_dp*a / (7.0_dp*a**2 - 2.0_dp*a + 4.0_dp)
    else
      a = 0.5_dp*(1.0_dp - sqrt(1.0_dp - CTo))
      ka = 1.0_dp / (1.0_dp-a)
    end if

    ! power coefficient
    H = 1.0_dp  ! per unit height
    Sref = 2.0_dp*r*H
    do i = 1,n
      Qp(i) = r*Tp(i)
    end do

    call pInt(n,thetavec,Qp,Pend)
    P = abs(Omega)*B/(2.0_dp*pi)*Pend
    CPo = P / (0.5_dp*rho*Vinf**3 * Sref)

end subroutine radialforce


! Residual function used for finding the root for solving the velocity field
! subroutine residual(nw,ntheta,nturb,naf,w,A,theta,k,af_data,cl_data,cd_data,r,chord,&
!   twist,delta,B,Omega,Vinf,rho,mu,q,ka,CTo,CPo,Rp,Tp,Zp)
!     implicit none
!
!     integer, parameter :: dp = kind(0.d0)
!
!     ! in
!     integer, intent(in) :: nw,ntheta,nturb,naf
!     real(dp), dimension(nw), intent(in) :: w
!     real(dp), dimension(naf), intent(in) :: af_data,cl_data,cd_data
!     real(dp), intent(in) :: r,chord,twist,delta,B,Omega,Vinf,rho,mu
!     integer, intent(in) :: m,n,int
!
!     ! out
!     real(dp), dimension(nw), intent(out) :: res
!
!     ! local
!     integer :: i,j,k,l
!     real(dp) :: pi,theta
!     real(dp), dimension(36) :: xd,yd,velx_int,vely_int,intex,intey
!     intrinsic sin
!     intrinsic cos
!     intrinsic sqrt
!     intrinsic abs
!     pi = 3.1415926535897932_dp
!
!     q = 0.0_dp
!     ka = 0.0_dp
!
!     do i = 1,nturb
!
!
!
!
!   function residual(w::Array{Float64,1}, A::Array{Float64,2}, theta::Array{Float64,1},
!       k::Array{Float64,1}, turbines::Array{Turbine,1}, env::Environment)
!
!       # setup
!       ntheta = length(theta)
!       nturbines = length(turbines)  #  int(length(w)/2/ntheta)
!       q = zeros(ntheta*nturbines)
!       ka = 0.0
!
!       for i in eachindex(turbines)
!           idx = (i-1)*ntheta+1:i*ntheta
!
!           u = w[idx]
!           v = w[ntheta*nturbines + idx]
!
!           q[idx], ka, _, _, _, _, _ = radialforce(u, v, theta, turbines[i], env)
!       end
!
!       if nturbines == 1  # if only one turbine use the k from the analysis
!           k = [ka]
!       end  # otherwise, use k that was input to this function
!
!       # reformat to multiply in correct locations
!       kmult = repeat(k, inner=[ntheta])
!       kmult = [kmult; kmult]
!
!       return (A*q).*kmult - w
!   end
!
! end subroutine residual


! Calculating effective velocities around given turbine due to wake interaction
subroutine overlap(t,xt,yt,diat,rot,chord,blades,x0,y0,dia,velf,loc1,loc2,loc3,spr1,spr2,&
  skw1,skw2,scl1,scl2,scl3,m,n,int,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t
    real(dp), dimension(t), intent(in) :: xt,yt,diat,rot
    real(dp), intent(in) :: x0,y0,dia,velf,chord,blades
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
    integer, intent(in) :: m,n,int

    ! out
    real(dp), dimension(36), intent(out) :: velx,vely

    ! local
    integer :: i,j,k,l
    real(dp) :: pi,theta
    real(dp), dimension(36) :: xd,yd,velx_int,vely_int,intex,intey
    intrinsic sin
    intrinsic cos
    intrinsic sqrt
    intrinsic abs
    pi = 3.1415926535897932_dp

    intex = 0.0_dp
    intey = 0.0_dp

    ! finding points around the flight path of the blades
    do i = 1,36
      theta = (2.0_dp*pi/36.0_dp)*i
      xd(i) = x0 + cos(theta)*(dia/2.0_dp)
      yd(i) = y0 + sin(theta)*(dia/2.0_dp)
    end do

    ! sum of squares of velocity deficits
    do j = 1,t
      do k = 1,36
        call vel_field(xt(j),yt(j),xd(k),yd(k),diat(j),rot(j),chord,blades,&
        loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,velf,m,n,int,&
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
    do l = 1,36
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
subroutine overlappoint(t,xt,yt,diat,rot,chord,blades,x0,y0,dia,velf,loc1,loc2,loc3,spr1,spr2,&
  skw1,skw2,scl1,scl2,scl3,m,n,int,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t
    real(dp), dimension(t), intent(in) :: xt,yt,diat,rot
    real(dp), intent(in) :: x0,y0,dia,velf,chord,blades
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
    integer, intent(in) :: m,n,int

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
      call vel_field(xt(i),yt(i),x0,y0,diat(i),rot(i),chord,blades,&
      loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,velf,m,n,int,&
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


! Creating velocity vectors of an isolated turbine (free stream wind speed)
subroutine isolate(x0,y0,dia,chord,blades,Omega,velf,loc1,loc2,loc3,&
  spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,int,velx,vely)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x0,y0,dia,velf,chord,blades,Omega
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
    integer, intent(in) :: m,n,int

    ! out
    real(dp), dimension(36), intent(out) :: velx,vely

    ! local
    integer :: i,j
    real(dp) :: pi
    real(dp), dimension(36) :: theta,xd,yd!,velx_int,vely_int
    intrinsic sin
    intrinsic cos
    pi = 3.1415926535897932_dp

    do i = 1,36
      theta(i) = (2.0_dp*pi/36.0_dp)*(i+1)-(2.0_dp*pi/36.0_dp)/2.0_dp
      xd(i) = x0 + cos(theta(i))*(dia/2.0_dp)
      yd(i) = y0 + sin(theta(i))*(dia/2.0_dp)
    end do

    do j = 1,36
      call vel_field(0.0_dp,0.0_dp,xd(j),yd(j),dia,Omega,chord,blades,&
      loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,velf,m,n,int,&
      velx(j),vely(j))
      velx(j) = velx(j)! + Omega*(dia/2.0_dp)*cos(theta(j))/velf
      vely(j) = vely(j)! + Omega*(dia/2.0_dp)*sin(theta(j))/velf
    end do

    ! do i = 1,36
    !   velx(i) = 1.0_dp
    !   vely(i) = 0.
    ! end do
    ! velx = (/0.459751,0.282597,0.0826669,-0.113647,-0.289426,-0.431141,-0.521327,&
    ! -0.560051,-0.537055,-0.482605,-0.420738,-0.356635,-0.292542,-0.225452,-0.133013,&
    ! -0.0208681,0.0851436,0.165534,0.0846263,-0.140613,-0.363099,-0.549279,-0.677296,&
    ! -0.76224,-0.837259,-0.910581,-0.974004,-1.03006,-1.04345,-0.988548,-0.866052,&
    ! -0.689424,-0.457304,-0.193809,0.0943889,0.377583/)
    ! vely = (/0.385218,0.524883,0.587173,0.583378,0.523866,0.418706,0.283568,0.135876,&
    ! -0.00149835,-0.109984,-0.19216,-0.256426,-0.309558,-0.358509,-0.397126,-0.407762,&
    ! -0.384917,-0.342492,-0.28303,-0.195373,-0.093928,-0.00860521,0.046295,0.0793963,&
    ! 0.099107,0.104145,0.0896153,0.0498272,-0.0163959,-0.0914874,-0.155433,-0.189781,&
    ! -0.179738,-0.119601,-0.00105869,0.183809/)

end subroutine isolate


! Calculating power and coefficient of power of single or multiple turbines
subroutine powercalc(t,f,xt,yt,diat,rot,x0,y0,dia,velf,loc1,loc2,loc3,spr1,spr2,&
  skw1,skw2,scl1,scl2,scl3,af_data,cl_data,cd_data,chord,twist,delta,blades,Omega,&
  H,rho,mu,power_tot,CPo)
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    integer, intent(in) :: t,f
    real(dp), dimension(t), intent(in) :: xt,yt,diat,rot
    real(dp), dimension(f), intent(in) :: af_data,cl_data,cd_data
    real(dp), intent(in) :: x0,y0,dia,velf
    real(dp), dimension(10), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
    real(dp), intent(in) :: chord,twist,delta,blades,Omega,H,rho,mu

    ! out
    real(dp), intent(out) :: power_tot,CPo

    ! local
    integer :: p,m,n,i,int
    real(dp) :: pi,ka,CTo
    real(dp), dimension(36) :: velx,vely,thetavec,q,Rp,Tp,Zp
    intrinsic sin
    intrinsic cos
    intrinsic sqrt
    intrinsic abs
    pi = 3.1415926535897932_dp

    p = 36
    m = 110
    n = 100
    int = 1 ! 2D Simpson's Rule
    ! int = 2 ! 2D Trapezoidal Rule

    do i = 1,36
      thetavec(i) = (2.0_dp*pi/36.0_dp)*i-(2.0_dp*pi/36.0_dp)/2.0_dp
    end do
    ! thetavec = (/0.0872665,0.261799,0.436332,0.610865,0.785398,0.959931,1.13446,&
    ! 1.309,1.48353,1.65806,1.8326,2.00713,2.18166,2.35619,2.53073,2.70526,2.87979,&
    ! 3.05433,3.22886,3.40339,3.57792,3.75246,3.92699,4.10152,4.27606,4.45059,4.62512,&
    ! 4.79966,4.97419,5.14872,5.32325,5.49779,5.67232,5.84685,6.02139,6.19592/)

    ! call overlap(t,xt,yt,diat,rot,chord,blades,x0,y0,dia,velf,loc1,loc2,loc3,spr1,spr2,&
    !   skw1,skw2,scl1,scl2,scl3,m,n,int,velx,vely)

    call isolate(0.0_dp,0.0_dp,dia,chord,blades,Omega,velf,loc1,loc2,loc3,&
      spr1,spr2,skw1,skw2,scl1,scl2,scl3,m,n,int,velx,vely)

    call radialforce(p,f,velx,vely,thetavec,af_data,cl_data,cd_data,dia/2.0_dp,chord,&
        twist,delta,blades,Omega,velf,rho,mu,q,ka,CTo,CPo,Rp,Tp,Zp)

    ! print *, thetavec*180.0_dp/pi
    ! print *, velx
    ! print *, vely

    power_tot = (0.5_dp*rho*velf**3)*(dia*H)*CPo

end subroutine powercalc

! ----------------------------------------
! ----------------------------------------
! ----------------------------------------
! ----------------------------------------
! ---------- helper methods --------------
! ----------------------------------------
! ----------------------------------------
! ----------------------------------------
! ----------------------------------------

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

    inte = 0.0
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
    inte = inte + dtheta * 0.5_dp*(f(1) + f(n))

    integral = inte

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
    do i = 1,n-1
      if (xval .lt. x(i)) then
        yval = y(i-1) + (xval - x(i-1))*((y(i)-y(i-1))/(x(i)-x(i-1)))
        exit
      else if (xval .eq. x(i)) then
        yval = y(i)
        exit
      end if

    end do

end subroutine interpolate


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


! Brent's Method for root finding of actuator cylinder velocity field
! www.fing.edu.uy/if/cursos/fiscomp/extras/numrec/book/f9.pdf
! subroutine brent(n,x1,x2,tol)
!     implicit none
!
!     integer, parameter :: dp = kind(0.d0)
!
!     ! in
!     integer, intent(in) :: n
!     real(dp), intent(in) :: tol
!     real(dp), dimension(10), intent(in) :: x1,x2
!
!     ! out
!     real(dp), intent(out) :: res
!
!     ! local
!     real(dp) :: eps,itmax,iter,a,b,c,d,e,fa,fb,fc,p,q,r,s,tol1,xm
!     intrinsic abs
!
!     a = x1
!     b = x2
!     call func(,fa)
!     call func(,fb)
!     if (((fa .gt. 0.0_dp) .and. (fb .gt. 0.0_dp)) .or. ((fa .lt. 0.0_dp) .and. &
!     (fb .lt. 0.0_dp))) then
!       pause "root must be bracketed for Brent's Method"
!
!     c = b
!     fc = fb
!       ! rename a, b, and c and adjust bounding interval d
!       c = a
!       fc = fa
!       d = b - a
!       e = d
!     end if
!     if (abs(fc) .lt. abs(fb)) then
!       a = b
!       b = c
!       c = a
!       fa = fb
!       fb = fc
!       fc = fa
!     end if
!
!
! end subroutine brent

! To build for Python interface: f2py -c  --opt=-O2 -m _vawtwake VAWT_Wake_Model.f90
