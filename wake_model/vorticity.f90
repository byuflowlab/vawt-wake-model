! To build for Python interface: f2py -c  --opt=-O2 -m _vortmodel vorticity.f90

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
    *(1.0_dp-erf((loc + skw*spr**2.0_dp - y)/(sqrt(2.0_dp)*spr)))

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


    if (loc2 < 0.0_dp) then
      loc = loc1*xd*xd + 0.01_dp*xd + loc3
    else
      loc = loc1*xd*xd + loc2*xd + loc3
    end if
    spr = spr1*xd + spr2
    skw = skw1*xd + skw2

    if (scl2 < 0.05_dp) then
      scl = scl1/(1.0_dp + exp(0.05_dp*(xd - scl3)))
    else
      scl = scl1/(1.0_dp + exp(scl2*(xd - scl3)))
    end if

    ! Limiting the fits to maximum values the Modified Gaussian Distribution can handle
    if (loc < 0.2_dp) then
      loc = 0.2_dp
    end if
    if (spr < -0.5_dp) then
      spr = -0.5_dp
    end if
    if (skw > 0.0_dp) then
      skw = 0.0_dp
    end if

    call dist(yd,loc,spr,skw,scl,g1)
    call dist(yd,-loc,-spr,-skw,-scl,g2)

    gam_lat = (g1 - g2)

    ! if (gam_lat /= gam_lat) then
    !   print *, gam_lat, loc, spr, skw, scl
    ! end if

end subroutine gamma

! Creating the integral for VAWT_Wake_Model.py
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
    call gamma(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gammav)

    inte = gammav*((y - y0)/((x - x0)*(x - x0) + (y - y0)*(y - y0)))

end subroutine integrandx

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
    call gamma(x,y,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,gammav)

    inte = gammav*((x0 - x)/((x - x0)*(x - x0) + (y - y0)*(y - y0)))

end subroutine integrandy

subroutine velcalc(x0,y0,dia,rot,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,&
  velf,a,b,c,d,m,n,vel)
    ! http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html
    implicit none

    integer, parameter :: dp = kind(0.d0)

    ! in
    real(dp), intent(in) :: x0,y0,dia,rot
    real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
    real(dp), intent(in) :: velf,a,b,c,d
    integer, intent(in) :: m,n

    ! out
    real(dp), intent(out) :: vel

    ! local
    integer :: i,j
    real(dp) :: h,k,velx,vely,xsum,ysum,intlim,pi2
    real(dp), dimension(2*m) :: xdiv
    real(dp), dimension(2*n) :: ydiv
    real(dp) :: xval1,xval2,xval3,xval4,yval1,yval2,yval3,yval4
    real(dp) :: xsum1,xsum2,xsum3,xsum4,xsum5,xsum6,xsum7,xsum8,xsum9,xsum10,xsum11,xsum12
    real(dp) :: ysum1,ysum2,ysum3,ysum4,ysum5,ysum6,ysum7,ysum8,ysum9,ysum10,ysum11,ysum12
    intrinsic sqrt
    pi2 = 6.28318530718_dp

    h = (b - a)/(2.0_dp*m)
    k = (d - c)/(2.0_dp*n)
    ! print *, h,k
    do i = 1,2*m
      xdiv(i) = a + i*h
    end do

    do j = 1,2*n
      ydiv(j) = c + j*k
    end do
    ! print *, 'xdiv',xdiv
    ! print *, 'ydiv',ydiv

    velx = 0.0_dp
    vely = 0.0_dp

    call integrandx(c,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval1)
    if (xval1 .ne. xval1) then
      print *, 'xval1'
    end if
    call integrandx(d,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval2)
    if (xval2 .ne. xval2) then
      print *, 'xval2'
    end if
    call integrandx(c,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval3)
    if (xval3 .ne. xval3) then
      print *, 'xval3'
    end if
    call integrandx(d,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xval4)
    if (xval4 .ne. xval4) then
      print *, 'xval4'
    end if

    call integrandy(c,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval1)
    if (yval1 .ne. yval1) then
      print *, 'yval1'
    end if
    call integrandy(d,a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval2)
    if (yval2 .ne. yval2) then
      print *, 'yval2'
    end if
    call integrandy(c,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval3)
    if (yval3 .ne. yval3) then
      print *, 'yval3'
    end if
    call integrandy(d,b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yval4)
    if (yval4 .ne. yval4) then
      print *, 'yval4'
    end if

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
        ! if (xsum .eq. xsum) then
          xsum1 = xsum1 + xsum
        ! end if
        call integrandy(ydiv(2*j-1),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ! if (ysum .eq. ysum) then
          ysum1 = ysum1 + ysum
        ! end if
      end if
    end do
    do j = 1,n-1
      if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
        call integrandx(ydiv(2*j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        ! if (xsum .eq. xsum) then
          xsum2 = xsum2 + xsum
        ! end if
        call integrandy(ydiv(2*j),a,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ! if (ysum .eq. ysum) then
          ysum2 = ysum2 + ysum
        ! end if
      end if
    end do
    do j = 1,n
      if ((ydiv(2*j-1) .gt. intlim) .or. (ydiv(2*j-1) .lt. -intlim)) then
        call integrandx(ydiv(2*j-1),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        ! if (xsum .eq. xsum) then
          xsum3 = xsum3 + xsum
        ! end if
        call integrandy(ydiv(2*j-1),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ! if (ysum .eq. ysum) then
          ysum3 = ysum3 + ysum
        ! end if
      end if
    end do
    do j = 1,n-1
      if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
        call integrandx(ydiv(2*j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
        ! if (xsum .eq. xsum) then
          xsum4 = xsum4 + xsum
        ! end if
        call integrandy(ydiv(2*j),b,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
        ! if (ysum .eq. ysum) then
          ysum4 = ysum4 + ysum
        ! end if
      end if
    end do
    do i = 1,m
      call integrandx(c,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
      ! if (xsum .eq. xsum) then
        xsum5 = xsum5 + xsum
      ! end if
      call integrandy(c,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
      ! if (ysum .eq. ysum) then
        ysum5 = ysum5 + ysum
      ! end if
    end do
    do i = 1,m
      call integrandx(d,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
      ! if (xsum .eq. xsum) then
        xsum6 = xsum6 + xsum
      ! end if
      call integrandy(d,xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
      ! if (ysum .eq. ysum) then
        ysum6 = ysum6 + ysum
      ! end if
    end do
    do i = 1,m-1
      call integrandx(c,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
      ! if (xsum .eq. xsum) then
        xsum7 = xsum7 + xsum
      ! end if
      call integrandy(c,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
      ! if (ysum .eq. ysum) then
        ysum7 = ysum7 + ysum
      ! end if
    end do
    do i = 1,m-1
      call integrandx(d,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
      ! if (xsum .eq. xsum) then
        xsum8 = xsum8 + xsum
      ! end if
      call integrandy(d,xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
      ! if (ysum .eq. ysum) then
        ysum8 = ysum8 + ysum
      ! end if
    end do

    do j = 1,n
      do i = 1,m
        if ((ydiv(2*j-1) .gt. intlim) .or. (ydiv(2*j-1) .lt. -intlim)) then
          call integrandx(ydiv(2*j-1),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          ! if (xsum .eq. xsum) then
            xsum9 = xsum9 + xsum
          ! end if
          call integrandy(ydiv(2*j-1),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ! if (ysum .eq. ysum) then
            ysum9 = ysum9 + ysum
          ! end if
        end if
      end do
    end do
    do j = 1,n-1
      do i = 1,m
        if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
          call integrandx(ydiv(2*j),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          ! if (xsum .eq. xsum) then
            xsum10 = xsum10 + xsum
          ! end if
          call integrandy(ydiv(2*j),xdiv(2*i-1),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ! if (ysum .eq. ysum) then
            ysum10 = ysum10 + ysum
          ! end if
        end if
      end do
    end do
    do j = 1,n
      do i = 1,m-1
        if ((ydiv(2*j-1) .gt. intlim) .or. (ydiv(2*j-1) .lt. -intlim)) then
          call integrandx(ydiv(2*j-1),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          ! if (xsum .eq. xsum) then
            xsum11 = xsum11 + xsum
          ! end if
          call integrandy(ydiv(2*j-1),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ! if (ysum .eq. ysum) then
            ysum11 = ysum11 + ysum
          ! end if
        end if
      end do
    end do
    do j = 1,n-1
      do i = 1,m-1
        if ((ydiv(2*j) .gt. intlim) .or. (ydiv(2*j) .lt. -intlim)) then
          call integrandx(ydiv(2*j),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xsum)
          ! if (xsum .eq. xsum) then
            xsum12 = xsum12 + xsum
          ! end if
          call integrandy(ydiv(2*j),xdiv(2*i),x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,ysum)
          ! if (ysum .eq. ysum) then
            ysum12 = ysum12 + ysum
          ! end if
        end if
      end do
    end do

    if (xsum1 .ne. xsum1) then
      print *, 'xsum1'
    end if
    if (xsum2 .ne. xsum2) then
      print *, 'xsum2'
    end if
    if (xsum3 .ne. xsum3) then
      print *, 'xsum3'
    end if
    if (xsum4 .ne. xsum4) then
      print *, 'xsum4'
    end if
    if (xsum5 .ne. xsum5) then
      print *, 'xsum5'
    end if
    if (xsum6 .ne. xsum6) then
      print *, 'xsum6'
    end if
    if (xsum7 .ne. xsum7) then
      print *, 'xsum7'
    end if
    if (xsum8 .ne. xsum8) then
      print *, 'xsum8'
    end if
    if (xsum9 .ne. xsum9) then
      print *, 'xsum9'
    end if
    if (xsum10 .ne. xsum10) then
      print *, 'xsum10'
    end if
    if (xsum11 .ne. xsum11) then
      print *, 'xsum11'
    end if
    if (xsum12 .ne. xsum12) then
      print *, 'xsum12'
    end if

    if (ysum1 .ne. ysum1) then
      print *, 'ysum1'
    end if
    if (ysum2 .ne. ysum2) then
      print *, 'ysum2'
    end if
    if (ysum3 .ne. ysum3) then
      print *, 'ysum3'
    end if
    if (ysum4 .ne. ysum4) then
      print *, 'ysum4'
    end if
    if (ysum5 .ne. ysum5) then
      print *, 'ysum5'
    end if
    if (ysum6 .ne. ysum6) then
      print *, 'ysum6'
    end if
    if (ysum7 .ne. ysum7) then
      print *, 'ysum7'
    end if
    if (ysum8 .ne. ysum8) then
      print *, 'ysum8'
    end if
    if (ysum9 .ne. ysum9) then
      print *, 'ysum9'
    end if
    if (ysum10 .ne. ysum10) then
      print *, 'ysum10'
    end if
    if (ysum11 .ne. ysum11) then
      print *, 'ysum11'
    end if
    if (ysum12 .ne. ysum12) then
      print *, 'ysum12'
    end if




    velx = (h*k*(xval1 + xval2 + xval3 + xval4 + 4.0_dp*xsum1 + 2.0_dp*xsum2 +&
    4.0_dp*xsum3 + 2.0_dp*xsum4 + 4.0_dp*xsum5 + 4.0_dp*xsum6 + 2.0_dp*xsum7 + 2.0_dp*xsum8 +&
    16.0_dp*xsum9 + 8.0_dp*xsum10 + 8.0_dp*xsum11 + 4.0_dp*xsum12)/9.0_dp)*(rot/pi2)

    vely = (h*k*(yval1 + yval2 + yval3 + yval4 + 4.0_dp*ysum1 + 2.0_dp*ysum2 +&
    4.0_dp*ysum3 + 2.0_dp*ysum4 + 4.0_dp*ysum5 + 4.0_dp*ysum6 + 2.0_dp*ysum7 + 2.0_dp*ysum8 +&
    16.0_dp*ysum9 + 8.0_dp*ysum10 + 8.0_dp*ysum11 + 4.0_dp*ysum12)/9.0_dp)*(rot/pi2)

    if (velx .ne. velx) then
      print *, 'velx'
    end if

    if (vely .ne. vely) then
      print *, 'vely'
    end if


    ! print *, xval1
    ! print *, (velx+velf)/velf
    ! print *, vely/velf

    vel = sqrt((velx + velf)*(velx + velf) + vely*vely)
    ! vel = (velx+velf)

end subroutine velcalc


subroutine radialforce(y,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,inte)
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

    inte = gammav*((x0 - x)/((x - x0)*(x - x0) + (y - y0)*(y - y0)))

end subroutine radialforce


! To build for Python interface: f2py -c  --opt=-O2 -m _vortmodel vorticity.f90
