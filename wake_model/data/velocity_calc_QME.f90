! module load python/2/7
! module load compiler_gnu/4.9.2
! f2py -c  --opt=-O2 -m _velcalcqme velocity_calc_QME.f90

subroutine overlay(xt,ys,coef,fun)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: xt,ys
  real(dp), dimension(10), intent(in) :: coef
  ! out
  real(dp), intent(out) :: fun
  ! local

  fun = coef(1) + coef(2)*xt + coef(3)*ys + coef(4)*xt**2 + coef(5)*xt*ys + &
  coef(6)*ys**2 + coef(7)*xt**3 + coef(8)*xt**2*ys + coef(9)*xt*ys**2 + &
  coef(10)*ys**3

end subroutine overlay


subroutine veldist(dn,lat,spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3,vel)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: dn,lat,spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3
  ! out
  real(dp), intent(out) :: vel
  ! local
  real(dp) :: pow_v,exp_v,quad_v,scl_v
  intrinsic exp
  intrinsic abs

  pow_v = pow1-pow2*dn**pow3
  exp_v = exp(-spr1*abs(lat)**pow_v)
  quad_v = spr2*(lat-skw)**2 - 1.0_dp
  scl_v = scl3*scl2*scl1*exp(scl2*dn)*exp(-scl1*exp(scl2*dn))

  vel = exp_v*quad_v*scl_v + 1.0_dp

end subroutine veldist


subroutine sheet(ndata,xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,&
  coef7,coef8,vel)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: ndata
  real(dp), dimension(10), intent(in) :: coef0,coef1,coef2,coef3,coef4,coef5
  real(dp), dimension(10), intent(in) :: coef6,coef7,coef8
  real(dp), dimension(ndata), intent(in) :: xttr,ystr,posdn,poslt
  ! out
  real(dp), dimension(ndata), intent(out) :: vel
  ! local
  integer :: i
  real(dp) :: spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3

  do i = 1,ndata
    call overlay(xttr(i),ystr(i),coef0,spr1)
    call overlay(xttr(i),ystr(i),coef1,pow1)
    call overlay(xttr(i),ystr(i),coef2,pow2)
    call overlay(xttr(i),ystr(i),coef3,pow3)
    call overlay(xttr(i),ystr(i),coef4,spr2)
    call overlay(xttr(i),ystr(i),coef5,skw)
    call overlay(xttr(i),ystr(i),coef6,scl1)
    call overlay(xttr(i),ystr(i),coef7,scl2)
    call overlay(xttr(i),ystr(i),coef8,scl3)

    if (spr1 < 0.0_dp) then
      spr1 = 0.0_dp
    end if
    if (pow1 < 0.0_dp) then
      pow1 = 0.0_dp
    end if
    if (pow2 < 0.0_dp) then
      pow2 = 0.0_dp
    end if
    if (pow3 < 0.0_dp) then
      pow3 = 0.0_dp
    end if
    if (spr2 < 0.0_dp) then
      spr2 = 0.0_dp
    end if
    if (scl1 < 0.0_dp) then
      scl1 = 0.0_dp
    end if
    if (scl2 < 0.0_dp) then
      scl2 = 0.0_dp
    end if
    if (scl3 < 0.0_dp) then
      scl3 = 0.0_dp
    end if

    call veldist(posdn(i),poslt(i),spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3,&
      vel(i))

  end do

end subroutine sheet
