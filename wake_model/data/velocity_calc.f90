!f2py -c  --opt=-O2 -m _velcalc velocity_calc.f90


subroutine veldist(dn,lat,spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3,vel)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: dn,lat,spr1,pow1,pow2,pow3,spr2,skw,scl1,scl2,scl3
  ! out
  real(dp), intent(out) :: vel
  ! local
  real(dp) :: sdv_v,spr_v,f1,f2,pi
  intrinsic exp
  intrinsic sqrt
  intrinsic abs
  intrinsic log
  intrinsic max
  ! constants
  pi = 3.1415926535897932_dp

  pow = pow1-pow2*dn

  scl_v = scl3*scl2*scl1*exp(scl2*dn)*exp(-scl1*exp(scl2*dn))


  sprlog = (log(0.01_dp)/-spr1)**(1.0_dp/pow)
  sprlat = lat/sprlog

  exp_v = (1.0_dp/(1.0_dp-skw))*(-exp(-spr1*fabs(lat)**pow)+skw)/(max(1.0_dp,fabs(sprlat)**pow3))+1.0_dp

  vel = scl_v*exp_v


end subroutine veldist


subroutine sheet(ndata,xttr,ystr,posdn,poslt,coef0,coef1,coef2,coef3,coef4,coef5,coef6,&
  coef7,coef8,coef9,coef10,coef11,coef12,coef13,coef14,coef15,coef16,coef17,coef18,&
  coef19,coef20,coef21,coef22,coef23,coef24,coef25,coef26,coef27,sdv_gom,spr_gom,vel)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  integer, intent(in) :: ndata
  real(dp), intent(in) :: sdv_gom,spr_gom
  real(dp), dimension(4), intent(in) :: coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8
  real(dp), dimension(4), intent(in) :: coef9,coef10,coef11,coef12,coef13,coef14,coef15,coef16
  real(dp), dimension(4), intent(in) :: coef17,coef18,coef19,coef20,coef21,coef22,coef23,coef24
  real(dp), dimension(4), intent(in) :: coef25,coef26,coef27
  real(dp), dimension(ndata), intent(in) :: xttr,ystr,posdn,poslt
  ! out
  real(dp), dimension(ndata), intent(out) :: vel
  ! local
  integer :: i
  real(dp) :: men,sdv1,sdv2,sdv3,sdv4,rat,tns
  real(dp) :: spr1,spr2,spr3,spr4,scl1,scl2,scl3

  do i = 1,ndata
    call overlay(xttr(i),ystr(i),coef0,coef14,men)
    call overlay(xttr(i),ystr(i),coef1,coef15,sdv1)
    call overlay(xttr(i),ystr(i),coef2,coef16,sdv2)
    call overlay(xttr(i),ystr(i),coef3,coef17,sdv3)
    call overlay(xttr(i),ystr(i),coef4,coef18,sdv4)
    call overlay(xttr(i),ystr(i),coef5,coef19,rat)
    call overlay(xttr(i),ystr(i),coef6,coef20,tns)
    call overlay(xttr(i),ystr(i),coef7,coef21,spr1)
    call overlay(xttr(i),ystr(i),coef8,coef22,spr2)
    call overlay(xttr(i),ystr(i),coef9,coef23,spr3)
    call overlay(xttr(i),ystr(i),coef10,coef24,spr4)
    call overlay(xttr(i),ystr(i),coef11,coef25,scl1)
    call overlay(xttr(i),ystr(i),coef12,coef26,scl2)
    call overlay(xttr(i),ystr(i),coef13,coef27,scl3)

    call veldist(posdn(i),poslt(i),men,sdv1,sdv2,sdv3,sdv4,rat,tns,spr1,spr2,spr3,spr4,&
      scl1,scl2,scl3,sdv_gom,spr_gom,vel(i))

  end do

end subroutine sheet
