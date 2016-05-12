module find_fit_module
use minpack, only: lmdif1
implicit none
private
public find_fit

contains

subroutine find_fit(data_x, data_y, data_z, expr, pars)
integer, parameter :: dp = kind(0.d0)
real(dp), intent(in) :: data_x(:), data_y(:), data_z(:)
interface
    function expr(x, pars) result(y)

    implicit none
    integer, parameter :: dp = kind(0.d0)
    real(dp), intent(in) :: x(:), pars(:)
    real(dp) :: y(size(x))
    end function
end interface
real(dp), intent(inout) :: pars(:)

real(dp) :: tol, fvec(size(data_x))
integer :: iwa(size(pars)), info, m, n
real(dp), allocatable :: wa(:)

tol = sqrt(epsilon(1._dp))
m = size(fvec)
n = size(pars)
allocate(wa(m*n + 5*n + m))
call lmdif1(fcn, m, n, pars, fvec, tol, info)!, iwa, wa, size(wa))
if (info /= 1) stop "failed to converge"

contains

subroutine fcn(m, n, x, fvec, iflag)
integer, parameter :: dp = kind(0.d0)
integer, intent(in) :: m, n, iflag
real(dp), intent(in) :: x(n)
real(dp), intent(out) :: fvec(m)
! Suppress compiler warning:
fvec(1) = iflag
fvec = data_z - expr(data_x, data_y, x)
end subroutine

end subroutine

end module
