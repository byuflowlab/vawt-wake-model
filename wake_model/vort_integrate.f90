! f2py -c  --opt=-O2 -m _vort_integrate vort_integrate.f90


subroutine qags ( a, b, epsabs, epsrel, result, &
  x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,x)

!*****************************************************************************80
!
!! QAGS estimates the integral of a function.
!
!  Discussion:
!
!    The routine calculates an approximation RESULT to a definite integral
!      I = integral of F over (A,B),
!    hopefully satisfying
!      || I - RESULT || <= max ( EPSABS, EPSREL * ||I|| ).
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Input, real ( kind = 8 ) EPSABS, EPSREL, the absolute and relative accuracy requested.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of || I - RESULT ||.
!
!    Output, integer ( kind = 8 ) NEVAL, the number of times the integral was evaluated.
!
!    Output, integer ( kind = 8 ) IER, error flag.
!                     ier = 0 normal and reliable termination of the
!                             routine. it is assumed that the requested
!                             accuracy has been achieved.
!                     ier > 0 abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. it is assumed that the
!                             requested accuracy has not been achieved.
!                         = 1 maximum number of subdivisions allowed
!                             has been achieved. one can allow more sub-
!                             divisions by increasing the data value of
!                             limit in qags (and taking the according
!                             dimension adjustments into account).
!                             however, if this yields no improvement
!                             it is advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. if the position of a
!                             local difficulty can be determined (e.g.
!                             singularity, discontinuity within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the sub-
!                             ranges. if possible, an appropriate
!                             special-purpose integrator should be used,
!                             which is designed for handling the type
!                             of difficulty involved.
!                         = 2 the occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             the error may be under-estimated.
!                         = 3 extremely bad integrand behavior occurs
!                             at some  points of the integration
!                             interval.
!                         = 4 the algorithm does not converge. roundoff
!                             error is detected in the extrapolation
!                             table. it is presumed that the requested
!                             tolerance cannot be achieved, and that the
!                             returned result is the best which can be
!                             obtained.
!                         = 5 the integral is probably divergent, or
!                             slowly convergent. it must be noted that
!                             divergence can occur with any other value
!                             of ier.
!                         = 6 the input is invalid, because
!                             epsabs < 0 and epsrel < 0,
!                             result, abserr and neval are set to zero.
!
!  Local Parameters:
!
!           alist     - list of left end points of all subintervals
!                       considered up to now
!           blist     - list of right end points of all subintervals
!                       considered up to now
!           rlist(i)  - approximation to the integral over
!                       (alist(i),blist(i))
!           rlist2    - array of dimension at least limexp+2 containing
!                       the part of the epsilon table which is still
!                       needed for further computations
!           elist(i)  - error estimate applying to rlist(i)
!           maxerr    - pointer to the interval with largest error
!                       estimate
!           errmax    - elist(maxerr)
!           erlast    - error on the interval currently subdivided
!                       (before that subdivision has taken place)
!           area      - sum of the integrals over the subintervals
!           errsum    - sum of the errors over the subintervals
!           errbnd    - requested accuracy max(epsabs,epsrel*
!                       abs(result))
!           *****1    - variable for the left interval
!           *****2    - variable for the right interval
!           last      - index for subdivision
!           nres      - number of calls to the extrapolation routine
!           numrl2    - number of elements currently in rlist2. if an
!                       appropriate approximation to the compounded
!                       integral has been obtained it is put in
!                       rlist2(numrl2) after numrl2 has been increased
!                       by one.
!           small     - length of the smallest interval considered
!                       up to now, multiplied by 1.5
!           erlarg    - sum of the errors over the intervals larger
!                       than the smallest interval considered up to now
!           extrap    - logical variable denoting that the routine is
!                       attempting to perform extrapolation i.e. before
!                       subdividing the smallest interval we try to
!                       decrease the value of erlarg.
!           noext     - logical variable denoting that extrapolation
!                       is no longer allowed (true value)
!
  implicit none

  integer ( kind = 8 ), parameter :: limit = 500

  real ( kind = 8 ) a
  real ( kind = 8 ) abseps
  real ( kind = 8 ) abserr
  real ( kind = 8 ) alist(limit)
  real ( kind = 8 ) area
  real ( kind = 8 ) area1
  real ( kind = 8 ) area12
  real ( kind = 8 ) area2
  real ( kind = 8 ) a1
  real ( kind = 8 ) a2
  real ( kind = 8 ) b
  real ( kind = 8 ) blist(limit)
  real ( kind = 8 ) b1
  real ( kind = 8 ) b2
  real ( kind = 8 ) correc
  real ( kind = 8 ) defabs
  real ( kind = 8 ) defab1
  real ( kind = 8 ) defab2
  real ( kind = 8 ) dres
  real ( kind = 8 ) elist(limit)
  real ( kind = 8 ) epsabs
  real ( kind = 8 ) epsrel
  real ( kind = 8 ) erlarg
  real ( kind = 8 ) erlast
  real ( kind = 8 ) errbnd
  real ( kind = 8 ) errmax
  real ( kind = 8 ) error1
  real ( kind = 8 ) error2
  real ( kind = 8 ) erro12
  real ( kind = 8 ) errsum
  real ( kind = 8 ) ertest
  logical extrap
  integer ( kind = 8 ) id
  integer ( kind = 8 ) ier
  integer ( kind = 8 ) ierro
  integer ( kind = 8 ) iord(limit)
  integer ( kind = 8 ) iroff1
  integer ( kind = 8 ) iroff2
  integer ( kind = 8 ) iroff3
  integer ( kind = 8 ) jupbnd
  integer ( kind = 8 ) k
  integer ( kind = 8 ) ksgn
  integer ( kind = 8 ) ktmin
  integer ( kind = 8 ) last
  logical noext
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) neval
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) nrmax
  integer ( kind = 8 ) numrl2
  real ( kind = 8 ) resabs
  real ( kind = 8 ) reseps
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) rlist(limit)
  real ( kind = 8 ) rlist2(52)
  real ( kind = 8 ) small

  real ( kind = 8 ) x0
  real ( kind = 8 ) y0
  real ( kind = 8 ) dia
  real ( kind = 8 ) loc1
  real ( kind = 8 ) loc2
  real ( kind = 8 ) loc3
  real ( kind = 8 ) spr1
  real ( kind = 8 ) spr2
  real ( kind = 8 ) skw1
  real ( kind = 8 ) skw2
  real ( kind = 8 ) scl1
  real ( kind = 8 ) scl2
  real ( kind = 8 ) scl3
  logical xyrun
  real ( kind = 8 ) x

!
!  The dimension of rlist2 is determined by the value of
!  limexp in QEXTR (rlist2 should be of dimension
!  (limexp+2) at least).
!
!  Test on validity of parameters.
!
  ier = 0
  neval = 0
  last = 0
  result = 0.0E+00
  abserr = 0.0E+00
  alist(1) = a
  blist(1) = b
  rlist(1) = 0.0E+00
  elist(1) = 0.0E+00

  if ( epsabs < 0.0E+00 .and. epsrel < 0.0E+00 ) then
    ier = 6
    return
  end if
!
!  First approximation to the integral.
!
  ! print *, 'first xyrun', xyrun
  ierro = 0
  call qk21 ( a, b, result, abserr, defabs, resabs, &
    x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,x)
!
!  Test on accuracy.
!
  dres = abs ( result )
  errbnd = max ( epsabs, epsrel * dres )
  last = 1
  rlist(1) = result
  elist(1) = abserr
  iord(1) = 1

  if ( abserr <= 1.0E+02 * epsilon ( defabs ) * defabs .and. &
    abserr > errbnd ) then
    ier = 2
  end if

  if ( limit == 1 ) then
    ier = 1
  end if

  if ( ier /= 0 .or. (abserr <= errbnd .and. abserr /= resabs ) .or. &
    abserr == 0.0E+00 ) go to 140
!
!  Initialization.
!
  rlist2(1) = result
  errmax = abserr
  maxerr = 1
  area = result
  errsum = abserr
  abserr = huge ( abserr )
  nrmax = 1
  nres = 0
  numrl2 = 2
  ktmin = 0
  extrap = .false.
  noext = .false.
  iroff1 = 0
  iroff2 = 0
  iroff3 = 0

  if ( dres >= (1.0E+00 - 5.0E+01* epsilon ( defabs ) ) * defabs ) then
    ksgn = 1
  else
    ksgn = -1
  end if

  do last = 2, limit
!
!  Bisect the subinterval with the nrmax-th largest error estimate.
!
    a1 = alist(maxerr)
    b1 = 5.0E-01 * ( alist(maxerr) + blist(maxerr) )
    a2 = b1
    b2 = blist(maxerr)
    erlast = errmax
    call qk21 ( a1, b1, area1, error1, resabs, defab1, &
      x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,x)
    call qk21 ( a2, b2, area2, error2, resabs, defab2, &
      x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,x)
!
!  Improve previous approximations to integral and error
!  and test for accuracy.
!
    area12 = area1+area2
    erro12 = error1+error2
    errsum = errsum+erro12-errmax
    area = area+area12-rlist(maxerr)

    if ( defab1 == error1 .or. defab2 == error2 ) go to 15

    if ( abs ( rlist(maxerr) - area12) > 1.0E-05 * abs(area12) &
      .or. erro12 < 9.9E-01 * errmax ) go to 10

    if ( extrap ) then
      iroff2 = iroff2+1
    else
      iroff1 = iroff1+1
    end if

10  continue

    if ( last > 10 .and. erro12 > errmax ) then
      iroff3 = iroff3+1
    end if

15  continue

    rlist(maxerr) = area1
    rlist(last) = area2
    errbnd = max ( epsabs, epsrel*abs(area) )
!
!  Test for roundoff error and eventually set error flag.
!
    if ( iroff1+iroff2 >= 10 .or. iroff3 >= 20 ) then
      ier = 2
    end if

    if ( iroff2 >= 5 ) then
      ierro = 3
    end if
!
!  Set error flag in the case that the number of subintervals
!  equals limit.
!
    if ( last == limit ) then
      ier = 1
    end if
!
!  Set error flag in the case of bad integrand behavior
!  at a point of the integration range.
!
    if ( max ( abs(a1),abs(b2)) <= (1.0E+00+1.0E+03* epsilon ( a1 ) )* &
      (abs(a2)+1.0E+03* tiny ( a2 ) ) ) then
      ier = 4
    end if
!
!  Append the newly-created intervals to the list.
!
    if ( error2 <= error1 ) then
      alist(last) = a2
      blist(maxerr) = b1
      blist(last) = b2
      elist(maxerr) = error1
      elist(last) = error2
    else
      alist(maxerr) = a2
      alist(last) = a1
      blist(last) = b1
      rlist(maxerr) = area2
      rlist(last) = area1
      elist(maxerr) = error2
      elist(last) = error1
    end if
!
!  Call QSORT to maintain the descending ordering
!  in the list of error estimates and select the subinterval
!  with nrmax-th largest error estimate (to be bisected next).
!
    call qsort ( limit, last, maxerr, errmax, elist, iord, nrmax )

    if ( errsum <= errbnd ) go to 115

    if ( ier /= 0 ) then
      exit
    end if

    if ( last == 2 ) go to 80
    if ( noext ) go to 90

    erlarg = erlarg-erlast

    if ( abs(b1-a1) > small ) then
      erlarg = erlarg+erro12
    end if
!
!  Test whether the interval to be bisected next is the
!  smallest interval.
!
    if ( .not. extrap ) then
      if ( abs(blist(maxerr)-alist(maxerr)) > small ) go to 90
      extrap = .true.
      nrmax = 2
    end if

!40  continue
!
!  The smallest interval has the largest error.
!  Before bisecting decrease the sum of the errors over the
!  larger intervals (erlarg) and perform extrapolation.
!
    if ( ierro /= 3 .and. erlarg > ertest ) then

      id = nrmax
      jupbnd = last

      if ( last > (2+limit/2) ) then
        jupbnd = limit+3-last
      end if

      do k = id, jupbnd
        maxerr = iord(nrmax)
        errmax = elist(maxerr)
        if ( abs(blist(maxerr)-alist(maxerr)) > small ) then
          go to 90
        end if
        nrmax = nrmax+1
      end do

    end if
!
!  Perform extrapolation.
!
!60  continue

    numrl2 = numrl2+1
    rlist2(numrl2) = area
    call qextr ( numrl2, rlist2, reseps, abseps, res3la, nres )
    ktmin = ktmin+1

    if ( ktmin > 5 .and. abserr < 1.0E-03 * errsum ) then
      ier = 5
    end if

    if ( abseps < abserr ) then

      ktmin = 0
      abserr = abseps
      result = reseps
      correc = erlarg
      ertest = max ( epsabs,epsrel*abs(reseps))

      if ( abserr <= ertest ) then
        exit
      end if

    end if
!
!  Prepare bisection of the smallest interval.
!
    if ( numrl2 == 1 ) then
      noext = .true.
    end if

    if ( ier == 5 ) then
      exit
    end if

    maxerr = iord(1)
    errmax = elist(maxerr)
    nrmax = 1
    extrap = .false.
    small = small * 5.0E-01
    erlarg = errsum
    go to 90

80  continue

    small = abs ( b - a ) * 3.75E-01
    erlarg = errsum
    ertest = errbnd
    rlist2(2) = area

90  continue

  end do
!
!  Set final result and error estimate.
!
  if ( abserr == huge ( abserr ) ) then
    go to 115
  end if

  if ( ier + ierro == 0 ) then
    go to 110
  end if

  if ( ierro == 3 ) then
    abserr = abserr + correc
  end if

  if ( ier == 0 ) then
    ier = 3
  end if

  if ( result /= 0.0E+00 .and. area /= 0.0E+00 ) then
    go to 105
  end if

  if ( abserr > errsum ) go to 115
  if ( area == 0.0E+00 ) go to 130
  go to 110

105 continue

  if ( abserr/abs(result) > errsum/abs(area) ) go to 115
!
!  Test on divergence.
!
110 continue

  if ( ksgn == (-1) .and. max ( abs(result),abs(area)) <=  &
   defabs*1.0E-02 ) go to 130

  if ( 1.0E-02 > (result/area) .or. (result/area) > 1.0E+02 &
   .or. errsum > abs(area) ) then
    ier = 6
  end if

  go to 130
!
!  Compute global integral sum.
!
115 continue

  result = sum ( rlist(1:last) )

  abserr = errsum

130 continue

  if ( 2 < ier ) then
    ier = ier - 1
  end if

140 continue

  neval = 42*last-21

  return
end

subroutine qk21 ( a, b, result, &
  x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,x)

!*****************************************************************************80
!
!! QK21 carries out a 21 point Gauss-Kronrod quadrature rule.
!
!  Discussion:
!
!    This routine approximates
!      I = integral ( A <= X <= B ) F(X) dx
!    with an error estimate, and
!      J = integral ( A <= X <= B ) | F(X) | dx
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, external real ( kind = 8 ) F, the name of the function routine, of the form
!      function f ( x )
!      real ( kind = 8 ) f
!      real ( kind = 8 ) x
!    which evaluates the integrand function.
!
!    Input, real ( kind = 8 ) A, B, the limits of integration.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!    RESULT is computed by applying the 21-point Kronrod rule (resk)
!    obtained by optimal addition of abscissae to the 10-point Gauss
!    rule (resg).
!
!    Output, real ( kind = 8 ) ABSERR, an estimate of | I - RESULT |.
!
!    Output, real ( kind = 8 ) RESABS, approximation to the integral of the absolute
!    value of F.
!
!    Output, real ( kind = 8 ) RESASC, approximation to the integral | F-I/(B-A) |
!    over [A,B].
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) absc
  real ( kind = 8 ) abserr
  real ( kind = 8 ) b
  real ( kind = 8 ) centr
  real ( kind = 8 ) dhlgth
  real ( kind = 8 ) fc
  real ( kind = 8 ) fsum
  real ( kind = 8 ) fval1
  real ( kind = 8 ) fval2
  real ( kind = 8 ) fv1(10)
  real ( kind = 8 ) fv2(10)
  real ( kind = 8 ) hlgth
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jtw
  integer ( kind = 8 ) jtwm1
  real ( kind = 8 ) resabs
  real ( kind = 8 ) resasc
  real ( kind = 8 ) resg
  real ( kind = 8 ) resk
  real ( kind = 8 ) reskh
  real ( kind = 8 ) result
  real ( kind = 8 ) wg(5)
  real ( kind = 8 ) wgk(11)
  real ( kind = 8 ) xgk(11)

  real ( kind = 8 ) x0
  real ( kind = 8 ) y0
  real ( kind = 8 ) dia
  real ( kind = 8 ) loc1
  real ( kind = 8 ) loc2
  real ( kind = 8 ) loc3
  real ( kind = 8 ) spr1
  real ( kind = 8 ) spr2
  real ( kind = 8 ) skw1
  real ( kind = 8 ) skw2
  real ( kind = 8 ) scl1
  real ( kind = 8 ) scl2
  real ( kind = 8 ) scl3
  logical xyrun
  real ( kind = 8 ) x
!
!           the abscissae and weights are given for the interval (-1,1).
!           because of symmetry only the positive abscissae and their
!           corresponding weights are given.
!
!           xgk    - abscissae of the 21-point Kronrod rule
!                    xgk(2), xgk(4), ...  abscissae of the 10-point
!                    Gauss rule
!                    xgk(1), xgk(3), ...  abscissae which are optimally
!                    added to the 10-point Gauss rule
!
!           wgk    - weights of the 21-point Kronrod rule
!
!           wg     - weights of the 10-point Gauss rule
!
  data xgk(1),xgk(2),xgk(3),xgk(4),xgk(5),xgk(6),xgk(7),xgk(8), &
    xgk(9),xgk(10),xgk(11)/ &
       9.956571630258081E-01,     9.739065285171717E-01, &
       9.301574913557082E-01,     8.650633666889845E-01, &
       7.808177265864169E-01,     6.794095682990244E-01, &
       5.627571346686047E-01,     4.333953941292472E-01, &
       2.943928627014602E-01,     1.488743389816312E-01, &
       0.000000000000000E+00/
!
  data wgk(1),wgk(2),wgk(3),wgk(4),wgk(5),wgk(6),wgk(7),wgk(8), &
    wgk(9),wgk(10),wgk(11)/ &
       1.169463886737187E-02,     3.255816230796473E-02, &
       5.475589657435200E-02,     7.503967481091995E-02, &
       9.312545458369761E-02,     1.093871588022976E-01, &
       1.234919762620659E-01,     1.347092173114733E-01, &
       1.427759385770601E-01,     1.477391049013385E-01, &
       1.494455540029169E-01/
!
  data wg(1),wg(2),wg(3),wg(4),wg(5)/ &
       6.667134430868814E-02,     1.494513491505806E-01, &
       2.190863625159820E-01,     2.692667193099964E-01, &
       2.955242247147529E-01/
!
!
!           list of major variables
!
!           centr  - mid point of the interval
!           hlgth  - half-length of the interval
!           absc   - abscissa
!           fval*  - function value
!           resg   - result of the 10-point Gauss formula
!           resk   - result of the 21-point Kronrod formula
!           reskh  - approximation to the mean value of f over (a,b),
!                    i.e. to i/(b-a)
!
  centr = 5.0E-01*(a+b)
  hlgth = 5.0E-01*(b-a)
  dhlgth = abs(hlgth)
!
!  Compute the 21-point Kronrod approximation to the
!  integral, and estimate the absolute error.
!
  xyrun = xyrun
  ! print *, 'xyrun', xyrun

  if (xyrun .eqv. .false.) then
    ! print *, 'inte'
    resg = 0.0E+00
    !fc = f(centr)
    call integrand(centr,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
    scl1,scl2,scl3,fc)

    resk = wgk(11)*fc
    resabs = abs(resk)

    do j = 1, 5
      jtw = 2*j
      absc = hlgth*xgk(jtw)
      ! fval1 = f(centr-absc)
      ! fval2 = f(centr+absc)
      call integrand(centr-absc,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval1)
      call integrand(centr+absc,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval2)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 5
      jtwm1 = 2*j-1
      absc = hlgth*xgk(jtwm1)
      ! fval1 = f(centr-absc)
      ! fval2 = f(centr+absc)
      call integrand(centr-absc,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval1)
      call integrand(centr+absc,x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval2)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do
  else if (xyrun .eqv. .true.) then
    ! print *, 'infunc'
    resg = 0.0E+00
    !fc = f(centr)
    call infunc(centr,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
    scl1,scl2,scl3,fc)

    resk = wgk(11)*fc
    resabs = abs(resk)

    do j = 1, 5
      jtw = 2*j
      absc = hlgth*xgk(jtw)
      ! fval1 = f(centr-absc)
      ! fval2 = f(centr+absc)
      call infunc(centr-absc,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval1)
      call infunc(centr+absc,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval2)
      fv1(jtw) = fval1
      fv2(jtw) = fval2
      fsum = fval1+fval2
      resg = resg+wg(j)*fsum
      resk = resk+wgk(jtw)*fsum
      resabs = resabs+wgk(jtw)*(abs(fval1)+abs(fval2))
    end do

    do j = 1, 5
      jtwm1 = 2*j-1
      absc = hlgth*xgk(jtwm1)
      ! fval1 = f(centr-absc)
      ! fval2 = f(centr+absc)
      call infunc(centr-absc,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval1)
      call infunc(centr+absc,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,&
      scl1,scl2,scl3,fval2)
      fv1(jtwm1) = fval1
      fv2(jtwm1) = fval2
      fsum = fval1+fval2
      resk = resk+wgk(jtwm1)*fsum
      resabs = resabs+wgk(jtwm1)*(abs(fval1)+abs(fval2))
    end do
  end if

  reskh = resk*5.0E-01
  resasc = wgk(11)*abs(fc-reskh)

  do j = 1, 10
    resasc = resasc+wgk(j)*(abs(fv1(j)-reskh)+abs(fv2(j)-reskh))
  end do

  result = resk*hlgth
  resabs = resabs*dhlgth
  resasc = resasc*dhlgth
  abserr = abs((resk-resg)*hlgth)

  if ( resasc /= 0.0E+00.and.abserr /= 0.0E+00) then
    abserr = resasc*min ( 1.0E+00,(2.0E+02*abserr/resasc)**1.5E+00)
  end if

  if ( resabs > tiny ( resabs ) /(5.0E+01* epsilon ( resabs ) )) then
    abserr = max (( epsilon ( resabs ) *5.0E+01)*resabs,abserr)
  end if

  return
end

subroutine qsort ( limit, last, maxerr, ermax, elist, iord, nrmax )

!*****************************************************************************80
!
!! QSORT maintains the order of a list of local error estimates.
!
!  Discussion:
!
!    This routine maintains the descending ordering in the list of the
!    local error estimates resulting from the interval subdivision process.
!    At each call two error estimates are inserted using the sequential
!    search top-down for the largest error estimate and bottom-up for the
!    smallest error estimate.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) LIMIT, the maximum number of error estimates the list can
!    contain.
!
!    Input, integer ( kind = 8 ) LAST, the current number of error estimates.
!
!    Input/output, integer ( kind = 8 ) MAXERR, the index in the list of the NRMAX-th
!    largest error.
!
!    Output, real ( kind = 8 ) ERMAX, the NRMAX-th largest error = ELIST(MAXERR).
!
!    Input, real ( kind = 8 ) ELIST(LIMIT), contains the error estimates.
!
!    Input/output, integer ( kind = 8 ) IORD(LAST).  The first K elements contain
!    pointers to the error estimates such that ELIST(IORD(1)) through
!    ELIST(IORD(K)) form a decreasing sequence, with
!      K = LAST
!    if
!      LAST <= (LIMIT/2+2),
!    and otherwise
!      K = LIMIT+1-LAST.
!
!    Input/output, integer ( kind = 8 ) NRMAX.
!
  implicit none

  integer ( kind = 8 ) last

  real ( kind = 8 ) elist(last)
  real ( kind = 8 ) ermax
  real ( kind = 8 ) errmax
  real ( kind = 8 ) errmin
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ibeg
  integer ( kind = 8 ) iord(last)
  integer ( kind = 8 ) isucc
  integer ( kind = 8 ) j
  integer ( kind = 8 ) jbnd
  integer ( kind = 8 ) jupbn
  integer ( kind = 8 ) k
  integer ( kind = 8 ) limit
  integer ( kind = 8 ) maxerr
  integer ( kind = 8 ) nrmax


!
!  Check whether the list contains more than two error estimates.
!
  if ( last <= 2 ) then
    iord(1) = 1
    iord(2) = 2
    go to 90
  end if
!
!  This part of the routine is only executed if, due to a
!  difficult integrand, subdivision increased the error
!  estimate. in the normal case the insert procedure should
!  start after the nrmax-th largest error estimate.
!
  errmax = elist(maxerr)

  do i = 1, nrmax-1

    isucc = iord(nrmax-1)

    if ( errmax <= elist(isucc) ) then
      exit
    end if

    iord(nrmax) = isucc
    nrmax = nrmax-1

  end do
!
!  Compute the number of elements in the list to be maintained
!  in descending order.  This number depends on the number of
!  subdivisions still allowed.
!
  jupbn = last

  if ( (limit/2+2) < last ) then
    jupbn = limit+3-last
  end if

  errmin = elist(last)
!
!  Insert errmax by traversing the list top-down, starting
!  comparison from the element elist(iord(nrmax+1)).
!
  jbnd = jupbn-1
  ibeg = nrmax+1

  do i = ibeg, jbnd
    isucc = iord(i)
    if ( elist(isucc) <= errmax ) then
      go to 60
    end if
    iord(i-1) = isucc
  end do

  iord(jbnd) = maxerr
  iord(jupbn) = last
  go to 90
!
!  Insert errmin by traversing the list bottom-up.
!
60 continue

  iord(i-1) = maxerr
  k = jbnd

  do j = i, jbnd
    isucc = iord(k)
    if ( errmin < elist(isucc) ) then
      go to 80
    end if
    iord(k+1) = isucc
    k = k-1
  end do

  iord(i) = last
  go to 90

80 continue

  iord(k+1) = last
!
!  Set maxerr and ermax.
!
90 continue

  maxerr = iord(nrmax)
  ermax = elist(maxerr)

  return
end

subroutine qextr ( n, epstab, result, abserr, res3la, nres )

!*****************************************************************************80
!
!! QEXTR carries out the Epsilon extrapolation algorithm.
!
!  Discussion:
!
!    The routine determines the limit of a given sequence of approximations,
!    by means of the epsilon algorithm of P. Wynn.  An estimate of the
!    absolute error is also given.  The condensed epsilon table is computed.
!    Only those elements needed for the computation of the next diagonal
!    are preserved.
!
!  Author:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner
!
!  Reference:
!
!    Robert Piessens, Elise de Doncker-Kapenger,
!    Christian Ueberhuber, David Kahaner,
!    QUADPACK, a Subroutine Package for Automatic Integration,
!    Springer Verlag, 1983
!
!  Parameters:
!
!    Input, integer ( kind = 8 ) N, indicates the entry of EPSTAB which contains
!    the new element in the first column of the epsilon table.
!
!    Input/output, real ( kind = 8 ) EPSTAB(52), the two lower diagonals of the triangular
!    epsilon table.  The elements are numbered starting at the right-hand
!    corner of the triangle.
!
!    Output, real ( kind = 8 ) RESULT, the estimated value of the integral.
!
!    Output, real ( kind = 8 ) ABSERR, estimate of the absolute error computed from
!    RESULT and the 3 previous results.
!
!    ?, real ( kind = 8 ) RES3LA(3), the last 3 results.
!
!    Input/output, integer ( kind = 8 ) NRES, the number of calls to the routine.  This
!    should be zero on the first call, and is automatically updated
!    before return.
!
!  Local Parameters:
!
!           e0     - the 4 elements on which the
!           e1       computation of a new element in
!           e2       the epsilon table is based
!           e3                 e0
!                        e3    e1    new
!                              e2
!           newelm - number of elements to be computed in the new
!                    diagonal
!           error  - error = abs(e1-e0)+abs(e2-e1)+abs(new-e2)
!           result - the element in the new diagonal with least value
!                    of error
!           limexp is the maximum number of elements the epsilon table
!           can contain. if this number is reached, the upper diagonal
!           of the epsilon table is deleted.
!
  implicit none

  real ( kind = 8 ) abserr
  real ( kind = 8 ) delta1
  real ( kind = 8 ) delta2
  real ( kind = 8 ) delta3
  real ( kind = 8 ) epsinf
  real ( kind = 8 ) epstab(52)
  real ( kind = 8 ) error
  real ( kind = 8 ) err1
  real ( kind = 8 ) err2
  real ( kind = 8 ) err3
  real ( kind = 8 ) e0
  real ( kind = 8 ) e1
  real ( kind = 8 ) e1abs
  real ( kind = 8 ) e2
  real ( kind = 8 ) e3
  integer ( kind = 8 ) i
  integer ( kind = 8 ) ib
  integer ( kind = 8 ) ib2
  integer ( kind = 8 ) ie
  integer ( kind = 8 ) indx
  integer ( kind = 8 ) k1
  integer ( kind = 8 ) k2
  integer ( kind = 8 ) k3
  integer ( kind = 8 ) limexp
  integer ( kind = 8 ) n
  integer ( kind = 8 ) newelm
  integer ( kind = 8 ) nres
  integer ( kind = 8 ) num
  real ( kind = 8 ) res
  real ( kind = 8 ) result
  real ( kind = 8 ) res3la(3)
  real ( kind = 8 ) ss
  real ( kind = 8 ) tol1
  real ( kind = 8 ) tol2
  real ( kind = 8 ) tol3



  nres = nres+1
  abserr = huge ( abserr )
  result = epstab(n)

  if ( n < 3 ) then
    abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))
    return
  end if

  limexp = 50
  epstab(n+2) = epstab(n)
  newelm = (n-1)/2
  epstab(n) = huge ( epstab(n) )
  num = n
  k1 = n

  do i = 1, newelm

    k2 = k1-1
    k3 = k1-2
    res = epstab(k1+2)
    e0 = epstab(k3)
    e1 = epstab(k2)
    e2 = res
    e1abs = abs(e1)
    delta2 = e2-e1
    err2 = abs(delta2)
    tol2 = max ( abs(e2),e1abs)* epsilon ( e2 )
    delta3 = e1-e0
    err3 = abs(delta3)
    tol3 = max ( e1abs,abs(e0))* epsilon ( e0 )
!
!  If e0, e1 and e2 are equal to within machine accuracy, convergence
!  is assumed.
!
    if ( err2 <= tol2 .and. err3 <= tol3 ) then
      result = res
      abserr = err2+err3
      abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))
      return
    end if

    e3 = epstab(k1)
    epstab(k1) = e1
    delta1 = e1-e3
    err1 = abs(delta1)
    tol1 = max ( e1abs,abs(e3))* epsilon ( e3 )
!
!  If two elements are very close to each other, omit a part
!  of the table by adjusting the value of N.
!
    if ( err1 <= tol1 .or. err2 <= tol2 .or. err3 <= tol3 ) go to 20

    ss = 1.0E+00/delta1+1.0E+00/delta2-1.0E+00/delta3
    epsinf = abs ( ss*e1 )
!
!  Test to detect irregular behavior in the table, and
!  eventually omit a part of the table adjusting the value of N.
!
    if ( epsinf > 1.0E-04 ) go to 30

20  continue

    n = i+i-1
    exit
!
!  Compute a new element and eventually adjust the value of RESULT.
!
30  continue

    res = e1+1.0E+00/ss
    epstab(k1) = res
    k1 = k1-2
    error = err2+abs(res-e2)+err3

    if ( error <= abserr ) then
      abserr = error
      result = res
    end if

  end do
!
!  Shift the table.
!
  if ( n == limexp ) then
    n = 2*(limexp/2)-1
  end if

  if ( (num/2)*2 == num ) then
    ib = 2
  else
    ib = 1
  end if

  ie = newelm+1

  do i = 1, ie
    ib2 = ib+2
    epstab(ib) = epstab(ib2)
    ib = ib2
  end do

  if ( num /= n ) then

    indx = num-n+1

    do i = 1, n
      epstab(i)= epstab(indx)
      indx = indx+1
    end do

  end if

  if ( nres < 4 ) then
    res3la(nres) = result
    abserr = huge ( abserr )
  else
    abserr = abs(result-res3la(3))+abs(result-res3la(2)) &
      +abs(result-res3la(1))
    res3la(1) = res3la(2)
    res3la(2) = res3la(3)
    res3la(3) = result
  end if

  abserr = max ( abserr,0.5E+00* epsilon ( result ) *abs(result))

  return
end

! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************
! *****************************************************************************

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

! Transforming double integral into single integral
subroutine infunc(x,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,yint)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: x,x0,y0,dia
  real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
  ! out
  real(dp), intent(out) :: yint
  ! local
  real(dp) :: a,b,epsabs,epsrel,xhold
  logical :: xyrun

  a = -4.0_dp*dia
  b = 4.0_dp*dia

  epsabs = 1.49e-8_dp
  epsrel = 1.49e-8_dp

  xyrun = .true.
  xhold = 1.0_dp
  ! print *, 'infunc func'

  call qags( a, b, epsabs, epsrel, yint, &
    x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,xhold)

  xyrun = .false.

  ! call integrand(x,xhold,x0,y0,dia,loc1,loc2,loc3,spr1,spr2,&
  !   skw1,skw2,scl1,scl2,scl3,yint)


end subroutine infunc

! Main run routine
subroutine runvort(x0,y0,dia,xbound,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,vel)
  implicit none
  integer, parameter :: dp = kind(0.d0)
  ! in
  real(dp), intent(in) :: x0,y0,dia,xbound
  real(dp), intent(in) :: loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3
  ! out
  real(dp), intent(out) :: vel
  ! local
  real(dp) :: a,b,epsabs,epsrel,x!,yint
  logical :: xyrun

  a = 0.0_dp
  b = xbound

  epsabs = 1.49e-8_dp
  epsrel = 1.49e-8_dp

  xyrun = .true.
  x = 5.0_dp

  call qags( a, b, epsabs, epsrel, vel, &
    x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,x)

  ! print *, yint
  !
  ! call qags( a, b, epsabs, epsrel, vel, &
  !   x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,xyrun,yint)

  ! print *, vel


end subroutine runvort
