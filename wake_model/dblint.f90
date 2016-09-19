! f2py -c  --opt=-O2 -m _dblint dblint.f90

MODULE dblint
  !
  !  Routines for double integration
  !  Peter Simon
  !  Converted to Fortran 90 on 11/2/98
  !  Modified: 2/20/2002 to increase nomax to 300000
  !
CONTAINS

  SUBROUTINE cdblin(fun,ax,bx,ay,by,aerr, RESULT,errest,nofun,flag)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    !  numerical integration of a complex function of two real variables   c
    !                                                                      c
    !  this routine provides an estimate of the integral of fun(x) over    c
    !  the region ax .le. x .le. bx,   ay(x) .le. y .le. by(x)             c
    !                                                                      c
    !  an automatic adaptive routine based on the 8-panel newton-cotes     c
    !  rule.  (double precision version)                                   c
    !                                                                      c
    !  input ...                                                           c
    !                                                                      c
    !    fun      the name of the complex*16 integrand function fun(x,y).  c
    !             fun must be declared external by the calling program.    c
    !                                                                      c
    !    ax, bx   the lower limit and upper limits of integration for the  c
    !             outer (over x) iterated integral.                        c
    !                                                                      c
    !    ay, by   double precision function subroutines provided by the    c
    !             user which take a single double precision argument (x)   c
    !             and return the lower and upper limits, respectively, of  c
    !             the inner iterated integral (over y).                    c
    !                                                                      c
    !    aerr   an absolute error tolerance (double precision)             c
    !             aerr is normally set to a nonnegative value which        c
    !             is the desired absolute accuracy of the estimate.        c
    !             if aerr is set to zero, the routine will return          c
    !             with an error indication, except for very simple         c
    !             integrands.  if aerr is set to a negative value          c
    !             the routine will estimate the value of the integral      c
    !             using the minimum possible number of integrand           c
    !             evaluations, 81.                                         c
    !                                                                      c
    !  output ...                                                          c
    !                                                                      c
    !    result   an approximation to the integral hopefully satisfying    c
    !             the least stringent of the two error tolerances.         c
    !             (complex*16)                                             c
    !                                                                      c
    !    errest    an estimate of the magnitude of the actual error.       c
    !             (double precision)                                       c
    !                                                                      c
    !    nofun    the number of function values used in calculation of     c
    !             result.  (integer)                                       c
    !                                                                      c
    !    flag     a reliability indicator.  if flag is zero, then result   c
    !             probably satisfies the error tolerance.  if flag is      c
    !             greater than zero, then its integer part is a total      c
    !             count (for both outer and inner integration routines)    c
    !             of the number of intervals which have not converged.     c
    !             if flag has a nonzero fractional part, flag = xxx.yyy,   c
    !             then 0.yyy is the fraction of the outer interval (ax,bx) c
    !             which was left to do when the limit on nofun for the     c
    !             outer integral was approached.  if, however, flag has    c
    !             a nonzero integer part and zero fractional part, then    c
    !             the limit on nofun was reached one or more times during  c
    !             calculation(s) of the inner integral.  in this case,     c
    !             no information on the location of the trouble spot is    c
    !             returned to the user.  in any case, errest is still      c
    !             a useful estimate of the error actually attained.        c
    !             since the size of the interval where convergence failed  c
    !             is very small, the result of this routine is often       c
    !             acceptable when flag is not too large.                   c
    !                                                                      c
    !  description of algorithm:                                           c
    !              the double integral is evaluated as an iterated         c
    !              integral.  both the outer and inner integrals are       c
    !              evaluated using a modified version of quanc8            c
    !              which was taken from reference 1.                       c
    !              two major modifications were applied to quanc8.         c
    !              the first modification involved changing the routine    c
    !              to handle complex*16 functions and to use double        c
    !              precision throughout.  the second modification was to   c
    !              use a four panel newton-cotes rule on the first         c
    !              interval (a,b) and its two halves, in order to get      c
    !              an answer quickly, with few function evaluations, for   c
    !              particularly smooth integrands.  if the accuracy test   c
    !              for the four panel calculation fails, the previously    c
    !              implemented eight panel rule is used for the remainder  c
    !              of the calculations.                                    c
    !              the absolute error tolerance is divvied up between      c
    !              outer and inner integrals using the formula suggested   c
    !              in reference 2.                                         c
    !                                                                      c
    !                                                                      c
    !  references: 1.  forsythe, malcolm, and moler, 'computer             c
    !                  methods for mathematical computations,              c
    !                  prentiss-hall, 1977.                                c
    !              2.  fritsch, kahaner, and lyness, 'double integration   c
    !                  using one dimensional adaptive quadrature routines: c
    !                  a software interface problem,' acm trans. math.     c
    !                  software, vol 7, no. 1, march 1981, pp 46-75.       c
    !                                                                      c
    !  level of control                                                    c
    !                                                                      c
    !    this routine directly calls routines                              c
    !      cquany                                                          c
    !                                                                      c
    !                                                                      c
    !  this routine was initially entered on the resd ibm 4381 by          c
    !  p. simon, dept. 9282, x3726, on 13 oct. 1987.                       c
    !                                                                      c
    !  language:  fortran 77, as implemented on the resd ibm 4381          c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    IMPLICIT DOUBLE PRECISION (a-h, o-z)

    COMPLEX*16 fun, RESULT, area, f0, cor11, cor7
    COMPLEX*16 qprev, qnow, qdiff, qleft
    COMPLEX*16 qright(31), f(16), fsave(8,30), czero

    REAL flag, flagi
    DOUBLE PRECISION  abserr, errest, aerr, errin, errini
    DOUBLE PRECISION ax, bx, ay, by
    DOUBLE PRECISION w0, w1, w2, w3, w4
    DOUBLE PRECISION x0, stone, step, temp
    DOUBLE PRECISION esterr, tolerr
    DOUBLE PRECISION x(16)
    DOUBLE PRECISION xsave(8,30), epsini, wf0, wf1, wf2, xxx

    INTEGER nofun, nfuni
    INTEGER levmin, levmax, levout, nomax, nofin, lev, nim, i, j

    EXTERNAL fun, ay, by

    COMMON /xvalue/ xxx

    PARAMETER ( czero = (0.d0,0.d0) )
    !
    !  ***  stage 1 ***  general initialization
    !  set constants
    !
    PARAMETER ( levmin = 0 )
    PARAMETER ( levout = 6 )
    PARAMETER ( nomax = 300000)
    !      parameter ( z0 = 1.5873015873016d-2 )
    !      parameter ( z1 = 9.7751710654936d-4 )
    PARAMETER ( z0 = 1.d0 / 63.d0 )
    PARAMETER ( z1 = 1.d0 / 1023.d0 )
    !
    !  trouble when nofun reaches nofin
    !
    !      parameter ( w0 =   3956.d0 / 14175.d0 )
    !      parameter ( w1 =  23552.d0 / 14175.d0 )
    !      parameter ( w2 =  -3712.d0 / 14175.d0 )
    !      parameter ( w3 =  41984.d0 / 14175.d0 )
    !      parameter ( w4 = -18160.d0 / 14175.d0 )
    PARAMETER ( w0 =   0.27908289241623d0 )
    PARAMETER ( w1 = 1.6615167548501d0 )
    PARAMETER ( w2 =  -0.26186948853616d0 )
    PARAMETER ( w3 =  2.9618342151675d0 )
    PARAMETER ( w4 = -1.2811287477954d0 )

    !      parameter ( wf0 =  7.d0 / 90.d0 )
    !      parameter ( wf1 = 32.d0 / 90.d0 )
    !      parameter ( wf2 = 12.d0 / 90.d0 )
    PARAMETER ( wf0 =  7.7777777777778d-2 )
    PARAMETER ( wf1 = 0.35555555555556d0 )
    PARAMETER ( wf2 = 0.13333333333333d0 )
    !
    !  first executable statement
    !
    levmax = 30
    nofin = nomax - 8*(levmax-levout+2**(levout+1))
    !
    !  initialize running sums to zero
    !
    flag = 0.0
    RESULT = czero
    cor11 = czero
    errest = 0.d0
    nofun = 0
    errin = 0.0
    IF ( ax .EQ. bx ) RETURN
    !
    !  choose error tolerances for outer and inner integrals
    !
    abserr = 0.5d0 * aerr
    epsini = aerr / (  2.90758d0 * (bx - ax) )
    !
    !  ***  stage 2 ***  initialization for first interval
    !
    lev = 0
    nim = 1
    x0 = ax
    x(16) = bx

    CALL cquany(fun,x0,ay,by,epsini,f0,errini,nfuni,flagi)


    nofun = nofun + nfuni
    errin = MAX(errin, errini)
    flag = flag + INT(flagi)
    stone = 0.0625d0 * (bx-ax)
    x(8) = 0.5d0 * (x0 + x(16))
    x(4) = 0.5d0 * (x0 + x(8))
    x(12) = 0.5d0 * (x(8) + x(16))
    x(2) = 0.5d0 * (x0 + x(4))
    x(6) = 0.5d0 * (x(4) + x(8))
    x(10) = 0.5d0 * (x(8) + x(12))
    x(14) = 0.5d0 * (x(12) + x(16))
    DO j = 2, 16, 2
       CALL cquany(fun,x(j),ay,by,epsini,f(j),errini,nfuni,flagi)
       nofun = nofun + nfuni
       errin = MAX(errin, errini)
       flag = flag + INT(flagi)
    END DO
    !
    !  estimate total integral using 4 panel rule
    !
    qprev = ( wf0 * (f0 + f(16)) &
         &          +  wf1 * (f(4) + f(12)) + wf2 * f(8) ) * (bx-ax)
    qnow = (  wf0 * ( f0 + 2.d0 * f(8) + f(16) ) &
         & + wf1 * ( f(2) + f(6) + f(10) + f(14) ) &
         & + wf2 * ( f(4) + f(12) ) ) * (0.5d0*(bx-ax))

    qdiff = qnow - qprev
    cor7 = z0 * qdiff
    esterr = ABS( cor7 )
    tolerr =  abserr
    IF ( esterr .LE. tolerr .OR. abserr .LT. 0.d0 ) THEN
       RESULT = qnow + cor7
       errest = esterr  + (bx - ax) * errin
       RETURN
    END IF
    !
    !  since 4 panel rule failed, calculate initial value using 8 panel rule
    !
    qprev = ( &
         & w0 * (f0   + f(16)) &
         & +   w1 * (f(2) + f(14)) &
         & +   w2 * (f(4) + f(12)) &
         & +   w3 * (f(6) + f(10)) &
         & +   w4 * f(8) &
         & ) * (0.125d0*(bx-ax))
    area = qprev
    !
    !  ***  stage 3 ***  central calculation
    !  requires qprev, x0, x2, x4, ..., x16, f0, f2, f4, ..., f16.
    !  calculates x1, x3, ..., x15, f1, f3, ..., f15, qleft, qnow, qdiff,
    !  area.
    !
30  CONTINUE
    x(1) = 0.5d0 * (x0 + x(2))
    !     f(1) = fun(x(1))
    CALL cquany(fun,x(1),ay,by,epsini,f(1),errini,nfuni,flagi)
    nofun = nofun + nfuni
    errin = MAX(errin, errini)
    flag = flag + INT(flagi)
    DO j = 3, 15, 2
       x(j) = 0.5d0 * (x(j-1) + x(j+1))
       !       f(j) = fun(x(j))
       CALL cquany(fun,x(j),ay,by,epsini,f(j),errini,nfuni,flagi)
       nofun = nofun + nfuni
       errin = MAX(errin, errini)
       flag = flag + INT(flagi)
    END DO

    step = 0.0625d0 * (x(16) - x0)
    qleft = ( &
         & w0 * (f0   + f(8)) &
         & +   w1 * (f(1) + f(7)) &
         & +   w2 * (f(2) + f(6)) &
         & +   w3 * (f(3) + f(5)) &
         & +   w4 * f(4) &
         & ) * step
    qright(lev+1) =  ( &
         & w0 * (f(8)  + f(16)) &
         & +   w1 * (f(9)  + f(15)) &
         & +   w2 * (f(10) + f(14)) &
         & +   w3 * (f(11) + f(13)) &
         & +   w4 *  f(12) &
         &   ) * step
    qnow = qleft + qright(lev+1)
    qdiff = qnow - qprev
    area = area + qdiff
    !
    !  ***  stage 4 ***  interval convergence test
    !
    esterr = z1 * ABS(qdiff)
    tolerr = abserr * (step /stone)

    IF ( lev .LT. levmin ) GOTO 50
    IF ( lev .GE. levmax ) GOTO 62
    IF ( nofun .GT. nofin ) GOTO 60
    IF (esterr .LE. tolerr ) GOTO 70
    !
    !  ***  stage 5 ***  no convergence
    !  locate next interval
    !
50  CONTINUE
    nim = 2 * nim
    lev = lev + 1
    !
    !  store right hand elements for future use
    !
    DO i = 1, 8
       fsave(i,lev) = f(i+8)
       xsave(i,lev) = x(i+8)
    END DO
    !
    !  assemble left hand elements for immediate use
    !
    qprev = qleft
    DO i = 1, 8
       j = -i
       f(2*j+18) = f(j+9)
       x(2*j+18) = x(j+9)
    END DO
    GOTO 30
    !
    !  ***  stage 6 ***  trouble section
    !  number of function values is about to exceed limit
    !
60  CONTINUE
    nofin = 2 * nofin
    levmax = levout
    flag = flag + (bx - x0) / (bx - ax)
    GOTO 70
    !
    !  current level is levmax
    !
62  CONTINUE
    flag = flag + 1.0
    !  ***  stage 7  *** interval converged
    !  add contributions into running sums
    !
70  CONTINUE
    RESULT = RESULT + qnow
    errest = errest + esterr
    cor11 = cor11 + z1 * qdiff
    !
    !  locate next interval
    !
72  CONTINUE
    IF (nim .EQ. 2*(nim/2)) GOTO 75
    nim = nim/2
    lev = lev - 1
    GOTO 72
75  CONTINUE
    nim = nim + 1
    IF ( lev .LE. 0) GOTO 80
    !
    !  assemble elements required for the next interval
    !
    qprev = qright(lev)
    x0 = x(16)
    f0 = f(16)
    DO i = 1, 8
       f(2*i) = fsave(i, lev)
       x(2*i) = xsave(i, lev)
    END DO
    GOTO 30
    !
    !  ***  stage 8  *** finalize and return
    !
80  CONTINUE
    RESULT = RESULT + cor11
    !
    !  make sure errest not less than roundoff level
    !
    errest = errest + (bx - ax) * errin
    IF ( errest .EQ. 0.d0 ) RETURN
82  CONTINUE
    temp = ABS(RESULT) + errest
    IF ( temp .NE. ABS(RESULT) ) RETURN
    errest = 2.d0 * errest
    GOTO 82
  END SUBROUTINE cdblin

  SUBROUTINE cquany(cfun,y,aa,bb,abserr,RESULT,errest,nofun,flag)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    !  this routine provides an estimate of the integral of cfun(y,x) from c
    !  x = a to x = b to a user provided tolerance, where cfun(y,x) is a   c
    !  user supplied complex*16 function of two double-precision           c
    !  arguments, y and x.                                                 c
    !                                                                      c
    !  an automatic adaptive routine based on the 8-panel newton-cotes     c
    !  rule.  (complex*16 version)                                         c
    !                                                                      c
    !  input ...                                                           c
    !                                                                      c
    !    cfun     the name of the complex*16 integrand function cfun(x).   c
    !             cfun must be declared external by the calling program.   c
    !                                                                      c
    !    y        the value of the second argument of cfun(x,y)            c
    !                                                                      c
    !    aa, bb   the names of user supplied double precision function     c
    !             subroutines which return the inner integration limits    c
    !             given the current value of the outer variable of         c
    !             integration. they must be declared external by the       c
    !             calling subroutine.                                      c
    !                                                                      c
    !    abserr   an absolute error tolerance (convergence criterion)      c
    !             should be nonnegative.   (double precision)              c
    !             if abserr is less than zero, the routine will compute    c
    !             an estimate of the integral using the minimum possible   c
    !             number of function evaluations (9), then return.         c
    !                                                                      c
    !  output ...                                                          c
    !                                                                      c
    !    result   an approximation to the integral hopefully satisfying    c
    !             the error tolerance. (complex*16)                        c
    !                                                                      c
    !    errest   an estimate of the magnitude of the actual error.        c
    !             (double precision)                                       c
    !                                                                      c
    !    nofun    the number of function values used in calculation of     c
    !             result.  (integer)                                       c
    !                                                                      c
    !    flag     a reliability indicator.  if flag is zero, then result   c
    !             probably satisfies the error tolerance.  if flag is      c
    !             xxx.yyy, then xxx = the number of intervals which have   c
    !             not converged and 0.yyy = the fraction of the interval   c
    !             left to do when the limit on nofun was approached.       c
    !             (double precision)                                       c
    !                                                                      c
    !  reference:  this routine is a modified version of quanc8            c
    !              which was taken from the textbook 'computer             c
    !              methods for mathematical computations,' by              c
    !              forsythe, malcolm, and moler,  prentiss-hall, 1977.     c
    !              the first modification involved changing the routine    c
    !              to handle complex*16 functions and to use double        c
    !              precision throughout.  the second modification is to    c
    !              use a four panel newton-cotes rule on the first         c
    !              interval (a,b) and its two halves, in order to get      c
    !              an answer quickly, with few function evaluations, for   c
    !              particularly smooth integrands.  if the accuracy test   c
    !              for the four panel calculation fails, the previously    c
    !              implemented eight panel rule is used for the remainder  c
    !              of the calculations.                                    c
    !                                                                      c
    !  language:  fortran 77, as implemented on the resd ibm 4381          c
    !                                                                      c
    !  this routine was initially entered on the resd ibm 4381 by          c
    !  p. simon, dept. 9282, x3726, on 13 oct. 1987.                       c
    !  Modified 4/7/94 to add external statement.                          c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    EXTERNAL   cfun, aa, bb
    COMPLEX*16 cfun, RESULT, area, f0, cor11, cor7, fun
    COMPLEX*16 qprev, qnow, qdiff, qleft
    COMPLEX*16 qright(31), f(16), fsave(8,30), czero

    REAL flag
    DOUBLE PRECISION a, b, abserr, errest, y
    DOUBLE PRECISION w0, w1, w2, w3, w4
    DOUBLE PRECISION x0, stone, step, temp
    DOUBLE PRECISION esterr, tolerr
    DOUBLE PRECISION x(16)
    DOUBLE PRECISION xsave(8,30), t, aa, bb, wf0, wf1, wf2

    INTEGER nofun
    INTEGER levmin, levmax, levout, nomax, nofin, lev, nim, i, j

    PARAMETER ( czero = (0.d0,0.d0) )
    !
    !  ***  stage 1 ***  general initialization
    !  set constants
    !
    PARAMETER ( levmin = 0 )
    PARAMETER ( levout = 6 )
    PARAMETER ( nomax = 300000 )
    !      parameter ( z0 = 1.d0 / 63.d0 )
    !      parameter ( z1 = 1.d0 / 1023.d0 )
    PARAMETER ( z0 = 1.5873015873016d-2 )
    PARAMETER ( z1 = 9.7751710654936d-4 )
    !
    !  trouble when nofun reaches nofin
    !
    !      parameter ( w0 =   3956.d0 / 14175.d0 )
    !      parameter ( w1 =  23552.d0 / 14175.d0 )
    !      parameter ( w2 =  -3712.d0 / 14175.d0 )
    !      parameter ( w3 =  41984.d0 / 14175.d0 )
    !      parameter ( w4 = -18160.d0 / 14175.d0 )
    PARAMETER ( w0 =   0.27908289241623d0 )
    PARAMETER ( w1 = 1.6615167548501d0 )
    PARAMETER ( w2 =  -0.26186948853616d0 )
    PARAMETER ( w3 =  2.9618342151675d0 )
    PARAMETER ( w4 = -1.2811287477954d0 )

    !      parameter ( wf0 =  7.d0 / 90.d0 )
    !      parameter ( wf1 = 32.d0 / 90.d0 )
    !      parameter ( wf2 = 12.d0 / 90.d0 )
    PARAMETER ( wf0 =  7.7777777777778d-2 )
    PARAMETER ( wf1 = 0.35555555555556d0 )
    PARAMETER ( wf2 = 0.13333333333333d0 )

    fun(t) = cfun(y,t)
    !
    !  first executable statement
    !
    levmax = 30
    nofin = nomax - 8*(levmax-levout+2**(levout+1))
    !
    !  initialize running sums to zero
    !
    flag = 0.0
    RESULT = czero
    cor11 = czero
    errest = 0.d0
    nofun = 0
    !
    !  get integration limits
    !
    a = aa(y)
    b = bb(y)
    IF ( a .EQ. b ) RETURN
    !
    !  ***  stage 2 ***  initialization for first interval
    !
    lev = 0
    nim = 1
    x0 = a
    x(16) = b
    f0 = fun(x0)
    stone = 0.0625d0 * (b-a)
    x(8) = 0.5d0 * (x0 + x(16))
    x(4) = 0.5d0 * (x0 + x(8))
    x(12) = 0.5d0 * (x(8) + x(16))
    x(2) = 0.5d0 * (x0 + x(4))
    x(6) = 0.5d0 * (x(4) + x(8))
    x(10) = 0.5d0 * (x(8) + x(12))
    x(14) = 0.5d0 * (x(12) + x(16))
    DO j = 2, 16, 2
       f(j) = fun(x(j))
    END DO
    nofun = 9
    !
    !  estimate total integral using 4 panel rule
    !
    qprev = ( wf0 * (f0 + f(16)) &
         &          +  wf1 * (f(4) + f(12)) + wf2 * f(8) ) * (b-a)
    qnow = (  wf0 * ( f0 + 2.d0 * f(8) + f(16) ) &
         &        + wf1 * ( f(2) + f(6) + f(10) + f(14) ) &
         &        + wf2 * ( f(4) + f(12) ) ) * (0.5d0*(b-a))

    qdiff = qnow - qprev
    cor7 = z0 * qdiff
    esterr = ABS( cor7 )
    tolerr =  abserr
    IF ( esterr .LE. tolerr .OR. abserr .LT. 0.d0 ) THEN
       RESULT = qnow + cor7
       errest = esterr
       RETURN
    END IF
    !
    !  since 4 panel rule failed, calculate initial value using 8 panel rule
    !
    qprev = ( &
         &          w0 * (f0   + f(16)) &
         &      +   w1 * (f(2) + f(14)) &
         &      +   w2 * (f(4) + f(12)) &
         &      +   w3 * (f(6) + f(10)) &
         &      +   w4 * f(8) &
         &        ) * (0.125d0*(b-a))
    area = qprev
    !
    !  ***  stage 3 ***  central calculation
    !  requires qprev, x0, x2, x4, ..., x16, f0, f2, f4, ..., f16.
    !  calculates x1, x3, ..., x15, f1, f3, ..., f15, qleft, qnow, qdiff,
    !  area.
    !
30  CONTINUE
    x(1) = 0.5d0 * (x0 + x(2))
    f(1) = fun(x(1))
    DO j = 3, 15, 2
       x(j) = 0.5d0 * (x(j-1) + x(j+1))
       f(j) = fun(x(j))
    END DO
    nofun = nofun + 8
    step = 0.0625d0 * (x(16) - x0)
    qleft = ( &
         &          w0 * (f0   + f(8)) &
         &      +   w1 * (f(1) + f(7)) &
         &      +   w2 * (f(2) + f(6)) &
         &      +   w3 * (f(3) + f(5)) &
         &      +   w4 * f(4) &
         &        ) * step
    qright(lev+1) =  ( &
         &          w0 * (f(8)  + f(16)) &
         &      +   w1 * (f(9)  + f(15)) &
         &      +   w2 * (f(10) + f(14)) &
         &      +   w3 * (f(11) + f(13)) &
         &      +   w4 *  f(12) &
         &        ) * step
    qnow = qleft + qright(lev+1)
    qdiff = qnow - qprev
    area = area + qdiff
    !
    !  ***  stage 4 ***  interval convergence test
    !
    esterr = z1 * ABS(qdiff)
    tolerr = abserr * (step /stone)

    IF ( lev .LT. levmin ) GOTO 50
    IF ( lev .GE. levmax ) GOTO 62
    IF ( nofun .GT. nofin ) GOTO 60
    IF (esterr .LE. tolerr ) GOTO 70
    !
    !  ***  stage 5 ***  no convergence
    !  locate next interval
    !
50  CONTINUE
    nim = 2 * nim
    lev = lev + 1
    !
    !  store right hand elements for future use
    !
    DO i = 1, 8
       fsave(i,lev) = f(i+8)
       xsave(i,lev) = x(i+8)
    END DO
    !
    !  assemble left hand elements for immediate use
    !
    qprev = qleft
    DO i = 1, 8
       j = -i
       f(2*j+18) = f(j+9)
       x(2*j+18) = x(j+9)
    END DO
    GOTO 30
    !
    !  ***  stage 6 ***  trouble section
    !  number of function values is about to exceed limit
    !
60  CONTINUE
    nofin = 2 * nofin
    levmax = levout
    flag = flag + (b - x0) / (b - a)
    GOTO 70
    !
    !  current level is levmax
    !
62  CONTINUE
    flag = flag + 1.0
    !  ***  stage 7  *** interval converged
    !  add contributions into running sums
    !
70  CONTINUE
    RESULT = RESULT + qnow
    errest = errest + esterr
    cor11 = cor11 + z1 * qdiff
    !
    !  locate next interval
    !
72  CONTINUE
    IF (nim .EQ. 2*(nim/2)) GOTO 75
    nim = nim/2
    lev = lev - 1
    GOTO 72
75  CONTINUE
    nim = nim + 1
    IF ( lev .LE. 0) GOTO 80
    !
    !  assemble elements required for the next interval
    !
    qprev = qright(lev)
    x0 = x(16)
    f0 = f(16)
    DO i = 1, 8
       f(2*i) = fsave(i, lev)
       x(2*i) = xsave(i, lev)
    END DO
    GOTO 30
    !
    !  ***  stage 8  *** finalize and return
    !
80  CONTINUE
    RESULT = RESULT + cor11
    !
    !  make sure errest not less than roundoff level
    !
    IF ( errest .EQ. 0.d0 ) RETURN
82  CONTINUE
    temp = ABS(RESULT) + errest
    IF ( temp .NE. ABS(RESULT) ) RETURN
    errest = 2.d0 * errest
    GOTO 82
  END SUBROUTINE cquany


  SUBROUTINE ddblin(fun,ax,bx,ay,by,aerr, RESULT,errest,nofun,flag)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    !  numerical integration of a real*8 function of two real variables    c
    !                                                                      c
    !  this routine provides an estimate of the integral of fun(x) over    c
    !  the region ax .le. x .le. bx,   ay(x) .le. y .le. by(x)             c
    !                                                                      c
    !  an automatic adaptive routine based on the 8-panel newton-cotes     c
    !  rule.  (double precision version)                                   c
    !                                                                      c
    !  input ...                                                           c
    !                                                                      c
    !    fun      the name of the real*8 integrand function fun(x,y).      c
    !             fun must be declared external by the calling program.    c
    !                                                                      c
    !    ax, bx   the lower limit and upper limits of integration for the  c
    !             outer (over x) iterated integral.                        c
    !                                                                      c
    !    ay, by   double precision function subroutines provided by the    c
    !             user which take a single double precision argument (x)   c
    !             and return the lower and upper limits, respectively, of  c
    !             the inner iterated integral (over y).                    c
    !                                                                      c
    !    aerr   an absolute error tolerance (double precision)             c
    !             aerr is normally set to a nonnegative value which        c
    !             is the desired absolute accuracy of the estimate.        c
    !             if aerr is set to zero, the routine will return          c
    !             with an error indication, except for very simple         c
    !             integrands.  if aerr is set to a negative value          c
    !             the routine will estimate the value of the integral      c
    !             using the minimum possible number of integrand           c
    !             evaluations, 81.                                         c
    !                                                                      c
    !  output ...                                                          c
    !                                                                      c
    !    result   an approximation to the integral hopefully satisfying    c
    !             the least stringent of the two error tolerances.         c
    !             (double precision)                                             c
    !                                                                      c
    !    errest    an estimate of the magnitude of the actual error.       c
    !             (double precision)                                       c
    !                                                                      c
    !    nofun    the number of function values used in calculation of     c
    !             result.  (integer)                                       c
    !                                                                      c
    !    flag     a reliability indicator.  if flag is zero, then result   c
    !             probably satisfies the error tolerance.  if flag is      c
    !             greater than zero, then its integer part is a total      c
    !             count (for both outer and inner integration routines)    c
    !             of the number of intervals which have not converged.     c
    !             if flag has a nonzero fractional part, flag = xxx.yyy,   c
    !             then 0.yyy is the fraction of the outer interval (ax,bx) c
    !             which was left to do when the limit on nofun for the     c
    !             outer integral was approached.  if, however, flag has    c
    !             a nonzero integer part and zero fractional part, then    c
    !             the limit on nofun was reached one or more times during  c
    !             calculation(s) of the inner integral.  in this case,     c
    !             no information on the location of the trouble spot is    c
    !             returned to the user.  in any case, errest is still      c
    !             a useful estimate of the error actually attained.        c
    !             since the size of the interval where convergence failed  c
    !             is very small, the result of this routine is often       c
    !             acceptable when flag is not too large.                   c
    !                                                                      c
    !  description of algorithm:                                           c
    !              the double integral is evaluated as an iterated         c
    !              integral.  both the outer and inner integrals are       c
    !              evaluated using a modified version of quanc8            c
    !              which was taken from reference 1.                       c
    !              two major modifications were applied to quanc8.         c
    !              the first modification involved changing the routine    c
    !              to handle double precision functions and to use double  c
    !              precision throughout.  the second modification was to   c
    !              use a four panel newton-cotes rule on the first         c
    !              interval (a,b) and its two halves, in order to get      c
    !              an answer quickly, with few function evaluations, for   c
    !              particularly smooth integrands.  if the accuracy test   c
    !              for the four panel calculation fails, the previously    c
    !              implemented eight panel rule is used for the remainder  c
    !              of the calculations.                                    c
    !              the absolute error tolerance is divvied up between      c
    !              outer and inner integrals using the formula suggested   c
    !              in reference 2.                                         c
    !                                                                      c
    !                                                                      c
    !  references: 1.  forsythe, malcolm, and moler, 'computer             c
    !                  methods for mathematical computations,              c
    !                  prentiss-hall, 1977.                                c
    !              2.  fritsch, kahaner, and lyness, 'double integration   c
    !                  using one dimensional adaptive quadrature routines: c
    !                  a software interface problem,' acm trans. math.     c
    !                  software, vol 7, no. 1, march 1981, pp 46-75.       c
    !                                                                      c
    !  level of control                                                    c
    !                                                                      c
    !    this routine directly calls routines                              c
    !      cquany                                                          c
    !                                                                      c
    !                                                                      c
    !  this routine was initially entered on the resd ibm 4381 by          c
    !  p. simon, dept. 9282, x3726, on 13 oct. 1987.                       c
    !                                                                      c
    !  language:  fortran 77, as implemented on the resd ibm 4381          c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    IMPLICIT DOUBLE PRECISION (a-h, o-z)

    DOUBLE PRECISION fun, RESULT, area, f0, cor11, cor7
    DOUBLE PRECISION qprev, qnow, qdiff, qleft
    DOUBLE PRECISION qright(31), f(16), fsave(8,30), czero

    REAL flag, flagi
    DOUBLE PRECISION  abserr, errest, aerr, errin, errini
    DOUBLE PRECISION ax, bx, ay, by
    DOUBLE PRECISION w0, w1, w2, w3, w4
    DOUBLE PRECISION x0, stone, step, temp
    DOUBLE PRECISION esterr, tolerr
    DOUBLE PRECISION x(16)
    DOUBLE PRECISION xsave(8,30), epsini, wf0, wf1, wf2, xxx

    INTEGER nofun, nfuni
    INTEGER levmin, levmax, levout, nomax, nofin, lev, nim, i, j

    EXTERNAL fun, ay, by

    COMMON /xvalue/ xxx

    PARAMETER ( czero = (0.d0,0.d0) )
    !
    !  ***  stage 1 ***  general initialization
    !  set constants
    !
    PARAMETER ( levmin = 0 )
    PARAMETER ( levout = 6 )
    PARAMETER ( nomax = 300000)
    !      parameter ( z0 = 1.5873015873016d-2 )
    !      parameter ( z1 = 9.7751710654936d-4 )
    PARAMETER ( z0 = 1.d0 / 63.d0 )
    PARAMETER ( z1 = 1.d0 / 1023.d0 )
    !
    !  trouble when nofun reaches nofin
    !
    !      parameter ( w0 =   3956.d0 / 14175.d0 )
    !      parameter ( w1 =  23552.d0 / 14175.d0 )
    !      parameter ( w2 =  -3712.d0 / 14175.d0 )
    !      parameter ( w3 =  41984.d0 / 14175.d0 )
    !      parameter ( w4 = -18160.d0 / 14175.d0 )
    PARAMETER ( w0 =   0.27908289241623d0 )
    PARAMETER ( w1 = 1.6615167548501d0 )
    PARAMETER ( w2 =  -0.26186948853616d0 )
    PARAMETER ( w3 =  2.9618342151675d0 )
    PARAMETER ( w4 = -1.2811287477954d0 )

    !      parameter ( wf0 =  7.d0 / 90.d0 )
    !      parameter ( wf1 = 32.d0 / 90.d0 )
    !      parameter ( wf2 = 12.d0 / 90.d0 )
    PARAMETER ( wf0 =  7.7777777777778d-2 )
    PARAMETER ( wf1 = 0.35555555555556d0 )
    PARAMETER ( wf2 = 0.13333333333333d0 )
    !
    !  first executable statement
    !
    levmax = 30
    nofin = nomax - 8*(levmax-levout+2**(levout+1))
    !
    !  initialize running sums to zero
    !
    flag = 0.0
    RESULT = czero
    cor11 = czero
    errest = 0.d0
    nofun = 0
    errin = 0.0
    IF ( ax .EQ. bx ) RETURN
    !
    !  choose error tolerances for outer and inner integrals
    !
    abserr = 0.5d0 * aerr
    epsini = aerr / (  2.90758d0 * (bx - ax) )
    !
    !  ***  stage 2 ***  initialization for first interval
    !
    lev = 0
    nim = 1
    x0 = ax
    x(16) = bx

    CALL dquany(fun,x0,ay,by,epsini,f0,errini,nfuni,flagi)


    nofun = nofun + nfuni
    errin = MAX(errin, errini)
    flag = flag + INT(flagi)
    stone = 0.0625d0 * (bx-ax)
    x(8) = 0.5d0 * (x0 + x(16))
    x(4) = 0.5d0 * (x0 + x(8))
    x(12) = 0.5d0 * (x(8) + x(16))
    x(2) = 0.5d0 * (x0 + x(4))
    x(6) = 0.5d0 * (x(4) + x(8))
    x(10) = 0.5d0 * (x(8) + x(12))
    x(14) = 0.5d0 * (x(12) + x(16))
    DO j = 2, 16, 2
       CALL dquany(fun,x(j),ay,by,epsini,f(j),errini,nfuni,flagi)
       nofun = nofun + nfuni
       errin = MAX(errin, errini)
       flag = flag + INT(flagi)
    END DO
    !
    !  estimate total integral using 4 panel rule
    !
    qprev = ( wf0 * (f0 + f(16)) &
         &          +  wf1 * (f(4) + f(12)) + wf2 * f(8) ) * (bx-ax)
    qnow = (  wf0 * ( f0 + 2.d0 * f(8) + f(16) ) &
         & + wf1 * ( f(2) + f(6) + f(10) + f(14) ) &
         & + wf2 * ( f(4) + f(12) ) ) * (0.5d0*(bx-ax))

    qdiff = qnow - qprev
    cor7 = z0 * qdiff
    esterr = ABS( cor7 )
    tolerr =  abserr
    IF ( esterr .LE. tolerr .OR. abserr .LT. 0.d0 ) THEN
       RESULT = qnow + cor7
       errest = esterr  + (bx - ax) * errin
       RETURN
    END IF
    !
    !  since 4 panel rule failed, calculate initial value using 8 panel rule
    !
    qprev = ( &
         & w0 * (f0   + f(16)) &
         & +   w1 * (f(2) + f(14)) &
         & +   w2 * (f(4) + f(12)) &
         & +   w3 * (f(6) + f(10)) &
         & +   w4 * f(8) &
         & ) * (0.125d0*(bx-ax))
    area = qprev
    !
    !  ***  stage 3 ***  central calculation
    !  requires qprev, x0, x2, x4, ..., x16, f0, f2, f4, ..., f16.
    !  calculates x1, x3, ..., x15, f1, f3, ..., f15, qleft, qnow, qdiff,
    !  area.
    !
30  CONTINUE
    x(1) = 0.5d0 * (x0 + x(2))
    !     f(1) = fun(x(1))
    CALL dquany(fun,x(1),ay,by,epsini,f(1),errini,nfuni,flagi)
    nofun = nofun + nfuni
    errin = MAX(errin, errini)
    flag = flag + INT(flagi)
    DO j = 3, 15, 2
       x(j) = 0.5d0 * (x(j-1) + x(j+1))
       !       f(j) = fun(x(j))
       CALL dquany(fun,x(j),ay,by,epsini,f(j),errini,nfuni,flagi)
       nofun = nofun + nfuni
       errin = MAX(errin, errini)
       flag = flag + INT(flagi)
    END DO

    step = 0.0625d0 * (x(16) - x0)
    qleft = ( &
         & w0 * (f0   + f(8)) &
         & +   w1 * (f(1) + f(7)) &
         & +   w2 * (f(2) + f(6)) &
         & +   w3 * (f(3) + f(5)) &
         & +   w4 * f(4) &
         & ) * step
    qright(lev+1) =  ( &
         & w0 * (f(8)  + f(16)) &
         & +   w1 * (f(9)  + f(15)) &
         & +   w2 * (f(10) + f(14)) &
         & +   w3 * (f(11) + f(13)) &
         & +   w4 *  f(12) &
         &   ) * step
    qnow = qleft + qright(lev+1)
    qdiff = qnow - qprev
    area = area + qdiff
    !
    !  ***  stage 4 ***  interval convergence test
    !
    esterr = z1 * ABS(qdiff)
    tolerr = abserr * (step /stone)

    IF ( lev .LT. levmin ) GOTO 50
    IF ( lev .GE. levmax ) GOTO 62
    IF ( nofun .GT. nofin ) GOTO 60
    IF (esterr .LE. tolerr ) GOTO 70
    !
    !  ***  stage 5 ***  no convergence
    !  locate next interval
    !
50  CONTINUE
    nim = 2 * nim
    lev = lev + 1
    !
    !  store right hand elements for future use
    !
    DO i = 1, 8
       fsave(i,lev) = f(i+8)
       xsave(i,lev) = x(i+8)
    END DO
    !
    !  assemble left hand elements for immediate use
    !
    qprev = qleft
    DO i = 1, 8
       j = -i
       f(2*j+18) = f(j+9)
       x(2*j+18) = x(j+9)
    END DO
    GOTO 30
    !
    !  ***  stage 6 ***  trouble section
    !  number of function values is about to exceed limit
    !
60  CONTINUE
    nofin = 2 * nofin
    levmax = levout
    flag = flag + (bx - x0) / (bx - ax)
    GOTO 70
    !
    !  current level is levmax
    !
62  CONTINUE
    flag = flag + 1.0
    !  ***  stage 7  *** interval converged
    !  add contributions into running sums
    !
70  CONTINUE
    RESULT = RESULT + qnow
    errest = errest + esterr
    cor11 = cor11 + z1 * qdiff
    !
    !  locate next interval
    !
72  CONTINUE
    IF (nim .EQ. 2*(nim/2)) GOTO 75
    nim = nim/2
    lev = lev - 1
    GOTO 72
75  CONTINUE
    nim = nim + 1
    IF ( lev .LE. 0) GOTO 80
    !
    !  assemble elements required for the next interval
    !
    qprev = qright(lev)
    x0 = x(16)
    f0 = f(16)
    DO i = 1, 8
       f(2*i) = fsave(i, lev)
       x(2*i) = xsave(i, lev)
    END DO
    GOTO 30
    !
    !  ***  stage 8  *** finalize and return
    !
80  CONTINUE
    RESULT = RESULT + cor11
    !
    !  make sure errest not less than roundoff level
    !
    errest = errest + (bx - ax) * errin
    IF ( errest .EQ. 0.d0 ) RETURN
82  CONTINUE
    temp = ABS(RESULT) + errest
    IF ( temp .NE. ABS(RESULT) ) RETURN
    errest = 2.d0 * errest
    GOTO 82
  END SUBROUTINE ddblin

  SUBROUTINE dquany(dfun,y,aa,bb,abserr,RESULT,errest,nofun,flag)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !                                                                      c
    !  this routine provides an estimate of the integral of dfun(y,x) from c
    !  x = a to x = b to a user provided tolerance, where dfun(y,x) is a   c
    !  user supplied double precision function of two double-precision           c
    !  arguments, y and x.                                                 c
    !                                                                      c
    !  an automatic adaptive routine based on the 8-panel newton-cotes     c
    !  rule.  (double precision version)                                         c
    !                                                                      c
    !  input ...                                                           c
    !                                                                      c
    !    dfun     the name of the double precision integrand function dfun(x).   c
    !             dfun must be declared external by the calling program.   c
    !                                                                      c
    !    y        the value of the second argument of dfun(x,y)            c
    !                                                                      c
    !    aa, bb   the names of user supplied double precision function     c
    !             subroutines which return the inner integration limits    c
    !             given the current value of the outer variable of         c
    !             integration. they must be declared external by the       c
    !             calling subroutine.                                      c
    !                                                                      c
    !    abserr   an absolute error tolerance (convergence criterion)      c
    !             should be nonnegative.   (double precision)              c
    !             if abserr is less than zero, the routine will compute    c
    !             an estimate of the integral using the minimum possible   c
    !             number of function evaluations (9), then return.         c
    !                                                                      c
    !  output ...                                                          c
    !                                                                      c
    !    result   an approximation to the integral hopefully satisfying    c
    !             the error tolerance. (double precision)                        c
    !                                                                      c
    !    errest   an estimate of the magnitude of the actual error.        c
    !             (double precision)                                       c
    !                                                                      c
    !    nofun    the number of function values used in calculation of     c
    !             result.  (integer)                                       c
    !                                                                      c
    !    flag     a reliability indicator.  if flag is zero, then result   c
    !             probably satisfies the error tolerance.  if flag is      c
    !             xxx.yyy, then xxx = the number of intervals which have   c
    !             not converged and 0.yyy = the fraction of the interval   c
    !             left to do when the limit on nofun was approached.       c
    !             (double precision)                                       c
    !                                                                      c
    !  reference:  this routine is a modified version of quanc8            c
    !              which was taken from the textbook 'computer             c
    !              methods for mathematical computations,' by              c
    !              forsythe, malcolm, and moler,  prentiss-hall, 1977.     c
    !              the first modification involved changing the routine    c
    !              to handle double precision functions and to use double        c
    !              precision throughout.  the second modification is to    c
    !              use a four panel newton-cotes rule on the first         c
    !              interval (a,b) and its two halves, in order to get      c
    !              an answer quickly, with few function evaluations, for   c
    !              particularly smooth integrands.  if the accuracy test   c
    !              for the four panel calculation fails, the previously    c
    !              implemented eight panel rule is used for the remainder  c
    !              of the calculations.                                    c
    !                                                                      c
    !  language:  fortran 77, as implemented on the resd ibm 4381          c
    !                                                                      c
    !  this routine was initially entered on the resd ibm 4381 by          c
    !  p. simon, dept. 9282, x3726, on 13 oct. 1987.                       c
    !  Modified 4/7/94 to add external statement.                          c
    !                                                                      c
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    EXTERNAL   dfun, aa, bb
    DOUBLE PRECISION dfun, RESULT, area, f0, cor11, cor7, fun
    DOUBLE PRECISION qprev, qnow, qdiff, qleft
    DOUBLE PRECISION qright(31), f(16), fsave(8,30), czero

    REAL flag
    DOUBLE PRECISION a, b, abserr, errest, y
    DOUBLE PRECISION w0, w1, w2, w3, w4
    DOUBLE PRECISION x0, stone, step, temp
    DOUBLE PRECISION esterr, tolerr
    DOUBLE PRECISION x(16)
    DOUBLE PRECISION xsave(8,30), t, aa, bb, wf0, wf1, wf2

    INTEGER nofun
    INTEGER levmin, levmax, levout, nomax, nofin, lev, nim, i, j

    PARAMETER ( czero = (0.d0,0.d0) )
    !
    !  ***  stage 1 ***  general initialization
    !  set constants
    !
    PARAMETER ( levmin = 0 )
    PARAMETER ( levout = 6 )
    PARAMETER ( nomax = 300000 )
    !      parameter ( z0 = 1.d0 / 63.d0 )
    !      parameter ( z1 = 1.d0 / 1023.d0 )
    PARAMETER ( z0 = 1.5873015873016d-2 )
    PARAMETER ( z1 = 9.7751710654936d-4 )
    !
    !  trouble when nofun reaches nofin
    !
    !      parameter ( w0 =   3956.d0 / 14175.d0 )
    !      parameter ( w1 =  23552.d0 / 14175.d0 )
    !      parameter ( w2 =  -3712.d0 / 14175.d0 )
    !      parameter ( w3 =  41984.d0 / 14175.d0 )
    !      parameter ( w4 = -18160.d0 / 14175.d0 )
    PARAMETER ( w0 =   0.27908289241623d0 )
    PARAMETER ( w1 = 1.6615167548501d0 )
    PARAMETER ( w2 =  -0.26186948853616d0 )
    PARAMETER ( w3 =  2.9618342151675d0 )
    PARAMETER ( w4 = -1.2811287477954d0 )

    !      parameter ( wf0 =  7.d0 / 90.d0 )
    !      parameter ( wf1 = 32.d0 / 90.d0 )
    !      parameter ( wf2 = 12.d0 / 90.d0 )
    PARAMETER ( wf0 =  7.7777777777778d-2 )
    PARAMETER ( wf1 = 0.35555555555556d0 )
    PARAMETER ( wf2 = 0.13333333333333d0 )

    fun(t) = dfun(y,t)
    !
    !  first executable statement
    !
    levmax = 30
    nofin = nomax - 8*(levmax-levout+2**(levout+1))
    !
    !  initialize running sums to zero
    !
    flag = 0.0
    RESULT = czero
    cor11 = czero
    errest = 0.d0
    nofun = 0
    !
    !  get integration limits
    !
    a = aa(y)
    b = bb(y)
    IF ( a .EQ. b ) RETURN
    !
    !  ***  stage 2 ***  initialization for first interval
    !
    lev = 0
    nim = 1
    x0 = a
    x(16) = b
    f0 = fun(x0)
    stone = 0.0625d0 * (b-a)
    x(8) = 0.5d0 * (x0 + x(16))
    x(4) = 0.5d0 * (x0 + x(8))
    x(12) = 0.5d0 * (x(8) + x(16))
    x(2) = 0.5d0 * (x0 + x(4))
    x(6) = 0.5d0 * (x(4) + x(8))
    x(10) = 0.5d0 * (x(8) + x(12))
    x(14) = 0.5d0 * (x(12) + x(16))
    DO j = 2, 16, 2
       f(j) = fun(x(j))
    END DO
    nofun = 9
    !
    !  estimate total integral using 4 panel rule
    !
    qprev = ( wf0 * (f0 + f(16)) &
         &          +  wf1 * (f(4) + f(12)) + wf2 * f(8) ) * (b-a)
    qnow = (  wf0 * ( f0 + 2.d0 * f(8) + f(16) ) &
         &        + wf1 * ( f(2) + f(6) + f(10) + f(14) ) &
         &        + wf2 * ( f(4) + f(12) ) ) * (0.5d0*(b-a))

    qdiff = qnow - qprev
    cor7 = z0 * qdiff
    esterr = ABS( cor7 )
    tolerr =  abserr
    IF ( esterr .LE. tolerr .OR. abserr .LT. 0.d0 ) THEN
       RESULT = qnow + cor7
       errest = esterr
       RETURN
    END IF
    !
    !  since 4 panel rule failed, calculate initial value using 8 panel rule
    !
    qprev = ( &
         &          w0 * (f0   + f(16)) &
         &      +   w1 * (f(2) + f(14)) &
         &      +   w2 * (f(4) + f(12)) &
         &      +   w3 * (f(6) + f(10)) &
         &      +   w4 * f(8) &
         &        ) * (0.125d0*(b-a))
    area = qprev
    !
    !  ***  stage 3 ***  central calculation
    !  requires qprev, x0, x2, x4, ..., x16, f0, f2, f4, ..., f16.
    !  calculates x1, x3, ..., x15, f1, f3, ..., f15, qleft, qnow, qdiff,
    !  area.
    !
30  CONTINUE
    x(1) = 0.5d0 * (x0 + x(2))
    f(1) = fun(x(1))
    DO j = 3, 15, 2
       x(j) = 0.5d0 * (x(j-1) + x(j+1))
       f(j) = fun(x(j))
    END DO
    nofun = nofun + 8
    step = 0.0625d0 * (x(16) - x0)
    qleft = ( &
         &          w0 * (f0   + f(8)) &
         &      +   w1 * (f(1) + f(7)) &
         &      +   w2 * (f(2) + f(6)) &
         &      +   w3 * (f(3) + f(5)) &
         &      +   w4 * f(4) &
         &        ) * step
    qright(lev+1) =  ( &
         &          w0 * (f(8)  + f(16)) &
         &      +   w1 * (f(9)  + f(15)) &
         &      +   w2 * (f(10) + f(14)) &
         &      +   w3 * (f(11) + f(13)) &
         &      +   w4 *  f(12) &
         &        ) * step
    qnow = qleft + qright(lev+1)
    qdiff = qnow - qprev
    area = area + qdiff
    !
    !  ***  stage 4 ***  interval convergence test
    !
    esterr = z1 * ABS(qdiff)
    tolerr = abserr * (step /stone)

    IF ( lev .LT. levmin ) GOTO 50
    IF ( lev .GE. levmax ) GOTO 62
    IF ( nofun .GT. nofin ) GOTO 60
    IF (esterr .LE. tolerr ) GOTO 70
    !
    !  ***  stage 5 ***  no convergence
    !  locate next interval
    !
50  CONTINUE
    nim = 2 * nim
    lev = lev + 1
    !
    !  store right hand elements for future use
    !
    DO i = 1, 8
       fsave(i,lev) = f(i+8)
       xsave(i,lev) = x(i+8)
    END DO
    !
    !  assemble left hand elements for immediate use
    !
    qprev = qleft
    DO i = 1, 8
       j = -i
       f(2*j+18) = f(j+9)
       x(2*j+18) = x(j+9)
    END DO
    GOTO 30
    !
    !  ***  stage 6 ***  trouble section
    !  number of function values is about to exceed limit
    !
60  CONTINUE
    nofin = 2 * nofin
    levmax = levout
    flag = flag + (b - x0) / (b - a)
    GOTO 70
    !
    !  current level is levmax
    !
62  CONTINUE
    flag = flag + 1.0
    !  ***  stage 7  *** interval converged
    !  add contributions into running sums
    !
70  CONTINUE
    RESULT = RESULT + qnow
    errest = errest + esterr
    cor11 = cor11 + z1 * qdiff
    !
    !  locate next interval
    !
72  CONTINUE
    IF (nim .EQ. 2*(nim/2)) GOTO 75
    nim = nim/2
    lev = lev - 1
    GOTO 72
75  CONTINUE
    nim = nim + 1
    IF ( lev .LE. 0) GOTO 80
    !
    !  assemble elements required for the next interval
    !
    qprev = qright(lev)
    x0 = x(16)
    f0 = f(16)
    DO i = 1, 8
       f(2*i) = fsave(i, lev)
       x(2*i) = xsave(i, lev)
    END DO
    GOTO 30
    !
    !  ***  stage 8  *** finalize and return
    !
80  CONTINUE
    RESULT = RESULT + cor11
    !
    !  make sure errest not less than roundoff level
    !
    IF ( errest .EQ. 0.d0 ) RETURN
82  CONTINUE
    temp = ABS(RESULT) + errest
    IF ( temp .NE. ABS(RESULT) ) RETURN
    errest = 2.d0 * errest
    GOTO 82
  END SUBROUTINE dquany

END MODULE dblint
