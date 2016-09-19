# Originally created by Travis Oliphant (2001) and Nathan Woods (2013) (nquad &c)

from scipy.integrate import _quadpack

def quad(func, a, b, args=(), full_output=0, epsabs=1.49e-8, epsrel=1.49e-8, limit=50):

    retval = _quadpack._qagse(func,a,b,args,full_output,epsabs,epsrel,limit)

    return retval[:-1]

def _infunc(x,func,gfun,hfun,more_args):
    a = gfun(x)
    b = hfun(x)
    myargs = (x,) + more_args

    return quad(func,a,b,args=myargs)[0]


def dblquad(func, a, b, gfun, hfun, args=(), epsabs=1.49e-8, epsrel=1.49e-8):

    return quad(_infunc, a, b, (func, gfun, hfun, args), epsabs=epsabs, epsrel=epsrel)


