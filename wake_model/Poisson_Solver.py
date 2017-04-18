from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from numpy import pi
from numpy.fft import fft,fft2,ifft2,ifft,irfft2,rfft2
from scipy.misc import derivative
import random as random
from mpl_toolkits.mplot3d import Axes3D
from VAWT_Wake_Model import velocity_field,coef_val
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

import _vawtwake

def _parameterval(tsr,sol,coef):
    """
    Creating polynomial surface based on given coefficients and calculating the point at a given TSR and solidity

    Parameters
    ----------
    tsr : float
        specified tip-speed ratio
    sol : float
        specified solidity
    coef : array
        the polynomial surface coefficients for a given EMG parameter

    Returns
    ----------
    surf : float
        the polynomial surface value for the EMG parameter based on tip-speed ratio and solidity
    """

    a = coef[0]
    b = coef[1]
    c = coef[2]
    d = coef[3]
    e = coef[4]
    f = coef[5]
    g = coef[6]
    h = coef[7]
    i = coef[8]
    j = coef[9]

    surf = a + b*tsr + c*sol + d*tsr**2 + e*tsr*sol + f*sol**2 + g*tsr**3 + h*tsr**2*sol + i*tsr*sol**2 + j*sol**3

    return surf

def xdiff(x,y,xt,yt,Vinf,dia,rot,c,B):
    return -velocity_field(xt,yt,x,y,Vinf,dia,rot,c,B,param=None,veltype='vort')*rot

def ydiff(y,x,xt,yt,Vinf,dia,rot,c,B):
    return velocity_field(xt,yt,x,y,Vinf,dia,rot,c,B,param=None,veltype='vort')*rot



def dst(x,axis=-1):
    """Discrete Sine Transform (DST-I)

    Implemented using 2(N+1)-point FFT
    xsym = r_[0,x,0,-x[::-1]]
    DST = (-imag(fft(xsym))/2)[1:(N+1)]

    adjusted to work over an arbitrary axis for entire n-dim array
    """
    n = len(x.shape)
    N = x.shape[axis]
    slices = [None]*3
    for k in range(3):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))
    newshape = list(x.shape)
    newshape[axis] = 2*(N+1)
    xsym = np.zeros(newshape,np.float)
    slices[0][axis] = slice(1,N+1)
    slices[1][axis] = slice(N+2,None)
    slices[2][axis] = slice(None,None,-1)
    for k in range(3):
        slices[k] = tuple(slices[k])
    xsym[slices[0]] = x
    xsym[slices[1]] = -x[slices[2]]
    DST = fft(xsym,axis=axis)
    #print xtilde
    return (-(DST.imag)/2)[slices[0]]

def dst2(x,axes=(-1,-2)):
    return dst(dst(x,axis=axes[0]),axis=axes[1])

def idst2(x,axes=(-1,-2)):
    return dst(dst(x,axis=axes[0]),axis=axes[1])

def fft_poisson(f,h):
    # Method taken from "2D Fast Poisson Solver" by arsenous
    # https://arsenous.wordpress.com/2013/03/22/221/

    m,n = f.shape
    f_bar = np.zeros([n,n])
    u_bar = f_bar			# make of the same shape
    u = u_bar

    f_bar = idst2(f)			# f_bar= fourier transform of f

    f_bar = f_bar*(2/n + 1)**2  #Normalize

    lam = np.arange(1,n+1)
    lam = -4/h**2 * (np.sin((lam*pi) / (2*(n + 1))))**2 #$compute $\lambda_x$
    #for rectangular domain add $lambda_y$
    for i in xrange(0,n):
        for j in xrange(0,n):
            u_bar[i,j] = (f_bar[i,j]) / (lam[i] + lam[j])

    u = dst2(u_bar)				#sine transform back
    u= u*(2/(n + 1))**2			#normalize ,change for rectangular domain
    return u

dia = 6.
#set bounds a,b,parameters
a = -5.*dia
b = 30.*dia
w = 5.*dia
alpha = 8				#alpha is grid points=2^alpha
n = 2**alpha
# n = 100
L = b-a					#length of system

xe = np.linspace(a,b,n)
ye = np.linspace(-w,w,n)
x,y = np.meshgrid(xe,ye)

h = L/(n)			#size
h2 = h**2			#h squared
hx2 = h2 			#started with a cube,hx2=hy2
hy2 = h2
f = np.zeros([n,n])


xt = 0.
yt = 0.
Vinf = 15.
rad = dia/2.
tsr = 4.
rot = tsr*Vinf/rad
B = 3
c = 0.25
solidity = B*c/rad

coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = coef_val()

loc1 = _parameterval(tsr,solidity,coef0)
loc2 = _parameterval(tsr,solidity,coef1)
loc3 = _parameterval(tsr,solidity,coef2)
spr1 = _parameterval(tsr,solidity,coef3)
spr2 = _parameterval(tsr,solidity,coef4)
skw1 = _parameterval(tsr,solidity,coef5)
skw2 = _parameterval(tsr,solidity,coef6)
scl1 = _parameterval(tsr,solidity,coef7)
scl2 = _parameterval(tsr,solidity,coef8)
scl3 = _parameterval(tsr,solidity,coef9)

vel_type = 'x'
vel_type = 'y'
k = 1
for i in range(n):
    for j in range(n):
        if vel_type == 'x':
            # x-velocity
            yargs = (x[i,j],xt,yt,Vinf,dia,rot,c,B)
            f[i,j] = -derivative(ydiff,y[i,j],dx=1e-6,n=1,args=yargs,order=3)
        elif vel_type == 'y':
            # y-velocity
            yargs = (y[i,j],xt,yt,Vinf,dia,rot,c,B)
            f[i,j] = -derivative(xdiff,x[i,j],dx=1e-6,n=1,args=yargs,order=3)

        print int((k/(n*n))*100), '%'
        k += 1



nx,ny = np.array(f.shape)-1		 #last index
U = np.zeros([n,n])

# BOUNDARY CONDITIONS
far = 0.
U[0,:] = far
U[nx,:] = far
U[:,0] = far
U[:,ny] = far

#homogenize boundary condition
f[1,:] = f[1,:] + U[0,:]/hx2;
f[nx-1,:] = f[nx-1,:] + U[n-1,:]/hx2;
f[:,1] = f[:,1] + U[:,0]/hy2;
f[:,ny-1] = f[:,ny-1] + U[:,n-1]/hy2;

U = fft_poisson(f,h)

fs = 25
xi = -3.*dia # starting point in downstream direction
xf = 17.0*dia # ending point in downstream direction
yd = -2.5*dia # lateral extent on down side
yu = 2.5*dia # lateral extent on up side
fig = plt.figure(2,figsize=(19,5))
fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
if vel_type == 'x':
    lb = 0.15 # lower bound on velocity to display
    ub = 1.15 # upper bound on velocity to display
elif vel_type == 'y':
    lb = -0.35 # lower bound on velocity to display
    ub = 0.35 # upper bound on velocity to display
ran = 32 # number of contours between the velocity bounds
bounds = np.linspace(lb,ub,ran)
v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
if vel_type == 'x':
    CS = plt.contourf(x/dia,y/dia,(U + Vinf)/Vinf,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
elif vel_type == 'y':
    CS = plt.contourf(x/dia,y/dia,U/Vinf,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
CB = plt.colorbar(CS, ticks=v) # creating colorbar
if vel_type == 'x':
    CB.ax.set_ylabel(r'$u/U_\infty$',fontsize=fs)
elif vel_type == 'y':
    CB.ax.set_ylabel(r'$v/U_\infty$',fontsize=fs)
CB.ax.tick_params(labelsize=fs)
CB.ax.set_aspect(20)
plt.xlabel('$x/D$',fontsize=fs)
plt.ylabel('$y/D$',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlim(xi/dia,xf/dia)
plt.ylim(yd/dia,yu/dia)
circ = plt.Circle((xt/dia,yt/dia),0.5,edgecolor='k',fill=False)
plt.gca().add_patch(circ)

plt.show()
