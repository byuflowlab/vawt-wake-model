from os import path
import numpy as np
from numpy import pi,sin,cos,arccos,fabs
from scipy.integrate import quad
from scipy.optimize import root
import h5py
import VAWT_Wake_Model as vwm

import _vawtwake

"""
applies for both Ay and Rx depending on which function ifunc(x, y, phi)
is passed in
"""
def panelIntegration(xvec,yvec,thetavec,ifunc):

    # initialize
    nx = np.size(xvec)
    ntheta = np.size(thetavec)
    dtheta = thetavec[1] - thetavec[0]  # assumes equally spaced
    A = np.zeros((nx,ntheta))

    for i in range(nx):
        # redefine function so it has one parameter for use in quad
        integrand = lambda phi: ifunc(xvec[i], yvec[i], phi)

        for j in range(ntheta):
            # an Adaptive Gauss-Kronrod quadrature integration.  Tried trapz but this was faster.
            A[i, j], error = quad(integrand, thetavec[j]-dtheta/2., thetavec[j]+dtheta/2., epsabs=1e-10)

    return A

"""
integrand used for computing Dx
"""
def Dxintegrand(x,y,phi):
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    # v1 and v2 should never both be zero b.c. we never integrate self.  RxII handles that case.
    return (v1*sin(phi) - v2*cos(phi))/(2.*pi*(v1*v1 + v2*v2))


"""
integrand used for computing Ay
"""
def Ayintegrand(x,y,phi):
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    if fabs(v1) < 1e-12 and abs(v2) < 1e-12:  # occurs when integrating self, function symmetric around singularity, should integrate to zero
        return 0.
    else:
        return (v1*cos(phi) + v2*sin(phi))/(2.*pi*(v1*v1 + v2*v2))

def AyIJ(xvec,yvec,thetavec):
    return panelIntegration(xvec, yvec, thetavec, Ayintegrand)

def DxIJ(xvec,yvec,thetavec):
    return panelIntegration(xvec, yvec, thetavec, Dxintegrand)

def WxIJ(xvec,yvec,thetavec):

    # initialize
    nx = np.size(xvec)
    ntheta = np.size(thetavec)
    dtheta = thetavec[1] - thetavec[0]  # assumes equally spaced
    Wx = np.zeros((nx,ntheta))

    for i in range(nx):
        if yvec[i] >= -1.0 and yvec[i] <= 1.0 and xvec[i] >= 0.0 and xvec[i]^2 + yvec[i]^2 >= 1.0:
        # if yvec[i] >= -1.0 and yvec[i] <= 1.0 and (xvec[i] >= 0.0 or (xvec[i] >= -1 and xvec[i]^2 + yvec[i]^2 <= 1.0))
            thetak = arccos(yvec[i])
            for j in range(ntheta):
                if thetavec[j] + dtheta/2. > thetak:
                    k = j
                    break
            Wx[i, k] = -1.0
            Wx[i, ntheta-k+1] = 1.0

    return Wx

def DxII(thetavec):

    # initialize
    ntheta = np.size(thetavec)
    dtheta = thetavec[1] - thetavec[0]  # assumes equally spaced
    Rx = dtheta/(4.*pi)*np.ones((ntheta,ntheta))

    for i in range(ntheta):
        if i <= ntheta/2.-1:
            Rx[i, i] = (-1. + 1./ntheta)/2.
        else:
            Rx[i, i] = (1. + 1./ntheta)/2.

    return Rx

def WxII(thetavec):

    # initialize
    ntheta = np.size(thetavec)
    Wx = np.zeros((ntheta,ntheta))

    for i in range(ntheta/2,ntheta):
        Wx[i, ntheta-1-i] = -1.

    return Wx


def precomputeMatrices(ntheta):

    # precompute self influence matrices

    # setup discretization (all the same, and uniformly spaced in theta)
    dtheta = 2.*pi/ntheta
    theta = np.arange(dtheta/2.,2.*pi,dtheta)

    Dxself = DxII(theta)
    Wxself = WxII(theta)
    Ayself = AyIJ(-sin(theta), cos(theta), theta)

    # write to file
    basepath = path.join(path.dirname(path.realpath('__file__')), 'data')
    fdata = basepath + path.sep + 'theta-%d.h5' %ntheta
    with h5py.File(fdata, 'w') as hf:
        hf.create_dataset('theta', data=theta)
        hf.create_dataset('Dx', data=Dxself)
        hf.create_dataset('Wx', data=Wxself)
        hf.create_dataset('Ay', data=Ayself)


def matrixAssemble(centerX,centerY,radii,ntheta):
    """
    centerX, centerY: array of x,y coordinates for centers of the VAWTs in the farm
    radii: corresponding array of their radii
    """

    try:
        basepath = path.join(path.dirname(path.realpath('__file__')), 'data')
        fdata = basepath + path.sep + 'theta-%d.h5' %ntheta
        with h5py.File(fdata,'r') as hf:
            theta = np.array(hf.get('theta'))
            Dxself = np.array(hf.get('Dx'))
            Wxself = np.array(hf.get('Wx'))
            Ayself = np.array(hf.get('Ay'))
    except:
        precomputeMatrices(ntheta)

        basepath = path.join(path.dirname(path.realpath('__file__')), 'data')
        fdata = basepath + path.sep + 'theta-%d.h5' %ntheta
        with h5py.File(fdata,'r') as hf:
            theta = np.array(hf.get('theta'))
            Dxself = np.array(hf.get('Dx'))
            Wxself = np.array(hf.get('Wx'))
            Ayself = np.array(hf.get('Ay'))

    # print Dxself


    # initialize global matrices
    nturbines = np.size(radii)
    Dx = np.zeros((nturbines*ntheta,nturbines*ntheta))
    Wx = np.zeros((nturbines*ntheta,nturbines*ntheta))
    Ay = np.zeros((nturbines*ntheta,nturbines*ntheta))

    # iterate through turbines
    for I in range(nturbines):
        for J in range(nturbines):

            # find normalized i locations relative to center of turbine J
            x = (centerX[I]-radii[I]*sin(theta) - centerX[J])/radii[J]
            y = (centerY[I]+radii[I]*cos(theta) - centerY[J])/radii[J]

            # self-influence is precomputed
            if I == J:
                Dxsub = Dxself
                Wxsub = Wxself
                Aysub = Ayself

            # pairs can be mapped for same radius
            elif J < I and radii[I] == radii[J]:

                # grab cross-diagonal I,J -> J,I matrix
                for K in range((J)*ntheta,(J+1)*ntheta):
                    for L in range((I)*ntheta,(I+1)*ntheta):
                        Dxsub[K,L] = Dx[K,L]
                        Aysub[K,L] = Ay[K,L]

                # mapping index for coefficients that are the same
                for M in range(ntheta/2,ntheta):
                    for N in range(0,ntheta/2):

                        # directly map over
                        Dxsub = Dxsub[M,N]
                        Aysub = Aysub[M,N]

                # wake term must be recomptued
                Wxsub = WxIJ(x, y, theta)

            # # if VAWTs are very far apart we can approximate some of the influence coefficients
            # elseif approxfar && sqrt((centerX[I]-centerX[J])^2 + (centerY[I]-centerY[J])^2) > 10*radii[I]
            #     println("far apart")
            #     xc = (centerX[I] - centerX[J])/radii[J]
            #     yc = (centerY[I] - centerY[J])/radii[J]

            #     Rxsub = RxIJFar(xc, yc, theta)
            #     Wxsub = zeros(ntheta, ntheta)  # should have negligible wake impact
            #     Aysub = AyIJFar(xc, yc, theta)

            else:
                Dxsub = DxIJ(x, y, theta)
                Wxsub = WxIJ(x, y, theta)
                Aysub = AyIJ(x, y, theta)

            # assemble into global matrix
            for O in range((I)*ntheta,(I+1)*ntheta):
                for P in range((J)*ntheta,(J+1)*ntheta):
                    Dx[O,P] = Dxsub[O,P]
                    Wx[O,P] = Wxsub[O,P]
                    Ay[O,P] = Aysub[O,P]

    Ax = Dx + Wx

    return Ax, Ay, theta


def residual(w,A,theta,af_data,cl_data,cd_data,r,chord,twist,delta,B,Omega,Vinf,Vinfx,Vinfy,rho,mu,interp):

    # setup
    ntheta = np.size(theta)

    uvec = np.zeros(ntheta)
    vvec = np.zeros(ntheta)

    for i in range(ntheta):
        uvec[i] = w[i]
        vvec[i] = w[ntheta + i]

    q,k,_,_,_,_,_ = _vawtwake.radialforce(uvec,vvec,theta,af_data,cl_data,cd_data,r,chord,twist,delta,B,Omega,Vinf,Vinfx,Vinfy,rho,mu,interp)

    # reformat to multiply in correct locations
    kmult = np.ones(2*ntheta)*k

    return np.dot(A,q)*kmult - w


def actuatorcylinder(centerX,centerY,radii,ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,Omega,Vinf,rho,mu,interp):

    # assemble global matrices
    Ax, Ay, theta = matrixAssemble(centerX, centerY, radii, ntheta)
    A = np.vstack((Ax,Ay))

    # setup
    ntheta = np.size(theta)
    Vinfx = np.zeros(ntheta)
    Vinfy = np.zeros(ntheta)
    tol_root = 1e-6

    w0 = np.zeros(ntheta*2)
    res = root(residual,w0,args=(A,theta,af_data,cl_data,cd_data,r,chord,twist,delta,B,Omega,Vinf,Vinfx,Vinfy,rho,mu,interp),method='hybr',tol=tol_root)
    w = res.x

    uvec = np.zeros(ntheta)
    vvec = np.zeros(ntheta)
    for i in range(ntheta):
        uvec[i] = w[i]
        vvec[i] = w[ntheta + i]

    _,_,Ct,Cp,Rp,Tp,Zp = _vawtwake.radialforce(uvec,vvec,theta,af_data,cl_data,cd_data,r,chord,twist,delta,B,Omega,Vinf,Vinfx,Vinfy,rho,mu,interp)

    return uvec,vvec,Ct,Cp,Rp,Tp,Zp


def induced_vel(r,af_data,cl_data,cd_data,chord,twist,delta,B,rot,velf,rho,mu,ntheta,interp):

    velx,vely,_,_,_,_,_ = actuatorcylinder(np.array([0.0]),np.array([0.0]),np.array([r]),ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp)

    return velx,vely

def CtCpcalc(r,af_data,cl_data,cd_data,chord,twist,delta,B,rot,velf,rho,mu,ntheta,interp):

    _,_,Ct,Cp,Rp,Tp,Zp = actuatorcylinder(np.array([0.0]),np.array([0.0]),np.array([r]),ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp)

    return Ct,Cp