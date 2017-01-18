from os import path
import numpy as np
from numpy import pi,sin,cos,arccos,fabs
from scipy.integrate import quad
from scipy.optimize import root
import h5py
import VAWT_Wake_Model as vwm

import _vawtwake

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

def Ayintegrand(x,y,phi):
    v1 = x + sin(phi)
    v2 = y - cos(phi)
    if fabs(v1) < 1e-12 and abs(v2) < 1e-12:  # occurs when integrating self, function symmetric around singularity, should integrate to zero
        return 0.
    else:
        return (v1*cos(phi) + v2*sin(phi))/(2.*pi*(v1*v1 + v2*v2))

def AyIJ(xvec,yvec,thetavec):
    return panelIntegration(xvec, yvec, thetavec, Ayintegrand)

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

def matrixAssemble(ntheta):
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

    # initialize global matrices
    Dx = np.zeros((ntheta,ntheta))
    Wx = np.zeros((ntheta,ntheta))
    Ay = np.zeros((ntheta,ntheta))

    # self-influence is precomputed
    Dxsub = Dxself
    Wxsub = Wxself
    Aysub = Ayself

    # assemble into global matrix
    for O in range(0,ntheta):
        for P in range(0,ntheta):
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

def actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,Omega,Vinf,rho,mu,interp,Vinfx,Vinfy):

    # assemble global matrices
    Ax, Ay, theta = matrixAssemble(ntheta)
    A = np.vstack((Ax,Ay))

    # setup
    ntheta = np.size(theta)
    # Vinfx = np.zeros(ntheta)
    # Vinfy = np.zeros(ntheta)
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