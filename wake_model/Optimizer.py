from pyoptsparse import Optimization, SNOPT, pyOpt_solution
from os import path
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
import VAWT_Wake_Model as vwm
from AC_InducedVel import induced_vel

import _vawtwake


def obj_func(xdict):
    global dia
    global rot
    global chord
    global twist
    global delta
    global B
    global H
    global rho
    global mu
    global velf

    global af_data
    global cl_data
    global cd_data
    global coef0
    global coef1
    global coef2
    global coef3
    global coef4
    global coef5
    global coef6
    global coef7
    global coef8
    global coef9

    global funcs

    ntheta = 36
    
    Vars = xdict['xvars']
    x = xdict['xvars']
    y = xdict['yvars']
    funcs = {}

    power_turb = np.zeros_like(x)

    for i in range(np.size(x)):
        xt = np.delete(x,i)
        yt = np.delete(y,i)
        diat = np.delete(dia,i)
        rott = np.delete(rot,i)
        xt = np.insert(xt,0,x[i])
        yt = np.insert(yt,0,y[i])
        diat = np.insert(diat,0,dia[i])
        rott = np.insert(rott,0,rot[i])

        velx,vely = induced_vel(dia[i]/2.,af_data,cl_data,cd_data,chord,twist,delta,B,rot[i],velf,rho,mu,ntheta)

        power_turb[i],_ = _vawtwake.powercalc(xt,yt,diat,rott,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely)
        print i+1
    
    power = np.sum(power_turb)
    print power 
    funcs['obj'] = (-1.0*power)
    
    sep = sep_func(np.append(x,y))
    funcs['sep'] = sep
    funcs['power_turb'] = power_turb
    
    fail = False

    return funcs, fail


def sep_func(loc):
    global turb_dia
    
    space = 2. # rotor diameters apart
    
    n = np.size(loc)/2
    x = loc[0:n]
    y = loc[n:]
    sep = np.zeros((n-1)*n/2)

    k = 0
    for i in range(0, n):
        for j in range(i+1, n):
            sep[k] = sqrt((x[j]-x[i])**2+(y[j]-y[i])**2)
            k += 1

    return sep - space*turb_dia


## Main
if __name__ == "__main__":
    
    # Option to plot the velocity field
    contour = True
    contour = False
    
    global dia
    global rot
    global chord
    global twist
    global delta
    global B
    global H
    global rho
    global mu
    global velf

    global af_data
    global cl_data
    global cd_data
    global coef0
    global coef1
    global coef2
    global coef3
    global coef4
    global coef5
    global coef6
    global coef7
    global coef8
    global coef9

    global funcs

    foildata = '/Users/ning1/Documents/FLOW Lab/VAWTAC/airfoils/NACA_0021.dat'

    af_data,cl_data,cd_data = vwm.airfoil_data(foildata)
    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()
    
    # define turbine size
    turb_dia = 6. # m
    turb_rot = 20. # rad/s
    
    # Scaling grid case
    nRows = 2     # number of rows and columns in grid
    spacing = 5     # turbine grid spacing in diameters
    
    # Set up position arrays
    points = np.linspace(spacing*turb_dia,nRows*spacing*turb_dia,nRows)
    xpoints, ypoints = np.meshgrid(points, points)
    xt = np.ndarray.flatten(xpoints)
    yt = np.ndarray.flatten(ypoints)
    dia = np.ones_like(xt)*turb_dia
    rot = np.ones_like(xt)*turb_rot

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.25
    H = 5.

    rho = 1.225
    mu = 1.7894e-5
    velf = 15.

    ntheta = 36
    velx,vely = induced_vel(dia[0]/2.,af_data,cl_data,cd_data,chord,twist,delta,B,rot[0],velf,rho,mu,ntheta)
    power_iso,_ = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([turb_dia]),np.array([turb_rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely)
    
    x0 = xt
    y0 = yt
    area = 2
    xlow = points[0]-spacing*area
    xupp = points[-1]+spacing*area
    ylow = points[0]-spacing*area
    yupp = points[-1]+spacing*area
    
    optProb = Optimization('VAWT_Power', obj_func)
    optProb.addObj('obj')
    
    n = np.size(x0)
    optProb.addVarGroup('xvars', n, 'c', lower=xlow, upper=xupp, value=x0)
    optProb.addVarGroup('yvars', n, 'c', lower=ylow, upper=yupp, value=y0)
    
    num_cons_sep = (n-1)*n/2
    optProb.addConGroup('sep', num_cons_sep, lower=0, upper=None)
    
    opt = SNOPT()
    opt.setOption('Scale option',0)
    res = opt(optProb, sens=None)
    print res
    
    pow = np.array(-1*res.fStar)
    xf = res.xStar['xvars']
    yf = res.xStar['yvars']
    power_turb = funcs['power_turb']
    
    print 'The power is:',pow,'W'
    print 'The isolated power is:',power_iso*np.size(xt),'W'
    print 'The x-locations:',xf
    print 'The y-locations:',yf
    print 'The power of each turbine (W):',power_turb
    print 'The isolated power of one turbine is:',power_iso,'W'
    
    
    plt.figure(1)
    for i in range(np.size(xt)):
        circ = plt.Circle((xt[i],yt[i]),dia[i]/2.,color='b',fill=True)
        plt.gca().add_patch(circ)
    for i in range(np.size(xt)):
        circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='r',fill=True)
        plt.gca().add_patch(circ)
    for i in range(np.size(xt)):
        plt.plot([xt[i], xf[i]], [yt[i], yf[i]], '--k')
    plt.xlim(points[0]-spacing*area,points[-1]+spacing*area)
    plt.ylim(points[0]-spacing*area,points[-1]+spacing*area)
    
    
    if contour == True:
    
        N = 10
        xplot = np.linspace(points[0]-spacing*area,points[-1]+spacing*area,N)
        yplot = np.linspace(points[0]-spacing*area,points[-1]+spacing*area,N)
        [X,Y] = np.meshgrid(xplot,yplot)
        Vel = np.zeros((N,N))
        
        k = 0
        for i in range(N):
            for j in range(N):
                Velx,Vely = _vawtwake.overlappoint(xt,yt,dia,rot,chord,B,X[i,j],Y[i,j],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,220,200,1)
                Vel[i,j] = sqrt((Velx+velf)**2 + (Vely)**2)
                k += 1
                print k
        
        plt.figure(2)
        lb = 0.0 # lower bound on velocity to display
        ub = velf # upper bound on velocity to display
        ran = 200 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,4) # setting the number of tick marks on colorbar
        CS = plt.contourf(X,Y,Vel,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.jet) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        for i in range(np.size(xt)):
            circ = plt.Circle((xt[i],yt[i]),dia[i]/2.,color='g',edgecolor='k',fill=True)
            plt.gca().add_patch(circ)
        
        plt.xlim(points[0]-spacing,points[-1]+spacing)
        plt.ylim(points[0]-spacing,points[-1]+spacing)
    
    plt.show()
