from pyoptsparse import Optimization, SNOPT, pyOpt_solution
import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt

from Wind_Farm_Layout import overlap,powerval


def obj_func(xdict):
    global dia
    global tsr
    global solidity
    global Cp
    global dens
    global velf
    global funcs
    
    Vars = xdict['xvars']
    x = xdict['xvars']
    y = xdict['yvars']
    funcs = {}
    
    veleff = np.zeros_like(x)
    power_turb = np.zeros_like(x)
    
    for i in range(np.size(x)):
        xp = np.delete(x,i)
        yp = np.delete(y,i)
        diap = np.delete(dia,i)
        tsrp = np.delete(tsr,i)
        solidityp = np.delete(solidity,i)
        veleff[i] = overlap(xp,yp,diap,tsrp,solidityp,x[i],y[i],dia[i],velf,True,False)
        power_turb[i] = powerval(Cp,dens,veleff[i],dia[i])
        # print i+1
    
    power = np.sum(power_turb)
    print power 
    funcs['obj'] = (-1.0*power)
    
    sep = sep_func(np.append(x,y))
    funcs['sep'] = sep
    funcs['veleff'] = veleff
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
    global tsr
    global solidity
    global Cp
    global dens
    global velf
    global turb_dia
    global funcs
    
    # define turbine size
    turb_dia = 6. # m
    
    # Scaling grid case
    nRows = 3     # number of rows and columns in grid
    spacing = 5     # turbine grid spacing in diameters
    
    # Set up position arrays
    points = np.linspace(spacing*turb_dia,nRows*spacing*turb_dia,nRows)
    xpoints, ypoints = np.meshgrid(points, points)
    xt = np.ndarray.flatten(xpoints)
    yt = np.ndarray.flatten(ypoints)
    dia = np.ones_like(xt)*turb_dia
    tsr = np.ones_like(xt)*4.0
    solidity = np.ones_like(xt)*0.25
    
    Cp = 0.4
    dens = 1.225
    velf = 15.
    
    power_iso = powerval(Cp,dens,velf,turb_dia)
    
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
    veleff = funcs['veleff']
    
    print 'The power is:',pow,'kJ'
    print 'The isolated power is:',power_iso*np.size(xt),'kJ'
    print 'The x-locations:',xf
    print 'The y-locations:',yf
    print 'The effective wind speeds:',veleff
    print 'The power of each turbine (kJ):',power_turb
    
    
    plt.figure(1)
    for i in range(np.size(xt)):
            circ = plt.Circle((xt[i],yt[i]),dia[i]/2.,color='b',fill=True)
            plt.gca().add_patch(circ)
    for i in range(np.size(xt)):
            circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='r',fill=True)
            plt.gca().add_patch(circ)
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
                Vel[i,j] = overlap(xf,yf,dia,tsr,solidity,X[i,j],Y[i,j],turb_dia,velf,True,False)
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
