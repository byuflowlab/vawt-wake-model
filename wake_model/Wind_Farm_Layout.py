import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt
from VAWT_Wake_Model import velocity_field


def overlap(xt,yt,diat,tsr,solidity,x0,y0,dia,velf,hub,out):
    """
    Calculating an effective velocity based on overlapping of wakes
    
    Parameters
    ----------
    xt : array
        downstream locations of turbines producing wakes (m)
    yt : array
        lateral locations of turbines producing wakes (m)
    diat : array
        diameters of each of the turbines producing wakes (m)
    tsr : array
        tip-speed ratio of each turbine
    solidity : array
        solidities of each turbine
    x0 : float
        downstream location of turbine being analyzed (m)
    y0 : float
        lateral location of turbine being analyzed (m)
    dia : float
        diameter of turbine being analyzed (m)
    velf : float
        free stream velocity (m/s)
    hub : bool
        option of using the hub velocity (True) or an average of
        velocities across the turbine's rotation (False)
    out : bool
        option of outputting the the information of turbines completed
        in the overlap calculation (True)
    
    Returns
    ----------
    veleff : float
        effective velocity in front of turbine (m/s)
    """
    
    N = np.size(xt) # number of turbines that overlap given turbine
    inte = 0. # initializing integration value
    r1 = y0 - dia/2. # lower bound of integration
    r2 = y0 + dia/2. # upper bound of integration
    n = 3
    if hub == True:
        y = np.array([y0]) # using hub velocity for effective velocity
    else:
        y = np.linspace(r1,r2,n) # using n velocity values for effective velocity
    
    vel = np.zeros_like(y)
    
    for i in range(N):
        for j in range(np.size(y)):
            vel[j] = 1. - velocity_field(xt[i],yt[i],x0,y[j],velf,diat[i],tsr[i],solidity[i])
        
        inte_n = np.average(vel)
        
        inte = inte + (inte_n)**2
        if out == True:
            print 'Turbine '+str(i+1)+' of '+str(N)+' completed'
    
    veleff = velf*(1. - sqrt(inte))
    
    return veleff
    

def powerval(Cp,dens,vel,dia):
    """
    Calculating turbine power of a 2D VAWT
    
    Parameters
    ----------
    Cp : float
        turbine power coefficient
    dens : float
        air density of wind farm (kg/m^3)
    vel : float
        effective velocity in front of turbine (m/s)
    dia : float
        turbine diameter (m)
    
    Returns
    ----------
    power : float
        turbine power (kJ)
    """
    
    power = 0.5*Cp*dens*dia*vel**3
    
    return power/1e3


if __name__ == "__main__":
    
    xt = np.array([0.])
    yt = np.array([0.])
    diat = np.array([6.])
    tsr = np.array([4.0])
    solidity = np.array([0.25])
    x0 = 11.
    y0 = 0.
    dia = 6.
    velf = 15.
    
    Cp = 0.4
    dens = 1.225
    
    power_iso = powerval(Cp,dens,velf,dia)
    print 'Isolated turbine power: ',power_iso,'kJ'
    
    veleff = overlap(xt,yt,diat,tsr,solidity,x0,y0,dia,velf,True,False)
    
    power = powerval(Cp,dens,veleff,dia)
    print 'Analyzed turbine power: ',power,'kJ'
    
    N = 10
    xplot = np.linspace(-4.*dia,4.*dia,N)
    yplot = np.linspace(-4.*dia,4.*dia,N)
    [X,Y] = np.meshgrid(xplot,yplot)
    P = np.zeros((N,N))
    
    k = 0
    for i in range(N):
        for j in range(N):
            veleffp = overlap(xt,yt,diat,tsr,solidity,X[i,j],Y[i,j],dia,velf,True,False)
            P[i,j] = powerval(Cp,dens,veleffp,dia)/power_iso
            k += 1
            print k
    
    plt.figure()
    lb = 0.0 # lower bound on velocity to display
    ub = 1.2 # upper bound on velocity to display
    ran = 200 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,7) # setting the number of tick marks on colorbar
    CS = plt.contourf(X,Y,P,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.parula) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    for i in range(np.size(xt)):
        circ = plt.Circle((xt[i],yt[i]),dia/2.,color='w',fill=True)
        plt.gca().add_patch(circ)
    plt.show()
    
    
    