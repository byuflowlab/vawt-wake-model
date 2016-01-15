import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt
from scipy.integrate import quad
from VAWT_Wake_Model import velocity_field


def Cr(y0,xt,yt,x0,velf,dia,tsr,solidity):
    return 1. - velocity_field(xt,yt,x0,y0,velf,dia,tsr,solidity)


def overlap(xt,yt,x0,y0,velf,dia,tsr,solidity):
    """
    Calculating an effective velocity based on overlapping of wakes
    
    Parameters
    ----------
    xt : array
        downstream locations of turbines producing wakes (m)
    yt : array
        lateral locations of turbines producing wakes (m)
    x0 : float
        downstream location of turbine being analyzed (m)
    y0 : float
        lateral location of turbine being analyzed (m)
    velf : float
        free stream velocity (m/s)
    dia : array
        turbine diameters of each turbine (m)
    tsr : array
        tip-speed ratios of each turbine
    solidity : array
        solidities of each turbine
    
    Returns
    ----------
    veleff : float
        effective velocity in front of turbine (m/s)
    """
    
    N = np.size(xt) # number of turbines that overlap given turbine
    inte = 0. # initializing integration value
    r1 = y0 - dia/2. # lower bound of integration
    r2 = y0 + dia/2. # upper bound of integration
    
    for i in range(N):
        inte_n = quad(Cr,r1,r2,args=(xt[i],yt[i],x0,velf,dia[i],tsr[i],solidity[i]))
        inte = inte + (inte_n[0])**2
    
    veleff = velf*(1. - sqrt(inte))
    
    return veleff
    

def powerval(Cp,dens,vel,dia):
    """
    Calculating an effective velocity based on overlapping of wakes
    
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
    x0 = 12.
    y0 = 6.
    velf = 15.
    dia = np.array([6.])
    tsr = np.array([3.0])
    solidity = np.array([0.25])
    
    Cp = 0.4
    dens = 1.225
    
    veleff = overlap(xt,yt,x0,y0,velf,dia,tsr,solidity)
    power_iso = powerval(Cp,dens,velf,dia[0])
    power = powerval(Cp,dens,veleff,dia[0])
    
    print 'Isolated turbine power: ',power_iso,'kJ'
    print 'Analyzed turbine power: ',power,'kJ'
    
    
    