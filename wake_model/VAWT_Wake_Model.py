"""
Parameterized VAWT Wake Model using CFD vorticity data
Developed by Eric Tingey at Brigham Young University

This code models the wake behind a vertical-axis wind turbine based on
parameters like tip-speed ratio, solidity and wind speed by converting the
vorticity of the wake into velocity information. The model uses CFD data
obtained from STAR-CCM+ of simulated turbines to make the wake model as
accurate as possible.

Only valid for tip-speed ratios between 1.5 and 7.0 and solidities between
0.15 and 1.0. Reynolds numbers should also be around the range of 200,000 to
6,000,000.

In this code, up and down are sides of the wake according to:

--------------->-------------------------------------------------
--------------->------<-ROTATION-<------------UP-----------------
--------------->---------=====--------#################----------
--------------->------//       \\#############################---
-FREE-STREAM--->-----|| TURBINE ||########## WAKE ###############
----WIND------->-----||         ||###############################
--------------->------\\       //#############################---
--------------->---------=====--------#################----------
--------------->------>-ROTATION->-----------DOWN----------------
--------------->-------------------------------------------------
"""


import csv
from os import path
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from scipy.integrate import dblquad
from scipy.interpolate import RectBivariateSpline

import _vortmodel

# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'


def vorticity(tsr,solidity):
    """
    Using EMG distribution parameters to define the vorticity strength and shape
    
    Parameters
    ----------
    tsr : float
        tip-speed ratio
    solidity : float
        turbine solidity
    
    Returns
    ----------
    loc : array
        array of the location parameter (3 values)
    spr : array
        array of the spread parameter (2 values)
    skw : array
        array of the skew parameter (2 values)
    scl : array
        array of the scale parameter (3 values)
    """
    
    # Reading in csv file (vorticity database)
    basepath = path.join(path.dirname(path.realpath(__file__)),'data')
    fdata = basepath + path.sep + 'vortdatabase.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw,0)
            vortdat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d,float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            vortdat = np.vstack([vortdat,raw]) # adding entry to vorticity database array
        i += 1
    f.close()
    
    vortdat = np.delete(vortdat,(0),axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    vortdat = vortdat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_loc1 = vortdat[i*10]\ns'+sol+'_loc2 = vortdat[i*10+1]\ns'+sol+'_loc3 = vortdat[i*10+2]\ns'+sol+'_spr1 = vortdat[i*10+3]\ns'+sol+'_spr2 = vortdat[i*10+4]\ns'+sol+'_skw1 = vortdat[i*10+5]\ns'+sol+'_skw2 = vortdat[i*10+6]\ns'+sol+'_scl1 = vortdat[i*10+7]\ns'+sol+'_scl2 = vortdat[i*10+8]\ns'+sol+'_scl3 = vortdat[i*10+9]\n')
    
    # BIVARIATE SPLINE FITTING
    
    iz = np.size(sol_d)
    jz = np.size(tsr_d)
    
    # Initializing rectangular matrices
    Z_loc1 = np.zeros((iz,jz))
    Z_loc2 = np.zeros((iz,jz))
    Z_loc3 = np.zeros((iz,jz))
    Z_spr1 = np.zeros((iz,jz))
    Z_spr2 = np.zeros((iz,jz))
    Z_skw1 = np.zeros((iz,jz))
    Z_skw2 = np.zeros((iz,jz))
    Z_scl1 = np.zeros((iz,jz))
    Z_scl2 = np.zeros((iz,jz))
    Z_scl3 = np.zeros((iz,jz))
    
    # Transferring raw data into rectangular matrices
    for i in range(iz):
        for j in range(jz):
            sol = str(i+1)
            exec('Z_loc1[i,j] = s'+sol+'_loc1[j]')
            exec('Z_loc2[i,j] = s'+sol+'_loc2[j]')
            exec('Z_loc3[i,j] = s'+sol+'_loc3[j]')
            exec('Z_spr1[i,j] = s'+sol+'_spr1[j]')
            exec('Z_spr2[i,j] = s'+sol+'_spr2[j]')
            exec('Z_skw1[i,j] = s'+sol+'_skw1[j]')
            exec('Z_skw2[i,j] = s'+sol+'_skw2[j]')
            exec('Z_scl1[i,j] = s'+sol+'_scl1[j]')
            exec('Z_scl2[i,j] = s'+sol+'_scl2[j]')
            exec('Z_scl3[i,j] = s'+sol+'_scl3[j]')
    
    # Creating a rectangular bivariate spline of the parameter data
    s_loc1 = RectBivariateSpline(sol_d,tsr_d,Z_loc1)
    s_loc2 = RectBivariateSpline(sol_d,tsr_d,Z_loc2)
    s_loc3 = RectBivariateSpline(sol_d,tsr_d,Z_loc3)
    s_spr1 = RectBivariateSpline(sol_d,tsr_d,Z_spr1)
    s_spr2 = RectBivariateSpline(sol_d,tsr_d,Z_spr2)
    s_skw1 = RectBivariateSpline(sol_d,tsr_d,Z_skw1)
    s_skw2 = RectBivariateSpline(sol_d,tsr_d,Z_skw2)
    s_scl1 = RectBivariateSpline(sol_d,tsr_d,Z_scl1)
    s_scl2 = RectBivariateSpline(sol_d,tsr_d,Z_scl2)
    s_scl3 = RectBivariateSpline(sol_d,tsr_d,Z_scl3)
    
    # Selecting the specific parameters to use for TSR and solidity
    loc1 = s_loc1(solidity,tsr)
    loc2 = s_loc2(solidity,tsr)
    loc3 = s_loc3(solidity,tsr)
    spr1 = s_spr1(solidity,tsr)
    spr2 = s_spr2(solidity,tsr)
    skw1 = s_skw1(solidity,tsr)
    skw2 = s_skw2(solidity,tsr)
    scl1 = s_scl1(solidity,tsr)
    scl2 = s_scl2(solidity,tsr)
    scl3 = s_scl3(solidity,tsr)
    
    # Creating arrays of the parameters
    loc = np.array([loc1[0,0],loc2[0,0],loc3[0,0]])
    spr = np.array([spr1[0,0],spr2[0,0]])
    skw = np.array([skw1[0,0],skw2[0,0]])
    scl = np.array([scl1[0,0],scl2[0,0],scl3[0,0]])
    
    return loc,spr,skw,scl


def velocity_field(xt,yt,x0,y0,velf,dia,tsr,solidity):
    """
    Calculating normalized velocity from the vorticity data at (x0,y0)
    
    Parameters
    ----------
    xt : float
        downstream position of turbine (m)
    yt : float
        lateral position of turbine (m)
    x0 : float
        downstream distance from turbine in flow domain (m)
    y0 : float
        lateral distance from turbine in flow domation (m)
    velf : float
        free stream velocity (m/s)
    dia : float
        turbine diameter (m)
    tsr : float
        tip-speed ratio; [rotation rate (rad/s)]*[turbine radius (m)]/[free stream velocity (m/s)]
    solidity : float
        turbine solidity; [number of turbine blades]*[blade chord length (m)]/[turbine radius (m)]
    
    Returns
    ----------
    vel : float
        final normalized velocity at (x0,y0) with respect to the free stream velocity (m/s)
    """
    rad = dia/2.
    rot = tsr*velf/rad

    # Calculating EMG distribution parameters
    loc,spr,skw,scl = vorticity(tsr,solidity)
    
    # Translating the turbine position
    x0t = x0 - xt
    y0t = y0 - yt
    
    # Integration of the vorticity profile using Fortran code (vorticity.f90; _vortrun.so)
    vel_vs = dblquad(_vortmodel.integrand,0.,35.*dia,lambda x: -4.*dia,lambda x: 4.*dia, args=(x0t,y0t,dia,loc[0],loc[1],loc[2],spr[0],spr[1],skw[0],skw[1],scl[0],scl[1],scl[2]))
    
    # Calculating velocity deficit
    vel = (vel_vs[0]*(rot))/(2.*pi)
    vel = (vel + velf)/velf # normalization of velocity
    
    return vel

