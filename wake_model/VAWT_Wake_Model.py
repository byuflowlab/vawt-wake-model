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

In this code, the x and y coordinates are made according to:

--------------->--------------------------------------------------------
--------------->--------------------------------------------------------
--------------->---------=====--------#################-----------Y-----
--------------->------//       \\#############################----|-----
-FREE-STREAM--->-----|| TURBINE ||########## WAKE ###############-|___X-
----WIND------->-----||         ||###############################-------
--------------->------\\       //#############################----------
--------------->---------=====--------#################-----------------
--------------->--------------------------------------------------------
--------------->--------------------------------------------------------

The imported vorticity data also assumes symmetry in the wake and therefore
rotation direction is irrelevant. The imported velocity data assumes the
rotation is counter-clockwise from the top.

"""


import csv
from os import path
import numpy as np
import matplotlib.pyplot as plt
from numpy import pi,exp,fabs,sqrt
from scipy.integrate import dblquad

import _vortmodel

# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'


def velocity_field(xt,yt,x0,y0,velf,dia,tsr,solidity,cfd_data,param):
    """
    Calculating normalized velocity from the vorticity data at (x0,y0) in global flow domain

    Parameters
    ----------
    xt : float
        downstream position of turbine domain (m)
    yt : float
        lateral position of turbine in flow domain (m)
    x0 : float
        downstream position in flow domain to be calculated (m)
    y0 : float
        lateral position in flow domain to be calculated (m)
    velf : float
        free stream velocity (m/s)
    dia : float
        turbine diameter (m)
    tsr : float
        tip-speed ratio; [rotation rate (rad/s)]*[turbine radius (m)]/[free stream velocity (m/s)]
    solidity : float
        turbine solidity; [number of turbine blades]*[blade chord length (m)]/[turbine radius (m)]
    cfd_data : string
        specifying to use CFD vorticity ('vort') or velocity ('velo') for the basis of the wake model
    param : array
        the coefficients used for the EMG or SMG distributions

    Returns
    ----------
    vel : float
        final normalized velocity at (x0,y0) with respect to the free stream velocity (m/s)
    """
    rad = dia/2.
    rot = tsr*velf/rad

    # Translating the turbine position
    x0t = x0 - xt
    y0t = y0 - yt

    if cfd_data == 'vort':
        # Calculating EMG distribution parameters
        loc = param[0]
        spr = param[1]
        skw = param[2]
        scl = param[3]

        # Integration of the vorticity profile using Fortran code (vorticity.f90)
        vel_vs = dblquad(_vortmodel.integrand,0.,35.*dia,lambda x: -4.*dia,lambda x: 4.*dia, args=(x0t,y0t,dia,loc[0],loc[1],loc[2],spr[0],spr[1],skw[0],skw[1],scl[0],scl[1],scl[2]))

        # Calculating velocity deficit
        vel = (vel_vs[0]*(rot))/(2.*pi)
        vel = (vel + velf)/velf # normalization of velocity

    elif cfd_data == 'velo':
        # Normalizing the downstream and lateral positions by the turbine diameter
        x0d = x0/dia
        y0d = y0/dia

        # Calculating SMG distribution parameters
        men = param[0]
        spr = param[1]
        scl = param[2]
        rat = param[3]
        tns = param[4]

        # men = np.array( [-0.0007448786610163438, 0.011700465818493566, -0.005332505770684337] )
        # spr = np.array( [6.462355161093711, 7.079901300173991, 12.102886237210939] )
        # scl = np.array( [8.509272717226171, 7.023483471068396, 27.707846411384697] )
        # rat = np.array( [-2.107186196351149, 44.93845180541949] )
        # tns = np.array( [-1.4660542829002265, 30.936653231840257] )
        #
        # men = np.array( [-0.00059737414699399, 0.009890587474506057, -0.0016721254639608882] )
        # spr = np.array( [-0.005652314031253564, 0.06923002880544946, 0.526304136118912] )
        # scl = np.array( [6.639808894608728, 5.477607580787858, 21.13678312202297] )
        # rat = np.array( [-2.0794010451530873, 44.798557035611] )
        # tns = np.array( [-1.43164706657537, 30.761785195818447] )


        men_v = men[0]*x0d**2 + men[1]*x0d + men[2]
        if men_v > 0.5:
            men_v = 0.5
        elif men_v < -0.5:
            men_v = -0.5

        # spr_v = spr[2]/(spr[1]*sqrt(2.*pi))*exp(-(x0d-spr[0])**2/(2.*spr[1]**2))
        spr_v = spr[0]*x0d**2 + spr[1]*x0d + spr[2]
        if spr_v < 0.35:
            spr_v = 0.35
        elif spr_v > 4.:
            spr_v = 4.

        scl_v = scl[2]/(scl[1]*sqrt(2.*pi))*exp(-(x0d-scl[0])**2/(2.*scl[1]**2))

        rat_v = rat[0]*x0d + rat[1]
        if rat_v < 0.:
            rat_v = 0.

        tns_v = tns[0]*x0d + tns[1]
        if tns_v < 0.:
            tns_v = 0.

        vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(y0d+men_v)**2/(2.*spr_v**2)))*(1./(1 + exp(rat_v*fabs(y0d)-tns_v))) + 1. # Normal distribution with sigmoid weighting

        if x0 < xt:
            vel = 1. # Velocity is free stream in front and to the sides of the turbine

    return vel
