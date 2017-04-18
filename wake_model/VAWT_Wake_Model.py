"""
Parameterized VAWT Wake Model Python Code
Developed by Eric Tingey at Brigham Young University

This code models the wake behind a vertical-axis wind turbine based on
tip-speed ratio, solidity and wind speed by converting the vorticity of
the wake into velocity information. The use CFD data obtained from
STAR-CCM+ turbine simulations serve as the basis of the initial wake model.

Only valid for tip-speed ratios between 1.5 and 7.0 and solidities between
0.15 and 1.0. Reynolds numbers should also be around the range of 600,000 to
6,000,000.

In this code, the x and y coordinates (looking down on the turbine) are
made according to:

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
rotation direction is irrelevant.
"""
import numpy as np
from numpy import pi,fabs,sqrt,sin,cos,argmin
from scipy.integrate import _quadpack
from scipy.interpolate import UnivariateSpline
import csv
from os import path

from joblib import Parallel,delayed

import _vawtwake

##########################################################################################
# Double Integration Method using necessary Quadpack (SciPy) code (_qagse)
# Originally created by Travis Oliphant (2001) and Nathan Woods (2013) (nquad &c)
def _quad(func, a, b, args=(), full_output=0, epsabs=1.49e-8, epsrel=1.49e-8, limit=50):
    # Calling the _qagse code from Quadpack to perform a single integration

    retval = _quadpack._qagse(func,a,b,args,full_output,epsabs,epsrel,limit)

    return retval[:-1]

def _infunc(x,func,gfun,hfun,more_args):
    # Arranging a double integral into two single integrations

    a = gfun(x)
    b = hfun(x)
    myargs = (x,) + more_args

    return _quad(func,a,b,args=myargs)[0]

def _dblquad(func, a, b, gfun, hfun, args=(), epsabs=1.49e-8, epsrel=1.49e-8):
    # Performing a double integration using _infunc and _quad

    return _quad(_infunc, a, b, (func, gfun, hfun, args), epsabs=epsabs, epsrel=epsrel)
##########################################################################################


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


def coef_val():
    """
    The polynomial surface coefficients used for the EMG parameters
    Published coefficients from paper (4 Reynolds numbers; 1 2 2 1 1 1 2 1 2 2) may be used

    Parameters
    ----------
    no parameters

    Returns
    ----------
    loc1 : array
        the first location parameter coefficients
    loc2 : array
        the second location parameter coefficients
    loc3 : array
        the third location parameter coefficients
    spr1 : array
        the first spread parameter coefficients
    spr1 : array
        the second spread parameter coefficients
    skw1 : array
        the first skew parameter coefficients
    skw1 : array
        the second skew parameter coefficients
    scl1 : array
        the first scale parameter coefficients
    scl2 : array
        the second scale parameter coefficients
    scl3 : array
        the third scale parameter coefficients
    """

    basepath = path.join(path.dirname(path.realpath('__file__')), 'data')
    fdata = basepath + path.sep + 'VAWTPolySurfaceCoef_pub.csv' # published coefficients from paper
    # fdata = basepath + path.sep + 'VAWTPolySurfaceCoef.csv' # polynomial surface fitting coefficients

    loc1 = np.zeros(10)
    loc2 = np.zeros(10)
    loc3 = np.zeros(10)
    spr1 = np.zeros(10)
    spr2 = np.zeros(10)
    skw1 = np.zeros(10)
    skw2 = np.zeros(10)
    scl1 = np.zeros(10)
    scl2 = np.zeros(10)
    scl3 = np.zeros(10)

    f = open(fdata)
    csv_f = csv.reader(f)

    i = 0
    for row in csv_f:
        if i != 0:
            loc1[i-1] = float(row[0])
            loc2[i-1] = float(row[1])
            loc3[i-1] = float(row[2])
            spr1[i-1] = float(row[3])
            spr2[i-1] = float(row[4])
            skw1[i-1] = float(row[5])
            skw2[i-1] = float(row[6])
            scl1[i-1] = float(row[7])
            scl2[i-1] = float(row[8])
            scl3[i-1] = float(row[9])
        i += 1

    f.close()

    return loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3


def airfoil_data(file):
    """
    Reading in an airfoil data file and fitting a spline to the data to smooth it out

    Parameters
    ----------
    file : string
        path to airfoil data file (.dat; typically from '/data/airfoils/')

    Returns
    ----------
    af_data_smooth : array
        the smoothed airfoil angles of attack
    cl_data_smooth : array
        the smoothed airfoil coefficients of lift
    cd_data_smooth : array
        the smoothed airfoil coefficients of drag
    """
    af_data = np.array([])
    cl_data = np.array([])
    cd_data = np.array([])

    f = open(file, 'r')
    for i in range(13):
        f.readline() # skipping preliminary lines
    for line in f:
        line = line.strip()
        columns = line.split()
        if columns[0] == 'EOT':
            break
        else:
            af_data = np.append(af_data,float(columns[0]))
            cl_data = np.append(cl_data,float(columns[1]))
            cd_data = np.append(cd_data,float(columns[2]))
    f.close()

    # Smoothing data with a univariate spline ('s' is the smoothing factor)
    clsmooth = UnivariateSpline(af_data,cl_data,s=0.1)
    cdsmooth = UnivariateSpline(af_data,cd_data,s=0.001)

    # Reassigning imported airfoil data with the smoothed airfoil data
    af_data_smooth = np.linspace(af_data[0],af_data[-1],4000)
    cl_data_smooth = clsmooth(af_data_smooth)
    cd_data_smooth = cdsmooth(af_data_smooth)

    return af_data_smooth,cl_data_smooth,cd_data_smooth


def velocity_field(xt,yt,x0,y0,velf,dia,rot,chord,B,param=None,veltype='all',integration='simp',m=220,n=200):
    """
    Calculating normalized velocity from the vorticity data at (x0,y0) in global flow domain

    Parameters
    ----------
    xt : float
        downstream position of surrounding turbine in flow domain (m)
    yt : float
        lateral position of surrounding turbine in flow domain (m)
    x0 : float
        downstream position in flow domain to be calculated (m)
    y0 : float
        lateral position in flow domain to be calculated (m)
    velf : float
        free stream velocity (m/s)
    dia : float
        turbine diameter (m)
    rot : float
        turbine rotation rate (rad/s)
    param : array
        the coefficients used for the EMG distributions ('None' will provide the coefficients using VAWTPolySurfaceCoef.csv)
        param should be an array of length 10 with each of the EMG parameters corresponding to loc, spr, skw, and scl
    veltype : string
        the type of velocity to calculate ('all': velocity magnitude, 'x': x-induced velocity, 'y': y-induced velocity,
        'ind': vector of both x- and y-induced velocities without free stream, 'vort': vorticity profile neglecting integration)
    integration : string
        the type of integration method used ('simp': Simpson's Rule, 'gskr': 21 Point Gauss-Kronrod Rule)
    m : int
        the number of downstream divisions requested for Simpson's Rule (must be divisible by 2); neglected otherwise
    n : int
        the number of downstream divisions requested for Simpson's Rule (must be divisible by 2); neglected otherwise

    Returns
    ----------
    vel : float
        final normalized velocity at (x0,y0) with respect to the free stream velocity (m/s)
    """
    rad = dia/2.
    tsr = rad*fabs(rot)/velf
    solidity = (chord*B)/rad

    # Translating the turbine position
    x0t = x0 - xt
    y0t = y0 - yt

    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = coef_val()

    # Calculating EMG distribution parameters (based on polynomial surface fitting)
    if param is None:
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

    else:
        # Reading in EMG distribution parameters
        loc1 = param[0]
        loc2 = param[1]
        loc3 = param[2]
        spr1 = param[3]
        spr2 = param[4]
        skw1 = param[5]
        skw2 = param[6]
        scl1 = param[7]
        scl2 = param[8]
        scl3 = param[9]

    ###################################
    if veltype == 'vort':
        # VORTICITY CALCULATION (NO INTEGRATION)
        if x0t < 0.:
            vel = 0.
        else:
            vel = _vawtwake.vorticitystrength(x0t,y0t,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3)/rot
    ###################################
    else:
        # Integration of the vorticity profile to calculate velocity
        if integration == 'simp':
            # SIMPSON'S RULE INTEGRATION (must use polynomial surface coefficients from VAWTPolySurfaceCoef.csv)
            inte = 1 # Simpson's Rule
            # inte = 2 # Trapezoidal Rule (optional ability of the code-- faster but less accurate)

            if param is not None:
                print "**** Using polynomial surface coefficients from VAWTPolySurfaceCoef.csv for Simpson's rule integration ****"

            vel_xs,vel_ys = _vawtwake.vel_field(xt,yt,x0,y0,dia,rot,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,inte)

            if veltype == 'all':
                vel = sqrt((vel_xs*velf + velf)**2 + (vel_ys*velf)**2)/velf
            elif veltype == 'x':
                vel = (vel_xs*velf + velf)/velf
            elif veltype == 'y':
                vel = vel_ys
            elif veltype == 'ind':
                vel = np.array([vel_xs,vel_ys])
    ###################################
        elif integration == 'gskr':
            # 21-POINT GAUSS-KRONROD RULE QUADRATURE INTEGRATION
            xbound = (scl3+5.)*dia
            argval = (x0t,y0t,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3)
            if veltype == 'all' or veltype == 'x' or veltype == 'ind':
                vel_x = _dblquad(_vawtwake.integrandx,0.,xbound,lambda x: -1.*dia,lambda x: 1.*dia,args=argval)
                vel_xs = (vel_x[0]*fabs(rot))/(2.*pi)
            if veltype == 'all' or veltype == 'y' or veltype == 'ind':
                vel_y = _dblquad(_vawtwake.integrandy,0.,xbound,lambda x: -1.*dia,lambda x: 1.*dia,args=argval)
                vel_ys = (vel_y[0]*fabs(rot))/(2.*pi)

            if veltype == 'all':
                vel = sqrt((vel_xs + velf)**2 + (vel_ys)**2)/velf
            elif veltype == 'x':
                vel = (vel_xs + velf)/velf
            elif veltype == 'y':
                vel = vel_ys/velf
            elif veltype == 'ind':
                vel = np.array([vel_xs,vel_ys])/velf
    ###################################

    return vel

def overlap(p,xt,yt,diat,rott,chord,B,x0,y0,dia,velf,pointcalc,param=None,veltype='ind',integration='gskr'):
    """
    Calculating wake velocities around a turbine based on wake overlap from surrounding turbines
    (using the 21-point Gauss-Kronrod rule quadrature integration; Simpson's rule integration can be used via VAWT_Wake_Model.f90)

    Parameters
    ----------
    p : int
        number of points to calculate the velocity around a turbine (typically 36)
    xt : array
        downstream positions of surrounding turbine(s) in flow domain (m)
    yt : array
        lateral position of surrounding turbine(s) in flow domain (m)
    diat : array
        diameters of surrounding turbines (m)
    rott : array
        rotation rates of surrounding turbines (rad/s)
    chord : float
        chord length of the turbines (m)
    B : int
        number of turbine blades
    x0 : float
        downstream position in flow domain of turbine to be calculated (m)
    y0 : float
        lateral position in flow domain of turbine to be calculated (m)
    dia : float
        diameter of turbine to be calculated (m)
    velf : float
        free stream velocity (m/s)
    pointcalc : bool
        calculate the overlap at a point (True) or at p points around the blade flight path (False)
    param : array
        the coefficients used for the EMG distributions ('None' will provide the published coefficients automatically)
    veltype : string
        the type of velocity to calculate ('all': velocity magnitude, 'x': x-induced velocity, 'y': y-induced velocity,
        'ind': vector of both x- and y-induced velocities without free stream, 'vort': vorticity profile neglecting integration)
    integration : string
        the type of integration method used ('simp': Simpson's Rule, 'gskr': 21 Point Gauss-Kronrod Rule)
    m : int
        the number of downstream divisions requested for Simpson's Rule (must be divisible by 2); neglected otherwise
    n : int
        the number of downstream divisions requested for Simpson's Rule (must be divisible by 2); neglected otherwise

    Returns
    ----------
    velx : array
        final induced x-velocity at each point around the turbine being calculated (m/s)
    vely : array
        final induced y-velocity at each point around the turbine being calculated (m/s)
    """
    # initializing local variables and arrays
    t = np.size(xt) # number of turbines
    xd = np.zeros(p)
    yd = np.zeros(p)
    velx = np.zeros(p)
    vely = np.zeros(p)
    velx_int = np.zeros(p)
    vely_int = np.zeros(p)

    # Use parallelization (with joblib)
    parallel = True
    # parallel = False

    # finding points around the flight path of the blades
    for i in range(p):
        if pointcalc == False:
            theta = (2.0*pi/p)*i-(2.0*pi/p)/2.0
            xd[i] = x0 - sin(theta)*(dia/2.0)
            yd[i] = y0 + cos(theta)*(dia/2.0)
        elif pointcalc == True:
            xd[0] = x0
            yd[0] = y0
    intex = np.zeros(p)
    intey = np.zeros(p)

    if (t == 1): # coupled configuration (only two VAWTs)
        if pointcalc == False:
            if parallel == True:
                wake = Parallel(n_jobs=-1)(delayed(velocity_field)(xt[0],yt[0],xd[j],yd[j],velf,diat[0],rott[0],chord,B,param,veltype,integration) for j in range(p) )
                for i in range(p):
                    velx[i] = wake[i][0]*velf
                    vely[i] = wake[i][1]*velf
            elif parallel == False:
                for j in range(p):
                    wake = velocity_field(xt[0],yt[0],xd[j],yd[j],velf,diat[0],rott[0],chord,B,param,veltype,integration)
                    velx[j] = wake[0]*velf
                    vely[j] = wake[1]*velf
        elif pointcalc == True:
            wake = velocity_field(xt[0],yt[0],xd[0],yd[0],velf,diat[0],rott[0],chord,B,param,veltype,integration)
            velx[0] = wake[0]*velf
            vely[0] = wake[1]*velf

    else: # multiple turbine wake overlap
        if pointcalc == False:
            if parallel == True:
                wake = Parallel(n_jobs=-1)(delayed(velocity_field)(xt[w],yt[w],xd[q],yd[q],velf,diat[w],rott[w],chord,B,param,veltype,integration) for w in range(t) for q in range(p) )
            for j in range(t):
                for k in range(p):
                    if parallel == True:
                        velx_int[k] = -wake[k+j*p][0]
                        vely_int[k] = wake[k+j*p][1]
                    elif parallel == False:
                        wake = velocity_field(xt[j],yt[j],xd[k],yd[k],velf,diat[j],rott[j],chord,B,param,veltype,integration)
                        velx_int[k] = -wake[0]
                        vely_int[k] = wake[1]

                    # sum of squares of velocity deficits
                    if (velx_int[k] >= 0.0):
                        intex[k] = intex[k] + (velx_int[k])**2
                    else:
                        intex[k] = intex[k] - (velx_int[k])**2

                    if (vely_int[k] >= 0.0):
                        intey[k] = intey[k] + (vely_int[k])**2
                    else:
                        intey[k] = intey[k] - (vely_int[k])**2
        elif pointcalc == True:
            for j in range(t):
                wake = velocity_field(xt[j],yt[j],xd[0],yd[0],velf,diat[j],rott[j],chord,B,param,veltype,integration)
                velx_int[0] = -wake[0]
                vely_int[0] = wake[1]

                # sum of squares of velocity deficits
                if (velx_int[0] >= 0.0):
                    intex[0] = intex[0] + (velx_int[0])**2
                else:
                    intex[0] = intex[0] - (velx_int[0])**2

                if (vely_int[0] >= 0.0):
                    intey[0] = intey[0] + (vely_int[0])**2
                else:
                    intey[0] = intey[0] - (vely_int[0])**2

        # square root of sum of squares
        for l in range(p):
            if (intex[l] >= 0.0):
                velx[l] = -velf*(sqrt(intex[l]))
            else:
                velx[l] = velf*(sqrt(fabs(intex[l])))

            if (intey[l] >= 0.0):
                vely[l] = velf*(sqrt(intey[l]))
            else:
                vely[l] = -velf*(sqrt(fabs(intey[l])))

    return velx,vely


def wake_order(x,y,dia,xt,yt,diat,rott):
    """
    Determining the turbine wakes to include in wake overlap calculation

    Parameters
    ----------
    x : float
        downstream position of given turbine (m)
    y : float
        lateral position of given turbine (m)
    dia : float
        diameter of given turbine (m)
    xt : array
        downstream positions of surrounding turbines (m)
    yt : array
        lateral positions of surrounding turbines (m)
    diat : array
        diameters of surrounding turbines (m)
    rott : array
        rotation rates of surrounding turbines (rad/s)
    pen1 : float
        penalty for downstream direction outside of acceptable boundaries
    pen2 : float
        penalty for lateral direction outside of acceptable boundaries

    Returns
    ----------
    xt/xo : array
        downstream positions of selected surrounding turbines (m)
    yt/yo : array
        lateral positions of selected surrounding turbines (m)
    diat/diao : array
        diameters of selected surrounding turbines (m)
    rott/roto : array
        rotation rates of selected surrounding turbines (rad/s)
    """
    n = np.size(xt) # number of surrounding turbines
    keep = 1 # minimum number of turbines to consider in final selection

    if n <= keep: # use all of the surrounding turbines
        return xt,yt,diat,rott

    else:
        pen1 = 10000.
        pen2 = 10000.
        order = sqrt((x-xt)**2 + (y-yt)**2) # distance between given and surrounding turbines
        down = xt - x # downstream distance between given and surrounding turbines
        lat = fabs(yt - y) # lateral distance between given and surrounding turbines

        for i in range(n):
            if order[i] <= 6.*dia:
                keep += 1
            else:
                if down[i] >= 0.:
                    order[i] = order[i] + pen1
                if lat[i] > 1.5*dia:
                    order[i] = order[i] + pen2

        # setting up arrays
        xo = np.zeros(keep)
        yo = np.zeros(keep)
        diao = np.zeros(keep)
        roto = np.zeros(keep)
        for j in range(keep):
            val = argmin(order)
            xo[j] = xt[val]
            yo[j] = yt[val]
            diao[j] = diat[val]
            roto[j] = rott[val]

            order[val] = order[val] + 1e10

        return xo,yo,diao,roto
