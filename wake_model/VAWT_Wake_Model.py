"""
Parameterized VAWT Wake Model using CFD vorticity data
Developed by Eric Tingey at Brigham Young University

This code models the wake behind a vertical-axis wind turbine based on
parameters like tip-speed ratio, solidity and wind speed by converting the
vorticity of the wake into velocity information. The model uses CFD data
obtained from STAR-CCM+ of simulated turbines to make the wake model as
accurate as possible.

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
from numpy import pi, exp, fabs, sqrt, log
# from scipy.integrate import dblquad # integration with complete module
from Integrate import dblquad # integration with simplification file (Integrate.py)

import _vawtwake

def _parameterval(tsr,sol,coef):
    # Creating polynomial surface based on given coefficients and calculating the point at a given TSR and solidity

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

    return a + b*tsr + c*sol + d*tsr**2 + e*tsr*sol + f*sol**2 + g*tsr**3 + h*tsr**2*sol + i*tsr*sol**2 + j*sol**3


def velocity_field(xt,yt,x0,y0,velf,dia,rot,chord,B,param=None,veltype='all',integration='simp',m=220,n=200):
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
    rot : float
        turbine rotation rate (rad/s)
    solidity : float
        turbine solidity; [number of turbine blades]*[blade chord length (m)]/[turbine radius (m)]
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
    vel : float
        final normalized velocity at (x0,y0) with respect to the free stream velocity (m/s)
    """
    rad = dia/2.
    tsr = rad*rot/velf
    solidity = (chord*B)/rad

    # Translating the turbine position
    x0t = x0 - xt
    y0t = y0 - yt

    # Published coefficients (read 4: 1 2 2 1 1 1 2 1 2 2)
    coef0 = np.array( [0.0025703809856661534, -0.0007386258659065129, 0.004595508188667984, 0.000380123563204793, -0.0005090098755683027, 0.005744581813281894, -4.103393770815313e-05, -0.0014146918534486358, -0.013975958482495927, 0.0] )
    coef1 = np.array( [-0.5047504670963536, 0.23477391362058556, 0.8414256436198028, -0.04252528528617351, -0.06962875967504166, -0.6566907653208429, 0.002839318332370807, 0.00571803958194812, 0.0070744372783060295, 0.22805286438890995] )
    coef2 = np.array( [0.2878345841026334, 0.11512552658662782, 0.7303949879914625, -0.007035517839387948, -0.18284850673545897, -0.5241921153256568, -0.0003704899921255296, 0.010972527139685873, 0.04380801537377295, 0.1724129349605399] )
    coef3 = np.array( [0.08234816067475287, -0.03530687906626052, -0.3662863944976986, 0.003240141344532779, 0.12172015102204112, 0.2993048183466721, 0.0, -0.009253185586804007, -0.057469126406649716, -0.07257633583877886] )
    coef4 = np.array( [-0.07083579909945328, 0.016182024377569406, 0.1985436342461859, 0.0017738254727425816, -0.09111094817943823, -0.06561408122153217, -0.0005115133402638633, 0.009434288536679505, 0.022392136905926813, 0.0] )
    coef5 = np.array( [-1.6712830849073221, 1.5625053380692426, -6.180392756736983, -0.20407668040293722, -4.6476103643607685, 29.380064536220306, 0.0, 0.7502978877582536, -0.16358232641365608, -19.937609244085568] )
    coef6 = np.array( [-3.423561091777921, -9.228795430171687, 86.95722105482042, 2.772872601988039, -11.968168333741515, -150.61261090270446, -0.24715316589674527, 0.5283723108899993, 4.537286811245538, 82.50581844010263] )
    coef7 = np.array( [-0.19815381951708524, 0.08438758133540872, 1.2650146439483734, -0.007606115512168328, -0.2747023984740461, -0.8844640101378567, 0.0, 0.01870057580949183, 0.0699898278743648, 0.2794360008051127] )
    coef8 = np.array( [2.3932787625531815, -2.020874419612962, -8.938221963838357, 0.576323845480877, 2.8782448498416944, 16.598492450314534, -0.04746016700352029, -0.197101203594028, -1.3860007472886064, -8.289767128060362] )
    coef9 = np.array( [104.40501489600803, -29.942999569370276, -174.42008279158216, 3.708514822202037, 25.14336546356742, 132.35546551746415, -0.16479555172343271, -1.351556690339512, -6.721810844025761, -40.39565289044579] )

    # Calculating EMG distribution parameters (looking up or reading in)
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
        # Calculating EMG distribution parameters
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
            # SIMPSON'S RULE INTEGRATION
            inte = 1 # Simpson's Rule
            # inte = 2 # Trapezoidal Rule (optional ability of the code-- faster but less accurate)

            vel_xs,vel_ys = _vawtwake.vel_field(xt,yt,x0,y0,dia,rot,chord,B,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,velf,m,n,inte)

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
            # 21 POINT GAUSS-KRONROD RULE QUADRATURE INTEGRATION
            xbound = (scl3+5.)*dia
            argval = (x0t,y0t,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3)
            if veltype == 'all' or veltype == 'x' or veltype == 'ind':
                vel_x = dblquad(_vawtwake.integrandx,0.,xbound,lambda x: -1.*dia,lambda x: 1.*dia,args=argval)
                vel_xs = (vel_x[0]*rot)/(2.*pi)
            if veltype == 'all' or veltype == 'y' or veltype == 'ind':
                vel_y = dblquad(_vawtwake.integrandy,0.,xbound,lambda x: -1.*dia,lambda x: 1.*dia,args=argval)
                vel_ys = (vel_y[0]*rot)/(2.*pi)

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