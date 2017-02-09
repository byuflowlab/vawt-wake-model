"""
Parameterized VAWT Wake Model using CFD vorticity data
Developed by Eric Tingey at Brigham Young University

This code models the wake behind a vertical-axis wind turbine based on
tip-speed ratio, solidity and wind speed by converting the vorticity of
the wake into velocity information. The model uses CFD data obtained
from STAR-CCM+ of simulated turbines to make the wake model as accurate
as possible.

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
from numpy import pi,fabs,sqrt,sin,cos
# from scipy.integrate import dblquad # integration with complete module
from Integrate import dblquad # integration with simplification file (Integrate.py)
from scipy.interpolate import UnivariateSpline

from joblib import Parallel,delayed

import _vawtwake


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
    Published values from paper (4 Reynolds numbers; 1 2 2 1 1 1 2 1 2 2) with TSR and solidity orders shown

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
    loc1 = np.array( [0.0025703809856661534, -0.0007386258659065129, 0.004595508188667984, 0.000380123563204793, -0.0005090098755683027, 0.005744581813281894, -4.103393770815313e-05, -0.0014146918534486358, -0.013975958482495927, 0.0] ) # TSR: 3; Solidity: 2
    loc2 = np.array( [-0.5047504670963536, 0.23477391362058556, 0.8414256436198028, -0.04252528528617351, -0.06962875967504166, -0.6566907653208429, 0.002839318332370807, 0.00571803958194812, 0.0070744372783060295, 0.22805286438890995] ) # TSR: 3; Solidity: 3
    loc3 = np.array( [0.2878345841026334, 0.11512552658662782, 0.7303949879914625, -0.007035517839387948, -0.18284850673545897, -0.5241921153256568, -0.0003704899921255296, 0.010972527139685873, 0.04380801537377295, 0.1724129349605399] ) # TSR: 3; Solidity: 3
    spr1 = np.array( [0.08234816067475287, -0.03530687906626052, -0.3662863944976986, 0.003240141344532779, 0.12172015102204112, 0.2993048183466721, 0.0, -0.009253185586804007, -0.057469126406649716, -0.07257633583877886] ) # TSR: 2; Solidity: 3
    spr2 = np.array( [-0.07083579909945328, 0.016182024377569406, 0.1985436342461859, 0.0017738254727425816, -0.09111094817943823, -0.06561408122153217, -0.0005115133402638633, 0.009434288536679505, 0.022392136905926813, 0.0] ) # TSR: 3; Solidity: 2
    skw1 = np.array( [-1.6712830849073221, 1.5625053380692426, -6.180392756736983, -0.20407668040293722, -4.6476103643607685, 29.380064536220306, 0.0, 0.7502978877582536, -0.16358232641365608, -19.937609244085568] ) # TSR: 2; Solidity: 3
    skw2 = np.array( [-3.423561091777921, -9.228795430171687, 86.95722105482042, 2.772872601988039, -11.968168333741515, -150.61261090270446, -0.24715316589674527, 0.5283723108899993, 4.537286811245538, 82.50581844010263] ) # TSR: 3; Solidity: 3
    scl1 = np.array( [-0.19815381951708524, 0.08438758133540872, 1.2650146439483734, -0.007606115512168328, -0.2747023984740461, -0.8844640101378567, 0.0, 0.01870057580949183, 0.0699898278743648, 0.2794360008051127] ) # TSR: 2; Solidity: 3
    scl2 = np.array( [2.3932787625531815, -2.020874419612962, -8.938221963838357, 0.576323845480877, 2.8782448498416944, 16.598492450314534, -0.04746016700352029, -0.197101203594028, -1.3860007472886064, -8.289767128060362] ) # TSR: 3; Solidity: 3
    scl3 = np.array( [104.40501489600803, -29.942999569370276, -174.42008279158216, 3.708514822202037, 25.14336546356742, 132.35546551746415, -0.16479555172343271, -1.351556690339512, -6.721810844025761, -40.39565289044579] ) # TSR: 3; Solidity: 3

    # loc1 = np.array( [0.0024881782234644594, -0.0007162124520772835, 0.004840905629292527, 0.00035940014241851975, -0.0004835618309358686, 0.0057312930423697225, -3.749132238596887e-05, -0.0014448678350805107, -0.013838697667670624, 0.0] )
    # loc2 = np.array( [-0.5148970510737115, 0.23664951901881176, 0.8508796336915305, -0.0422155314186519, -0.06965958435002093, -0.6590584244238057, 0.002780441298933052, 0.005527629993589326, 0.0064771325789259974, 0.2262447334696142] )
    # loc3 = np.array( [0.31053085704135697, 0.11264440683388684, 0.7085100386881281, -0.007855277883163883, -0.18261141342231252, -0.5252607546488571, -0.0002593237635706134, 0.011353734667742903, 0.04461611956627692, 0.17847967604672085] )
    # spr1 = np.array( [0.08160365473705117, -0.035209146348268265, -0.3657089321462334, 0.0032095158671521656, 0.12215088328116319, 0.29924788767922195, 0.0, -0.009337895738912123, -0.05728786434141602, -0.0734184612863306] )
    # spr2 = np.array( [-0.07046108538731136, 0.016662058569062665, 0.20178424265599015, 0.0016797319754888462, -0.09315936851491065, -0.07057840461085574, -0.00048637093956067695, 0.009366304029411365, 0.024572217032978416, 0.0] )
    # skw1 = np.array( [-1.6617923988070358, 1.5456256500296082, -6.15462425653742, -0.20731723544238276, -4.636823166329576, 29.322817794328987, 0.0, 0.7439115823342768, -0.16196462404243628, -20.075376228055998] )
    # skw2 = np.array( [-3.42356282753647, -9.209290631300641, 87.90211972541844, 2.773568799442361, -11.695774984376406, -151.05284030686784, -0.24689398793699247, 0.5033111918389676, 4.197594495917425, 82.25216982468561] )
    # scl1 = np.array( [-0.26931621108474796, 0.10907476831328872, 1.4442635090507794, -0.00968517038054994, -0.3248122217820121, -0.9934068714474691, 0.0, 0.021876975741408893, 0.08792346901752275, 0.2905428329957239] )
    # scl2 = np.array( [-1.065367001488832, 0.11114929355575819, -3.3141966280358215, 0.08516024237491295, 2.6804608016839513, 6.030096841771481, -0.011684187520652582, -0.16727191293057872, -1.5716817442222772, -1.7538392539592593] )
    # scl3 = np.array( [103.58755782882461, -30.004921914818564, -171.8165540416405, 3.757680076380713, 24.89320887829779, 130.45773327654743, -0.167933443248442, -1.3695260144709467, -6.503525472627104, -40.20025748470305] )

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
    tsr = rad*fabs(rot)/velf
    solidity = (chord*B)/rad

    # Translating the turbine position
    x0t = x0 - xt
    y0t = y0 - yt

    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = coef_val()

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
            # 21 POINT GAUSS-KRONROD RULE QUADRATURE INTEGRATION
            xbound = (scl3+5.)*dia
            argval = (x0t,y0t,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3)
            if veltype == 'all' or veltype == 'x' or veltype == 'ind':
                vel_x = dblquad(_vawtwake.integrandx,0.,xbound,lambda x: -1.*dia,lambda x: 1.*dia,args=argval)
                vel_xs = (vel_x[0]*fabs(rot))/(2.*pi)
            if veltype == 'all' or veltype == 'y' or veltype == 'ind':
                vel_y = dblquad(_vawtwake.integrandy,0.,xbound,lambda x: -1.*dia,lambda x: 1.*dia,args=argval)
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

def overlap(p,xt,yt,diat,rott,chord,B,x0,y0,dia,velf,param=None,veltype='ind',integration='gskr'):
    """
    Calculating wake velocities around a turbine based on wake overlap from surrounding turbines
    (using the 21 point Gauss-Kronrod rule quadrature integration; Simpson's rule integration can be used via VAWT_Wake_Model.f90)

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

    # finding points around the flight path of the blades
    for i in range(p):
      theta = (2.0*pi/p)*i-(2.0*pi/p)/2.0
      xd[i] = x0 - sin(theta)*(dia/2.0)
      yd[i] = y0 + cos(theta)*(dia/2.0)

    intex = np.zeros(p)
    intey = np.zeros(p)

    if (t == 1): # coupled configuration (only two VAWTs)
        wake = Parallel(n_jobs=-1)(delayed(velocity_field)(xt[0],yt[0],xd[j],yd[j],velf,diat[0],rott[0],chord,B,param,veltype,integration) for j in range(p) )
        for i in range(p):
            velx[i] = wake[i][0]*velf
            vely[i] = wake[i][1]*velf
        # for j in range(p):
        #     wake = velocity_field(xt[0],yt[0],xd[j],yd[j],velf,diat[0],rott[0],chord,B,param,veltype,integration)
        #     velx[j] = wake[0]*velf
        #     vely[j] = wake[1]*velf

    else: # multiple turbine wake overlap
        for j in range(t):
            for k in range(p):
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