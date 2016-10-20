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
rotation direction is irrelevant. The imported velocity data assumes the
rotation is counter-clockwise from the top.

"""
import numpy as np
from numpy import pi, exp, fabs, sqrt, log
# from scipy.integrate import dblquad
from Integrate import dblquad
import matplotlib.pyplot as plt

import _vortmodel
import _vort_integrate
import _vawtwake

# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'

def parameterval(tsr,sol,coef):
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


def velocity_field(xt, yt, x0, y0, velf, dia, tsr, solidity, cfd_data='vort2', param=None, veltype='all', m=1, n=1):
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
    if veltype != 'velfort':
        x0t = x0 - xt
        y0t = y0 - yt

    if cfd_data == 'vort':
        # Calculating EMG distribution parameters
        loc = param[0]
        spr = param[1]
        skw = param[2]
        scl = param[3]
        
        # Integration of the vorticity profile using Fortran code (vorticity.f90)
        xbound = 45.*dia/max(1, np.int(tsr*solidity))
        vel_vs = dblquad(_vortmodel.integrandx, 0., xbound, lambda x: -4.*dia, lambda x: 4.*dia,
                         args=(x0t, y0t, dia, loc[0], loc[1], loc[2], spr[0], spr[1],
                               skw[0], skw[1], scl[0], scl[1], scl[2]))
        
        # Calculating velocity deficit
        vel = (vel_vs[0]*rot)/(2.*pi)
        vel = (vel + velf)/velf  # normalization of velocity

    elif cfd_data == 'velo':
        # Normalizing the downstream and lateral positions by the turbine diameter
        x0d = x0/dia
        y0d = y0/dia

        # # Calculating SMG distribution parameters
        # men = param[0]
        # spr = param[1]
        # scl = param[2]
        # rat = param[3]
        # wdt = param[4]
        #
        # men_v = men[0]*x0d**2 + men[1]*x0d + men[2]
        # if men_v > 1.5:
        #     men_v = 1.5
        # elif men_v < -1.5:
        #     men_v = -1.5
        #
        # spr_v = spr[2]*spr[1]*spr[0]*exp(spr[1]*x0d)*exp(spr[0])*exp(-spr[0]*exp(spr[1]*x0d)) + spr[3]
        #
        # scl_v = scl[2]*scl[1]*scl[0]*exp(scl[1]*x0d)*exp(scl[0])*exp(-scl[0]*exp(scl[1]*x0d))
        #
        # rat_v = rat[0]*x0d + rat[1]
        # if rat_v < 0.:
        #     rat_v = 0.
        #
        # wdt_v = wdt[0]*x0d + wdt[1]
        # if wdt_v < 0.:
        #     wdt_v = 0.
        #
        # vel = (-scl_v/(spr_v*sqrt(2.*pi))*exp(-(y0d+men_v)**2/(2.*spr_v**2))) * \
        #       (1./(1 + exp(rat_v*fabs(y0d)-wdt_v))) + 1.  # Normal distribution with sigmoid weighting

        # Calculating SMG distribution parameters
        men1 = param[0]
        sdv1 = param[1]
        rat1 = param[2]
        wdt1 = param[3]
        spr1 = param[4]
        scl1 = param[5]
        men2 = param[6]
        sdv2 = param[7]
        rat2 = param[8]
        wdt2 = param[9]
        spr2 = param[10]
        scl2 = param[11]
        men3 = param[12]
        sdv3 = param[13]
        rat3 = param[14]
        wdt3 = param[15]
        spr3 = param[16]
        scl3 = param[17]
        men4 = param[18]
        sdv4 = param[19]
        rat4 = param[20]
        wdt4 = param[21]
        spr4 = param[22]
        scl4 = param[23]
        p = param[24]
        q = param[25]

        #
        # f1 = -1./(sdv*sqrt(2.*pi))*exp(-(y0d-men)**2/(2.*sdv**2))*(1./(1.+exp(rat*fabs(y0d)-spr)))
        # f2 = scl3*scl2*scl1*exp(scl2*x0d)*exp(scl1)*exp(-scl1*exp(scl2*x0d))
        # vel = f1*f2 + 1.

        # #s3,150
        # men[0] = 0.0386330700933
        # sdv[0] = 0.0059298340497
        # sdv[1] = 0.0844659700584
        # sdv[2] = 18.870370513
        # sdv[3] = 0.162004601856
        # rat[0] = 20.9494045568
        # wdt[0] = 8.4990599166
        # spr[0] = 0.836092366775
        # spr[1] = 0.916710794393
        # spr[2] = 0.0
        # spr[3] = 1.91314592022
        # scl[0] = 1.0
        # scl[1] = 0.0345051895021
        # scl[2] = 6.10731305639
        # #s3,175
        # men[0] = 0.41083553787
        # sdv[0] = 3.52518060451
        # sdv[1] = 1.22591011537e-05
        # sdv[2] = 19.6488391747
        # sdv[3] = 0.233826597495
        # rat[0] = 67.9357459734
        # wdt[0] = 4.01469967395
        # spr[0] = 1.26316250777e-09
        # spr[1] = 0.86004212226
        # spr[2] = 0.0466368647533
        # spr[3] = 7.70202114915
        # scl[0] = 0.663204353332
        # scl[1] = 0.0663868879855
        # scl[2] = 28.5130674608
        # #s2,450
        # men = -0.0654893352396
        # sdv1 = 0.167326286784
        # sdv2 = 0.289238791814
        # sdv3 = 4.92522066656
        # sdv4 = 0.439504215756
        # rat = 19.9541588088
        # wdt = 16.2640455645
        # spr1 = 0.127551157116
        # spr2 = 0.352282606588
        # spr3 = 1.14681685429
        # spr4 = 0.780097482758
        # scl1 = 0.27649709454
        # scl2 = 0.203510889145
        # scl3 = 21.3948777262

        # men = -0.0711937539023
        # sdv1 = 0.206060841675
        # sdv2 = 0.272499100971
        # sdv3 = 5.99939629668
        # sdv4 = 0.537639817765
        # rat = 18.0195576223
        # wdt = 16.9505222198
        # spr1 = 0.552010770521
        # spr2 = 0.10205357438
        # spr3 = 21.5938168453
        # spr4 = 0.0
        # scl1 = 0.325063169077
        # scl2 = 0.186084737017
        # scl3 = 33.4285641474

        # sdv_v1 = sdv[2]*sdv[1]*sdv[0]*exp(sdv[1]*x0d)*exp(sdv[0])*exp(-sdv[0]*exp(sdv[1]*x0d))+sdv[3]
        # # sdv_v = 1.
        #
        # spr_v = spr[2]*spr[1]*spr[0]*exp(spr[1]*x0d)*exp(spr[0])*exp(-spr[0]*exp(spr[1]*x0d))+spr[3]
        # # spr_v = 1.
        #
        # f1 = -1./(sdv_v*sqrt(2.*pi))*exp(-((y0d/spr_v)-men[0])**2/(2.*sdv_v**2))*(1./(1.+exp(rat[0]*fabs((y0d/spr_v))-wdt[0])))
        # f2 = scl[2]*scl[1]*scl[0]*exp(scl[1]*x0d)*exp(scl[0])*exp(-scl[0]*exp(scl[1]*x0d))
        #
        # vel = f1*f2 + 1.

        if p == 0 and q == 0:
            # TSR and solidity are both in CFD data set
            sdv_v1 = sdv1[2]*sdv1[1]*sdv1[0]*exp(sdv1[1]*x0d)*exp(-sdv1[0]*exp(sdv1[1]*x0d))+sdv1[3]
            spr_v1 = spr1[2]*spr1[1]*spr1[0]*exp(spr1[1]*x0d)*exp(-spr1[0]*exp(spr1[1]*x0d))#+spr1[3]

            f11 = -1./(sdv_v1*sqrt(2.*pi))*exp(-((y0d/spr_v1)-men1[0])**2/(2.*sdv_v1**2))*(1./(1.+exp(rat1[0]*fabs((y0d/spr_v1))-wdt1[0])))
            f21 = scl1[2]*scl1[1]*scl1[0]*exp(scl1[1]*x0d)*exp(-scl1[0]*exp(scl1[1]*x0d))

            vel = f11*f21 + 1.

        elif p != 0 and q == 0:
            # solidity is in CFD data set and TSR needs interpolation
            sdv_v1 = sdv1[2]*sdv1[1]*sdv1[0]*exp(sdv1[1]*x0d)*exp(-sdv1[0]*exp(sdv1[1]*x0d))+sdv1[3]
            spr_v1 = spr1[2]*spr1[1]*spr1[0]*exp(spr1[1]*x0d)*exp(-spr1[0]*exp(spr1[1]*x0d))+spr1[3]

            f11 = -1./(sdv_v1*sqrt(2.*pi))*exp(-((y0d/spr_v1)-men1[0])**2/(2.*sdv_v1**2))*(1./(1.+exp(rat1[0]*fabs((y0d/spr_v1))-wdt1[0])))
            f21 = scl1[2]*scl1[1]*scl1[0]*exp(scl1[1]*x0d)*exp(-scl1[0]*exp(scl1[1]*x0d))

            vel1 = f11*f21 + 1.

            sdv_v2 = sdv2[2]*sdv2[1]*sdv2[0]*exp(sdv2[1]*x0d)*exp(-sdv2[0]*exp(sdv2[1]*x0d))+sdv2[3]
            spr_v2 = spr2[2]*spr2[1]*spr2[0]*exp(spr2[1]*x0d)*exp(-spr2[0]*exp(spr2[1]*x0d))+spr2[3]

            f12 = -1./(sdv_v2*sqrt(2.*pi))*exp(-((y0d/spr_v2)-men2[0])**2/(2.*sdv_v2**2))*(1./(1.+exp(rat2[0]*fabs((y0d/spr_v2))-wdt2[0])))
            f22 = scl2[2]*scl2[1]*scl2[0]*exp(scl2[1]*x0d)*exp(-scl2[0]*exp(scl2[1]*x0d))

            vel2 = f12*f22 + 1.

            vel = (1.-p)*vel1 + p*vel2

        elif p == 0 and q != 0:
            # TSR is in CFD data set and solidity needs interpolation
            sdv_v3 = sdv3[2]*sdv3[1]*sdv3[0]*exp(sdv3[1]*x0d)*exp(-sdv3[0]*exp(sdv3[1]*x0d))+sdv3[3]
            spr_v3 = spr3[2]*spr3[1]*spr3[0]*exp(spr3[1]*x0d)*exp(-spr3[0]*exp(spr3[1]*x0d))+spr3[3]

            f13 = -1./(sdv_v3*sqrt(2.*pi))*exp(-((y0d/spr_v3)-men3[0])**2/(2.*sdv_v3**2))*(1./(1.+exp(rat3[0]*fabs((y0d/spr_v3))-wdt3[0])))
            f23 = scl3[2]*scl3[1]*scl3[0]*exp(scl3[1]*x0d)*exp(-scl3[0]*exp(scl3[1]*x0d))

            vel3 = f13*f23 + 1.

            sdv_v4 = sdv4[2]*sdv4[1]*sdv4[0]*exp(sdv4[1]*x0d)*exp(-sdv4[0]*exp(sdv4[1]*x0d))+sdv4[3]
            spr_v4 = spr4[2]*spr4[1]*spr4[0]*exp(spr4[1]*x0d)*exp(-spr4[0]*exp(spr4[1]*x0d))+spr4[3]

            f14 = -1./(sdv_v4*sqrt(2.*pi))*exp(-((y0d/spr_v4)-men4[0])**2/(2.*sdv_v4**2))*(1./(1.+exp(rat4[0]*fabs((y0d/spr_v4))-wdt4[0])))
            f24 = scl4[2]*scl4[1]*scl4[0]*exp(scl4[1]*x0d)*exp(-scl4[0]*exp(scl4[1]*x0d))

            vel4 = f14*f24 + 1.

            vel = (1.-q)*vel3 + q*vel4

        else:
            # both TSR and solidity need interpolation
            sdv_v1 = sdv1[2]*sdv1[1]*sdv1[0]*exp(sdv1[1]*x0d)*exp(-sdv1[0]*exp(sdv1[1]*x0d))+sdv1[3]
            spr_v1 = spr1[2]*spr1[1]*spr1[0]*exp(spr1[1]*x0d)*exp(-spr1[0]*exp(spr1[1]*x0d))+spr1[3]

            f11 = -1./(sdv_v1*sqrt(2.*pi))*exp(-((y0d/spr_v1)-men1[0])**2/(2.*sdv_v1**2))*(1./(1.+exp(rat1[0]*fabs((y0d/spr_v1))-wdt1[0])))
            f21 = scl1[2]*scl1[1]*scl1[0]*exp(scl1[1]*x0d)*exp(-scl1[0]*exp(scl1[1]*x0d))

            vel1 = f11*f21 + 1.

            sdv_v2 = sdv2[2]*sdv2[1]*sdv2[0]*exp(sdv2[1]*x0d)*exp(-sdv2[0]*exp(sdv2[1]*x0d))+sdv2[3]
            spr_v2 = spr2[2]*spr2[1]*spr2[0]*exp(spr2[1]*x0d)*exp(-spr2[0]*exp(spr2[1]*x0d))+spr2[3]

            f12 = -1./(sdv_v2*sqrt(2.*pi))*exp(-((y0d/spr_v2)-men2[0])**2/(2.*sdv_v2**2))*(1./(1.+exp(rat2[0]*fabs((y0d/spr_v2))-wdt2[0])))
            f22 = scl2[2]*scl2[1]*scl2[0]*exp(scl2[1]*x0d)*exp(-scl2[0]*exp(scl2[1]*x0d))

            vel2 = f12*f22 + 1.

            sdv_v3 = sdv3[2]*sdv3[1]*sdv3[0]*exp(sdv3[1]*x0d)*exp(-sdv3[0]*exp(sdv3[1]*x0d))+sdv3[3]
            spr_v3 = spr3[2]*spr3[1]*spr3[0]*exp(spr3[1]*x0d)*exp(-spr3[0]*exp(spr3[1]*x0d))+spr3[3]

            f13 = -1./(sdv_v3*sqrt(2.*pi))*exp(-((y0d/spr_v3)-men3[0])**2/(2.*sdv_v3**2))*(1./(1.+exp(rat3[0]*fabs((y0d/spr_v3))-wdt3[0])))
            f23 = scl3[2]*scl3[1]*scl3[0]*exp(scl3[1]*x0d)*exp(-scl3[0]*exp(scl3[1]*x0d))

            vel3 = f13*f23 + 1.

            sdv_v4 = sdv4[2]*sdv4[1]*sdv4[0]*exp(sdv4[1]*x0d)*exp(-sdv4[0]*exp(sdv4[1]*x0d))+sdv4[3]
            spr_v4 = spr4[2]*spr4[1]*spr4[0]*exp(spr4[1]*x0d)*exp(-spr4[0]*exp(spr4[1]*x0d))+spr4[3]

            f14 = -1./(sdv_v4*sqrt(2.*pi))*exp(-((y0d/spr_v4)-men4[0])**2/(2.*sdv_v4**2))*(1./(1.+exp(rat4[0]*fabs((y0d/spr_v4))-wdt4[0])))
            f24 = scl4[2]*scl4[1]*scl4[0]*exp(scl4[1]*x0d)*exp(-scl4[0]*exp(scl4[1]*x0d))

            vel4 = f14*f24 + 1.

            vel = (((1.-p)*vel1 + p*vel2) + ((1.-q)*vel3 + q*vel4))/2.

        if x0 < xt:
            vel = 1.  # Velocity is free stream in front and to the sides of the turbine

    elif cfd_data == 'vort2':
        # Calculating EMG distribution parameters (looking up or reading in)
        if param is None:
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

            loc1 = parameterval(tsr,solidity,coef0)
            loc2 = parameterval(tsr,solidity,coef1)
            loc3 = parameterval(tsr,solidity,coef2)
            spr1 = parameterval(tsr,solidity,coef3)
            spr2 = parameterval(tsr,solidity,coef4)
            skw1 = parameterval(tsr,solidity,coef5)
            skw2 = parameterval(tsr,solidity,coef6)
            scl1 = parameterval(tsr,solidity,coef7)
            scl2 = parameterval(tsr,solidity,coef8)
            scl3 = parameterval(tsr,solidity,coef9)


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


        # Integration of the vorticity profile using Fortran code (vorticity.f90)
        # xbound = 35.*dia/max(1, int(tsr*solidity))
        xbound = (scl3+5.)*dia
        if veltype == 'all' or veltype == 'x' or veltype == 'ind':
            vel_x = dblquad(_vortmodel.integrandx, 0., xbound, lambda x: -1.*dia, lambda x: 1.*dia,
                             args=(x0t, y0t, dia, loc1, loc2, loc3, spr1, spr2,
                                   skw1, skw2, scl1, scl2, scl3))
            vel_xs = (vel_x[0]*rot)/(2.*pi)
        if veltype == 'all' or veltype == 'y' or veltype == 'ind':
            vel_y = dblquad(_vortmodel.integrandy, 0., xbound, lambda x: -1.*dia, lambda x: 1.*dia,
                             args=(x0t, y0t, dia, loc1, loc2, loc3, spr1, spr2,
                                   skw1, skw2, scl1, scl2, scl3))
            vel_ys = (vel_y[0]*rot)/(2.*pi)
        if veltype == 'vort':
            if x0t < 0.:
                vel = 0.
            else:
                vel = _vawtwake.vorticitystrength(x0t,y0t,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3)

        # vel_vs = _vort_integrate.runvort(x0t,y0t,dia,xbound,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3)
        if veltype == 'velfort':
            # m = xbound/0.2
            # n = 100.
            m = 110.
            n = 100.
            int = 1 #Simpson's Rule
            # int = 2 #Trapezoidal Rule

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

            # vel = _vortmodel.velcalc(x0t,y0t,dia,rot,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,velf,0.,xbound,-1.*dia,1.*dia,m,n)
            vel_xs,vel_ys = _vawtwake.vel_field(xt,yt,x0,y0,dia,rot,solidity*rad/3.,3.,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,velf,m,n,int)

        # Calculating velocity deficit
        if veltype == 'all':
            vel = sqrt((vel_xs + velf)**2 + (vel_ys)**2)
        elif veltype == 'x':
            vel = vel_xs + velf
        elif veltype == 'y':
            vel = vel_ys
        elif veltype == 'ind':
            vel = np.array([vel_xs,vel_ys])
        elif veltype == 'velfort':
            # vel = sqrt((vel_xs*velf + velf)**2 + (vel_ys*velf)**2)
            vel = vel_xs*velf + velf

        vel = vel/velf


    elif cfd_data == 'velo2':
        x0d = x0/dia
        y0d = y0/dia

        # Calculating SMG distribution parameters
        spr1 = param[0]
        pow1 = param[1]
        pow2 = param[2]
        spr2 = param[3]
        skw = param[4]
        scl1 = param[5]
        scl2 = param[6]
        scl3 = param[7]

        pow = pow1-pow2*x0d**2
        if pow < 1.5:
            pow = 1.5
        # spr2_fix = 1./((log(0.0001)/-spr1)**(1./pow)+fabs(skw))**2
        # if spr2 > spr2_fix:
        #     spr2 = spr2_fix
        if spr2 < 0.:
            spr2 = 0.
        exp_v = exp(-spr1*fabs(y0d)**pow)
        quad = spr2*(y0d-skw)**2-1.
        scl_v = scl3*scl2*scl1*exp(scl2*x0d)*exp(-scl1*exp(scl2*x0d))

        vel = exp_v*quad*scl_v + 1.

        if x0 < xt:
            vel = 1.  # Velocity is free stream in front and to the sides of the turbine

    return vel,m,n


