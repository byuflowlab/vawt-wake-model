import numpy as np
from numpy import cos,sin,pi
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
from numpy import sqrt
import VAWT_Wake_Model as vwm
from ACsingle import actuatorcylinder
from os import path

import time,sys
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

import _vawtwake

# Progress bar for plotting
def progress_bar(percent,tasks,runtime):
    bar_long = 40
    timeleft = (runtime)*(tasks*(1.-percent))
    if timeleft < 120.:
        status = 'Working... '+str(int(timeleft))+' seconds left'
    else:
        timeleft = timeleft/60.
        status = 'Working... '+str(int(timeleft))+' minutes left'
    if percent == 1:
        status = 'Complete\n'
    bar_seg = int(round(bar_long*percent))
    text = '\rStatus: [{0}] {1}% {2}'.format( '='*bar_seg + ' '*(bar_long - bar_seg), percent*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()


turbine = 'uppsala'
turbine = 'windspire'
# turbine = 'delft'

method = 'power_curve'
method = 'coupled_vawt_power'
# method = 'polar_power'
method = 'velocity_overlap'

int_type = 'simp'
# int_type = 'gskr'

rotdir = 'corot'
rotdir = 'counter'


ntheta = 36
interp = 1 # linear airfoil interpolation
# interp = 2 # cubic spline interpolation

m = 220#1000
n = 200#1000
d_space = 2.0

basepath = path.join(path.dirname(path.realpath('__file__')), 'data/airfoils')

if turbine == 'uppsala':
    xt = np.array([0.])
    yt = np.array([0.])
    dia = 6.
    diat = np.ones_like(xt)*dia
    velf = 15.
    r = dia/2.
    tsrd = 4.
    rot = tsrd*velf/r

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.25
    H = 5.
    foildata = basepath + path.sep + 'NACA_0021.dat'

    rho = 1.225
    mu = 1.7894e-5

elif turbine == 'windspire':
    xt = np.array([0.,2.])
    yt = np.array([0.,-2.])
    dia = 1.2
    diat = np.ones_like(xt)*dia
    velf = 1.
    r = dia/2.
    tsrd = 2.625
    rot = tsrd*velf/r

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.128
    H = 6.1
    foildata = basepath + path.sep + 'du06w200.dat'

    rho = 1.225
    mu = 1.7894e-5

elif turbine == 'delft':
    xt = np.array([0.])
    yt = np.array([0.])
    dia = 1.
    diat = np.ones_like(xt)*dia
    velf = 10.
    r = dia/2.
    tsrd = 4.
    rot = tsrd*velf/r

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.1*r
    H = 1.
    foildata = basepath + path.sep + 'NACA_0015.dat'

    rho = 1.225
    mu = 1.7894e-5

af_data,cl_data,cd_data = vwm.airfoil_data(foildata)
coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if method == 'power_curve':

    N = 100
    low = 1.5
    high = 4.5

    tsr = np.linspace(low,high,N)
    cp_plot = np.zeros_like(tsr)

    time0 = time.time()
    for j in range(np.size(tsr)):
        rot = np.ones_like(xt)*(tsr[j]*velf/r)

        _,_,_,Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot[0],velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))

        cp_plot[j] = Cp
        runtime = time.time()-time0
        progress_bar(float(j+1)/N,N,runtime)
        time0 = time.time()

    if turbine == 'windspire':
        juliatsr = np.linspace(1.5,4.5,100)
        juliacp = np.array([0.0716233,0.0762639,0.0814103,0.0871416,0.0922382,0.0967217,0.101195,0.106454,0.111312,0.115768,0.12053,0.126162,0.130722,0.135023,0.140747,0.147514,0.154634,0.162332,0.171205,0.179536,0.190313,0.202898,0.214106,0.222564,0.229452,0.235479,0.240708,0.24525,0.249453,0.253586,0.25762,0.261289,0.264816,0.268426,0.271909,0.275398,0.27914,0.282687,0.286247,0.289792,0.293371,0.296686,0.299124,0.299987,0.299863,0.299155,0.29802,0.296533,0.294743,0.292687,0.290393,0.287883,0.285171,0.282271,0.27919,0.275937,0.272518,0.268937,0.265199,0.261305,0.257259,0.253063,0.248716,0.24422,0.239576,0.234785,0.229846,0.224761,0.21953,0.214153,0.20863,0.202961,0.197146,0.191185,0.185077,0.178823,0.172421,0.165872,0.159176,0.152333,0.145343,0.138206,0.130922,0.12349,0.115911,0.108185,0.100312,0.0922939,0.0841295,0.0758191,0.0673629,0.0587608,0.0500128,0.0411189,0.0320789,0.0228926,0.0135599,0.00408066,-0.00554549,-0.0153188])

    plt.figure()
    if turbine == 'windspire':
        plt.plot(tsr,cp_plot,'bo-',label='Python')
        plt.plot(juliatsr,juliacp,'r-',label='Julia')
        plt.legend(loc=1)
    else:
        plt.plot(tsr,cp_plot,'bo-')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$C_p$')

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

elif method == 'coupled_vawt_power':

    nT = 72
    nR = 21
    Theta = np.linspace(-180, 180, nT)*pi/180.0
    Radii = np.linspace(2.0, 12.0, nR)
    X = np.zeros((nR,nT))
    Y = np.zeros((nR,nT))
    P = np.zeros((nR,nT))
    P1 = np.zeros((nR,nT))
    P2 = np.zeros((nR,nT))

    rot1 = rot
    if rotdir == 'corot':
        rot2 = rot
    elif rotdir == 'counter':
        rot2 = -rot

    _,_,_,Cp_iso,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))

    iter = 0
    time0 = time.time()
    for i in range(nR):
        for j in range(nT):

            centerX = Radii[i]*r*cos(Theta[j])
            centerY = Radii[i]*r*sin(Theta[j])

            X[i,j] = centerX
            Y[i,j] = centerY

            if int_type == 'simp':
                wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([centerX]),np.array([centerY]),np.array([dia]),np.array([rot2]),chord,B,0.,0.,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
                wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot1]),chord,B,centerX,centerY,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
            elif int_type == 'gskr':
                wakex1,wakey1 = vwm.overlap(ntheta,np.array([centerX]),np.array([centerY]),np.array([dia]),np.array([rot2]),chord,B,0.,0.,dia,velf,False)
                wakex2,wakey2 = vwm.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot1]),chord,B,centerX,centerY,dia,velf,False)

            _,_,_,Cp1,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,velf,rho,mu,interp,wakex1,wakey1)
            _,_,_,Cp2,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot2,velf,rho,mu,interp,wakex2,wakey2)


            P[i,j] = (Cp1+Cp2)/(2*Cp_iso)
            P1[i,j] = Cp1/Cp_iso
            P2[i,j] = Cp2/Cp_iso
            iter += 1
            runtime = time.time()-time0
            progress_bar(float(iter)/(nR*nT),nR*nT,runtime)
            time0 = time.time()

    fs = 20
    plt.figure()
    lb = 0.9 # lower bound on velocity to display
    ub = 1.1 # upper bound on velocity to display
    ran = 100 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,5) # setting the number of tick marks on colorbar
    CS = plt.contourf(X/dia,Y/dia,P,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel('Normalized Power',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    plt.xlabel('$x/D$',fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)
    if rotdir == 'counter':
        plt.title('Combined Power: Counter-Rotating')
    elif rotdir == 'corot':
        plt.title('Combined Power: Co-Rotating')
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlim(-6.1,6.1)
    plt.ylim(-6.1,6.1)
    for i in range(np.size(xt)):
        circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
        plt.gca().add_patch(circ1)
        circ = mpatches.Circle((xt[i]/dia,yt[i]/dia),r,color='k',linestyle=':',fill=False,linewidth=1)
        ell1 = mpatches.Ellipse((xt[i]/dia,r),chord,chord*0.25,color='k',fill=True)
        ell2 = mpatches.Ellipse((r*cos(210.*pi/180.),r*sin(210.*pi/180.)),chord,chord*0.25,angle=120.,color='k',fill=True)
        ell3 = mpatches.Ellipse((r*cos(330.*pi/180.),r*sin(330.*pi/180.)),chord,chord*0.25,angle=240.,color='k',fill=True)
        plt.gca().add_patch(circ)
        plt.gca().add_patch(ell1)
        plt.gca().add_patch(ell2)
        plt.gca().add_patch(ell3)
        plt.plot((xt[i]/dia,xt[i]/dia),(yt[i]/dia,r),'k',linewidth=1.)
        plt.plot((xt[i]/dia,r*cos(210.*pi/180.)),(yt[i]/dia,r*sin(210.*pi/180.)),'k',linewidth=1.)
        plt.plot((xt[i]/dia,r*cos(330.*pi/180.)),(yt[i]/dia,r*sin(330.*pi/180.)),'k',linewidth=1.)

    plt.figure()
    lb = -0.2 # lower bound on velocity to display
    ub = 1.2 # upper bound on velocity to display
    ran = 8 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,8) # setting the number of tick marks on colorbar
    CS = plt.contourf(X/dia,Y/dia,P1,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel('Normalized Power',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    plt.xlabel('$x/D$',fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)
    if rotdir == 'counter':
        plt.title('Center Turbine Power: Counter-Rotating')
    elif rotdir == 'corot':
        plt.title('Center Turbine Power: Co-Rotating')
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlim(-6.1,6.1)
    plt.ylim(-6.1,6.1)
    for i in range(np.size(xt)):
        circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
        plt.gca().add_patch(circ1)
        circ = mpatches.Circle((xt[i]/dia,yt[i]/dia),r,color='k',linestyle=':',fill=False,linewidth=1)
        ell1 = mpatches.Ellipse((xt[i]/dia,r),chord,chord*0.25,color='k',fill=True)
        ell2 = mpatches.Ellipse((r*cos(210.*pi/180.),r*sin(210.*pi/180.)),chord,chord*0.25,angle=120.,color='k',fill=True)
        ell3 = mpatches.Ellipse((r*cos(330.*pi/180.),r*sin(330.*pi/180.)),chord,chord*0.25,angle=240.,color='k',fill=True)
        plt.gca().add_patch(circ)
        plt.gca().add_patch(ell1)
        plt.gca().add_patch(ell2)
        plt.gca().add_patch(ell3)
        plt.plot((xt[i]/dia,xt[i]/dia),(yt[i]/dia,r),'k',linewidth=1.)
        plt.plot((xt[i]/dia,r*cos(210.*pi/180.)),(yt[i]/dia,r*sin(210.*pi/180.)),'k',linewidth=1.)
        plt.plot((xt[i]/dia,r*cos(330.*pi/180.)),(yt[i]/dia,r*sin(330.*pi/180.)),'k',linewidth=1.)

    plt.figure()
    lb = -0.2 # lower bound on velocity to display
    ub = 1.2 # upper bound on velocity to display
    ran = 8 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,8) # setting the number of tick marks on colorbar
    CS = plt.contourf(X/dia,Y/dia,P2,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel('Normalized Power',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    plt.xlabel('$x/D$',fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)
    if rotdir == 'counter':
        plt.title('Swept Turbine Power: Counter-Rotating')
    elif rotdir == 'corot':
        plt.title('Swept Turbine Power: Co-Rotating')
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlim(-6.1,6.1)
    plt.ylim(-6.1,6.1)
    for i in range(np.size(xt)):
        circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
        plt.gca().add_patch(circ1)
        circ = mpatches.Circle((xt[i]/dia,yt[i]/dia),r,color='k',linestyle=':',fill=False,linewidth=1)
        ell1 = mpatches.Ellipse((xt[i]/dia,r),chord,chord*0.25,color='k',fill=True)
        ell2 = mpatches.Ellipse((r*cos(210.*pi/180.),r*sin(210.*pi/180.)),chord,chord*0.25,angle=120.,color='k',fill=True)
        ell3 = mpatches.Ellipse((r*cos(330.*pi/180.),r*sin(330.*pi/180.)),chord,chord*0.25,angle=240.,color='k',fill=True)
        plt.gca().add_patch(circ)
        plt.gca().add_patch(ell1)
        plt.gca().add_patch(ell2)
        plt.gca().add_patch(ell3)
        plt.plot((xt[i]/dia,xt[i]/dia),(yt[i]/dia,r),'k',linewidth=1.)
        plt.plot((xt[i]/dia,r*cos(210.*pi/180.)),(yt[i]/dia,r*sin(210.*pi/180.)),'k',linewidth=1.)
        plt.plot((xt[i]/dia,r*cos(330.*pi/180.)),(yt[i]/dia,r*sin(330.*pi/180.)),'k',linewidth=1.)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

elif method == 'polar_power':

    rot1 = rot
    if rotdir == 'corot':
        rot2 = rot
    elif rotdir == 'counter':
        rot2 = -rot

    theta_mod = np.linspace(-pi,pi,72)
    cp_mod = np.zeros_like(theta_mod)

    _,_,_,Cp_iso,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))

    if turbine == 'windspire':
        theta_cfd = np.array([270.1114379016532, 246.7988999869089, 225.19062704244018, 201.86037312349106, 179.88586513201912, 516.7008614759567, 494.80957252795724, 471.6647172260273, 449.6585021210285, 426.6668191463685, 404.64546467066145, 381.7268571885939, 359.8858651320191, 336.61717949856603, 314.96577551188346, 292.03610306036995, 270.1114379016532])
        cp_cfd = np.array([0.5137292963669421, 0.9407530012258202, 1.0583347246905448, 1.100368119669435, 1.127448529615236, 1.1014079312070053, 1.0569344196425345, 0.940520231361647, 0.5098077736897938, 0.7368409217354737, 1.0070289977142592, 1.0697784260663252, 1.0882329965142843, 1.0679385516564928, 1.0097941954631366, 0.7333565294554432, 0.5137292963669421])
        theta_ac = np.linspace(-pi,pi,72)
        cp_ac = np.array([0.445982, 0.51201, 0.615172, 0.712552, 0.776262, 0.814086, 0.881157, 0.93065, 0.995502, 1.00445, 1.01043, 1.01692, 1.02396, 1.03157, 1.03922, 1.04553, 1.05007, 1.05269, 1.05322, 1.05163, 1.04803, 1.04259, 1.03549, 1.02769, 1.02037, 1.01361, 1.00738, 1.00164, 0.971863, 0.898155, 0.849936, 0.792628, 0.738928, 0.658594, 0.557209, 0.469759, 0.493928, 0.573753, 0.617474, 0.682425, 0.741807, 0.816748, 0.907066, 0.978855, 0.991511, 0.996756, 1.00322, 1.01079, 1.01907, 1.02736, 1.03491, 1.04111, 1.04541, 1.0474, 1.04672, 1.04353, 1.03821, 1.03127, 1.02326, 1.01488, 1.00689, 0.99984, 0.993982, 0.991411, 0.94704, 0.869073, 0.790213, 0.714075, 0.653963, 0.605099, 0.532343, 0.445982])

    # for i in range(np.size(theta_mod)):
    #     x1 = 0.
    #     y1 = 0.
    #     x2 = 2.*cos(theta_mod[i])
    #     y2 = 2.*sin(theta_mod[i])
    #
    #     if int_type == 'simp':
    #         wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([rot2]),chord,B,x1,y1,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
    #         wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot1]),chord,B,x2,y2,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
    #     elif int_type == 'gskr':
    #         wakex1,wakey1 = vwm.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([rot2]),chord,B,x1,y1,dia,velf,False)
    #         wakex2,wakey2 = vwm.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot1]),chord,B,x2,y2,dia,velf,False)
    #
    #     _,_,_,Cp1,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,velf,rho,mu,interp,wakex1,wakey1)
    #     _,_,_,Cp2,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot2,velf,rho,mu,interp,wakex2,wakey2)
    #
    #     cp_mod[i] = (Cp1+Cp2)/(2*Cp_iso)
    #     print i
    #
    # print cp_mod.tolist()

    cp_mod = np.array([0.3440348622914387, 0.3486250945961999, 0.39894157748818737, 0.5863430562737842, 0.7129970832708875, 0.7855921597699241, 0.8470482538334438, 0.9013953505819924, 0.9601455893263708, 1.0194288950952934, 1.033693122128879, 1.0413605061884483, 1.0506381788443249, 1.060058261050523, 1.0689654414310268, 1.0774951768607937, 1.0851606029911507, 1.089650292761621, 1.0905920596163592, 1.0878281653806119, 1.0816003402932872, 1.0732751972167345, 1.0645617181210305, 1.0554416781535598, 1.0457958116850572, 1.03744720692213, 1.0289739659522332, 0.979338853504734, 0.9328965276766208, 0.8687096863492634, 0.8054556953066844, 0.7613318568157301, 0.6397558676110174, 0.4771682732669759, 0.3632896339495982, 0.3434522229055874, 0.34628206102316905, 0.3878422594267405, 0.5172162997465365, 0.6005684340690277, 0.6994852713991858, 0.7595821027782338, 0.8536741565040549, 0.9534604104561609, 0.9971821817839681, 1.0184854079232033, 1.0262396891829249, 1.0331792074982886, 1.0403024456086578, 1.047134816841171, 1.0531219782184522, 1.0580477306590146, 1.0611763655557784, 1.0624636404921124, 1.0620357786391121, 1.0598530877808174, 1.0557547705362558, 1.0502530689232352, 1.0437918778935287, 1.0367521209323047, 1.0296423808559845, 1.0225430345381388, 1.021304931623804, 0.998891912589102, 0.9049627047290153, 0.8271948349522203, 0.7381963424190326, 0.6500036047549428, 0.5814230730821028, 0.4473189803940036, 0.3574685280060074, 0.34403486229203556])

    fs = 20
    fig = plt.figure(figsize=(8,5))
    fig.subplots_adjust(left=0.08,right=0.65)
    ax = plt.subplot(111, projection='polar')
    ax.set_theta_zero_location("W")
    plt.plot(theta_cfd*pi/180, cp_cfd, 'b--', linewidth=2, label='CFD')
    plt.plot(theta_ac-pi/2., cp_ac, 'g-.', linewidth=2, label='AC')
    plt.plot(theta_mod-pi/2., cp_mod, 'r-', linewidth=2, label='Wake Model')
    plt.xticks(fontsize=fs)
    plt.yticks([.5,1],fontsize=fs)
    ax.set_rlabel_position(240)
    thetaticks = np.arange(0,360,45)
    ax.set_thetagrids(thetaticks, frac=1.13)
    ax.plot([0],[0],marker=r'$\colon$',ms=65,color='k')
    plt.ylim(0,1.3)

    plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

elif method == 'velocity_overlap':
    N = 100
    xp = np.linspace(-3,10,N)
    yp = np.linspace(-5,5,N)
    X,Y = np.meshgrid(xp,yp)
    P = np.zeros((N,N))
    iter = 0
    time0 = time.time()
    for i in range(N):
        for j in range(N):
            if int_type == 'simp':
                veleffpx,veleffpy = _vawtwake.overlap(ntheta,xt,yt,diat,np.ones_like(xt)*rot,chord,B,X[i,j],Y[i,j],dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,220,200,1,0)
            elif int_type == 'gskr':
                veleffpx,veleffpy = vwm.overlap(ntheta,xt,yt,diat,np.ones_like(xt)*rot,chord,B,X[i,j],Y[i,j],dia,velf,True)

            P[i,j] = sqrt((veleffpx[0]+velf)**2 + (veleffpy[0])**2)/velf
            # P[i,j] = (veleffpx[0]+velf)/velf
            # P[i,j] = veleffpy[0]/velf
            iter += 1
            runtime = time.time()-time0
            progress_bar(float(iter)/(N*N),N*N,runtime)
            time0 = time.time()

    fs = 20
    plt.figure()
    lb = 0.15 # lower bound on velocity to display
    ub = 1.15 # upper bound on velocity to display
    ran = 32 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
    CS = plt.contourf(X/dia,Y/dia,P/velf,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel(r'$velocity mag/U_\infty$',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    plt.xlabel('$x/D$',fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    for i in range(np.size(xt)):
        circ = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='w',fill=True)
        plt.gca().add_patch(circ)

plt.show()