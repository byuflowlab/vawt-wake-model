import numpy as np
from numpy import cos,sin,pi,fabs
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

# Option to choose the type of turbine
turbine = 'windspire'   # 1.2 kW Mariah Windspire
# turbine = 'uppsala'   # 12 kW H-Rotor Uppsala University VAWT
# turbine = 'delft'     # Delft University of Technology H-Rotor VAWT

# Option to choose the type of calculation to run
method = 'power_curve'              # create a TSR/Cp power curve
# method = 'coupled_vawt_power'     # create power plot of two coupled turbines
# method = 'polar_power'            # create a power plot of two coupled turbines with respect to wind direction
# method = 'velocity_overlap'       # create a velocity plot of overlapping turbine wakes

# Option to choose the method of integration
int_type = 'simp'       # use Simpson's Rule integration (Fortran code)
m = 220                 # number of divisions in the downstream direction (for Simpson's Rule)
n = 200                 # number of divisions in the lateral direction (for Simpson's Rule)
# int_type = 'gskr'     # use 21 Point Gauss-Kronrod Rule Quadrature integration

# Option to choose rotation direction of two turbines
rotdir = 'corot'        # co-rotating (both counter-clockwise from above) 
rotdir = 'counter'      # counter-rotating (+ counter-clockwise, - clockwise from above)

ntheta = 36     # number of actuator cylinder points around the flight path of the blades

# Option to choose airfoil data (lift and drag) interpolation method
interp = 1      # linear airfoil interpolation
interp = 2      # cubic spline interpolation

airfoilpath = path.join(path.dirname(path.realpath('__file__')), 'data/airfoils') # airfoil directory

if turbine == 'windspire':
    xt = np.array([0.])
    yt = np.array([0.])
    dia = 1.2
    diat = np.ones_like(xt)*dia
    Vinf = 8.#1.
    r = dia/2.
    tsrd = 2.625
    rot = tsrd*Vinf/r

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.128
    H = 6.1
    foildata = airfoilpath + path.sep + 'du06w200.dat'

    rho = 1.225
    mu = 1.7894e-5

elif turbine == 'uppsala':
    xt = np.array([0.])
    yt = np.array([0.])
    dia = 6.
    diat = np.ones_like(xt)*dia
    Vinf = 15.
    r = dia/2.
    tsrd = 4.
    rot = tsrd*Vinf/r

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.25
    H = 5.
    foildata = airfoilpath + path.sep + 'NACA_0021.dat'

    rho = 1.225
    mu = 1.7894e-5

elif turbine == 'delft':
    xt = np.array([0.])
    yt = np.array([0.])
    dia = 1.
    diat = np.ones_like(xt)*dia
    Vinf = 9.3
    r = dia/2.
    tsrd = 4.5
    rot = tsrd*Vinf/r

    twist = 0.0
    delta = 0.0
    B = 2
    chord = 0.06
    H = 1.
    foildata = airfoilpath + path.sep + 'NACA_0018.dat'

    rho = 1.225
    mu = 1.7894e-5

# Reading in airfoil data
af_data,cl_data,cd_data = vwm.airfoil_data(foildata)

# Reading in VAWT wake model coefficients
coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()

if method != 'power_curve':
    Cp_iso,Tpp,Vnp,Vtp = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,fabs(rot),Vinf,rho,interp,np.zeros(ntheta),np.zeros(ntheta))
    Cpp = (fabs(rot)*B/(2.*pi*rho*Vinf**3))*Tpp
    _,Tpn,Vnn,Vtn = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-fabs(rot),Vinf,rho,interp,np.zeros(ntheta),np.zeros(ntheta))
    Cpn = (fabs(rot)*B/(2.*pi*rho*Vinf**3))*Tpn

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

if method == 'power_curve':
    print 'Power Curve (TSR/Cp) Using Actuator Cylinder Model'

    N = 100
    low = 1.#1.5
    high = 4.5

    tsr = np.linspace(low,high,N)
    cp_plot = np.zeros_like(tsr)

    time0 = time.time()
    for j in range(np.size(tsr)):
        rot = np.ones_like(xt)*(tsr[j]*Vinf/r)

        _,_,Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot[0],Vinf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))

        cp_plot[j] = Cp
        runtime = time.time()-time0
        progress_bar(float(j+1)/N,N,runtime)
        time0 = time.time()

    if turbine == 'windspire':
        julia_tsr = np.linspace(1.5,4.5,100)
        julia_cp = np.array([0.0716233,0.0762639,0.0814103,0.0871416,0.0922382,0.0967217,0.101195,0.106454,0.111312,0.115768,0.12053,0.126162,0.130722,0.135023,0.140747,0.147514,0.154634,0.162332,0.171205,0.179536,0.190313,0.202898,0.214106,0.222564,0.229452,0.235479,0.240708,0.24525,0.249453,0.253586,0.25762,0.261289,0.264816,0.268426,0.271909,0.275398,0.27914,0.282687,0.286247,0.289792,0.293371,0.296686,0.299124,0.299987,0.299863,0.299155,0.29802,0.296533,0.294743,0.292687,0.290393,0.287883,0.285171,0.282271,0.27919,0.275937,0.272518,0.268937,0.265199,0.261305,0.257259,0.253063,0.248716,0.24422,0.239576,0.234785,0.229846,0.224761,0.21953,0.214153,0.20863,0.202961,0.197146,0.191185,0.185077,0.178823,0.172421,0.165872,0.159176,0.152333,0.145343,0.138206,0.130922,0.12349,0.115911,0.108185,0.100312,0.0922939,0.0841295,0.0758191,0.0673629,0.0587608,0.0500128,0.0411189,0.0320789,0.0228926,0.0135599,0.00408066,-0.00554549,-0.0153188])

    plt.figure()
    if turbine == 'windspire':
        plt.plot(tsr,cp_plot,'bo-')#,label='Python')
        # plt.plot(julia_tsr,julia_cp,'r-',label='Julia (for validation)')
        # plt.legend(loc=1)
    else:
        plt.plot(tsr,cp_plot,'bo-')
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$C_p$')

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

elif method == 'coupled_vawt_power':
    print 'Coupled VAWT Power Plot'

    nT = 72
    nR = 21
    Theta = np.linspace(-180,180,nT)*pi/180.0
    Radii = np.linspace(2.0,12.0,nR)
    X = np.zeros((nR,nT))
    Y = np.zeros((nR,nT))
    P = np.zeros((nR,nT))
    P1 = np.zeros((nR,nT))
    P2 = np.zeros((nR,nT))
    xt = np.array([0.]) # adjusting turbine x-position for plot
    yt = np.array([0.]) # adjusting turbine y-position for plot

    thetavec = np.zeros(ntheta)
    for i in range(ntheta):
        thetavec[i] = (2.*pi/ntheta)*(i+1)-(2.*pi/ntheta)/2.

    rot1 = rot
    if rotdir == 'corot':
        rot2 = rot
    elif rotdir == 'counter':
        rot2 = -rot

    iter = 0
    time0 = time.time()
    for i in range(nR):
        for j in range(nT):

            centerX = Radii[i]*r*cos(Theta[j])
            centerY = Radii[i]*r*sin(Theta[j])

            X[i,j] = centerX
            Y[i,j] = centerY

            if int_type == 'simp':
                wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([centerX]),np.array([centerY]),np.array([dia]),np.array([rot2]),chord,B,0.,0.,dia,Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
                wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot1]),chord,B,centerX,centerY,dia,Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
            elif int_type == 'gskr':
                wakex1,wakey1 = vwm.overlap(ntheta,np.array([centerX]),np.array([centerY]),np.array([dia]),np.array([rot2]),chord,B,0.,0.,dia,Vinf,False)
                wakex2,wakey2 = vwm.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot1]),chord,B,centerX,centerY,dia,Vinf,False)

            _,Cp1 = _vawtwake.powercalc(thetavec,Vinf,wakex1,wakey1,Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,rot1,r,1.,af_data,cl_data,cd_data,twist,rho,interp)
            _,Cp2 = _vawtwake.powercalc(thetavec,Vinf,wakex2,wakey2,Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,rot2,r,1.,af_data,cl_data,cd_data,twist,rho,interp)

            P[i,j] = (Cp1+Cp2)/(2*Cp_iso)
            P1[i,j] = Cp1/Cp_iso
            P2[i,j] = Cp2/Cp_iso
            iter += 1
            runtime = time.time()-time0
            progress_bar(float(iter)/(nR*nT),nR*nT,runtime)
            time0 = time.time()

    fs = 15

    # Combined Power
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

    # Center Turbine Power
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

    # Swept Turbine Power
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
    print 'Coupled VAWT Power Plot with Repsect to Wind Direction'
    # Option to recalculate the Cp values (True) or use the values from a previous calculation (False)
    calculate = True
    # calculate = False

    thetavec = np.zeros(ntheta)
    for i in range(ntheta):
        thetavec[i] = (2.*pi/ntheta)*(i+1)-(2.*pi/ntheta)/2.

    rot1 = rot
    if rotdir == 'corot':
        rot2 = rot
    elif rotdir == 'counter':
        rot2 = -rot

    theta_mod = np.linspace(-pi,pi,72)
    cp_mod = np.zeros_like(theta_mod)

    if turbine == 'windspire':
        theta_cfd = np.array([270.1114379016532, 246.7988999869089, 225.19062704244018, 201.86037312349106, 179.88586513201912, 516.7008614759567, 494.80957252795724, 471.6647172260273, 449.6585021210285, 426.6668191463685, 404.64546467066145, 381.7268571885939, 359.8858651320191, 336.61717949856603, 314.96577551188346, 292.03610306036995, 270.1114379016532])
        cp_cfd = np.array([0.5137292963669421, 0.9407530012258202, 1.0583347246905448, 1.100368119669435, 1.127448529615236, 1.1014079312070053, 1.0569344196425345, 0.940520231361647, 0.5098077736897938, 0.7368409217354737, 1.0070289977142592, 1.0697784260663252, 1.0882329965142843, 1.0679385516564928, 1.0097941954631366, 0.7333565294554432, 0.5137292963669421])
        theta_ac = np.linspace(-pi,pi,72)
        cp_ac = np.array([0.445982, 0.51201, 0.615172, 0.712552, 0.776262, 0.814086, 0.881157, 0.93065, 0.995502, 1.00445, 1.01043, 1.01692, 1.02396, 1.03157, 1.03922, 1.04553, 1.05007, 1.05269, 1.05322, 1.05163, 1.04803, 1.04259, 1.03549, 1.02769, 1.02037, 1.01361, 1.00738, 1.00164, 0.971863, 0.898155, 0.849936, 0.792628, 0.738928, 0.658594, 0.557209, 0.469759, 0.493928, 0.573753, 0.617474, 0.682425, 0.741807, 0.816748, 0.907066, 0.978855, 0.991511, 0.996756, 1.00322, 1.01079, 1.01907, 1.02736, 1.03491, 1.04111, 1.04541, 1.0474, 1.04672, 1.04353, 1.03821, 1.03127, 1.02326, 1.01488, 1.00689, 0.99984, 0.993982, 0.991411, 0.94704, 0.869073, 0.790213, 0.714075, 0.653963, 0.605099, 0.532343, 0.445982])

    if calculate == True:
        print '---Calculating values for polar plot---'
        time0 = time.time()
        for i in range(np.size(theta_mod)):
            dia_sep = 2.0 # diameter separation
            x1 = 0.
            y1 = 0.
            x2 = dia_sep*cos(theta_mod[i])
            y2 = dia_sep*sin(theta_mod[i])

            if int_type == 'simp':
                wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([rot2]),chord,B,x1,y1,dia,Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
                wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot1]),chord,B,x2,y2,dia,Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
            elif int_type == 'gskr':
                wakex1,wakey1 = vwm.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([rot2]),chord,B,x1,y1,dia,Vinf,False)
                wakex2,wakey2 = vwm.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot1]),chord,B,x2,y2,dia,Vinf,False)

            # Cp1,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,Vinf,rho,interp,wakex1,wakey1)
            # Cp2,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot2,Vinf,rho,interp,wakex2,wakey2)

            _,Cp1 = _vawtwake.powercalc(thetavec,Vinf,wakex1,wakey1,Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,rot1,r,1.,af_data,cl_data,cd_data,twist,rho,interp)
            _,Cp2 = _vawtwake.powercalc(thetavec,Vinf,wakex2,wakey2,Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,rot2,r,1.,af_data,cl_data,cd_data,twist,rho,interp)


            cp_mod[i] = (Cp1+Cp2)/(2*Cp_iso)

            runtime = time.time()-time0
            progress_bar(float(i+1)/(np.size(theta_mod)),np.size(theta_mod),runtime)
            time0 = time.time()

        print '\nCopy to "elif calculate == False:" section for faster creation of same plot:'
        print 'cp_mod = np.array(',cp_mod.tolist(),')'

    elif calculate == False:
        cp_mod = np.array( [0.6719408687671133, 0.6576406044815359, 0.6123650054132985, 0.6958393975693847, 0.7844590781339361, 0.8097784706905345, 0.8397044941851077, 0.8823882305819531, 0.9532489314083461, 1.0268736094392068, 1.035289020619549, 1.0376661222064318, 1.0412420471255437, 1.046410969308487, 1.0530315965029144, 1.060796881204869, 1.0689012788277716, 1.075811137885089, 1.0812423449575215, 1.085147001182814, 1.0874704998611848, 1.0881197121340338, 1.0870474361845492, 1.0843688114346157, 1.0803666510381609, 1.0754484569776432, 1.0669429622941666, 1.019921112753567, 0.9638102786194865, 0.9211272828173458, 0.9008038004972125, 0.9158954363033139, 0.8579093650961004, 0.7810008106205151, 0.685916333889993, 0.6750458098133678, 0.6661972372520877, 0.6383931923110998, 0.6283230087142456, 0.7372745794154836, 0.8172175319473947, 0.8102673548654481, 0.8450644642886022, 0.9229158163976937, 0.9724897493257025, 1.0330501237167238, 1.0364808824362037, 1.0392531278801003, 1.0436299279418968, 1.0495594778796045, 1.0567931537478852, 1.0649242288002347, 1.0725399468069692, 1.0787119551176314, 1.083389909027947, 1.0865160019194, 1.088010244440272, 1.0877964331110452, 1.085894063987666, 1.0825124884046584, 1.077982016869052, 1.072423071223065, 1.0722473384278324, 1.0227242537817838, 0.946049548847798, 0.9401457992773887, 0.9156684335963895, 0.8918150032376986, 0.8569452687425566, 0.7131310834895461, 0.6795926153464056, 0.6719408687671135] ) # Windspire turbine, counter-rotating, 2D separation, Simpson integration, 36 actuator cylinder points, cubic spline airfoil interpolation

    fs = 15
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
    print 'Velocity Plot with Multiple VAWT Wake Overlap'

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
                veleffpx,veleffpy = _vawtwake.overlap(ntheta,xt,yt,diat,np.ones_like(xt)*rot,chord,B,X[i,j],Y[i,j],dia,Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,220,200,1,0)
            elif int_type == 'gskr':
                veleffpx,veleffpy = vwm.overlap(ntheta,xt,yt,diat,np.ones_like(xt)*rot,chord,B,X[i,j],Y[i,j],dia,Vinf,True)

            P[i,j] = sqrt((veleffpx[0]+Vinf)**2 + (veleffpy[0])**2)/Vinf
            # P[i,j] = (veleffpx[0]+Vinf)/Vinf
            # P[i,j] = veleffpy[0]/Vinf
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
    CS = plt.contourf(X/dia,Y/dia,P/Vinf,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
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

plt.show() # display all created figures
