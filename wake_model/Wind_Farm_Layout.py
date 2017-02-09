import numpy as np
from numpy import cos,sin,pi
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rcParams
from numpy import sqrt
import VAWT_Wake_Model as vwm
from ACsingle import actuatorcylinder,actuatorcylinder2
from os import path
from scipy.interpolate import interp1d,UnivariateSpline,Akima1DInterpolator
import sys
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

from Wind_Farm_Layout_data import airfoil,ac_data

from akima import Akima, akima_interp

import _vawtwake

# Progress bar for plotting
def progress_bar(percent):
    bar_long = 40
    status = 'Working...'
    if percent == 1:
        status = 'Complete\n'
    bar_seg = int(round(bar_long*percent))
    text = '\rStatus: [{0}] {1}% {2}'.format( '='*bar_seg + ' '*(bar_long - bar_seg), percent*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

method = 'power'
turbine = 'windspire'
# turbine = 'delft'

method = 'overlap'
# turbine = 'uppsala'

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
    diat = np.ones_like(xt)*6.
    velf = 15.
    dia = 6.
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
    xt = np.array([0.])
    yt = np.array([0.])
    diat = np.ones_like(xt)*1.2
    velf = 15.
    dia = 1.2
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
    diat = np.ones_like(xt)*1.2
    velf = 10.
    dia = 1.
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
# af_data,cl_data,cd_data = airfoil()
coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()


if method == 'power':
    plot = 'induced'
    # plot = 'Cptsr'

    if plot == 'induced':
        tsr = np.array([tsrd])
    elif plot == 'Cptsr':
        tsr = np.linspace(1.5,4.5,100)
        cp_plot = np.zeros_like(tsr)

    for j in range(np.size(tsr)):
        rot = np.ones_like(xt)*(tsr[j]*velf/r)

        theta = np.zeros(ntheta)
        thetaj = np.zeros(36)
        # thetaj = np.zeros(6)
        for i in range(ntheta):
            theta[i] = (2.*np.pi/ntheta)*(i+1)-(2.*np.pi/ntheta)/2.
        for i in range(36):
            thetaj[i] = (2.*np.pi/36)*(i+1)-(2.*np.pi/36)/2.
        # for i in range(6):
        #     thetaj[i] = (2.*np.pi/6)*(i+1)-(2.*np.pi/6)/2.

        velx,vely,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot[0],velf,rho,mu,interp)
        velx2,vely2,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot[0],velf,rho,mu,interp)

        if plot == 'induced':
            uvec = np.array([0.122831,0.0106121,-0.10353,-0.204446,-0.286066,-0.327597,-0.325175,-0.306816,-0.281749,-0.185508,-0.168015,-0.14789,-0.127679,-0.152689,-0.124647,-0.0876765,-0.0270031,0.0291537,-0.0135476,-0.162562,-0.309975,-0.400652,-0.467356,-0.437018,-0.480841,-0.520683,-0.552635,-0.697772,-0.735698,-0.758577,-0.757071,-0.689664,-0.560565,-0.391287,-0.188057,0.0422926])
            vvec = np.array([0.307522,0.344698,0.33951,0.299252,0.229196,0.140454,0.0555828,-0.0154241,-0.0840132,-0.114607,-0.113872,-0.122619,-0.124973,-0.146116,-0.185361,-0.215959,-0.231804,-0.226199,-0.208966,-0.185256,-0.156949,-0.128431,-0.103292,-0.074472,-0.0411212,-0.00933489,0.0277341,0.0483948,0.0500565,0.0507668,0.0445597,0.0395081,0.0510002,0.0855453,0.145189,0.229797])

            thetaj = np.array([5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0])

            plt.figure(figsize=(15, 6))
            plt.subplot(1,2,1)
            plt.plot(theta*180./np.pi,velx,'b',label='Model')
            plt.plot(thetaj,uvec,'r--',label='AC')
            plt.xlabel('Rotation Degree')
            plt.ylabel('X-Velocity/Vinf')
            plt.xlim(0,360)
            plt.legend(loc=4)

            # plt.figure()
            plt.subplot(1,2,2)
            plt.plot(theta*180./np.pi,vely,'b',label='Model')
            plt.plot(thetaj,vvec,'r--',label='AC')
            # plt.plot((180,180),(-.3,.39),'r-')
            plt.xlabel('Rotation Degree')
            plt.ylabel('Y-Velocity/Vinf')
            plt.xlim(0,360)
            plt.legend(loc=4)
            # plt.savefig('/Users/ning1/Documents/FLOW Lab/induced_compare1_D'+str(d_space)+'.png')

        power,Cp = _vawtwake.powercalc(xt,yt,diat,rot,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
        if plot == 'induced':
            print 'Power:',power,'W'
            print 'Coefficient of Power:',Cp
        elif plot == 'Cptsr':
            # _,_,_,CP,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot[0],velf,rho,mu,interp)

            cp_plot[j] = Cp
            print j+1
        juliatsr = np.linspace(1.5,4.5,100)
        juliacp = np.array([0.0716233,0.0762639,0.0814103,0.0871416,0.0922382,0.0967217,0.101195,0.106454,0.111312,0.115768,0.12053,0.126162,0.130722,0.135023,0.140747,0.147514,0.154634,0.162332,0.171205,0.179536,0.190313,0.202898,0.214106,0.222564,0.229452,0.235479,0.240708,0.24525,0.249453,0.253586,0.25762,0.261289,0.264816,0.268426,0.271909,0.275398,0.27914,0.282687,0.286247,0.289792,0.293371,0.296686,0.299124,0.299987,0.299863,0.299155,0.29802,0.296533,0.294743,0.292687,0.290393,0.287883,0.285171,0.282271,0.27919,0.275937,0.272518,0.268937,0.265199,0.261305,0.257259,0.253063,0.248716,0.24422,0.239576,0.234785,0.229846,0.224761,0.21953,0.214153,0.20863,0.202961,0.197146,0.191185,0.185077,0.178823,0.172421,0.165872,0.159176,0.152333,0.145343,0.138206,0.130922,0.12349,0.115911,0.108185,0.100312,0.0922939,0.0841295,0.0758191,0.0673629,0.0587608,0.0500128,0.0411189,0.0320789,0.0228926,0.0135599,0.00408066,-0.00554549,-0.0153188])
    if plot == 'Cptsr':
        plt.figure()
        plt.plot(tsr,cp_plot,'bo-',label='Python')
        plt.plot(juliatsr,juliacp,'r-',label='Julia')
        plt.xlabel(r'$\lambda$')
        plt.ylabel(r'$C_p$')
        plt.legend(loc=1)

    plt.show()
































elif method == 'overlap':
    #Plotting
    plot = 'power'
    # plot = 'powerpoint'
    # plot = 'polar'
    # plot = 'sweep'
    # plot = 'velocity'
    N = 50
    nT = 72
    nR = 21
    # nT = 18
    # nR = 11
    dist = 6.
    # xplot = np.linspace(-dist*dia,dist*dia,N)
    # yplot = np.linspace(-dist*dia,dist*dia,N)
    Theta = np.linspace(-180, 180, nT)*pi/180.0
    Radii = np.linspace(2.0, 12.0, nR)
    # [X,Y] = np.meshgrid(xplot,yplot)
    # P = np.zeros((N,N))
    # P1 = np.zeros((N,N))
    # P2 = np.zeros((N,N))
    X = np.zeros((nR,nT))
    Y = np.zeros((nR,nT))
    P = np.zeros((nR,nT))
    P1 = np.zeros((nR,nT))
    P2 = np.zeros((nR,nT))

    rotplot = 'corot'
    rotplot = 'counter'

    if plot == 'power':
        velx,vely,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
        velx2,vely2,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
        # q,k,A,velx,vely = actuatorcylinder2(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
        # q2,k2,A2,velx2,vely2 = actuatorcylinder2(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))

        print q,k
        print q2,k2

        power_iso,cp_iso = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)

        rot1 = rot
        if rotplot == 'corot':
            rot2 = rot
        elif rotplot == 'counter':
            rot2 = -rot

        x3 = 10.
        y3 = 10.

        k = 0
        for i in range(nR):
            for j in range(nT):

                centerX = Radii[i]*r*cos(Theta[j])
                centerY = Radii[i]*r*sin(Theta[j])

                X[i,j] = centerX
                Y[i,j] = centerY

                # xd = np.insert(xt,1,X[i,j])
                # yd = np.insert(yt,1,Y[i,j])
                # power1,_ = _vawtwake.powercalc(xd,yd,np.array([dia,dia]),np.array([rot,-rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
                # xd = np.insert(xt,0,X[i,j])
                # yd = np.insert(yt,0,Y[i,j])
                # power2,_ = _vawtwake.powercalc(xd,yd,np.array([dia,dia]),np.array([-rot,rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx2,vely2,interp)


                # if np.fabs(centerY)/dia >= 3.:
                wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([centerX]),np.array([centerY]),np.array([dia]),np.array([rot2]),chord,B,0.,0.,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
                wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot1]),chord,B,centerX,centerY,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)

                # wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([centerX,x3]),np.array([centerY,y3]),np.array([dia,dia]),np.array([rot2,rot1]),chord,B,0.,0.,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
                # wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([0.,x3]),np.array([0.,y3]),np.array([dia,dia]),np.array([rot1,rot1]),chord,B,centerX,centerY,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
                # else:
                #     wakex1,wakey1 = vwm.overlap(ntheta,np.array([centerX]),np.array([centerY]),np.array([dia]),np.array([rot2]),chord,B,0.,0.,dia,velf)
                #     wakex2,wakey2 = vwm.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot1]),chord,B,centerX,centerY,dia,velf)

                # thetaj = np.array([5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0])
                #
                # _,_,_,Cp1,_,_,_ = _vawtwake.radialforce(velx,vely,thetaj,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,velf,wakex1,wakey1,rho,mu,interp)
                # _,_,_,Cp2,_,_,_ = _vawtwake.radialforce(velx2,vely2,thetaj,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot2,velf,wakex2,wakey2,rho,mu,interp)

                # fac = 0.6
                # wakex1 = wakex1*fac
                # wakey1 = wakey1*fac
                # wakex2 = wakex2*fac
                # wakey2 = wakey2*fac

                _,_,_,Cp1,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,velf,rho,mu,interp,wakex1,wakey1)
                _,_,_,Cp2,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot2,velf,rho,mu,interp,wakex2,wakey2)
                # q,k11,A,uvec,vvec = actuatorcylinder2(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot1,velf,rho,mu,interp,wakex1,wakey1)
                # q,k22,A,uvec,vvec = actuatorcylinder2(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot2,velf,rho,mu,interp,wakex2,wakey2)

                # P[i,j] = (power1+power2)/(2*power_iso)
                # P1[i,j] = power1/power_iso
                # P2[i,j] = power2/power_iso
                P[i,j] = (Cp1+Cp2)/(2*cp_iso)
                P1[i,j] = Cp1/cp_iso
                P2[i,j] = Cp2/cp_iso
                k += 1
                progress_bar(float(k)/(nR*nT))
                # print k,'of',N*N

        plt.figure()
        lb = 0.9 # lower bound on velocity to display
        ub = 1.1 # upper bound on velocity to display
        ran = 100 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,5) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,P,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        for i in range(np.size(xt)):
            circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
            plt.gca().add_patch(circ1)
            circ2 = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='k',fill=True)
            plt.gca().add_patch(circ2)
        if rotplot == 'counter':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_windspire_counterrot_tot.png')
        elif rotplot == 'corot':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_windspire_corot_tot.png')
        # plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_uppsala_counterrot_tot.png')

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
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xlim(-6.1,6.1)
        plt.ylim(-6.1,6.1)
        for i in range(np.size(xt)):
            circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
            plt.gca().add_patch(circ1)
            # circ2 = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='k',fill=True)
            # plt.gca().add_patch(circ2)
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
        if rotplot == 'counter':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Thesis/BYU_ME_Thesis_Template/figures/chp_power/windspire_power.pdf')
        elif rotplot == 'corot':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Thesis/BYU_ME_Thesis_Template/figures/chp_power/windspire_corotpower.pdf')

        plt.figure()
        lb = 0.9 # lower bound on velocity to display
        ub = 1.1 # upper bound on velocity to display
        ran = 100 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,5) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,P,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.parula) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        CB.ax.set_ylabel('Normalized Power',fontsize=fs)
        CB.ax.tick_params(labelsize=fs)
        plt.xlabel('$x/D$',fontsize=fs)
        plt.ylabel('$y/D$',fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        for i in range(np.size(xt)):
            circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
            plt.gca().add_patch(circ1)
            # circ2 = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='k',fill=True)
            # plt.gca().add_patch(circ2)
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
        if rotplot == 'counter':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Thesis/BYU_ME_Thesis_Template/figures/chp_power/windspire_power_p.pdf')
        elif rotplot == 'corot':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Thesis/BYU_ME_Thesis_Template/figures/chp_power/windspire_corotpower_p.pdf')

        plt.figure()
        lb = -0.2 # lower bound on velocity to display
        ub = 1.2 # upper bound on velocity to display
        ran = 8 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,8) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,P1,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.jet) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        for i in range(np.size(xt)):
            circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
            plt.gca().add_patch(circ1)
            circ2 = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='k',fill=True)
            plt.gca().add_patch(circ2)
        if rotplot == 'counter':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_windspire_counterrot_cp1.png')
        elif rotplot == 'corot':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_windspire_corot_cp1.png')
        # plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_uppsala_counterrot_cp1.png')

        plt.figure()
        lb = -0.2 # lower bound on velocity to display
        ub = 1.2 # upper bound on velocity to display
        ran = 8 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,8) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,P2,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.jet) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        for i in range(np.size(xt)):
            circ1 = plt.Circle((xt[i]/dia,yt[i]/dia),1.0,color='w',fill=True)
            plt.gca().add_patch(circ1)
            circ2 = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='k',fill=True)
            plt.gca().add_patch(circ2)
        if rotplot == 'counter':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_windspire_counterrot_cp2.png')
        elif rotplot == 'corot':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_windspire_corot_cp2.png')
        # plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/overlap_power_uppsala_counterrot_cp2.png')


        # plt.savefig('/Users/ning1/Documents/FLOW Lab/overlap_power_windspire_counterrot5.png')

    elif plot == 'powerpoint':

        x1 = 0.0
        x2 = 3.
        y1 = 0.
        y2 = 0.0
        x3 = 15.
        y3 = 15.

        # # wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2,x3]),np.array([y2,y3]),np.array([dia,dia]),np.array([-rot,rot]),chord,B,x1,y1,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        # # wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([x1,x3]),np.array([y1,y3]),np.array([dia,dia]),np.array([rot,rot]),chord,B,x2,y2,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        #
        # wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([-rot]),chord,B,x1,y1,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        # wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot]),chord,B,x2,y2,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        # print wakex1,wakey1,wakex2,wakey2

        # N = 50
        # x1 = np.linspace(-10,10.,N)
        # y1 = np.linspace(-10.0,10.,N)
        # [X,Y] = np.meshgrid(x1,y1)
        # P = np.zeros((N,N))
        # P2 = np.zeros((N,N))
        # P3 = np.zeros((N,N))

        cp_iso = 0.283118078

        # wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([-rot]),chord,B,x1,y1,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        # wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot]),chord,B,x2,y2,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)

        wakex1,wakey1 = vwm.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([-rot]),chord,B,x1,y1,dia,velf)
        wakex2,wakey2 = vwm.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot]),chord,B,x2,y2,dia,velf)

        _,_,_,Cp1,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,wakex1,wakey1)
        _,_,_,Cp2,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,wakex2,wakey2)

        print wakex1,wakey1
        print wakex2,wakey2

        print Cp1,Cp2,(Cp1+Cp2)/(2.*cp_iso)

        # for i in range(N):
        #     for j in range(N):
        #         wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2,x3]),np.array([y2,y3]),np.array([dia,dia]),np.array([-rot,rot]),chord,B,X[i,j],Y[i,j],dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        #         wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([X[i,j],x3]),np.array([Y[i,j],y3]),np.array([dia,dia]),np.array([rot,rot]),chord,B,x2,y2,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        #         # print wakex1,wakey1,wakex2,wakey2
        #
        #         uvec1,vvec1,_,cp1w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,wakex1,wakey1)
        #         uvec2,vvec2,_,cp2w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,wakex2,wakey2)
        #         P[i,j] = cp1w
        #         P2[i,j] = cp2w
        #         P3[i,j] = (cp1w+cp2w)/(2*cp_iso)
        #
        #         print i,j
        # plt.figure()
        # CS = plt.contourf(X,Y,P,100,cmap=plt.cm.parula)
        # CB = plt.colorbar(CS)
        #
        # plt.figure()
        # CS = plt.contourf(X,Y,P2,100,cmap=plt.cm.parula)
        # CB = plt.colorbar(CS)
        #
        # plt.figure()
        # CS = plt.contourf(X,Y,P3,100,cmap=plt.cm.parula)
        # CB = plt.colorbar(CS)
        # plt.show()

        # wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([-rot]),chord,B,x1,y1,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        # wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot]),chord,B,x2,y2,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        #
        # uvec1,vvec1,_,cp1w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,wakex1,wakey1)
        # uvec2,vvec2,_,cp2w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,wakex2,wakey2)
        #
        # print cp1w
        # print cp2w


        # alpha = np.linspace(-180.,180.,4000)
        # clp = np.zeros_like(alpha)
        # cdp = np.zeros_like(alpha)
        #
        # # print af_data
        # # print cl_data
        # # print cd_data
        #
        # # cl = Akima(af_data,cl_data)
        # # cd = Akima(af_data,cd_data)
        # # print cl
        # for i in range(4000):
        #     # clp[i] = _vawtwake.interpolate(af_data,cl_data,alpha[i])
        #     # cdp[i] = _vawtwake.interpolate(af_data,cd_data,alpha[i])
        #
        #     clp[i] = _vawtwake.splineint(af_data,cl_data,alpha[i])
        #     cdp[i] = _vawtwake.splineint(af_data,cd_data,alpha[i])
        #
        #
        # plt.figure()
        # plt.plot(alpha,clp)
        # plt.plot(af_data,cl_data,'b-')
        #
        # plt.figure()
        # plt.plot(alpha,cdp)
        # plt.plot(af_data,cd_data,'b-')
        # plt.show()
        # print velf,rot
        # velx,vely,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
        # velx2,vely2,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
        # power_iso,cp_iso = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
        # print velx.tolist()
        # print vely.tolist()
        #
        # print 'Cp',cp_iso
        # u [0.130053,0.0172343,-0.100554,-0.206186,-0.285872,-0.331198,-0.339091,-0.311855,-0.257649,-0.198857,-0.17342,-0.153704,-0.144239,-0.15129,-0.130448,-0.0854709,-0.0249459,0.0319706,-0.0110843,-0.163447,-0.310304,-0.414919,-0.467859,-0.468693,-0.493591,-0.531481,-0.574983,-0.664414,-0.742709,-0.779144,-0.763469,-0.691039,-0.562558,-0.385705,-0.175121,0.0539504]
        # v [0.307092,0.348335,0.345646,0.304658,0.234351,0.146861,0.0558436,-0.0246643,-0.0807229,-0.104103,-0.111845,-0.119991,-0.129019,-0.153226,-0.191164,-0.222886,-0.236586,-0.230104,-0.212217,-0.187211,-0.156823,-0.127404,-0.101221,-0.0729178,-0.0406995,-0.00845896,0.0243131,0.0494383,0.057264,0.0517809,0.0422214,0.0387825,0.0504856,0.0839942,0.143082,0.227117]
        # Power Calculations Here
        # Vn [-0.207433,-0.0731861,0.0668608,0.205752,0.339254,0.463615,0.575387,0.671081,0.746561,0.789021,0.769467,0.716294,0.626996,0.491781,0.342162,0.184493,0.0238377,-0.139286,-0.297599,-0.397348,-0.433608,-0.439952,-0.447855,-0.477045,-0.476163,-0.454744,-0.421281,-0.338617,-0.263345,-0.222048,-0.217972,-0.245892,-0.292262,-0.335737,-0.351701,-0.31811]
        # Vt [3.77752,3.69773,3.58625,3.45,3.29568,3.12891,2.95492,2.77928,2.60928,2.45147,2.30303,2.15859,2.02847,1.91652,1.80306,1.70196,1.62194,1.5769,1.65834,1.86541,2.0662,2.21881,2.32029,2.37999,2.44787,2.51191,2.56374,2.605,2.63628,2.67141,2.72608,2.81605,2.95437,3.14624,3.38474,3.65515]
        # W [3.78321,3.69845,3.58687,3.45613,3.31309,3.16307,3.01042,2.85915,2.71399,2.57532,2.42817,2.27433,2.12316,1.97861,1.83524,1.71193,1.62211,1.58304,1.68483,1.90726,2.11121,2.262,2.36312,2.42732,2.49375,2.55274,2.59812,2.62691,2.6494,2.68062,2.73478,2.82676,2.96879,3.16411,3.40296,3.66896]
        # phi [-0.0548573,-0.0197896,0.0186415,0.0595677,0.102578,0.147101,0.192315,0.236923,0.278672,0.311386,0.32245,0.320401,0.299783,0.251181,0.187538,0.107979,0.014696,-0.0881003,-0.177566,-0.209872,-0.206856,-0.195744,-0.190672,-0.197819,-0.192122,-0.179095,-0.162867,-0.129263,-0.0995623,-0.0829294,-0.0797883,-0.0870973,-0.0986043,-0.106308,-0.103536,-0.0868121]
        # alpha [-0.0548573,-0.0197896,0.0186415,0.0595677,0.102578,0.147101,0.192315,0.236923,0.278672,0.311386,0.32245,0.320401,0.299783,0.251181,0.187538,0.107979,0.014696,-0.0881003,-0.177566,-0.209872,-0.206856,-0.195744,-0.190672,-0.197819,-0.192122,-0.179095,-0.162867,-0.129263,-0.0995623,-0.0829294,-0.0797883,-0.0870973,-0.0986043,-0.106308,-0.103536,-0.0868121]
        # cn [-0.243912,-0.0475339,0.179997,0.416552,0.636101,0.811945,0.918305,0.934152,0.851448,0.722397,0.703965,0.70644,0.75775,0.917858,0.911058,0.660586,0.156551,-0.407525,-0.741689,-0.845499,-0.835826,-0.80025,-0.783989,-0.806892,-0.788643,-0.746658,-0.693301,-0.575459,-0.457693,-0.383887,-0.369205,-0.40299,-0.453612,-0.485887,-0.474415,-0.401696]
        # ct [-0.00414974,-0.0158909,-0.0143025,0.00496253,0.0413601,0.0842859,0.118537,0.129982,0.105678,0.0554206,0.0436964,0.0456307,0.0716141,0.126453,0.115743,0.046591,-0.0152082,0.0161143,0.0830747,0.107179,0.104857,0.0964645,0.0926981,0.0980149,0.0937719,0.0841931,0.0723306,0.0471744,0.0244407,0.0125274,0.0104123,0.015409,0.0237285,0.0295173,0.0274193,0.0152093]
        # q [-0.177797,-0.0331142,0.117942,0.253408,0.355601,0.413728,0.42385,0.388922,0.319406,0.24401,0.211388,0.186103,0.173965,0.183007,0.156279,0.0985987,0.0209791,-0.0520126,-0.107228,-0.156639,-0.189735,-0.208536,-0.222973,-0.242126,-0.249779,-0.247802,-0.238347,-0.202244,-0.163621,-0.14049,-0.140632,-0.163999,-0.203617,-0.247746,-0.279797,-0.275394]
        # qdyn [8.76651,8.37811,7.88022,7.31621,6.72315,6.12807,5.55087,5.00704,4.5115,4.06226,3.61132,3.16821,2.76104,2.39788,2.06296,1.79506,1.61164,1.53494,1.73868,2.22804,2.73003,3.13395,3.42041,3.60879,3.80901,3.99134,4.13451,4.22666,4.29933,4.40126,4.58091,4.89423,5.39842,6.13208,7.09285,8.24503]
        # Rp [0.273697,0.0509753,-0.181557,-0.390091,-0.547405,-0.636885,-0.652466,-0.598699,-0.491688,-0.375624,-0.325407,-0.286484,-0.267798,-0.281717,-0.240573,-0.151781,-0.0322948,0.0800671,0.165064,0.241128,0.292075,0.321017,0.34324,0.372724,0.384506,0.381462,0.366907,0.311331,0.251875,0.216267,0.216486,0.252457,0.313445,0.381376,0.430714,0.423936]
        # Tp [-0.00465647,-0.0170414,-0.0144265,0.00464729,0.035593,0.0661133,0.0842216,0.0833057,0.0610262,0.028817,0.0201986,0.0185047,0.0253093,0.0388121,0.0305628,0.0107051,-0.00313729,0.003166,0.0184884,0.0305663,0.0366416,0.0386963,0.0405843,0.0452756,0.0457188,0.0430136,0.0382786,0.025522,0.0134501,0.00705743,0.00610531,0.00965314,0.0163964,0.0231683,0.0248935,0.0160513]
        # Zp [0.0,0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
        # Cp 0.2831196606299424

        # x0 = 0.*dia
        # y0 = 2.*dia
        #
        # xd = np.insert(xt,1,x0)
        # yd = np.insert(yt,1,y0)
        # power1,cp1 = _vawtwake.powercalc(xd,yd,np.array([dia,dia]),np.array([rot,-rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx2,vely2,interp)
        # xd = np.insert(xt,0,x0)
        # yd = np.insert(yt,0,y0)
        # power2,cp2 = _vawtwake.powercalc(xd,yd,np.array([dia,dia]),np.array([-rot,rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
        #
        # power_norm = (power1+power2)/(2*power_iso)
        #
        # print 'Power 1:',power1,'W'
        # print 'Power 2:',power2,'W'
        # print 'Power Iso:',power_iso,'W'
        # print 'Total Power',power1+power2,'W'
        # print 'Double Iso',2*power_iso,'W'
        # print 'Normalized Power:',power_norm
        # print cp_iso,cp1,cp2
        # thetaj = np.array([5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0])
        #
        # velxac = np.array([0.0679113,-0.0506645,-0.165054,-0.259275,-0.323055,-0.351683,-0.34497,-0.307675,-0.246814,-0.20387,-0.179822,-0.161136,-0.161364,-0.155433,-0.125498,-0.0739787,-0.0113531,0.0418549,0.00165605,-0.140511,-0.288316,-0.400703,-0.46464,-0.482198,-0.486472,-0.517742,-0.554531,-0.616727,-0.701399,-0.749291,-0.752127,-0.703099,-0.599818,-0.444353,-0.244792,-0.0117645])
        # velyac = np.array([0.230144,0.255665,0.239088,0.189377,0.117213,0.0345924,-0.0460268,-0.112563,-0.150632,-0.16313,-0.169162,-0.174522,-0.187745,-0.217349,-0.252865,-0.27929,-0.286706,-0.277196,-0.259792,-0.234537,-0.202028,-0.16985,-0.141694,-0.114553,-0.0835294,-0.0505825,-0.0183377,0.0106717,0.0259373,0.0253131,0.0176296,0.0118827,0.0172345,0.0414225,0.0899208,0.163296])
        #
        # velxac2 = np.array([0.0418549,-0.0113531,-0.0739787,-0.125498,-0.155433,-0.161364,-0.161136,-0.179822,-0.20387,-0.246814,-0.307676,-0.34497,-0.351683,-0.323055,-0.259275,-0.165054,-0.0506645,0.0679113,-0.0117645,-0.244792,-0.444353,-0.599818,-0.703099,-0.752127,-0.749291,-0.701399,-0.616727,-0.554531,-0.517742,-0.486472,-0.482198,-0.46464,-0.400703,-0.288316,-0.140511,0.00165605])
        # velyac2 = np.array([0.277196,0.286706,0.27929,0.252865,0.217349,0.187745,0.174522,0.169162,0.16313,0.150632,0.112563,0.0460268,-0.0345924,-0.117213,-0.189377,-0.239088,-0.255665,-0.230144,-0.163296,-0.0899208,-0.0414225,-0.0172345,-0.0118827,-0.0176296,-0.0253131,-0.0259374,-0.0106718,0.0183378,0.0505825,0.0835294,0.114553,0.141694,0.16985,0.202028,0.234537,0.259792])
        #
        #
        # wakex,wakey = _vawtwake.overlap(ntheta,np.array([x0]),np.array([y0]),np.array([dia]),np.array([-rot]),chord,B,0.,0.,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        # wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot]),chord,B,x0,y0,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        #
        # velxac_wind = np.array([0.00674861,0.00385562,0.00100442,-0.00179478,-0.00452194,-0.00714264,-0.00960157,-0.0118134,-0.0136509,-0.0149294,-0.0153915,-0.0146978,-0.0124394,-0.00820023,-0.00170269,0.00696293,0.0171138,0.0274666,0.036516,0.0431439,0.047,0.0484014,0.0479645,0.0462974,0.0438687,0.0409984,0.0378915,0.0346747,0.0314247,0.0281858,0.0249822,0.0218253,0.0187187,0.0156613,0.0126496,0.00967974])
        # velyac_wind = np.array([0.079387,0.08046,0.0818534,0.0835976,0.0857331,0.088311,0.0913922,0.0950455,0.0993411,0.104337,0.110055,0.116432,0.123256,0.130082,0.136161,0.140487,0.142031,0.140189,0.135166,0.127954,0.119858,0.111955,0.104864,0.0988119,0.0937938,0.0897071,0.0864237,0.0838219,0.0817974,0.0802649,0.0791571,0.078422,0.0780214,0.077929,0.0781286,0.0786137])
        #
        #
        # # print wakex/velf,wakey/velf,wakex2/velf,wakey2/velf
        #
        # plt.figure()
        # plt.subplot(1,2,1)
        # plt.plot(thetaj,wakex2)
        # plt.plot(thetaj,velxac_wind)
        # plt.subplot(1,2,2)
        # plt.plot(thetaj,wakey2)
        # plt.plot(thetaj,velyac_wind)

    elif plot == 'polar':
        cpdata = np.array([0.5137292963669421, 0.9407530012258202, 1.0583347246905448, 1.100368119669435, 1.127448529615236, 1.1014079312070053, 1.0569344196425345, 0.940520231361647, 0.5098077736897938, 0.7368409217354737, 1.0070289977142592, 1.0697784260663252, 1.0882329965142843, 1.0679385516564928, 1.0097941954631366, 0.7333565294554432, 0.5137292963669421])
        thetadata = np.array([270.1114379016532, 246.7988999869089, 225.19062704244018, 201.86037312349106, 179.88586513201912, 516.7008614759567, 494.80957252795724, 471.6647172260273, 449.6585021210285, 426.6668191463685, 404.64546467066145, 381.7268571885939, 359.8858651320191, 336.61717949856603, 314.96577551188346, 292.03610306036995, 270.1114379016532])
        theta = np.linspace(-pi,pi,72)
        cptheta = np.array([0.445982, 0.51201, 0.615172, 0.712552, 0.776262, 0.814086, 0.881157, 0.93065, 0.995502, 1.00445, 1.01043, 1.01692, 1.02396, 1.03157, 1.03922, 1.04553, 1.05007, 1.05269, 1.05322, 1.05163, 1.04803, 1.04259, 1.03549, 1.02769, 1.02037, 1.01361, 1.00738, 1.00164, 0.971863, 0.898155, 0.849936, 0.792628, 0.738928, 0.658594, 0.557209, 0.469759, 0.493928, 0.573753, 0.617474, 0.682425, 0.741807, 0.816748, 0.907066, 0.978855, 0.991511, 0.996756, 1.00322, 1.01079, 1.01907, 1.02736, 1.03491, 1.04111, 1.04541, 1.0474, 1.04672, 1.04353, 1.03821, 1.03127, 1.02326, 1.01488, 1.00689, 0.99984, 0.993982, 0.991411, 0.94704, 0.869073, 0.790213, 0.714075, 0.653963, 0.605099, 0.532343, 0.445982])

        thetamod = np.linspace(-pi,pi,72)

        cpmod = np.zeros_like(thetamod)
        cp_iso = 0.283118078

        # for i in range(np.size(thetamod)):
        #     x1 = 0.
        #     y1 = 0.
        #     x2 = 2.*cos(thetamod[i])
        #     y2 = 2.*sin(thetamod[i])
        #
        #     wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([x2]),np.array([y2]),np.array([dia]),np.array([-rot]),chord,B,x1,y1,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        #     wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([x1]),np.array([y1]),np.array([dia]),np.array([rot]),chord,B,x2,y2,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        #
        #     fac = 0.8
        #     wakex1 = wakex1*fac
        #     wakey1 = wakey1*fac
        #     wakex2 = wakex2*fac
        #     wakey2 = wakey2*fac
        #
        #     uvec1,vvec1,_,cp1w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,wakex1,wakey1)
        #     uvec2,vvec2,_,cp2w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,wakex2,wakey2)
        #     cpmod[i] = (cp1w+cp2w)/(2*cp_iso)
        #     print i
        #
        # print cpmod.tolist()

        cpmod = np.array([0.34403486229934377, 0.34862509462216756, 0.39894157768009036, 0.5863430562386762, 0.7129970832610787, 0.7855921597522427, 0.8470482539630853, 0.9013953506018821, 0.9601455893616664, 1.0194288948705703, 1.033693122129043, 1.0413605061884208, 1.0506381788442232, 1.060058261050193, 1.0689654414305627, 1.0774951768601828, 1.08516060298998, 1.0896502927606777, 1.090592059615426, 1.087828165379764, 1.0816003402926102, 1.073275197216215, 1.0645617181205838, 1.0554416781533134, 1.0457958116850024, 1.0374472069223197, 1.0289739659759387, 0.9793388535624352, 0.9328965276589404, 0.8687096863324896, 0.8054556953197295, 0.7613318567938111, 0.6397558675986322, 0.4771682732414631, 0.36328963392066166, 0.3434522228186334, 0.3462820609189695, 0.3878422593824303, 0.5172162997342789, 0.6005684340530796, 0.6994852713795089, 0.7595821027780201, 0.8536741564896726, 0.9534604104216948, 0.9971821817553583, 1.0184854079417303, 1.0262396891831798, 1.0331792074981805, 1.0403024456084002, 1.0471348168407635, 1.053121978219379, 1.0580477306582485, 1.0611763655548085, 1.0624636404913026, 1.0620357786382415, 1.0598530722064783, 1.0557547705350605, 1.0502530689227327, 1.0437918778931803, 1.0367521209320387, 1.0296423808559239, 1.0225430345387512, 1.0213049314743428, 0.9988919125859831, 0.9049627047469962, 0.8271948352315188, 0.7381963424220676, 0.6500036047389534, 0.5814230730140153, 0.4473189805021868, 0.35746852803756585, 0.3440348622992079])

        fs = 23
        fig = plt.figure(figsize=(8,5))
        fig.subplots_adjust(left=0.08,right=.59)
        ax = plt.subplot(111, projection='polar')
        ax.set_theta_zero_location("W")
        plt.plot(thetadata*pi/180, cpdata, 'b--', linewidth=2, label='CFD')
        plt.plot(theta-pi/2., cptheta, 'g-.', linewidth=2, label='AC')
        plt.plot(thetamod-pi/2., cpmod, 'r-', linewidth=2, label='Wake Model')
        plt.xticks(fontsize=fs)
        plt.yticks([.5,1],fontsize=fs)
        ax.set_rlabel_position(240)
        # ax.set_rlabel_
        thetaticks = np.arange(0,360,45)
        ax.set_thetagrids(thetaticks, frac=1.13)
        ax.plot([0],[0],marker=r'$\colon$',ms=65,color='k')
        plt.ylim(0,1.3)
        # plt.hold(True)
        # ax = plt.subplot(111)
        # arc = mpatches.Arc([0,0],dia/2.,dia/2.,angle=90,theta1=0,theta2=350,capstyle='round',linestyle='-',lw=10,color='k')
        # ax.add_patch(arc)
        # ax.add_line(arc)

        plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)



    elif plot == 'sweep':
        velx,vely,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
        velx2,vely2,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
        powerw_iso,cpw_iso = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)

        thetavec = np.zeros(ntheta)
        for i in range(ntheta):
            thetavec[i] = (2.0*np.pi/ntheta)*(i+1)-(2.0*np.pi/ntheta)/2.0

        Cp_iso = 0.2831196606299424
        down = np.array([-6,-4,-2,0,2,4,6])
        lat = np.array([6,4,2,0,-2,-4,-6])
        for i in range(7):
            for j in range(7):
                velxac1,velyac1,Cp1,velxac2,velyac2,Cp2,velxac_wind1,velyac_wind1,velxac_wind2,velyac_wind2 = ac_data(down[i],lat[j])

                wakex1,wakey1 = _vawtwake.overlap(ntheta,np.array([dia*down[i]]),np.array([dia*lat[j]]),np.array([dia]),np.array([-rot]),chord,B,0.,0.,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
                wakex2,wakey2 = _vawtwake.overlap(ntheta,np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot]),chord,B,dia*down[i],dia*lat[j],dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
                x0 = dia*down[i]
                y0 = dia*lat[j]

                # xd = np.insert(xt,0,x0)
                # yd = np.insert(yt,0,y0)
                # power2,cp2w = _vawtwake.powercalc(xd,yd,np.array([dia,dia]),np.array([-rot,rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx2,vely2,interp)
                # xd = np.insert(xt,1,x0)
                # yd = np.insert(yt,1,y0)
                # power1,cp1w = _vawtwake.powercalc(xd,yd,np.array([dia,dia]),np.array([rot,-rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
                #
                # power_norm = (power1+power2)/(2*powerw_iso)
                #
                # print 'Cp1', cp1w
                # print 'Cp2', cp2w

                # _,_,_,cp1w,_,_,_ = _vawtwake.radialforce(velx+wakex1/velf,vely+wakey1/velf,thetavec,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,np.zeros(ntheta),np.zeros(ntheta),rho,mu,interp)
                # _,_,_,cp2w,_,_,_ = _vawtwake.radialforce(velx2+wakex2/velf,vely2+wakey2/velf,thetavec,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,np.zeros(ntheta),np.zeros(ntheta),rho,mu,interp)

                _,_,_,CpAC1,_,_,_ = _vawtwake.radialforce(velxac1,velyac1,thetavec,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,np.zeros(ntheta),np.zeros(ntheta),rho,mu,interp)
                _,_,_,CpAC2,_,_,_ = _vawtwake.radialforce(velxac2,velyac2,thetavec,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,np.zeros(ntheta),np.zeros(ntheta),rho,mu,interp)

                _,_,_,CpACw1,_,_,_ = _vawtwake.radialforce(velx+velxac_wind1,vely+velyac_wind1,thetavec,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,np.zeros(ntheta),np.zeros(ntheta),rho,mu,interp)
                _,_,_,CpACw2,_,_,_ = _vawtwake.radialforce(velx2+velxac_wind2,vely2+velyac_wind2,thetavec,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,np.zeros(ntheta),np.zeros(ntheta),rho,mu,interp)

                # fac = 0.8
                # wakex1 = wakex1*fac
                # wakey1 = wakey1*fac
                # wakex2 = wakex2*fac
                # wakey2 = wakey2*fac

                uvec1,vvec1,_,cp1w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,rot,velf,rho,mu,interp,wakex1,wakey1)
                uvec2,vvec2,_,cp2w,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,r,chord,twist,delta,B,-rot,velf,rho,mu,interp,wakex2,wakey2)

                thetaj = np.array([5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0])

                plt.figure((i+1))
                plt.subplot(7,4,(1)+(4*j))
                # plt.plot(thetaj,velx+wakex1/velf,'b-',label='Model')
                plt.plot(thetaj,uvec1+wakex1/velf,'b-',label='Model')
                plt.plot(thetaj,velxac1,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(2)+(4*j))
                # plt.plot(thetaj,vely+wakey1/velf,'b-',label='Model')
                plt.plot(thetaj,vvec1+wakey1/velf,'b-',label='Model')
                plt.plot(thetaj,velyac1,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(3)+(4*j))
                # plt.plot(thetaj,velx2+wakex2/velf,'b-',label='Model')
                plt.plot(thetaj,uvec2+wakex2/velf,'b-',label='Model')
                plt.plot(thetaj,velxac2,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(4)+(4*j))
                # plt.plot(thetaj,vely2+wakey2/velf,'b-',label='Model')
                plt.plot(thetaj,vvec2+wakey2/velf,'b-',label='Model')
                plt.plot(thetaj,velyac2,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcomp_'+str(down[i])+'.png')

                plt.figure((400+i))
                plt.subplot(7,4,(1)+(4*j))
                plt.plot(thetaj,wakex1/velf,'b-',label='Model')
                plt.plot(thetaj,velxac_wind1,'r--',label='AC')
                plt.xlim(0,360)
                if j == 3 and i < 3:
                    plt.ylim(-1,-.5)
                elif j == 3 and i == 3:
                    plt.ylim(-1,1)
                else:
                    plt.ylim(-.15,.15)
                # plt.legend(loc=3)
                plt.subplot(7,4,(2)+(4*j))
                plt.plot(thetaj,wakey1/velf,'b-',label='Model')
                plt.plot(thetaj,velyac_wind1,'r--',label='AC')
                plt.xlim(0,360)
                if j == 3 and i == 3:
                    plt.ylim(-1,1)
                else:
                    plt.ylim(-.15,.15)
                # plt.legend(loc=3)
                plt.subplot(7,4,(3)+(4*j))
                plt.plot(thetaj,wakex2/velf,'b-',label='Model')
                plt.plot(thetaj,velxac_wind2,'r--',label='AC')
                plt.xlim(0,360)
                if j == 3 and i > 3:
                    plt.ylim(-1,-.5)
                elif j == 3 and i == 3:
                    plt.ylim(-1,1)
                else:
                    plt.ylim(-.15,.15)
                # plt.legend(loc=3)
                plt.subplot(7,4,(4)+(4*j))
                plt.plot(thetaj,wakey2/velf,'b-',label='Model')
                plt.plot(thetaj,velyac_wind2,'r--',label='AC')
                plt.xlim(0,360)
                if j == 3 and i == 3:
                    plt.ylim(-1,1)
                else:
                    plt.ylim(-.15,.15)
                # plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcompwake_'+str(down[i])+'.png')

                plt.figure((500+i))
                plt.subplot(7,4,(1)+(4*j))
                plt.plot(thetaj,velx,'b-',label='Model')
                plt.plot(thetaj,velxac1-velxac_wind1,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(2)+(4*j))
                plt.plot(thetaj,vely,'b-',label='Model')
                plt.plot(thetaj,velyac1-velyac_wind1,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(3)+(4*j))
                plt.plot(thetaj,velx2,'b-',label='Model')
                plt.plot(thetaj,velxac2-velxac_wind2,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(4)+(4*j))
                plt.plot(thetaj,vely2,'b-',label='Model')
                plt.plot(thetaj,velyac2-velyac_wind2,'r--',label='AC')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcompind_'+str(down[i])+'.png')

                plt.figure((600+i))
                plt.subplot(7,4,(1)+(4*j))
                plt.plot(thetaj,wakex1/velf-velxac_wind1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(2)+(4*j))
                plt.plot(thetaj,wakey1/velf-velyac_wind1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(3)+(4*j))
                plt.plot(thetaj,wakex2/velf-velxac_wind2,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(4)+(4*j))
                plt.plot(thetaj,wakey2/velf-velyac_wind2,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcompdiffwake_'+str(down[i])+'.png')

                plt.figure((700+i))
                plt.subplot(7,4,(1)+(4*j))
                plt.plot(thetaj,velx+wakex1/velf-velxac1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(2)+(4*j))
                plt.plot(thetaj,vely+wakey1/velf-velyac1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(3)+(4*j))
                plt.plot(thetaj,velx2+wakex2/velf-velxac2,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(4)+(4*j))
                plt.plot(thetaj,vely2+wakey2/velf-velyac2,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcompdiffvel_'+str(down[i])+'.png')

                plt.figure((800+i))
                plt.subplot(7,4,(1)+(4*j))
                plt.plot(thetaj,velx-(velxac1-velxac_wind1),'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(2)+(4*j))
                plt.plot(thetaj,vely-(velyac1-velyac_wind1),'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(3)+(4*j))
                plt.plot(thetaj,velx2-(velxac2-velxac_wind2),'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(4)+(4*j))
                plt.plot(thetaj,vely2-(velyac2-velyac_wind2),'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcompdiffind_'+str(down[i])+'.png')

                plt.figure((900+i))
                plt.subplot(7,4,(1)+(4*j))
                plt.plot(thetaj,velx-velxac1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(2)+(4*j))
                plt.plot(thetaj,vely-velxac1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(3)+(4*j))
                plt.plot(thetaj,velx2-velxac1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.subplot(7,4,(4)+(4*j))
                plt.plot(thetaj,vely2-velxac1,'b-',label='Model')
                plt.xlim(0,360)
                # plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcompdiffind2_'+str(down[i])+'.png')

                plt.figure(100+i)
                plt.plot(lat[j],cp1w,'bo',label='Model')
                plt.plot(lat[j],Cp1,'ro',label='AC')
                plt.plot(lat[j],CpAC1,'.',color='lime',label='ACvel')
                plt.plot(lat[j],CpACw1,'.',color='c',label='ACvelw')
                plt.plot((-6,6),(cpw_iso,cpw_iso),'b--')
                plt.plot((-6,6),(Cp_iso,Cp_iso),'r:')
                plt.ylim(-0.05,0.35)
                if j == 0:
                    plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcomp_cp1_'+str(down[i])+'.png')
                plt.figure(200+i)
                plt.plot(lat[j],cp2w,'bo',label='Model')
                plt.plot(lat[j],Cp2,'ro',label='AC')
                plt.plot(lat[j],CpAC2,'.',color='lime',label='ACvel')
                plt.plot(lat[j],CpACw2,'.',color='c',label='ACvelw')
                plt.plot((-6,6),(cpw_iso,cpw_iso),'b--')
                plt.plot((-6,6),(Cp_iso,Cp_iso),'r:')
                plt.ylim(-0.05,0.35)
                if j == 0:
                    plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcomp_cp2_'+str(down[i])+'.png')
                plt.figure(300+i)
                plt.plot(lat[j],(cp1w+cp2w)/(2*cpw_iso),'bo',label='Model')
                plt.plot(lat[j],(Cp1+Cp2)/(2*Cp_iso),'ro',label='AC')
                plt.plot(lat[j],(CpAC1+CpAC2)/(2*Cp_iso),'.',color='lime',label='ACvel')
                plt.plot(lat[j],(CpACw1+CpACw2)/(2*Cp_iso),'.',color='c',label='ACvelw')
                # plt.plot(lat[j],(cp1w+cp2w),'bo',label='Model')
                # plt.plot(lat[j],(Cp1+Cp2),'ro',label='AC')
                # plt.plot(lat[j],(CpAC1+CpAC2),'.',color='lime',label='ACvel')
                plt.plot((-6,6),(1,1),'k:')
                plt.ylim(0.4,1.1)
                # plt.ylim(0,0.7)
                if j == 0:
                    plt.legend(loc=3)
                plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Sweep/fullcomp_cptot_'+str(down[i])+'.png')

                print i, j

        plt.close()


    elif plot == 'velocity':
        k = 0
        for i in range(N):
            for j in range(N):
                veleffpx,veleffpy = _vawtwake.overlappoint(xt,yt,diat,np.ones_like(xt)*rot,chord,B,X[i,j],Y[i,j],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,220,200,1)

                P[i,j] = sqrt((veleffpx+velf)**2 + (veleffpy)**2)
                # P[i,j] = veleffpx
                k += 1
                print k,'of',N*N

        plt.figure()
        lb = 0.15 # lower bound on velocity to display
        ub = 1.15 # upper bound on velocity to display
        ran = 32 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,7) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,P/velf,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
        CB = plt.colorbar(CS) # creating colorbar
        CB.ax.set_ylabel(r'$vel_{mag}/U_\infty$')
        plt.xlabel('x/D')
        plt.ylabel('y/D')
        for i in range(np.size(xt)):
            circ = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='w',fill=True)
            plt.gca().add_patch(circ)

        plt.savefig('/Users/ning1/Documents/FLOW Lab/overlap_velocity.png')

plt.show()


