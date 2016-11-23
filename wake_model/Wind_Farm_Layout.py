import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt
import VAWT_Wake_Model as vwm
from AC_InducedVel import induced_vel
from os import path

import _vawtwake


method = 'power'
turbine = 'windspire'

method = 'overlap'
# turbine = 'uppsala'

ntheta = 36
interp = 1 # linear airfoil interpolation
# interp = 2 # cubic spline interpolation

basepath = path.join(path.dirname(path.realpath('__file__')), 'data/airfoils')

if turbine == 'uppsala':
    xt = np.array([0.])
    yt = np.array([0.])
    diat = np.ones_like(xt)*6.
    velf = 15.
    dia = 6.
    r = dia/2.
    rot = 4.*velf/r

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
    velf = 10.
    dia = 1.2
    r = dia/2.
    rot = 2.6*velf/r

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.128
    H = 6.1
    foildata = basepath + path.sep + 'du06w200.dat'

    rho = 1.225
    mu = 1.7894e-5



af_data,cl_data,cd_data = vwm.airfoil_data(foildata)
coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()


if method == 'power':
    plot = 'induced'
    # plot = 'Cptsr'

    if plot == 'induced':
        tsr = np.array([2.6])
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

        velx,vely = induced_vel(r,af_data,cl_data,cd_data,chord,twist,delta,B,rot[0],velf,rho,mu,ntheta,interp)

        if plot == 'induced':
            # uvec = np.array([0.459751,0.282597,0.0826669,-0.113647,-0.289426,-0.431141,-0.521327,-0.560051,-0.537055,-0.482605,-0.420738,-0.356635,-0.292542,-0.225452,-0.133013,-0.0208681,0.0851436,0.165534,0.0846263,-0.140613,-0.363099,-0.549279,-0.677296,-0.76224,-0.837259,-0.910581,-0.974004,-1.03006,-1.04345,-0.988548,-0.866052,-0.689424,-0.457304,-0.193809,0.0943889,0.377583])
            # vvec = np.array([0.385218,0.524883,0.587173,0.583378,0.523866,0.418706,0.283568,0.135876,-0.00149835,-0.109984,-0.19216,-0.256426,-0.309558,-0.358509,-0.397126,-0.407762,-0.384917,-0.342492,-0.28303,-0.195373,-0.093928,-0.00860521,0.046295,0.0793963,0.099107,0.104145,0.0896153,0.0498272,-0.0163959,-0.0914874,-0.155433,-0.189781,-0.179738,-0.119601,-0.00105869,0.183809])
            uvec = np.array([0.122831,0.0106121,-0.10353,-0.204446,-0.286066,-0.327597,-0.325175,-0.306816,-0.281749,-0.185508,-0.168015,-0.14789,-0.127679,-0.152689,-0.124647,-0.0876765,-0.0270031,0.0291537,-0.0135476,-0.162562,-0.309975,-0.400652,-0.467356,-0.437018,-0.480841,-0.520683,-0.552635,-0.697772,-0.735698,-0.758577,-0.757071,-0.689664,-0.560565,-0.391287,-0.188057,0.0422926])
            vvec = np.array([0.307522,0.344698,0.33951,0.299252,0.229196,0.140454,0.0555828,-0.0154241,-0.0840132,-0.114607,-0.113872,-0.122619,-0.124973,-0.146116,-0.185361,-0.215959,-0.231804,-0.226199,-0.208966,-0.185256,-0.156949,-0.128431,-0.103292,-0.074472,-0.0411212,-0.00933489,0.0277341,0.0483948,0.0500565,0.0507668,0.0445597,0.0395081,0.0510002,0.0855453,0.145189,0.229797])
            # uvec = np.array([-0.233381,-0.199025,-0.111714,-0.387137,-0.5955,-0.649722])
            # vvec = np.array([0.152899,-0.0370119,-0.158819,-0.136832,0.00592026,0.173844])

            plt.figure()
            plt.subplot(1,2,1)
            plt.plot(theta*180./np.pi,velx,'b',label='Model')
            plt.plot(thetaj*180./np.pi,uvec,'r--',label='AC')
            plt.xlabel('Rotation Degree')
            plt.ylabel('X-Velocity (m/s)')
            plt.xlim(0,360)
            plt.legend(loc=3)

            # plt.figure()
            plt.subplot(1,2,2)
            plt.plot(theta*180./np.pi,vely,'b',label='Model')
            plt.plot(thetaj*180./np.pi,vvec,'r--',label='AC')
            plt.xlabel('Rotation Degree')
            plt.ylabel('Y-Velocity (m/s)')
            plt.xlim(0,360)
            plt.legend(loc=3)

        power,Cp = _vawtwake.powercalc(xt,yt,diat,rot,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
        if plot == 'induced':
            print 'Power:',power,'W'
            print 'Coefficient of Power:',Cp
        elif plot == 'Cptsr':
            cp_plot[j] = Cp
            print j+1
        juliatsr = np.linspace(1.5,4.5,100)
        juliacp = np.array([0.0721066,0.0758743,0.0801562,0.0894377,0.0939292,0.0968861,0.0994119,0.106429,0.108902,0.111458,0.119002,0.12232,0.131418,0.133951,0.136469,0.138613,0.154034,0.159047,0.170879,0.175563,0.185431,0.198986,0.218908,0.225036,0.230523,0.23483,0.243059,0.245855,0.257497,0.250424,0.252098,0.261542,0.262178,0.262951,0.268604,0.277377,0.278353,0.279409,0.28745,0.295071,0.307547,0.28698,0.287694,0.30161,0.29939,0.29684,0.28662,0.292104,0.289912,0.287794,0.285776,0.283298,0.280279,0.277104,0.273811,0.270607,0.267509,0.264453,0.261436,0.258224,0.253995,0.249774,0.245508,0.241207,0.236912,0.231929,0.226299,0.220169,0.21376,0.207267,0.200711,0.194078,0.187375,0.180629,0.17372,0.166689,0.159409,0.151967,0.144442,0.136821,0.129099,0.121218,0.113105,0.104808,0.0962765,0.0875658,0.078711,0.0697048,0.0605721,0.0513054,0.0419286,0.0324134,0.0227404,0.0129355,0.00300431,-0.00705294,-0.0172319,-0.0275424,-0.0380206,-0.0486291])
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
    # plot = 'velocity'
    N = 20
    dist = 6.
    xplot = np.linspace(-dist*dia,dist*dia,N)
    yplot = np.linspace(-dist*dia,dist*dia,N)
    [X,Y] = np.meshgrid(xplot,yplot)
    P = np.zeros((N,N))

    if plot == 'power':
        velx,vely = induced_vel(r,af_data,cl_data,cd_data,chord,twist,delta,B,rot,velf,rho,mu,ntheta,interp)

        power_iso,_ = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)

        k = 0
        for i in range(N):
            for j in range(N):
                xd = np.insert(xt,0,X[i,j])
                yd = np.insert(yt,0,Y[i,j])
                power1,_ = _vawtwake.powercalc(xd,yd,np.ones_like(xd)*dia,np.ones_like(xd)*-rot,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
                xd = np.insert(xt,1,X[i,j])
                yd = np.insert(yt,1,Y[i,j])
                power2,_ = _vawtwake.powercalc(xd,yd,np.ones_like(xd)*dia,np.ones_like(xd)*rot,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)

                P[i,j] = (power1+power2)/(2*power_iso)
                k += 1
                print k,'of',N*N

        plt.figure()
        lb = 0.7 # lower bound on velocity to display
        ub = 1.3 # upper bound on velocity to display
        ran = 100 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,5) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,P,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        for i in range(np.size(xt)):
            circ = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='k',fill=True)
            plt.gca().add_patch(circ)

        plt.savefig('/Users/ning1/Documents/FLOW Lab/overlap_power_windspire_counterrot3.png')

    elif plot == 'powerpoint':
        velx,vely = induced_vel(r,af_data,cl_data,cd_data,chord,twist,delta,B,rot,velf,rho,mu,ntheta)

        power_iso,cp_iso = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([dia]),np.array([rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)

        x0 = 12.
        y0 = 18.

        xd = np.insert(xt,0,x0)
        yd = np.insert(yt,0,y0)
        power1,cp1 = _vawtwake.powercalc(xd,yd,np.ones_like(xd)*dia,np.ones_like(xd)*rot,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
        xd = np.insert(xt,1,x0)
        yd = np.insert(yt,1,y0)
        power2,cp2 = _vawtwake.powercalc(xd,yd,np.ones_like(xd)*dia,np.ones_like(xd)*rot,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)

        power_norm = (power1+power2)/(2*power_iso)

        print 'Power 1:',power1,'W'
        print 'Power 2:',power2,'W'
        print 'Power Iso:',power_iso,'W'
        print 'Total Power',power1+power2,'W'
        print 'Double Iso',2*power_iso,'W'
        print 'Normalized Power:',power_norm
        print cp_iso,cp1,cp2


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


