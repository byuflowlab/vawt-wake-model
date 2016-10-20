import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt
from VAWT_Wake_Model import velocity_field

import _vawtwake


def overlap(xt,yt,diat,tsr,solidity,x0,y0,dia,velf,hub,out):
    """
    Calculating an effective velocity based on overlapping of wakes
    
    Parameters
    ----------
    xt : array
        downstream locations of turbines producing wakes (m)
    yt : array
        lateral locations of turbines producing wakes (m)
    diat : array
        diameters of each of the turbines producing wakes (m)
    tsr : array
        tip-speed ratio of each turbine
    solidity : array
        solidities of each turbine
    x0 : float
        downstream location of turbine being analyzed (m)
    y0 : float
        lateral location of turbine being analyzed (m)
    dia : float
        diameter of turbine being analyzed (m)
    velf : float
        free stream velocity (m/s)
    hub : bool
        option of using the hub velocity (True) or an average of
        velocities across the turbine's rotation (False)
    out : bool
        option of outputting the the information of turbines completed
        in the overlap calculation (True)
    
    Returns
    ----------
    veleff : float
        effective velocity in front of turbine (m/s)
    """
    
    N = np.size(xt) # number of turbines that overlap given turbine
    inte = 0. # initializing integration value
    r1 = y0 - dia/2. # lower bound of integration
    r2 = y0 + dia/2. # upper bound of integration
    n = 3
    if hub == True:
        y = np.array([y0]) # using hub velocity for effective velocity
    else:
        y = np.linspace(r1,r2,n) # using n velocity values for effective velocity
    
    vel = np.zeros_like(y)
    
    for i in range(N):
        for j in range(np.size(y)):
            vel[j] = velocity_field(xt[i],yt[i],x0,y[j],velf,diat[i],tsr[i],solidity[i])
        
        inte_n = np.average(vel)
        
        inte = inte + (inte_n)**2
        if out == True:
            print 'Turbine '+str(i+1)+' of '+str(N)+' completed'
    
    veleff = velf*(sqrt(inte))
    
    return veleff
    

def powerval(Cp,dens,vel,dia):
    """
    Calculating turbine power of a 2D VAWT
    
    Parameters
    ----------
    Cp : float
        turbine power coefficient
    dens : float
        air density of wind farm (kg/m^3)
    vel : float
        effective velocity in front of turbine (m/s)
    dia : float
        turbine diameter (m)
    
    Returns
    ----------
    power : float
        turbine power (kJ)
    """
    
    power = 0.5*Cp*dens*dia*vel**3
    
    return power/1e3


if __name__ == "__main__":
    
    xt = np.array([0.,0.])
    yt = np.array([0.,7.])
    diat = np.array([6.,6.])
    velf = 15.
    dia = 6.
    rot = np.array([4.*velf/(dia/2.),4.*velf/(dia/2.)])
    solidity = np.array([0.25,0.25])
    x0 = 12.
    y0 = 0.


    Omega = 4.*velf/(dia/2.)

    r = dia/2.
    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.25
    H = 5.

    plot = 'power'
    plot = 'velocity'

    rho = 1.225
    mu = 1.7894e-5

    af_data = np.array([])
    cl_data = np.array([])
    cd_data = np.array([])

    # f = open('/Users/ning1/Documents/FLOW Lab/VAWTAC/airfoils/NACA_0021.dat', 'r')
    f = open('/Users/ning1/Documents/FLOW Lab/VAWTAC/airfoils/du06w200.dat', 'r')
    for i in range(13):
        f.readline()
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

    # print af_data
    # print cl_data
    # print cd_data

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


    xt = np.array([0.])
    yt = np.array([0.])
    diat = np.array([0.6])
    velf = 1.
    dia = 1.2
    r = dia/2.
    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.128
    H = 6.1
    solidity = np.array([B*chord/r])

    tsr = np.linspace(2.6,2.6,1)
    cp_plot = np.zeros_like(tsr)

    for j in range(1):
        rot = np.array([tsr[j]*velf/(dia/2.)])

        x0 = 0.
        y0 = 0.

        Omega = tsr[j]*velf/(dia/2.)

        rho = 1.225
        mu = 1.7894e-5

        theta = np.zeros(36)
        xd = np.zeros(36)
        yd = np.zeros(36)
        for i in range(36):
            theta[i] = (2.*np.pi/36.)*(i+1)-(2.*np.pi/36.)/2.
            xd[i] = x0 + np.cos(theta[i])*r
            yd[i] = y0 + np.sin(theta[i])*r

        velx = np.zeros_like(xd)
        vely = np.zeros_like(yd)

        for i in range(36):
            velx[i],vely[i] = _vawtwake.vel_field(xt,yt,xd[i],yd[i],dia,rot,chord,B,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,velf,110.,100.,1)
            velx[i] = velx[i]# + Omega*(dia/2.)*np.cos(theta[i])/velf
            vely[i] = vely[i]# + Omega*(dia/2.)*np.sin(theta[i])/velf


        # print velx
        # print vely
        # print theta*180./np.pi

        uvec = np.array([0.459751,0.282597,0.0826669,-0.113647,-0.289426,-0.431141,-0.521327,-0.560051,-0.537055,-0.482605,-0.420738,-0.356635,-0.292542,-0.225452,-0.133013,-0.0208681,0.0851436,0.165534,0.0846263,-0.140613,-0.363099,-0.549279,-0.677296,-0.76224,-0.837259,-0.910581,-0.974004,-1.03006,-1.04345,-0.988548,-0.866052,-0.689424,-0.457304,-0.193809,0.0943889,0.377583])
        vvec = np.array([0.385218,0.524883,0.587173,0.583378,0.523866,0.418706,0.283568,0.135876,-0.00149835,-0.109984,-0.19216,-0.256426,-0.309558,-0.358509,-0.397126,-0.407762,-0.384917,-0.342492,-0.28303,-0.195373,-0.093928,-0.00860521,0.046295,0.0793963,0.099107,0.104145,0.0896153,0.0498272,-0.0163959,-0.0914874,-0.155433,-0.189781,-0.179738,-0.119601,-0.00105869,0.183809])

        plt.figure()
        plt.plot(theta*180./np.pi,velx,'b',label='Model')
        plt.plot(theta*180./np.pi,uvec,'r',label='AC')
        plt.xlabel('Rotation Degree')
        plt.ylabel('X-Velocity (m/s)')
        plt.xlim(0,360)
        plt.legend(loc=3)

        plt.figure()
        plt.plot(theta*180./np.pi,vely,'b',label='Model')
        plt.plot(theta*180./np.pi,vvec,'r',label='AC')
        plt.xlabel('Rotation Degree')
        plt.ylabel('Y-Velocity (m/s)')
        plt.xlim(0,360)
        plt.legend(loc=3)




        power,Cp = _vawtwake.powercalc(xt,yt,diat,rot,x0,y0,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,Omega,H,rho,mu)
        cp_plot[j] = Cp
        print j
        # print power,Cp

    # q,ka,CTo,CPo,Rp,Tp,Zp = _vawtwake.radialforce(uvec,vvec,thetavec,af_data,cl_data,cd_data,r,chord,twist,delta,B,Omega,velf,rho,mu)


    # plt.figure()
    # plt.plot(tsr,cp_plot,'bo-')
    # plt.xlabel('TSR')
    # plt.ylabel('Cp')
    plt.show()


    # power_iso = powerval(Cp,dens,velf,dia)
    # print 'Isolated turbine power: ',power_iso,'kJ'
    #
    # veleff = overlap(xt,yt,diat,tsr,solidity,x0,y0,dia,velf,True,False)
    #
    # power = powerval(Cp,dens,veleff,dia)
    # print 'Analyzed turbine power: ',power,'kJ'




    # #Plotting
    # plot = 'power'
    # plot = 'velocity'
    # N = 100
    # xplot = np.linspace(-2.*dia,4.*dia,N)
    # yplot = np.linspace(-2.*dia,4.*dia,N)
    # [X,Y] = np.meshgrid(xplot,yplot)
    # P = np.zeros((N,N))
    #
    # k = 0
    # for i in range(N):
    #     for j in range(N):
    #         veleffpx,veleffpy = _vawtwake.overlappoint(xt,yt,diat,rot,chord,B,X[i,j],Y[i,j],dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,110,100,1)
    #         if plot == 'power':
    #             P[i,j] = (Cp,dens,veleffp,dia)/power_iso
    #         elif plot == 'velocity':
    #             P[i,j] = sqrt((veleffpx+velf)**2 + (veleffpy)**2)
    #             # P[i,j] = veleffpx
    #         k += 1
    #         print k,'of',N*N
    #
    # plt.figure()
    # lb = 0.15 # lower bound on velocity to display
    # ub = 1.15 # upper bound on velocity to display
    # ran = 32 # number of contours between the velocity bounds
    # bounds = np.linspace(lb,ub,ran)
    # v = np.linspace(lb,ub,7) # setting the number of tick marks on colorbar
    # if plot == 'power':
    #     CS = plt.contourf(X,Y,P,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.parula) # plotting the contour plot
    #     CB = plt.colorbar(CS, ticks=v) # creating colorbar
    # elif plot == 'velocity':
    #     CS = plt.contourf(X/dia,Y/dia,P/velf,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
    #     CB = plt.colorbar(CS) # creating colorbar
    #     CB.ax.set_ylabel(r'$vel_{mag}/U_\infty$')
    #     plt.xlabel('x/D')
    #     plt.ylabel('y/D')
    # for i in range(np.size(xt)):
    #     circ = plt.Circle((xt[i]/dia,yt[i]/dia),0.5,color='w',fill=True)
    #     plt.gca().add_patch(circ)
    #
    # plt.savefig('/Users/ning1/Documents/FLOW Lab/overlap5.png')
    # plt.show()


