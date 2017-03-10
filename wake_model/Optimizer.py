from pyoptsparse import Optimization, SNOPT, pyOpt_solution
from os import path
import numpy as np
from numpy import sqrt,pi,sin,cos,fabs
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import VAWT_Wake_Model as vwm
from ACsingle import actuatorcylinder
from sys import argv

from joblib import Parallel, delayed

import _vawtwake
import _bpmvawtacoustic


def obj_func(xdict):
    global dia
    global rot
    global chord
    global twist
    global delta
    global B
    global H
    global rho
    global mu
    global velf
    global windroseDirections
    global windFrequencies

    global af_data
    global cl_data
    global cd_data
    global coef0
    global coef1
    global coef2
    global coef3
    global coef4
    global coef5
    global coef6
    global coef7
    global coef8
    global coef9

    global funcs
    global power_iso_tot
    global ntheta

    # ntheta = 72
    interp = 1 # linear airfoil interpolation
    # interp = 2 # cubic spline interpolation
    m = 220
    n = 200

    x = xdict['xvars']
    y = xdict['yvars']
    funcs = {}

    nturb = np.size(x)
    nwind = np.size(windroseDirections)

    power_turb = np.zeros(nturb)
    power_dir = np.zeros(nwind)

    rotw = np.zeros((nwind,nturb))
    k = 0
    for i in range(nwind):
        for j in range(nturb):
            rotw[i,j] = rot[k]
            k += 1

    winddir_turb = np.zeros_like(windroseDirections)
    for d in range(0, nwind):
        # Adjusting coordinate system
        winddir_turb[d] = 270. - windroseDirections[d]
        if winddir_turb[d] < 0.:
            winddir_turb[d] += 360.
        winddir_turb_rad = pi*winddir_turb[d]/180.0
        xw = x*cos(-winddir_turb_rad) - y*sin(-winddir_turb_rad)
        yw = x*sin(-winddir_turb_rad) + y*cos(-winddir_turb_rad)

        wakex,wakey = vawt_wake(xw,yw,dia,rotw[d],ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)
        res = Parallel(n_jobs=-1)(delayed(vawt_power2)(i,dia,rotw[d],ntheta,chord,B,velf,af_data,cl_data,cd_data,twist,delta,rho,mu,interp,wakex,wakey) for i in range(nturb) )
        for i in range(nturb):
            power_turb[i] = res[i][0]
            if i == 0:
                velx = res[i][1]
                vely = res[i][2]
            else:
                velx = np.append(velx,res[i][1])
                vely = np.append(vely,res[i][2])

        power_dir[d] = np.sum(power_turb)*windFrequencies[d]

        SPL_d = bpm_noise(x,y,windroseDirections[d],rotw[d],velx,vely,wakex,wakey)
        SPL_dir = np.array(SPL_d)
        if d == 0:
            SPL = SPL_dir
        else:
            SPL = np.append(SPL,SPL_dir)

    power = np.sum(power_dir)

    funcs['obj'] = (-1.*power/1e3)

    funcs['SPL'] = (SPL)/10.

    print 'Power:',power,'W (Isolated: '+str(power_iso_tot)+' W; '+str(power/power_iso_tot)+')   Max SPL:',max(SPL),'dB'

    sep = sep_func(np.append(x,y))
    funcs['sep'] = sep

    fail = False

    return funcs, fail


def vawt_wake(xw,yw,dia,rotw,ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n):

    t = np.size(xw)

    for i in range(t):
        xt = np.delete(xw,i)
        yt = np.delete(yw,i)
        diat = np.delete(dia,i)
        rott = np.delete(rotw,i)

        # wakexd,wakeyd = _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
        # wakexd,wakeyd = gskrint(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)
        wakexd,wakeyd = vwm.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,False)

        if i == 0:
            wakex = wakexd
            wakey = wakeyd
        else:
            wakex = np.append(wakex,wakexd)
            wakey = np.append(wakey,wakeyd)

    return wakex,wakey

def vawt_power(i,xw,yw,dia,rotw,ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,af_data,cl_data,cd_data,twist,delta,rho,mu,interp):

    # print i
    xt = np.delete(xw,i)
    yt = np.delete(yw,i)
    diat = np.delete(dia,i)
    rott = np.delete(rotw,i)

    # wakex,wakey = _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
    # wakex,wakey = gskrint(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)
    wakex,wakey = vwm.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,False)

    uvec,vvec,_,Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,rho,mu,interp,wakex,wakey)

    power_turb = (0.5*rho*velf**3)*(dia[i]*H)*Cp

    return power_turb,uvec,vvec,wakex,wakey

def vawt_power2(i,dia,rotw,ntheta,chord,B,velf,af_data,cl_data,cd_data,twist,delta,rho,mu,interp,wakext,wakeyt):

    # print i
    wakex = np.zeros(ntheta)
    wakey = np.zeros(ntheta)
    for j in range(ntheta):
        wakex[j] = wakext[j+ntheta*i]
        wakey[j] = wakeyt[j+ntheta*i]

    uvec,vvec,_,Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,rho,mu,interp,wakex,wakey)

    power_turb = (0.5*rho*velf**3)*(dia[i]*H)*Cp

    return power_turb,uvec,vvec


# SPL CALCULATION BASED ON BPM ACOUSTIC MODEL
def bpm_noise(turbineX,turbineY,winddir,rot,velx,vely,wakex,wakey):
    global turb_dia
    global obs
    global B
    global Hub
    global H
    global chord
    global velf
    global ntheta

    # ntheta = 72
    nobs = np.size(obs[:,0])

    noise_corr = 1.#0.8442597725840473
    nu = 1.78e-5
    c0 = 343.2
    psi = 14.0
    AR = 5.
    rad = turb_dia/2.

    div = 5

    c = np.ones(div)*chord
    c1 = c*0.5
    alpha = np.ones(div)*0.0
    high = np.linspace(0,H,div+1)

    SPL = Parallel(n_jobs=-1)(delayed(bpmnoise)(ntheta,turbineX,turbineY,obs[i],winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,velf,velx,vely,wakex,wakey) for i in range(nobs) )

    # SPL = np.zeros(nobs) #setting up vector for the constraints
    # k = 0
    # for i in range(nobs):
    #     SPL[k] = _bpmvawt.turbinepos(ntheta,turbineX,turbineY,obs[i],winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,velf,velx,vely,wakex,wakey)
    #     k += 1

    return SPL


def bpmnoise(ntheta,turbineX,turbineY,obs,winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,velf,velx,vely,wakex,wakey):
    return _bpmvawt.turbinepos(ntheta,turbineX,turbineY,obs,winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,velf,velx,vely,wakex,wakey)


def sep_func(loc):
    global turb_dia

    space = 1.65 # rotor diameters apart

    n = np.size(loc)/2
    x = loc[0:n]
    y = loc[n:]
    sep = np.zeros((n-1)*n/2)

    k = 0
    for i in range(0, n):
        for j in range(i+1, n):
            sep[k] = (x[j]-x[i])**2+(y[j]-y[i])**2
            k += 1

    return sep - (space*turb_dia)**2


## Main
if __name__ == "__main__":


    optimize = True
    # optimize = False
    plot = True
    # plot = False
    # Option to plot the velocity field
    # contour = True
    # contour = False

    global turb_dia
    global dia
    global rot
    global chord
    global twist
    global delta
    global B
    global H
    global rho
    global mu
    global velf
    global windroseDirections
    global windFrequencies

    global obs
    global Hub

    global af_data
    global cl_data
    global cd_data
    global coef0
    global coef1
    global coef2
    global coef3
    global coef4
    global coef5
    global coef6
    global coef7
    global coef8
    global coef9

    global funcs
    global power_iso_tot

    global ntheta

    # SPLlim = float(argv[1])
    # rotdir_spec = argv[2]

    SPLlim = 100.
    rotdir_spec = 'cn'


    print '\nSPL Limit:',SPLlim
    if rotdir_spec == 'cn':
        print 'Rotation Specification: Counter-rotating\n'
    elif rotdir_spec == 'co':
        print 'Rotation Specification: Co-rotating\n'

    foildata = '/Users/ning1/Documents/FLOW Lab/VAWTAC/airfoils/du06w200.dat'

    af_data,cl_data,cd_data = vwm.airfoil_data(foildata)
    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()

    # windroseDirections = np.array([270.])
    windroseDirections = np.array([270.])
    windFrequencies = np.array([1.])

    # define turbine size
    velf = 8.
    turb_dia = 1.2 # m
    tsrd = 2.625
    turb_rot = tsrd*velf/(turb_dia/2.) # rad/sec

    grid = 2.
    gridt = 24.
    nCols = 2
    add_space = 3.

    points1 = np.linspace(grid,gridt,nCols)
    points2 = np.linspace(grid+add_space,gridt+add_space,nCols)

    ypoints1, xpoints1 = np.meshgrid(points1, points1)
    ypoints2, xpoints2 = np.meshgrid(points1, points2)
    x01 = np.ndarray.flatten(xpoints1)
    y01 = np.ndarray.flatten(ypoints1)
    x02 = np.ndarray.flatten(xpoints2)
    y02 = np.ndarray.flatten(ypoints2)

    x0 = np.concatenate((x01,x02))
    y0 = np.concatenate((y01,y02))

    # x0 = np.array([0.,5.,0.,5.])
    # y0 = np.array([0.,0.,5.,5.])
    print 'x0:',x0.tolist()
    print 'y0:',y0.tolist()
    dia = np.ones_like(x0)*turb_dia

    # x0 = np.array([2.,5.,5.,8.])
    # y0 = np.array([2.,2.,8.,8.])
    # dia = np.ones_like(x0)*turb_dia
    # rot = np.array([turb_rot,-turb_rot,turb_rot,-turb_rot])

    if rotdir_spec == 'cn':
        # rot = np.zeros_like(x0)
        # for i in range(np.size(rot)):
        #     if i % 2 == 0:
        #         rot[i] = turb_rot
        #     else:
        #         rot[i] = -turb_rot
        rot1 = np.ones_like(x01)*turb_rot
        rot2 = np.ones_like(x02)*-turb_rot
        rot = np.concatenate((rot1,rot2))
    elif rotdir_spec == 'co':
        rot = np.ones_like(x0)*turb_rot
    print 'rot:',rot.tolist(),'\n'
    # rot = np.array([turb_rot,turb_rot,-turb_rot,-turb_rot])

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.128
    H = 6.1
    Hub = 2.

    rho = 1.225
    mu = 1.7894e-5

    ntheta = 72
    interp = 1 # linear airfoil interpolation
    # interp = 2 # cubic spline interpolation

    velx,vely,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[0]/2.,chord,twist,delta,B,rot[0],velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
    power_iso,cp_iso = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([turb_dia]),np.array([turb_rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)

    power_iso_tot = power_iso*np.size(x0)

    spaceval = 2.
    xlow = 0.#x0[0]-spaceval
    xupp = 30.#x0[-1]+spaceval
    ylow = 0.#y0[0]-spaceval
    yupp = 30.#y0[-1]+spaceval
    # xlow = 0.
    # xupp = 10.
    # ylow = 0.
    # yupp = 10.

    # grid_radius = 35.#sqrt((grid/2.)**2 + (grid/2.)**2) + 1.
    grid_radius = 25.
    grid_x = xupp/2.
    grid_y = yupp/2.

    nobs = 8
    obs_theta = np.linspace(-pi,pi,nobs+1)
    obs = np.zeros((nobs,3))
    for i in range(nobs):
        obs[i,0] = grid_x + grid_radius*cos(obs_theta[i])
        obs[i,1] = grid_y + grid_radius*sin(obs_theta[i])
        obs[i,2] = 2.

    # obs = np.array([[xlow-1,ylow-1,2],[xupp+1,ylow-1,2],[xupp+1,yupp+1,2],[xlow-1,yupp+1,2]])*1.


    if optimize == True:
        optProb = Optimization('VAWT_Power', obj_func)
        optProb.addObj('obj')

        n = np.size(x0)
        optProb.addVarGroup('xvars', n, 'c', lower=xlow, upper=xupp, value=x0)
        optProb.addVarGroup('yvars', n, 'c', lower=ylow, upper=yupp, value=y0)

        num_cons_sep = (n-1)*n/2
        optProb.addConGroup('sep', num_cons_sep, lower=0, upper=None)
        # num_cons_obs = np.size(obs[:,0])
        # optProb.addConGroup('SPL', num_cons_obs, lower=0, upper=SPLlim)

        opt = SNOPT()
        opt.setOption('Scale option',0)
        # if rotdir_spec == 'cn':
        #     opt.setOption('Print file','/fslhome/ebtingey/compute/optVAWT/data/SNOPT_print_SPL'+str(SPLlim)+'_turb'+str(n)+'_counterrot.out')
        #     opt.setOption('Summary file','/fslhome/ebtingey/compute/optVAWT/data/SNOPT_summary_SPL'+str(SPLlim)+'_turb'+str(n)+'_counterrot.out')
        # elif rotdir_spec == 'co':
        #     opt.setOption('Print file','/fslhome/ebtingey/compute/optVAWT/data/SNOPT_print_SPL'+str(SPLlim)+'_turb'+str(n)+'_corot.out')
        #     opt.setOption('Summary file','/fslhome/ebtingey/compute/optVAWT/data/SNOPT_summary_SPL'+str(SPLlim)+'_turb'+str(n)+'_corot.out')
        res = opt(optProb)#, sens=None)
        print res

        pow = np.array(-1*res.fStar)*1e4
        xf = res.xStar['xvars']
        yf = res.xStar['yvars']

        SPLd = funcs['SPL']
        SPLw = np.zeros((np.size(windroseDirections),np.size(obs[:,0])))
        k = 0
        for i in range(np.size(windroseDirections)):
            for j in range(np.size(obs[:,0])):
                SPLw[i,j] = SPLd[k]
                k += 1
    else:
        xf = x0
        yf = y0
        # xf = np.array([2.04299218,  10.69711973,   0.62759177,  -0.03202973])
        # yf = np.array([ -0.30868907,  5.82290592,  2.25557342,  0.38809135])

        # xf = np.array([ 2.23703622,  7.03994721,  4.07916559,  8.24416228])
        # yf = np.array([ 2.6406578,   4.31184911,  8.51934267,  5.89849381])

        # The power is: 2709.75262982 W ( 1.04235323524 )
        # The x-locations: [ 0.47307964  1.29739913  0.82146607  2.06454097]
        # The y-locations: [-0.40000015  2.2999387   6.73322378  4.12528618]
        # The power of each turbine (W): [ 666.68124805  671.4849696   671.89120384  699.69520833]

        input = {'xvars':xf,'yvars':yf}
        funcs,_ = obj_func(input)
        pow = -1*funcs['obj']*1e3
        SPLd = funcs['SPL']
        SPLw = np.zeros((np.size(windroseDirections),np.size(obs[:,0])))
        k = 0
        for i in range(np.size(windroseDirections)):
            for j in range(np.size(obs[:,0])):
                SPLw[i,j] = SPLd[k]
                k += 1

    print 'Wind Directions:',windroseDirections
    print 'The power is:',pow,'W (',pow/(power_iso*np.size(x0)),')'
    print 'The isolated power is:',power_iso*np.size(x0),'W'
    print 'The x-locations:',xf
    print 'The y-locations:',yf
    # print 'The power of each turbine (W):',power_turb
    # print 'The isolated power of one turbine is:',power_iso,'W'
    print 'SPL:',SPLw

    if rotdir_spec == 'cn':
        filename = '/fslhome/ebtingey/compute/optVAWT/data/Optrun_SPL'+str(SPLlim)+'_turb'+str(np.size(x0))+'_counterrot.txt'
    elif rotdir_spec == 'co':
        filename = '/fslhome/ebtingey/compute/optVAWT/data/Optrun_SPL'+str(SPLlim)+'_turb'+str(np.size(x0))+'_corot.txt'
    target = open(filename,'w')

    target.write('\nThe power is: '+str(pow)+' W ('+str(pow/(power_iso*np.size(x0)))+')\n')
    target.write('The isolated power is: '+str(power_iso*np.size(x0))+' W\n')
    target.write('Max SPL: '+str(np.max(SPLw))+' dB\n')
    target.write('\nWind Directions: '+str(windroseDirections.tolist())+' degrees\n')
    target.write('X-locations (initial): '+str(x0.tolist())+' m\n')
    target.write('X-locations (final): '+str(xf.tolist())+' m\n')
    target.write('Y-locations (initial): '+str(y0.tolist())+' m\n')
    target.write('Y-locations (final): '+str(yf.tolist())+' m\n')
    target.write('\nSPL: '+str(SPLw.tolist())+' dB\n')
    # target.write('The isolated power of one turbine is:',power_iso,'W')
    # target.write('RPM: '+str(rpmfw)+'\n')
    # target.write('Effective wind speeds: '+str(veleff)+' m/s\n')
    # target.write('\n-----------SNOPT-----------\n')
    # target.write(res)
    target.close()

    if plot == True:
        fs = 21
        ms = 10
        plt.figure(1,figsize=(11.5,8))
        plt.subplots_adjust(right=0.68)
        for i in range(np.size(x0)):
            if rot[i] > 0.:
                if i == 0:
                    plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='none',label='Original Turbines')
                else:
                    plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='none')
                # circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='k',fill=False)
                # plt.gca().add_patch(circ)
            elif rot[i] < 0.:
                plt.plot(x0[i],y0[i],'o',color='gray',markersize=ms,fillstyle='none')
                # circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='silver',fill=False)
                # plt.gca().add_patch(circ)
        for i in range(np.size(x0)):
            plt.plot([x0[i], xf[i]], [y0[i], yf[i]], '--k',zorder=-1)
        for i in range(np.size(xf)):
            if rot[i] > 0.:
                if i == 0:
                    plt.plot(xf[i],yf[i],'o',color='k',markersize=ms,fillstyle='full',label='Optimized Turbines')
                else:
                    plt.plot(xf[i],yf[i],'o',color='k',markersize=ms,fillstyle='full')
                # circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='k',fill=True)
                # plt.gca().add_patch(circ)
            elif rot[i] < 0.:
                plt.plot(xf[i],yf[i],'o',color='silver',markersize=ms,fillstyle='full')
                # circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='silver',fill=True)
                # plt.gca().add_patch(circ)
        plt.plot(obs[:,0],obs[:,1],'^',color='lime',markersize=ms,label='Observers')
        rect = plt.Rectangle((xlow,ylow), xupp-xlow,yupp-ylow, linestyle='dashed',linewidth=2,facecolor="#ffffff",fill=False,label='Boundaries')
        plt.gca().add_patch(rect)
        plt.annotate('N', xy=(-12, 56), xycoords='data', xytext=(-12, 49), textcoords='data', size=fs, va="center", ha="center", arrowprops=dict(arrowstyle='simple', facecolor='k'), horizontalalignment='right', verticalalignment='top',color='k')
        plt.annotate('Wind', xy=(15,-6),  xycoords='data', xytext=(-50*np.tan(45.*np.pi/180),-50), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=15), horizontalalignment='right', verticalalignment='top', fontsize=fs,color='k')
        plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xlabel('X-Position (m)',fontsize=fs)
        plt.ylabel('Y-Position (m)',fontsize=fs)
        if rotdir_spec == 'cn':
            # plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/SPLrunlayout_'+str(SPLlim)+'_'+str(nCols*nCols)+'_counterrot.pdf')
            plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/Optrunlayout_SPL'+str(SPLlim)+'_turb'+str(np.size(x0))+'_counterrot.png')
        elif rotdir_spec == 'co':
            # plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/SPLrunlayout_'+str(SPLlim)+'_'+str(nCols*nCols)+'_corot.pdf')
            plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/Optrunlayout_SPL'+str(SPLlim)+'_turb'+str(np.size(x0))+'_corot.png')

        # plt.figure(2,figsize=(11.5,8))
        # plt.subplots_adjust(right=0.68)
        # # for i in range(np.size(x0)):
        # #     if rot[i] > 0.:
        # #         if i == 0:
        # #             plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='none',label='Original Turbines')
        # #         else:
        # #             plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='none')
        # #         # circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='k',fill=False)
        # #         # plt.gca().add_patch(circ)
        # #     elif rot[i] < 0.:
        # #         plt.plot(x0[i],y0[i],'o',color='gray',markersize=ms,fillstyle='none')
        # #         # circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='silver',fill=False)
        # #         # plt.gca().add_patch(circ)
        # # for i in range(np.size(x0)):
        # #     plt.plot([x0[i], xf[i]], [y0[i], yf[i]], '--k',zorder=-1)
        # for i in range(np.size(xf)):
        #     if rot[i] > 0.:
        #         if i == 0:
        #             plt.plot(xf[i],yf[i],'o',color='k',markersize=ms,fillstyle='full',label='Optimized Turbines')
        #         else:
        #             plt.plot(xf[i],yf[i],'o',color='k',markersize=ms,fillstyle='full')
        #         # circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='k',fill=True)
        #         # plt.gca().add_patch(circ)
        #     elif rot[i] < 0.:
        #         plt.plot(xf[i],yf[i],'o',color='silver',markersize=ms,fillstyle='full')
        #         # circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='silver',fill=True)
        #         # plt.gca().add_patch(circ)
        # plt.plot(obs[:,0],obs[:,1],'^',color='lime',markersize=ms,label='Observers')
        # rect = plt.Rectangle((xlow,ylow), xupp-xlow,yupp-ylow, linestyle='dashed',linewidth=2,facecolor="#ffffff",fill=False,label='Boundaries')
        # plt.gca().add_patch(rect)
        # plt.annotate('N', xy=(-12, 56), xycoords='data', xytext=(-12, 49), textcoords='data', size=fs, va="center", ha="center", arrowprops=dict(arrowstyle='simple', facecolor='k'), horizontalalignment='right', verticalalignment='top',color='k')
        # plt.annotate('Wind', xy=(15,-6),  xycoords='data', xytext=(-50*np.tan(45.*np.pi/180),-50), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=15), horizontalalignment='right', verticalalignment='top', fontsize=fs,color='k')
        # plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
        # plt.xticks(fontsize=fs)
        # plt.yticks(fontsize=fs)
        # plt.xlabel('X-Position (m)',fontsize=fs)
        # plt.ylabel('Y-Position (m)',fontsize=fs)
        # if rotdir_spec == 'cn':
        #     plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/SPLrunlayout_'+str(SPLlim)+'_'+str(nCols*nCols)+'_counterrot2.pdf')
        #     plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/SPLrunlayout_'+str(SPLlim)+'_'+str(nCols*nCols)+'_counterrot2.png')
        # elif rotdir_spec == 'co':
        #     plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/SPLrunlayout_'+str(SPLlim)+'_'+str(nCols*nCols)+'_corot2.pdf')
        #     plt.savefig('/fslhome/ebtingey/compute/optVAWT/data/SPLrunlayout_'+str(SPLlim)+'_'+str(nCols*nCols)+'_corot2.png')


        # if contour == True:
        #
        #     N = 10
        #     xplot = np.linspace(x0[0]-spaceval-1,x0[-1]+spaceval+1,N)
        #     yplot = np.linspace(y0[0]-spaceval-1,y0[-1]+spaceval,N)
        #     [X,Y] = np.meshgrid(xplot,yplot)
        #     Vel = np.zeros((N,N))
        #
        #     k = 0
        #     for i in range(N):
        #         for j in range(N):
        #             Velx,Vely = _vawtwake.overlappoint(x0,y0,dia,rot,chord,B,X[i,j],Y[i,j],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,220,200,1,1)
        #             Vel[i,j] = sqrt((Velx+velf)**2 + (Vely)**2)
        #             k += 1
        #             print k
        #
        #     plt.figure(2)
        #     lb = 0.0 # lower bound on velocity to display
        #     ub = velf # upper bound on velocity to display
        #     ran = 200 # number of contours between the velocity bounds
        #     bounds = np.linspace(lb,ub,ran)
        #     v = np.linspace(lb,ub,4) # setting the number of tick marks on colorbar
        #     CS = plt.contourf(X,Y,Vel,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.jet) # plotting the contour plot
        #     CB = plt.colorbar(CS, ticks=v) # creating colorbar
        #     for i in range(np.size(x0)):
        #         circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='g',edgecolor='k',fill=True)
        #         plt.gca().add_patch(circ)
        #
        #     plt.xlim(x0[0]-spaceval-1,x0[-1]+spaceval+1)
        #     plt.ylim(y0[0]-spaceval-1,y0[-1]+spaceval+1)

        # plt.show()
