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

    ntheta = 36
    interp = 1 # linear airfoil interpolation
    # interp = 2 # cubic spline interpolation
    m = 220
    n = 200

    x = xdict['xvars']
    y = xdict['yvars']
    funcs = {}

    nturb = np.size(x)
    nwind = np.size(windroseDirections)

    # power_turb = np.zeros((nwind,nturb))
    power_turb = np.zeros(nturb)
    power_dir = np.zeros(nwind)

    # uvec = np.zeros(nwind*nturb)
    # vvec = np.zeros(nwind*nturb)
    # wakex = np.zeros(nwind*nturb)
    # wakey = np.zeros(nwind*nturb)

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

        # res = Parallel(n_jobs=-1)(delayed(vawt_power)(i,xw,yw,dia,rotw[d],ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,af_data,cl_data,cd_data,twist,delta,rho,mu,interp) for i in range(nturb) )
        #
        # for i in range(nturb):
        #     power_turb[i] = res[i][0]
        #     if i == 0:
        #         velx = res[i][1]
        #         vely = res[i][2]
        #         wakex = res[i][3]
        #         wakey = res[i][4]
        #     else:
        #         velx = np.append(velx,res[i][1])
        #         vely = np.append(vely,res[i][2])
        #         wakex = np.append(wakex,res[i][3])
        #         wakey = np.append(wakey,res[i][4])

        # for i in range(nturb):
        #     power_turb[i],uvec,vvec,wakexx,wakeyy = vawt_power(i,xw,yw,dia,rotw[d],ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,af_data,cl_data,cd_data,twist,delta,rho,mu,interp)
        #     if i == 0:
        #         velx = uvec
        #         vely = vvec
        #         wakex = wakexx
        #         wakey = wakeyy
        #     else:
        #         velx = np.append(velx,uvec)
        #         vely = np.append(vely,vvec)
        #         wakex = np.append(wakex,wakexx)
        #         wakey = np.append(wakey,wakeyy)
        #     # print i

        wakex,wakey = vawt_wake(xw,yw,dia,rotw[d],ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)
        res = res = Parallel(n_jobs=-1)(delayed(vawt_power2)(i,dia,rotw[d],ntheta,chord,B,velf,af_data,cl_data,cd_data,twist,delta,rho,mu,interp,wakex,wakey) for i in range(nturb) )
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

    print 'Power:',power,'W (Isolated: '+str(power_iso_tot)+' W; '+str(power/power_iso_tot)+')','Max SPL:',max(SPL)

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

        # wakexd,wakeyd = _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
        wakexd,wakeyd = gskrint(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)

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

    # wakex,wakey = _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
    wakex,wakey = gskrint(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)

    uvec,vvec,_,Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,rho,mu,interp,wakex,wakey)

    power_turb = (0.5*rho*velf**3)*(dia[i]*H)*Cp

    return power_turb,uvec,vvec,wakex,wakey

def vawt_power2(i,dia,rotw,ntheta,chord,B,velf,af_data,cl_data,cd_data,twist,delta,rho,mu,interp,wakext,wakeyt):
    global uvec_iso
    global vvec_iso

    # print i
    wakex = np.zeros(ntheta)
    wakey = np.zeros(ntheta)
    for j in range(ntheta):
        wakex[j] = wakext[j+ntheta*i]
        wakey[j] = wakeyt[j+ntheta*i]

    # uvec,vvec,_,Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,rho,mu,interp,wakex,wakey)
    # uvec,vvec,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
    theta = np.array([5.0,15.0,25.0,35.0,45.0,55.0,65.0,75.0,85.0,95.0,105.0,115.0,125.0,135.0,145.0,155.0,165.0,175.0,185.0,195.0,205.0,215.0,225.0,235.0,245.0,255.0,265.0,275.0,285.0,295.0,305.0,315.0,325.0,335.0,345.0,355.0])*np.pi/180.
    # _,_,_,Cp,_,_,_ = _vawtwake.radialforce(uvec,vvec,theta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,wakex,wakey,rho,mu,interp)
    print 'here'
    print 'uvec = np.array(',uvec_iso.tolist(),')'
    print 'vvec = np.array(',vvec_iso.tolist(),')'
    print 'wakex = np.array(',wakex.tolist(),')'
    print 'wakey = np.array(',wakey.tolist(),')'
    _,_,Ct,Cp,Rp,Tp,Zp = _vawtwake.radialforce(uvec_iso,vvec_iso,theta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,wakex,wakey,rho,mu,interp)

    # print 'ct',Ct
    # print 'cp',Cp
    # print 'rp',Rp
    # print 'tp = np.array(',Tp.tolist(),')',rotw[i]
    # print 'zp',Zp

    power_turb = (0.5*rho*velf**3)*(dia[i]*H)*Cp

    return power_turb,uvec_iso,vvec_iso



def overlap(ntheta,xt,yt,diat,rott,chord,B,xw,yw,dia,rotw,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n):
    return _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw,yw,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)

def gskrint(p,xt,yt,diat,rott,chord,B,x0,y0,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n):

    split = True
    split = False

    if split == True:
        t = np.size(xt)
        wakex = np.zeros(t*p)
        wakey = np.zeros(t*p)
        xd = np.zeros(p)
        yd = np.zeros(p)

        for i in range(p):
            theta = (2.0*pi/p)*i-(2.0*pi/p)/2.0
            xd[i] = x0 - sin(theta)*(dia/2.0)
            yd[i] = y0 + cos(theta)*(dia/2.0)

        for i in range(t):
            if fabs(yt[i] - y0)/diat[i] >= 3.:
                for j in range(p):
                    wakex[j+i*p],wakey[j+i*p] = _vawtwake.vel_field(xt[i],yt[i],xd[j],yd[j],diat[i],rott[i],chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
            else:
                for j in range(p):
                    wake = vwm.velocity_field(xt[i],yt[i],xd[j],yd[j],velf,diat[i],rott[i],chord,B,param=None,veltype='ind',integration='gskr')
                    wakex[j+i*p] = wake[0]
                    wakey[j+i*p] = wake[1]

        velx,vely = _vawtwake.overlap2(t,p,velf,wakex,wakey)


    elif split == False:
        velx,vely = vwm.overlap(p,xt,yt,diat,rott,chord,B,x0,y0,dia,velf)

    return velx,vely

def gskrintpoint(xt,yt,diat,rott,chord,B,x0,y0,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n):

    split = True
    split = False

    p = 1

    if split == True:
        t = np.size(xt)
        wakex = np.zeros(t*p)
        wakey = np.zeros(t*p)
        xd = np.zeros(p)
        yd = np.zeros(p)

        for i in range(p):
            # theta = (2.0*pi/p)*i-(2.0*pi/p)/2.0
            # xd[i] = x0 - sin(theta)*(dia/2.0)
            # yd[i] = y0 + cos(theta)*(dia/2.0)
            xd[i] = x0
            yd[i] = y0

        for i in range(t):
            # if fabs(yt[i] - y0)/diat[i] >= 3.:
            for j in range(p):
                wakex[j+i*p],wakey[j+i*p] = _vawtwake.vel_field(xt[i],yt[i],xd[j],yd[j],diat[i],rott[i],chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
            # else:
            #     for j in range(p):
            #         wake = vwm.velocity_field(xt[i],yt[i],xd[j],yd[j],velf,diat[i],rott[i],chord,B,param=None,veltype='ind',integration='gskr')
            #         wakex[j+i*p] = wake[0]
            #         wakey[j+i*p] = wake[1]

        velx,vely = _vawtwake.overlap2(t,p,velf,wakex,wakey)


    elif split == False:
        velx,vely = vwm.overlappoint(xt,yt,diat,rott,chord,B,x0,y0,dia,velf,-1)

    return velx,vely



# SPL CALCULATION BASED ON BPM ACOUSTIC MODEL
def bpm_noise(turbineX,turbineY,winddir,rot,velx,vely,wakex,wakey):
    global turb_dia
    global obs
    global B
    global Hub
    global H
    global chord
    global velf

    ntheta = 36
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
            sep[k] = sqrt((x[j]-x[i])**2+(y[j]-y[i])**2)
            k += 1

    return sep - space*turb_dia


## Main
if __name__ == "__main__":


    plot = True
    # plot = False
    splcon = False
    # splcon = True
    # Option to plot the velocity field
    contour = True
    contour = False

    rotdir_spec = 'cn'
    # rotdir_spec = 'co'

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
    global uvec_iso
    global vvec_iso

    # xdn = float(argv[1])
    # plotpoint = int(argv[2])

    foildata = '/Users/ning1/Documents/FLOW Lab/VAWTWakeModel/wake_model/data/airfoils/du06w200.dat'

    af_data,cl_data,cd_data = vwm.airfoil_data(foildata)
    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()

    windroseDirections = np.array([270.])
    windFrequencies = np.array([1.])

    # define turbine size
    velf = 8.
    turb_dia = 1.2 # m
    tsrd = 2.625
    turb_rot = tsrd*velf/(turb_dia/2.) # rad/sec
    print turb_rot


    grid = 2.
    gridt = 35.
    nCols = 3
    add_space = 3.

    # points1 = np.linspace(grid,gridt,nCols)
    # points2 = np.linspace(grid+add_space,gridt+add_space,nCols)
    #
    # ypoints1, xpoints1 = np.meshgrid(points1, points1)
    # ypoints2, xpoints2 = np.meshgrid(points1, points2)
    # x01 = np.ndarray.flatten(xpoints1)
    # y01 = np.ndarray.flatten(ypoints1)
    # x02 = np.ndarray.flatten(xpoints2)
    # y02 = np.ndarray.flatten(ypoints2)
    #
    # x0 = np.concatenate((x01,x02))
    # y0 = np.concatenate((y01,y02))

    # x0 = np.array([0.,5.,0.,5.])
    # y0 = np.array([0.,0.,5.,5.])
    # x0 = np.array([0.0, 0.0, 0.0, 15.0, 15.0, 15.0, 30.0, 30.0, 30.0, 3.0, 3.0, 3.0, 18.0, 18.0, 18.0, 33.0, 33.0, 33.0])
    # y0 = np.array([0.0, 15.0, 30.0, 0.0, 15.0, 30.0, 0.0, 15.0, 30.0, 0.0, 15.0, 30.0, 0.0, 15.0, 30.0, 0.0, 15.0, 30.0])

    x0 = np.array([0.,5.])
    y0 = np.array([0.,0.])
    dia = np.ones_like(x0)*turb_dia

    # print 'x0:',x0.tolist()
    # print 'y0:',y0.tolist()
    # dia = np.ones_like(x0)*turb_dia
    #
    # if rotdir_spec == 'cn':
    #     # rot = np.zeros_like(x0)
    #     # for i in range(np.size(rot)):
    #     #     if i % 2 == 0:
    #     #         rot[i] = turb_rot
    #     #     else:
    #     #         rot[i] = -turb_rot
    #     rot1 = np.ones_like(x01)*turb_rot
    #     rot2 = np.ones_like(x02)*-turb_rot
    #     rot = np.concatenate((rot1,rot2))
    # elif rotdir_spec == 'co':
    #     rot = np.ones_like(x0)*turb_rot
    # print 'rot:',rot.tolist(),'\n'
    # rot = np.array([turb_rot,turb_rot,-turb_rot,-turb_rot])
    rot = np.array([turb_rot,-turb_rot])

    twist = 0.0
    delta = 0.0
    B = 3
    chord = 0.128
    H = 6.1
    Hub = 2.

    rho = 1.225
    mu = 1.7894e-5

    ntheta = 36
    interp = 1 # linear airfoil interpolation
    # interp = 2 # cubic spline interpolation

    velx,vely,_,_,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[0]/2.,chord,twist,delta,B,rot[0],velf,rho,mu,interp,np.zeros(ntheta),np.zeros(ntheta))
    power_iso,cp_iso = _vawtwake.powercalc(np.array([0.]),np.array([0.]),np.array([turb_dia]),np.array([turb_rot]),velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,af_data,cl_data,cd_data,chord,twist,delta,B,H,rho,mu,velx,vely,interp)
    uvec_iso = velx
    vvec_iso = vely

    power_iso_tot = power_iso*np.size(x0)

    spaceval = 2.
    # xlow = 0.#x0[0]-spaceval
    # xupp = 40.#x0[-1]+spaceval
    # ylow = 0.#y0[0]-spaceval
    # yupp = 40.#y0[-1]+spaceval
    xlow = 0.
    xupp = 10.
    ylow = 0.
    yupp = 10.

    grid_radius = 10.#35.#sqrt((grid/2.)**2 + (grid/2.)**2) + 1.
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

    xdn = 2.

    xf = np.array([0.,xdn])
    yf = np.array([0.,0.])

    dec = 3
    power = np.zeros(dec)
    power1 = np.zeros(dec)
    power2 = np.zeros(dec)
    for pt in range(36):
        exec('indx'+str(pt+1)+' = np.zeros(dec)')
        exec('indy'+str(pt+1)+' = np.zeros(dec)')
    wakexp = np.zeros(dec)
    wakeyp = np.zeros(dec)
    adjust = np.linspace(1.1,1.104,dec)
    for q in range(dec):
        # xf[0] = xf[0] + adjust[q]
        # xf[1] = xf[1] + adjust[q]
        # yf[0] = yf[0] + adjust[q]
        yf[1] = yf[1] + adjust[q]
        m = 220
        n = 200
        interp = 2
        power_turb = np.zeros_like(xf)
        for j in range(np.size(xf)):
            rotw = np.zeros((1,np.size(xf)))
        k = 0
        for f in range(1):
            for j in range(np.size(xf)):
                rotw[f,j] = rot[k]
                k += 1

        winddir_turb = np.zeros_like(windroseDirections)
        for d in range(0, 1):
            # Adjusting coordinate system
            winddir_turb[d] = 270. - windroseDirections[d]
            if winddir_turb[d] < 0.:
                winddir_turb[d] += 360.
            winddir_turb_rad = pi*winddir_turb[d]/180.0
            xw = xf*cos(-winddir_turb_rad) - yf*sin(-winddir_turb_rad)
            yw = xf*sin(-winddir_turb_rad) + yf*cos(-winddir_turb_rad)
            # res = Parallel(n_jobs=-1)(delayed(vawt_power)(i,xw,yw,dia,rotw[d],ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,af_data,cl_data,cd_data,twist,delta,rho,mu,interp) for i in range(np.size(xf)) )
            #
            # for i in range(np.size(xf)):
            #     power_turb[i] = res[i][0]
            #     if i == 0:
            #         velx = res[i][1]
            #         vely = res[i][2]
            #         wakex = res[i][3]
            #         wakey = res[i][4]
            #     else:
            #         velx = np.append(velx,res[i][1])
            #         vely = np.append(vely,res[i][2])
            #         wakex = np.append(wakex,res[i][3])
            #         wakey = np.append(wakey,res[i][4])
            #
            # power[q] = np.sum(power_turb)
            wakex,wakey = vawt_wake(xw,yw,dia,rotw[d],ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)
            # xt = np.array([xw[0]])
            # yt = np.array([yw[0]])
            # diat = np.array([dia[0]])
            # rott = np.array([rotw[d,0]])
            # wakexp[q],wakeyp[q] = gskrintpoint(xt,yt,diat,rott,chord,B,xw[1],yw[1],dia[1],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)
            res = Parallel(n_jobs=1)(delayed(vawt_power2)(i,dia,rotw[d],ntheta,chord,B,velf,af_data,cl_data,cd_data,twist,delta,rho,mu,interp,wakex,wakey) for i in range(np.size(xf)) )
            for i in range(np.size(xf)):
                power_turb[i] = res[i][0]
                if i == 1:
                    for j in range(36):
                        exec('indx'+str(j+1)+'[q] = res[i][1][j]')
                        exec('indy'+str(j+1)+'[q] = res[i][2][j]')
                # if i == 0:
                #     velx = res[i][1]
                #     vely = res[i][2]
                # else:
                #     velx = np.append(velx,res[i][1])
                #     vely = np.append(vely,res[i][2])

        power[q] = np.sum(power_turb)
        power1[q] = power_turb[0]
        power2[q] = power_turb[1]

        # xf[0] = xf[0] - adjust[q]
        # xf[1] = xf[1] - adjust[q]
        # yf[0] = yf[0] - adjust[q]
        yf[1] = yf[1] - adjust[q]
        print q

    plt.figure()
    plt.plot(adjust,power)
    # plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_simp_'+str(int(xdn))+'_power.png')
    # plt.figure()
    # plt.plot(adjust,wakexp)
    # plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_simp2_'+str(int(xdn))+'_wakex.png')
    # plt.figure()
    # plt.plot(adjust,wakeyp)
    # plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_simp2_'+str(int(xdn))+'_wakey.png')

    # plt.figure()
    # plt.plot(adjust,power1)
    # plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_iso2_'+str(int(xdn))+'_power1.png')
    # plt.figure()
    # plt.plot(adjust,power2)
    # plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_iso2_'+str(int(xdn))+'_power2.png')

    # for i in range(36):
    #     plt.figure()
    #     exec('plt.plot(adjust,indx'+str(i+1)+')')
    #     plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_iso2_'+str(int(xdn))+'_'+str(i+1)+'_indx.png')
    #     plt.figure()
    #     exec('plt.plot(adjust,indy'+str(i+1)+')')
    #     plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_iso2_'+str(int(xdn))+'_'+str(i+1)+'_indy.png')

    # plt.figure()
    # plt.plot(adjust,wakexp)
    # plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_gskr3_'+str(plotpoint)+'_wakex.png')
    # plt.figure()
    # plt.plot(adjust,wakeyp)
    # plt.savefig('/fslhome/ebtingey/compute/optVAWT/wake_check_gskr3_'+str(plotpoint)+'_wakey.png')

    plt.show()

    # input = {'xvars':xf,'yvars':yf}
    # funcs,_ = obj_func(input)
    # pow = -1*funcs['obj']*1e3
    # SPLd = funcs['SPL']*10.
    # SPLw = np.zeros((np.size(windroseDirections),np.size(obs[:,0])))
    # k = 0
    # for i in range(np.size(windroseDirections)):
    #     for j in range(np.size(obs[:,0])):
    #         SPLw[i,j] = SPLd[k]
    #         k += 1

    # if rotdir_spec == 'cn':
    #     print 'Counter-rotating Turbines'
    # elif rotdir_spec == 'co':
    #     print 'Co-rotating Turbines'
    # print 'Wind Directions:',windroseDirections
    # print 'The power is:',pow,'W (',pow/(power_iso*np.size(x0)),')'
    # print 'The isolated power is:',power_iso*np.size(x0),'W'
    # print 'The x-locations:',xf
    # print 'The y-locations:',yf
    # # print 'The power of each turbine (W):',power_turb
    # # print 'The isolated power of one turbine is:',power_iso,'W'
    # print 'SPL:',SPLw
