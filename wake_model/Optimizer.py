from pyoptsparse import Optimization, SNOPT, pyOpt_solution
from os import path
import numpy as np
from numpy import sqrt,pi,sin,cos,fabs
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import VAWT_Wake_Model as vwm
from ACsingle import actuatorcylinder

from joblib import Parallel, delayed

import _vawtwake
import _bpmvawt
from spl75data import spldata


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

        # wakex,wakey = _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)
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

    # print i
    wakex = np.zeros(ntheta)
    wakey = np.zeros(ntheta)
    for j in range(ntheta):
        wakex[j] = wakext[j+ntheta*i]
        wakey[j] = wakeyt[j+ntheta*i]

    uvec,vvec,_,Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],velf,rho,mu,interp,wakex,wakey)

    power_turb = (0.5*rho*velf**3)*(dia[i]*H)*Cp

    return power_turb,uvec,vvec



def overlap(ntheta,xt,yt,diat,rott,chord,B,xw,yw,dia,rotw,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n):
    return _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw,yw,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1)

def gskrint(p,xt,yt,diat,rott,chord,B,x0,y0,dia,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n):

    split = True
    # split = False

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
    

    optimize = True
    optimize = False
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

    foildata = '/Users/ning1/Documents/FLOW Lab/VAWTAC/airfoils/du06w200.dat'

    af_data,cl_data,cd_data = vwm.airfoil_data(foildata)
    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()

    windroseDirections = np.array([270.])
    windFrequencies = np.array([1.])

    # define turbine size
    velf = 8.
    turb_dia = 1.2 # m
    tsrd = 2.625
    turb_rot = tsrd*velf/(turb_dia/2.) # rad/sec


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

    x0 = np.array([2.,8.])
    y0 = np.array([5.,5.])
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


    if optimize == True:
        optProb = Optimization('VAWT_Power', obj_func)
        optProb.addObj('obj')

        n = np.size(x0)
        optProb.addVarGroup('xvars', n, 'c', lower=xlow, upper=xupp, value=x0)
        optProb.addVarGroup('yvars', n, 'c', lower=ylow, upper=yupp, value=y0)

        num_cons_sep = (n-1)*n/2
        optProb.addConGroup('sep', num_cons_sep, lower=0, upper=None)
        # num_cons_obs = np.size(obs[:,0])
        # optProb.addConGroup('SPL', num_cons_obs, lower=0, upper=100./10.)

        opt = SNOPT()
        opt.setOption('Scale option',0)
        res = opt(optProb, sens=None)
        print res
    
        pow = np.array(-1*res.fStar)*1e3
        xf = res.xStar['xvars']
        yf = res.xStar['yvars']

        SPLd = funcs['SPL']*10.
        SPLw = np.zeros((np.size(windroseDirections),np.size(obs[:,0])))
        k = 0
        for i in range(np.size(windroseDirections)):
            for j in range(np.size(obs[:,0])):
                SPLw[i,j] = SPLd[k]
                k += 1
    else:
        # xf = np.array([2.04299218,  10.69711973,   0.62759177,  -0.03202973])
        # yf = np.array([ -0.30868907,  5.82290592,  2.25557342,  0.38809135])

        # xf = np.array([ 2.23703622,  7.03994721,  4.07916559,  8.24416228])
        # yf = np.array([ 2.6406578,   4.31184911,  8.51934267,  5.89849381])

        # The power is: 2709.75262982 W ( 1.04235323524 )
        # The x-locations: [ 0.47307964  1.29739913  0.82146607  2.06454097]
        # The y-locations: [-0.40000015  2.2999387   6.73322378  4.12528618]
        # The power of each turbine (W): [ 666.68124805  671.4849696   671.89120384  699.69520833]

        # The power is: 2807.44392777 W ( 1.07993188332 )
        # xf = np.array([4.757338,5.114716,4.529348,4.935540])
        # yf = np.array([-1.159434,12.000000,0.807396,10.028124])

        # xf = x0
        # yf = y0

        # if splcon == False:
        #     # no SPL
        #     xf = np.array([2.801866358429759, 2.848513388835996, 2.7840232552345183, 19.413051540423652, 19.294346000080274, 19.37538244358126, 35.71158929637296, 35.74003509165456, 35.916194210704674, 4.202034249425052, 4.248606365289341, 4.183953474229456, 20.81306522925305, 20.694564690314255, 20.775301580690723, 37.11181018164184, 37.14049970831674, 37.31625026257014])
        #     yf = np.array([2.759603512205314, 19.243522850507482, 35.55694553748992, 2.777604587365613, 19.19948769971568, 35.56523339362047, 2.769059798842588, 19.209334818367946, 35.61120185119654, 1.3596285353919102, 17.843472942144167, 34.15673288017082, 1.3774753720919206, 17.79956354668723, 34.165009629581604, 1.3691378355607284, 17.809656398916616, 34.21111461426202])
        #
        # elif splcon == True:
        #     # SPL 75
        #     xf = np.array([4.304407710922205, 2.742599728579465, 2.679006956311623, 19.2916124058625, 19.22321050813174, 19.182612281070345, 35.23224306969283, 35.67686946379533, 35.68759490207986, 5.87253134630314, 4.310328582887386, 4.2468315787080915, 20.859304340936088, 20.79089403217838, 20.75056006419225, 36.80032195857313, 37.24467754502207, 37.25548719139211])
        #     yf = np.array([4.08173378052814, 19.28439342060845, 35.62355805713165, 2.662644769852271, 19.163621171867387, 35.609260235839464, 3.0149675533919438, 19.172711904604267, 35.599708897111526, 2.872854494508013, 18.075002561662405, 34.41429135636939, 1.4532060505406381, 17.954171556929147, 34.4001532207326, 1.8060305318752128, 17.9634237584284, 34.39052992513346])


        # Co-rot
        # Wind Directions: [ 270.]
        # The power is: 2669.46810939 W ( 1.02685709742 )
        # The isolated power is: 2599.64908076 W
        # The x-locations: [-1.64436126  3.37422599  8.6434073   9.49867472]
        # The y-locations: [  3.20403137   4.9115479    7.98224744  10.01672136]
        # SPL: [[ 82.675871    80.78388182  78.45728358  78.42972317  80.50865692
        #    84.85919787  80.73753673  79.76310028]]

        xf = np.array([1.365011,8.741231])
        yf = np.array([4.126831,5.816073])

        dec = 500
        power = np.zeros(dec)
        adjust = np.linspace(-1.,1.,dec)
        for q in range(dec):
            # xf[0] = xf[0] + adjust[q]
            # xf[1] = xf[1] + adjust[q]
            # yf[0] = yf[0] + adjust[q]
            yf[1] = yf[1] + adjust[q]
            m = 220
            n = 200
            interp = 1
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
            res = Parallel(n_jobs=-1)(delayed(vawt_power)(i,xw,yw,dia,rotw[d],ntheta,chord,B,velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,af_data,cl_data,cd_data,twist,delta,rho,mu,interp) for i in range(np.size(xf)) )

            for i in range(np.size(xf)):
                power_turb[i] = res[i][0]
                if i == 0:
                    velx = res[i][1]
                    vely = res[i][2]
                    wakex = res[i][3]
                    wakey = res[i][4]
                else:
                    velx = np.append(velx,res[i][1])
                    vely = np.append(vely,res[i][2])
                    wakex = np.append(wakex,res[i][3])
                    wakey = np.append(wakey,res[i][4])

            power[q] = np.sum(power_turb)

            # xf[0] = xf[0] - adjust[q]
            # xf[1] = xf[1] - adjust[q]
            # yf[0] = yf[0] - adjust[q]
            yf[1] = yf[1] - adjust[q]
            print q

        plt.plot(adjust,power)
        # plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/optsweep1.png')
        # plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/optsweep2.png')
        # plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/optsweep3.png')
        plt.savefig('/Users/ning1/Documents/FLOW Lab/Wind Farm Layout Power/optsweep4.png')
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

    if rotdir_spec == 'cn':
        print 'Counter-rotating Turbines'
    elif rotdir_spec == 'co':
        print 'Co-rotating Turbines'
    print 'Wind Directions:',windroseDirections
    print 'The power is:',pow,'W (',pow/(power_iso*np.size(x0)),')'
    print 'The isolated power is:',power_iso*np.size(x0),'W'
    print 'The x-locations:',xf
    print 'The y-locations:',yf
    # print 'The power of each turbine (W):',power_turb
    # print 'The isolated power of one turbine is:',power_iso,'W'
    print 'SPL:',SPLw
    
    if plot == True:
        plottype = 'initial'
        plottype = 'final'

        fs = 21
        ms = 10
        if splcon == False:
            plt.figure(1,figsize=(11.5,8))
            plt.subplots_adjust(right=0.68)
        elif splcon == True:
            plt.figure(1,figsize=(11.5,9.7))
            plt.subplots_adjust(right=0.68,bottom=0.0)

        if plottype == 'initial':
            for i in range(np.size(x0)):
                if rot[i] > 0.:
                    if i == 0:
                        plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='full',label='Turbine Locations')
                    else:
                        plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='full')
                    # circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='k',fill=True)
                    # plt.gca().add_patch(circ)
                elif rot[i] < 0.:
                    plt.plot(x0[i],y0[i],'o',color='silver',markersize=ms,fillstyle='full')
                    # circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,color='silver',fill=True)
                    # plt.gca().add_patch(circ)
        elif plottype == 'final':
            # for i in range(np.size(x0)):
            #     if rot[i] > 0.:
            #         if i == 0:
            #             plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='none',label='Original Turbines')
            #         else:
            #             plt.plot(x0[i],y0[i],'o',color='k',markersize=ms,fillstyle='none')
            #         # circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='k',fill=False)
            #         # plt.gca().add_patch(circ)
            #     elif rot[i] < 0.:
            #         plt.plot(x0[i],y0[i],'o',color='gray',markersize=ms,fillstyle='none')
            #         # circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='silver',fill=False)
            #         # plt.gca().add_patch(circ)
            # for i in range(np.size(x0)):
            #     plt.plot([x0[i], xf[i]], [y0[i], yf[i]], '--k',zorder=-1)
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

        if splcon == True:
            # SPL 75
            velx,vely,wakex,wakey = spldata()

            n = 100
            xp = np.linspace(-20.0, 60.0, n)
            yp = np.linspace(-20.0, 60.0, n)
            [X, Y] = np.meshgrid(xp, yp)
            F = np.zeros((n, n))

            ntheta = 36
            # nobs = np.size(obs[:,0])

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

            point = 0
            for i in range(n):
                for j in range(n):
                    F[i,j] = _bpmvawt.turbinepos(ntheta,xf,yf,np.array([X[i,j],Y[i,j],2.]),windroseDirections[0],B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,velf,velx,vely,wakex,wakey)
                    # print F[i,j]
                    point += 1
                    print point, 'of', n**2

            lb = 68.0 #lower bound on velocity to display
            ub = 88.0 #upper bound on velocity to display
            ran = 100 #number of contours between the velocity bounds
            cbtix = np.linspace(lb,ub,5)
            bounds = np.linspace(lb,ub,ran)
            CS = plt.contourf(xp, yp, F, ran, cmap=plt.cm.parula, vmax=ub, vmin=lb, levels=bounds)
            # CSn = plt.contour(xp, yp, F, cmap=plt.cm.parula)
            CB = plt.colorbar(CS,ticks=cbtix,orientation='horizontal',pad=0.13)
            CB.ax.tick_params(labelsize=fs)
            CB.ax.set_xlabel('Sound Pressure Level (dB)',fontsize=fs)

        rect = plt.Rectangle((xlow,ylow), xupp-xlow,yupp-ylow, linestyle='dashed',linewidth=2,facecolor="#ffffff",fill=False,label='Boundaries')
        plt.gca().add_patch(rect)

        if splcon == False:
            plt.annotate('N', xy=(-12, 56), xycoords='data', xytext=(-12, 49), textcoords='data', size=fs, va="center", ha="center", arrowprops=dict(arrowstyle='simple', facecolor='k'), horizontalalignment='right', verticalalignment='top',color='k')
            plt.annotate('Wind', xy=(15,-6),  xycoords='data', xytext=(-50*np.tan(45.*np.pi/180),-50), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=15), horizontalalignment='right', verticalalignment='top', fontsize=fs,color='k')
        elif splcon == True:
            plt.annotate('N', xy=(-12, 56), xycoords='data', xytext=(-12, 49), textcoords='data', size=fs, va="center", ha="center", arrowprops=dict(arrowstyle='simple', facecolor='w'), horizontalalignment='right', verticalalignment='top',color='w')
            plt.annotate('Wind', xy=(15,-6),  xycoords='data', xytext=(-50*np.tan(45.*np.pi/180),-50), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=15), horizontalalignment='right', verticalalignment='top', fontsize=fs,color='w')

        plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
        plt.xlabel('X-Position (m)',fontsize=fs)
        plt.ylabel('Y-Position (m)',fontsize=fs)


        # plt.xlim(x0[0]-spaceval-1,x0[-1]+spaceval+1)
        # plt.ylim(y0[0]-spaceval-1,y0[-1]+spaceval+1)

        if splcon == True and plottype == 'final':
            plt.savefig('/Users/ning1/Documents/FLOW Lab/Thesis/BYU_ME_Thesis_Template/figures/chp_power/optimized_SPL2.pdf')


        if contour == True:

            N = 10
            xplot = np.linspace(x0[0]-spaceval-1,x0[-1]+spaceval+1,N)
            yplot = np.linspace(y0[0]-spaceval-1,y0[-1]+spaceval,N)
            [X,Y] = np.meshgrid(xplot,yplot)
            Vel = np.zeros((N,N))

            k = 0
            for i in range(N):
                for j in range(N):
                    Velx,Vely = _vawtwake.overlappoint(x0,y0,dia,rot,chord,B,X[i,j],Y[i,j],velf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,220,200,1)
                    Vel[i,j] = sqrt((Velx+velf)**2 + (Vely)**2)
                    k += 1
                    print k

            plt.figure(2)
            lb = 0.0 # lower bound on velocity to display
            ub = velf # upper bound on velocity to display
            ran = 200 # number of contours between the velocity bounds
            bounds = np.linspace(lb,ub,ran)
            v = np.linspace(lb,ub,4) # setting the number of tick marks on colorbar
            CS = plt.contourf(X,Y,Vel,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.jet) # plotting the contour plot
            CB = plt.colorbar(CS, ticks=v) # creating colorbar
            for i in range(np.size(x0)):
                circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,color='g',edgecolor='k',fill=True)
                plt.gca().add_patch(circ)

            plt.xlim(x0[0]-spaceval-1,x0[-1]+spaceval+1)
            plt.ylim(y0[0]-spaceval-1,y0[-1]+spaceval+1)

        plt.show()
