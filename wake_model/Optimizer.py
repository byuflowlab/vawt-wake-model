from pyoptsparse import Optimization, SNOPT, pyOpt_solution
from os import path
import numpy as np
from numpy import sqrt,pi,sin,cos,fabs
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import time
import VAWT_Wake_Model as vwm
from ACsingle import actuatorcylinder
from sys import argv

#from joblib import Parallel, delayed
from mpi4py import MPI

import _vawtwake
import _bpmvawtacoustic
import pdb

try:
    from mpi4py import MPI
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    status=MPI.Status()
    size = comm.Get_size()
except:
    raise ImportError('mpi4py is required for parallelization')

def enum(*sequential, **named):
    """Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)

def obj_func(xdict):
    if rank ==0:
        global dia
        global rot
        global chord
        global twist
        global delta
        global B
        global H
        global rho
        global mu
        global Vinf
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
        global interp

        #BPM Globals (Parameters)
        global turb_dia
        global obs
        global B
        global Hub
        global H
        global chord
        global Vinf
        global ntheta

        # Simpson's rule integration division
        m = 220
        n = 200

        x = xdict['xvars'] # turbine x-positions
        y = xdict['yvars'] # turbine y-positions
        funcs = {}

        nturb = np.size(x) # number of turbines
        nwind = np.size(windroseDirections) # number of wind directions

        power_turb = np.zeros(nturb)
        power_dir = np.zeros(nwind)
        # reordering rotation directions to a matrix of wind directions
        rotw = np.zeros((nwind,nturb))
        k = 0
        for i in range(nwind):
            for j in range(nturb):
                rotw[i,j] = rot[k]
                k += 1

        winddir_turb = np.zeros_like(windroseDirections)
        for d in range(0, nwind):
            # adjusting coordinate system for wind direction
            winddir_turb[d] = 270. - windroseDirections[d]
            if winddir_turb[d] < 0.:
                winddir_turb[d] += 360.
            winddir_turb_rad = pi*winddir_turb[d]/180.0
            xw = x*cos(-winddir_turb_rad) - y*sin(-winddir_turb_rad)
            yw = x*sin(-winddir_turb_rad) + y*cos(-winddir_turb_rad)
            print 'Precalculate Wake Components'
            # calculating wake velocity components
            wakex,wakey = vawt_wake(xw,yw,dia,rotw[d],ntheta,chord,B,Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n)
            print 'Wake Component Precalculations Complete'
            # power parameters
            #           0   1       2       3  4 5 6    7       8       9       10      11  12  13      14      15
            constsPWR = [dia,rotw[d],ntheta,chord,H,B,Vinf,af_data,cl_data,cd_data,twist,delta,rho,interp,wakex,wakey]
            #BPM Precalculations
            rot=rotw[d]
            winddir=windroseDirections[d]
            #BPM parameters
            #         0 1 2                         3   4       5   6
            constsBPM=[x,y,windroseDirections[d],rotw[d],wakex,wakey,i]
            print 'BPM Precalculations Complete'
            #While Loop for calculations (State Machine)
            calcneeded=nturb+nobs
            casig=0
            calccompleted=0
            tlocs=0
            power_turb=np.empty(nturb)
            SPL_d=np.empty(nwind*len(obs))
            data=7
            pdb.set_trace()
            print 'Master Slave Loop'
            while calccompleted<calcneeded:
                tag = None
                print 'Complete: ',calccompleted,'/',calcneeded
                data = comm.recv(source=MPI.ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
                print 'recieve'
                source=status.Get_source()
                tag=status.Get_tag()
                print 'Source:',source,'Tag:',tag
                if tag == tags.READY:
                    if tlocs<nturb:
                        ConstsPWR=np.append(constsPWR,tlocs)
                        comm.ssend(ConstsPWR,dest=source,tag=tags.PWR)
                        print 'PWR sent'
                        tlocs+=1
                    elif tlocs>=nturb:
                        obsn=tlocs-nturb
                        ConstsBPM=np.append(constsBPM,obsn)
                        print 'BPM sent'
                        comm.ssend(ConstsBPM,dest=source,tag=tags.BPM)
                        tlocs+=1
                    else:
                        comm.send(None,dest=source,tag=tags.SLEEP)
                elif tag==tags.SPWR:
                    loc=data[1]
                    pwr=data[0]
                    power_turb[loc]=pwr
                    calccompleted+=1
                elif tag==tags.SBPM:
                    loc=data[1]
                    SPLt=data[0]
                    SPL_d[loc]=SPLt
                    calccompleted+=1

            for i in range(nturb):
                power_turb[i] = res[i]
            power_dir[d] = np.sum(power_turb)*windFrequencies[d]

            # calculating noise (dB)
            #SPL_d = bpm_noise(x,y,windroseDirections[d],rotw[d],wakex,wakey) -function pasted below

            #For i in range(nobs) - Parrellization
            #SPL=bpmnoise(ntheta,turbineX,turbineY,obs[i],winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,Vinf,wakex,wakey)


            SPL_dir = np.array(SPL_d)
            if d == 0:
                SPL = SPL_dir
            else:
                SPL = np.append(SPL,SPL_dir)

        power = np.sum(power_dir)

        funcs['obj'] = (-1.*power/1e3)

        funcs['SPL'] = (SPL)/10.

        print 'Power:',power,'W (Isolated: '+str(power_iso_tot)+' W; '+str(power/power_iso_tot)+')   Max SPL:',max(SPL),'dB'

        # calculating separation between turbines
        sep = sep_func(np.append(x,y))
        funcs['sep'] = sep

    fail = False
    return funcs, fail


def vawt_wake(xw,yw,dia,rotw,ntheta,chord,B,Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n):
    global wake_method

    t = np.size(xw) # number of turbines

    for i in range(t):
        xt = np.delete(xw,i)
        yt = np.delete(yw,i)
        diat = np.delete(dia,i)
        rott = np.delete(rotw,i)

        # xtd = np.delete(xw,i)
        # ytd = np.delete(yw,i)
        # diatd = np.delete(dia,i)
        # rottd = np.delete(rotw,i)
        # xt,yt,diat,rott = vwm.wake_order(xw[i],yw[i],dia[i],xtd,ytd,diatd,rottd)

        if wake_method == 'simp':
            wakexd,wakeyd = _vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],Vinf,coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9,m,n,1,1)
        elif wake_method == 'gskr':
            wakexd,wakeyd = vwm.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],dia[i],Vinf,False)

        if i == 0:
            wakex = wakexd
            wakey = wakeyd
        else:
            wakex = np.append(wakex,wakexd)
            wakey = np.append(wakey,wakeyd)

    return wakex,wakey


def vawt_power(i,dia,rotw,ntheta,chord,H,B,Vinf,af_data,cl_data,cd_data,twist,delta,rho,interp,wakext,wakeyt):
    global thetavec
    global Vnp
    global Vnn
    global Vtp
    global Vtn
    global Cpp
    global Cpn

    global useAC

    wakex = np.zeros(ntheta)
    wakey = np.zeros(ntheta)
    for j in range(ntheta):
        wakex[j] = wakext[j+ntheta*i]
        wakey[j] = wakeyt[j+ntheta*i]

    if useAC == True:
        Cp,_,_,_ = actuatorcylinder(ntheta,af_data,cl_data,cd_data,dia[i]/2.,chord,twist,delta,B,rotw[i],Vinf,rho,interp,wakex,wakey)

        power_turb = (0.5*rho*Vinf**3)*(dia[i]*H)*Cp

    elif useAC == False:
        power_turb,Cp = _vawtwake.powercalc(thetavec,Vinf,wakex,wakey,Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,rotw[i],dia[i]/2.,H,af_data,cl_data,cd_data,twist,rho,interp)

    return power_turb


# SPL CALCULATION BASED ON BPM ACOUSTIC MODEL
def bpm_noise(turbineX,turbineY,winddir,rot,wakex,wakey,i):
    global turb_dia
    global obs
    global B
    global Hub
    global H
    global chord
    global Vinf
    global ntheta

    nobs = np.size(obs[:,0])

    noise_corr = 1.
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
    #SPL = Parallel(n_jobs=-1)(delayed(bpmnoise)(ntheta,turbineX,turbineY,obs[i],winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,Vinf,wakex,wakey) for i in range(nobs) )
    SPL=bpmnoise(ntheta,turbineX,turbineY,obs[i],winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,Vinf,wakex,wakey)
    return SPL


def bpmnoise(ntheta,turbineX,turbineY,obs,winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,Vinf,wakex,wakey):
    return _bpmvawtacoustic.turbinepos(ntheta,turbineX,turbineY,obs,winddir,B,Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,rot,Vinf,wakex,wakey)


def sep_func(loc):
    global turb_dia

    space = 1.65 # turbine diameters apart

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

    # RUN OPTIMIZATION
    optimize = True
    #optimize = False

    #MPI States
    #            recieve  EXIT   calc calc   sendP  sendB  sleep
    tags = enum('READY','EXIT', 'PWR','BPM','SPWR','SBPM','SLEEP')
    # PLOT RESULTS
    plot = True
    # plot = False
    plot_type = 'start'
    plot_type = 'start-finish'
    # plot_type = 'finish'

    # SAVE RESULTS
    saveresult = True
    # saveresult = False

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
    global Vinf
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
    global Vnp
    global Vnn
    global Vtp
    global Vtn
    global Cpp
    global Cpn

    global ntheta
    global thetavec
    global useAC
    global wake_method

    # SPLlim = float(argv[1])
    # rotdir_spec = argv[2]
    # ntheta = int(argv[3])
    # wake_method = argv[4]
    # nRows = int(argv[5])
    # nCols = int(argv[6])

    SPLlim = 100.           # sound pressure level limit of observers
    rotdir_spec = 'cn'      # rotation direction (cn- counter-rotating, co- co-rotating)
    ntheta = 72             # number of points around blade flight path
    wake_method = 'simp'    # wake model calculation using Simpson's rule
    wake_method = 'gskr'    # wake model calculation using 21-point Gauss-Kronrod
    nRows = 1               # number of paired group rows
    nCols = 2               # number of paired group columns
    if rank==0:
        print '\nSPL Limit:',SPLlim
    if rotdir_spec == 'cn':
        if rank==0:
            print 'Rotation Specification: Counter-rotating'
    elif rotdir_spec == 'co':
        if rank==0:
            print 'Rotation Specification: Co-rotating'
    if rank==0:
        print 'Points around VAWT:',ntheta
    if wake_method == 'simp':
        if rank==0:
            print "Using Simpson's Rule for Wake Calculation"
    elif wake_method == 'gskr':
        if rank==0:
            print "Using 21-Point Gauss-Kronrod Quadrature for Wake Calculation"
    if rank==0:
        print 'Rows of Paired Groups:',nRows
        print 'Columns of Paired Groups:',nCols,'\n'

    basepath = path.join(path.dirname(path.realpath('__file__')), 'data')
    foildata = basepath + path.sep + 'airfoils/du06w200.dat'
    foildata='/Users/ning3/OneDrive - BYU Office 365/Research/WES/Optimization/Combinedgit/vawt-wake-model/wake_model/data/airfoils/du06w200.dat'
    af_data,cl_data,cd_data = vwm.airfoil_data(foildata)
    coef0,coef1,coef2,coef3,coef4,coef5,coef6,coef7,coef8,coef9 = vwm.coef_val()

    # define wind specifications
    windroseDirections = np.array([205.,225.,245.])
    windFrequencies = np.array([0.25,0.50,0.25])
    nwind = np.size(windroseDirections)
    if rank==0:
        print 'wind:',windroseDirections

    # define turbine specifications
    Vinf = 8.                           # free stream velocity (m/s)
    turb_dia = 1.2                      # turbine diameter (m)
    tsrd = 2.625                        # tip-speed ratio
    turb_rot = tsrd*Vinf/(turb_dia/2.)  # turbine rotation rate (rad/sec)
    twist = 0.0                         # blade twist angle (rad)
    delta = 0.0                         # blade curvature angle (rad)
    B = 3                               # number of turbine blades
    chord = 0.128                       # chord length (m)
    H = 6.1                             # turbine blade height (m)
    Hub = 2.                            # turbine hub height (m)

    rho = 1.225                         # air density (kg/m^3)
    mu = 1.7894e-5                      # fluid viscosity (kg/ms)

    # interp = 1 # linear airfoil interpolation
    interp = 2 # cubic spline interpolation

    thetavec = np.zeros(ntheta)
    for i in range(ntheta):
        thetavec[i] = (2.*pi/ntheta)*(i+1)-(2.*pi/ntheta)/2.

    # setup initial turbine positions
    grid_start = 2.             # location of starting corner (m)
    pair_sep = 3.               # separation distance between pairs (m)
    group_sep = pair_sep + 10.  # separation distance between paired groups (m)

    x01 = np.zeros(nRows*nCols)
    y01 = np.zeros(nRows*nCols)
    x02 = np.zeros(nRows*nCols)
    y02 = np.zeros(nRows*nCols)

    x1 = np.linspace(grid_start,grid_start+group_sep*(nCols-1),nCols)
    y1 = np.linspace(grid_start,grid_start+group_sep*(nRows-1),nRows)
    x2 = np.linspace(grid_start+pair_sep,grid_start+pair_sep+group_sep*(nCols-1),nCols)
    # y2 = np.linspace(grid_start+pair_sep,grid_start+pair_sep+group_sep*(nRows-1),nRows) # use to angle the pairs

    k = 0
    for i in range(nRows):
        for j in range(nCols):
            x01[k] = x1[j]
            y01[k] = y1[i]
            x02[k] = x2[j]
            y02[k] = y1[i]
            # y02[k] = y2[i]
            k += 1

    x0 = np.concatenate((x01,x02))
    y0 = np.concatenate((y01,y02))
    nturb = np.size(x0)
    if rank==0:
        print 'x0:',x0.tolist()
        print 'y0:',y0.tolist()

    dia = np.ones(nturb)*turb_dia

    # specify turbine rotation direction
    if rotdir_spec == 'cn':
        rot1 = np.ones_like(x01)*turb_rot
        rot2 = np.ones_like(x02)*-turb_rot
        rot = np.concatenate((rot1,rot2))
    elif rotdir_spec == 'co':
        rot = np.ones(nturb)*turb_rot
    if rank==0:
        print 'rot:',rot.tolist(),'\n'
    for i in range(nwind-1):
        rot = np.append(rot,rot)

    # specifying boundary locations
    spaceval = 2.
    xlow = 0.
    xupp = max(x2[-1],y1[-1]) + 2.
    # xupp = max(x2[-1],y2[-1]) + 2.
    ylow = 0.
    yupp = max(x2[-1],y1[-1]) + 2.
    # yupp = max(x2[-1],y2[-1]) + 2.

    # specifying observer locations
    grid_x = xupp/2.
    grid_y = yupp/2.
    grid_radius = int(sqrt((grid_x)**2 + (grid_y)**2)) + 4.

    nobs = 8
    obs_theta = np.linspace(-pi,pi,nobs+1)
    obs = np.zeros((nobs,3))
    for i in range(nobs):
        obs[i,0] = grid_x + grid_radius*cos(obs_theta[i])
        obs[i,1] = grid_y + grid_radius*sin(obs_theta[i])
        obs[i,2] = 2.
    if rank==0:
        print nobs,'observers around a radius of ',grid_radius,'\n'

    # power value precompute (for CCW and CW directions)
    Cp_iso,Tpp,Vnp,Vtp = actuatorcylinder(ntheta,af_data,cl_data,cd_data,turb_dia/2.,chord,twist,delta,B,fabs(turb_rot),Vinf,rho,interp,np.zeros(ntheta),np.zeros(ntheta)) # CCW
    Cpp = (fabs(turb_rot)*B/(2.*pi*rho*Vinf**3))*Tpp
    _,Tpn,Vnn,Vtn = actuatorcylinder(ntheta,af_data,cl_data,cd_data,turb_dia/2.,chord,twist,delta,B,-fabs(turb_rot),Vinf,rho,interp,np.zeros(ntheta),np.zeros(ntheta)) # CW
    Cpn = (fabs(turb_rot)*B/(2.*pi*rho*Vinf**3))*Tpn

    power_iso = (0.5*rho*Vinf**3)*(dia[0]*H)*Cp_iso # isolated power of a single turbine (W)
    power_iso_tot = power_iso*nturb # total power of isolated turbines (W)

    # option to use actuator cylinder or not (use a correction factor method)
    useAC = True
    useAC = False
    toOpt = None
    if optimize == True and rank==0:
        num_workers = size - 1
        closed_workers = 0
        print("Master starting with %d workers" % num_workers)
        # optimization setup
        optProb = Optimization('VAWT_Power', obj_func)
        optProb.addObj('obj')

        n = np.size(x0)
        optProb.addVarGroup('xvars', n, 'c', lower=xlow, upper=xupp, value=x0)
        optProb.addVarGroup('yvars', n, 'c', lower=ylow, upper=yupp, value=y0)

        num_cons_sep = (n-1)*n/2
        optProb.addConGroup('sep', num_cons_sep, lower=0, upper=None)
        num_cons_obs = nobs*nwind
        optProb.addConGroup('SPL', num_cons_obs, lower=0, upper=SPLlim/10.)

        opt = SNOPT()
        opt.setOption('Scale option',0)
        if rotdir_spec == 'cn':
            opt.setOption('Print file',basepath + path.sep + 'optimization_results/SNOPT_print_SPL'+str(SPLlim)+'_turb'+str(n)+'_counterrot.out')
            opt.setOption('Summary file',basepath + path.sep + 'optimization_results/SNOPT_summary_SPL'+str(SPLlim)+'_turb'+str(n)+'_counterrot.out')
        elif rotdir_spec == 'co':
            opt.setOption('Print file',basepath + path.sep + 'optimization_results/SNOPT_print_SPL'+str(SPLlim)+'_turb'+str(n)+'_corot.out')
            opt.setOption('Summary file',basepath + path.sep + 'optimization_results/SNOPT_summary_SPL'+str(SPLlim)+'_turb'+str(n)+'_corot.out')
            toOpt = True
        # run optimization
    comm.Barrier()

    if optimize==True and rank!=0:
        name = MPI.Get_processor_name()
        print 'I am a worker with rank %d on %s.' % (rank,name)
        #pyOpt dummy mpi calls (pyOPT mpi calls that we do not use but must include)
        lists = ['xvars','yvars']
        other = ['SPL','sep']
        comm.send(None,dest=0)
        comm.gather(lists,root=0)
        dummy = comm.bcast(lists,root=0)
        variables = comm.bcast(lists,root=0)
        comm.gather(other,root=0)
        variables = comm.bcast(other,root=0)
        comm.gather(other,root=0)
        while True:
            comm.send([None,MPI.INT],dest=0,tag=tags.READY)
            print 'Worker',rank,'is Ready'
            task=comm.recv(source=0,tag=MPI.ANY_TAG,status=status)
            tag= status.Get_tag()
            source=status.Get_source()
            if tag==tags.PWR:
                print 'PWR assigned to %i'%rank
                dia = task[0]
                rotwl = task[1]
                ntheta=task[2]
                chord=task[3]
                H=task[4]
                B=task[5]
                Vinf = task[6]
                af_data=task[7]
                cl_data=task[8]
                cd_data=task[9]
                twist=task[10]
                delta=task[11]
                rho=task[12]
                interp=task[13]
                wakex=task[14]
                wakey=task[15]
                i=task[16]
                #Calculate Power
                res=vawt_power(i,dia,rotwl,ntheta,chord,H,B,Vinf,af_data,cl_data,cd_data,twist,delta,rho,interp,wakex,wakey)
                result=np.array([res,i])
                print 'calculate'
                comm.ssend(result,dest=0,tag=tags.SPWR)
                print 'sent'
            elif tag==tags.BPM:
                print 'BPM assigned to %i'%rank
                turbineX=task[0]
                turbineY=task[1]
                winddir=task[2]
                rot=task[3]
                wakex=task[4]
                wakey=task[5]
                i=task[6]
                #Calculate BPM
                SPL=bpm_noise(turbineX,turbineY,winddir,rot,wakex,wakey,i)
                result=np.append(SPL,i)
                comm.ssend(result,dest=0,tag=tags.SPWR)
            elif tag==tags.SLEEP:
                time.sleep(2)
            elif tag==tags.EXIT:
                break
        comm.send(None,dest=0,tag=tags.EXIT)


    if optimize==True and rank==0:
        res = opt(optProb)
        print 'Optimization Complete, Closing Workers'
        #'Close' workers and print Result
        if rank==0:
            while closed_workers<nworkers:
                data=comm.recv(souce=MPI>ANY_SOURCE,tag=MPI.ANY_TAG,status=status)
                source=status.Get_source()
                tag=status.Get_tag()
                if tag==tags.READY:
                    comm.send(None, dest=source, tag=tags.EXIT)
                elif tag==tags.EXIT:
                    print("Worker %d exited." % source)
                    closed_workers += 1
            #Print the result after 'closing' processes
            print res
        toOpt=False


    if optimize==True and rank==0:
        pow = np.array(-1*res.fStar)*1e3
        xf = res.xStar['xvars']
        yf = res.xStar['yvars']

        SPLd = funcs['SPL']*10.
        SPLw = np.zeros((nwind,nobs))
        k = 0
        for i in range(nwind):
            for j in range(nobs):
                SPLw[i,j] = SPLd[k]
                k += 1
    elif optimize==False:
        xf = x0
        yf = y0

        # xf = np.array([])
        # yf = np.array([])

        input = {'xvars':xf,'yvars':yf}
        funcs,_ = obj_func(input)
        pow = -1*funcs['obj']*1e3
        SPLd = funcs['SPL']*10.
        SPLw = np.zeros((nwind,nobs))
        k = 0
        for i in range(nwind):
            for j in range(nobs):
                SPLw[i,j] = SPLd[k]
                k += 1
    if rank==0:
        print 'Wind Directions:',windroseDirections
        print 'The power is:',pow,'W (',pow/(power_iso_tot),')'
        print 'The isolated power is:',power_iso_tot,'W'
        print 'The x-locations:',xf
        print 'The y-locations:',yf
        print 'SPL:',SPLw

        if saveresult == True:
            # writing results to text file
            if rotdir_spec == 'cn':
                filename = basepath + path.sep + 'optimization_results/Optrun_SPL'+str(SPLlim)+'_r'+str(nRows)+'_c'+str(nCols)+'_counterrot.txt'
        elif rotdir_spec == 'co':
            filename = basepath + path.sep + 'optimization_results/Optrun_SPL'+str(SPLlim)+'_r'+str(nRows)+'_c'+str(nCols)+'_corot.txt'
        target = open(filename,'w')

        target.write('\nThe power is: '+str(pow)+' W ('+str(pow/(power_iso_tot))+')\n')
        target.write('The isolated power is: '+str(power_iso_tot)+' W\n')
        target.write('Max SPL: '+str(np.max(SPLw))+' dB\n')
        target.write('\nWind Directions: '+str(windroseDirections.tolist())+' degrees\n')
        target.write('X-locations (initial): '+str(x0.tolist())+' m\n')
        target.write('X-locations (final): '+str(xf.tolist())+' m\n')
        target.write('Y-locations (initial): '+str(y0.tolist())+' m\n')
        target.write('Y-locations (final): '+str(yf.tolist())+' m\n')
        target.write('\nSPL: '+str(SPLw.tolist())+' dB\n')
        target.close()
        if plot == True:
            fs = 20
            ms = 10
            plt.figure(1,figsize=(11.5,8))
            plt.subplots_adjust(right=0.68)
            if plot_type == 'start':
                plt.plot(1e1000,1e1000,'o',color='k',markersize=ms,fillstyle='full',label='Turbines (CCW)')
                plt.plot(1e1000,1e1000,'o',color='silver',markersize=ms,fillstyle='full',label='Turbines (CW)')
            elif plot_type == 'start-finish':
                plt.plot(1e1000,1e1000,'o',color='k',markersize=ms,fillstyle='none',label='Original Turbines')
                plt.plot(1e1000,1e1000,'o',color='k',markersize=ms,fillstyle='full',label='Optimized (CCW)')
                plt.plot(1e1000,1e1000,'o',color='silver',markersize=ms,fillstyle='full',label='Optimized (CW)')
            elif plot_type == 'finish':
                plt.plot(1e1000,1e1000,'o',color='k',markersize=ms,fillstyle='full',label='Turbines (CCW)')
                plt.plot(1e1000,1e1000,'o',color='silver',markersize=ms,fillstyle='full',label='Turbines (CW)')
            for i in range(np.size(x0)):
                if rot[i] > 0.:
                    if plot_type == 'start':
                        circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,facecolor='k',edgecolor='k')
                        plt.gca().add_patch(circ)
                    elif plot_type == 'start-finish':
                        circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,facecolor='w',edgecolor='k')
                        plt.gca().add_patch(circ)
                elif rot[i] < 0.:
                    if plot_type == 'start':
                        circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,facecolor='silver',edgecolor='k')
                        plt.gca().add_patch(circ)
                    elif plot_type == 'start-finish':
                        circ = plt.Circle((x0[i],y0[i]),dia[i]/2.,facecolor='w',edgecolor='gray')
                        plt.gca().add_patch(circ)
            if plot_type == 'start-finish':
                for i in range(np.size(x0)):
                    plt.plot([x0[i], xf[i]], [y0[i], yf[i]], '--k',zorder=-1)
            for i in range(np.size(xf)):
                if rot[i] > 0.:
                    if plot_type == 'start-finish' or plot_type == 'finish':
                        circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,facecolor='k',edgecolor='k')
                        plt.gca().add_patch(circ)
                elif rot[i] < 0.:
                    if plot_type == 'start-finish' or plot_type == 'finish':
                        circ = plt.Circle((xf[i],yf[i]),dia[i]/2.,facecolor='silver',edgecolor='k')
                        plt.gca().add_patch(circ)
            plt.plot(obs[:,0],obs[:,1],'^',color='lime',markersize=ms,label='Observers')
            rect = plt.Rectangle((xlow,ylow), xupp-xlow,yupp-ylow, linestyle='dashed',linewidth=2,facecolor="#ffffff",fill=False,label='Boundaries')
            plt.gca().add_patch(rect)
            plt.annotate('N', xy=((min(obs[:,0])-5)*(2./3.), (max(obs[:,0])+5)*43./45.), xycoords='data', xytext=(0,-40), textcoords='offset points', size=fs, va="center", ha="center", arrowprops=dict(arrowstyle='simple', facecolor='k'), horizontalalignment='right', verticalalignment='top',color='k')
            plt.annotate('Wind', xy=(((max(obs[:,0])+5)-(min(obs[:,0])-5))/4.5,(min(obs[:,0])-5)*1./15.),  xycoords='data', xytext=(-50*np.tan(45.*np.pi/180),-50), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=15), horizontalalignment='right', verticalalignment='top', fontsize=fs,color='k') # wind at 225 degrees
            # plt.annotate('Wind', xy=(((max(obs[:,0])+5)-(min(obs[:,0])-5))/4.,(min(obs[:,0])-5)*4./15.),  xycoords='data', xytext=(-50,0), textcoords='offset points', arrowprops=dict(facecolor='skyblue',width=5,headwidth=15), horizontalalignment='right', verticalalignment='center', fontsize=fs,color='k') # wind at 90 degrees
            plt.legend(loc="upper left", bbox_to_anchor=(1,1),fontsize=fs)
            plt.xticks(fontsize=fs)
            plt.yticks(fontsize=fs)
            plt.xlabel('X-Position (m)',fontsize=fs)
            plt.ylabel('Y-Position (m)',fontsize=fs)
            plt.xlim(min(obs[:,0])-5,max(obs[:,0])+5)
            plt.ylim(min(obs[:,0])-5,max(obs[:,0])+5)
            if saveresult == True:
                if rotdir_spec == 'cn':
                    plt.savefig(basepath + path.sep + 'optimization_results/SPLrunlayout_'+str(SPLlim)+'_r'+str(nRows)+'_c'+str(nCols)+'_counterrot.png')
                elif rotdir_spec == 'co':
                    plt.savefig(basepath + path.sep + 'optimization_results/SPLrunlayout_'+str(SPLlim)+'_r'+str(nRows)+'_c'+str(nCols)+'_corot.png')

plt.show()
