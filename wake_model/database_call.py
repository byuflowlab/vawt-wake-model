import csv
from os import path
import numpy as np
from numpy import fabs
from scipy.interpolate import RectBivariateSpline,griddata


def vorticity(tsr, solidity):
    """
    Using EMG distribution parameters to define the vorticity strength and shape
    
    Parameters
    ----------
    tsr : float
        tip-speed ratio
    solidity : float
        turbine solidity
    
    Returns
    ----------
    loc : array
        array of the location parameter (3 values)
    spr : array
        array of the spread parameter (2 values)
    skw : array
        array of the skew parameter (2 values)
    scl : array
        array of the scale parameter (3 values)
    """
    
    # Reading in csv file (vorticity database)
    basepath = path.join(path.dirname(path.realpath(__file__)), 'data')
    fdata = basepath + path.sep + 'vortdatabase.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw,0)
            vortdat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d, float(row[1]))  # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            vortdat = np.vstack([vortdat, raw])  # adding entry to vorticity database array
        i += 1
    f.close()
    
    vortdat = np.delete(vortdat, 0, axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float)  # converting tip-speed ratio entries into floats
    vortdat = vortdat.astype(np.float)  # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_loc1 = vortdat[i*10]\ns'+sol+'_loc2 = vortdat[i*10+1]\ns'+sol+'_loc3 = vortdat[i*10+2]\ns'+sol+'_spr1 = vortdat[i*10+3]\ns'+sol+'_spr2 = vortdat[i*10+4]\ns'+sol+'_skw1 = vortdat[i*10+5]\ns'+sol+'_skw2 = vortdat[i*10+6]\ns'+sol+'_scl1 = vortdat[i*10+7]\ns'+sol+'_scl2 = vortdat[i*10+8]\ns'+sol+'_scl3 = vortdat[i*10+9]\n')
    
    # BIVARIATE SPLINE INTERPOLATION
    
    iz = np.size(sol_d)
    jz = np.size(tsr_d)
    
    # Initializing rectangular matrices
    z_loc1 = np.zeros((iz, jz))
    z_loc2 = np.zeros((iz, jz))
    z_loc3 = np.zeros((iz, jz))
    z_spr1 = np.zeros((iz, jz))
    z_spr2 = np.zeros((iz, jz))
    z_skw1 = np.zeros((iz, jz))
    z_skw2 = np.zeros((iz, jz))
    z_scl1 = np.zeros((iz, jz))
    z_scl2 = np.zeros((iz, jz))
    z_scl3 = np.zeros((iz, jz))
    
    # Transferring raw data into rectangular matrices
    for i in range(iz):
        for j in range(jz):
            sol = str(i+1)
            exec('z_loc1[i, j] = s'+sol+'_loc1[j]')
            exec('z_loc2[i, j] = s'+sol+'_loc2[j]')
            exec('z_loc3[i, j] = s'+sol+'_loc3[j]')
            exec('z_spr1[i, j] = s'+sol+'_spr1[j]')
            exec('z_spr2[i, j] = s'+sol+'_spr2[j]')
            exec('z_skw1[i, j] = s'+sol+'_skw1[j]')
            exec('z_skw2[i, j] = s'+sol+'_skw2[j]')
            exec('z_scl1[i, j] = s'+sol+'_scl1[j]')
            exec('z_scl2[i, j] = s'+sol+'_scl2[j]')
            exec('z_scl3[i, j] = s'+sol+'_scl3[j]')
    
    # Creating a rectangular bivariate spline of the parameter data
    s_loc1 = RectBivariateSpline(sol_d, tsr_d, z_loc1)
    s_loc2 = RectBivariateSpline(sol_d, tsr_d, z_loc2)
    s_loc3 = RectBivariateSpline(sol_d, tsr_d, z_loc3)
    s_spr1 = RectBivariateSpline(sol_d, tsr_d, z_spr1)
    s_spr2 = RectBivariateSpline(sol_d, tsr_d, z_spr2)
    s_skw1 = RectBivariateSpline(sol_d, tsr_d, z_skw1)
    s_skw2 = RectBivariateSpline(sol_d, tsr_d, z_skw2)
    s_scl1 = RectBivariateSpline(sol_d, tsr_d, z_scl1)
    s_scl2 = RectBivariateSpline(sol_d, tsr_d, z_scl2)
    s_scl3 = RectBivariateSpline(sol_d, tsr_d, z_scl3)
    
    # Selecting the specific parameters to use for TSR and solidity
    loc1 = s_loc1(solidity, tsr)
    loc2 = s_loc2(solidity, tsr)
    loc3 = s_loc3(solidity, tsr)
    spr1 = s_spr1(solidity, tsr)
    spr2 = s_spr2(solidity, tsr)
    skw1 = s_skw1(solidity, tsr)
    skw2 = s_skw2(solidity, tsr)
    scl1 = s_scl1(solidity, tsr)
    scl2 = s_scl2(solidity, tsr)
    scl3 = s_scl3(solidity, tsr)
    
    # Creating arrays of the parameters
    loc = np.array([loc1[0, 0], loc2[0, 0], loc3[0, 0]])
    spr = np.array([spr1[0, 0], spr2[0, 0]])
    skw = np.array([skw1[0, 0], skw2[0, 0]])
    scl = np.array([scl1[0, 0], scl2[0, 0], scl3[0, 0]])
    
    return loc, spr, skw, scl


def velocity(tsr, solidity):
    """
    Using SMG distribution parameters to define the velocity strength and shape
    
    Parameters
    ----------
    tsr : float
        tip-speed ratio
    solidity : float
        turbine solidity
    
    Returns
    ----------
    men : array
        array of the mean parameter (1 values)
    sdv : array
        array of the standard deviation parameter (4 values)
    rat : array
        array of the rate parameter (1 values)
    wdt : array
        array of the translation parameter (1 values)
    spr : array
        array of the spread parameter (4 values)
    scl : array
        array of the scale parameter (3 values)
    """
    # Reading in csv file (vorticity database)
    basepath = path.join(path.dirname(path.realpath(__file__)), 'data')
    # fdata = basepath + path.sep + 'velodatabase_SMG_surf2_edit.csv'
    fdata = basepath + path.sep + 'velodatabase_SMG_surf4.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw, 0)
            velodat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d, float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw, 0)
            velodat = np.vstack([velodat, raw]) # adding entry to vorticity database array
        i += 1
    f.close()

    velodat = np.delete(velodat, 0, axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    velodat = velodat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each SMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_men = velodat[i*14]\ns'+sol+'_sdv1 = velodat[i*14+1]\ns'+sol+'_sdv2 = velodat[i*14+2]\ns'+sol+'_sdv3 = velodat[i*14+3]\ns'+sol+'_sdv4 = velodat[i*14+4]\ns'+sol+'_rat = velodat[i*14+5]\ns'+sol+'_wdt = velodat[i*14+6]\ns'+sol+'_spr1 = velodat[i*14+7]\ns'+sol+'_spr2 = velodat[i*14+8]\ns'+sol+'_spr3 = velodat[i*14+9]\ns'+sol+'_spr4 = velodat[i*14+10]\ns'+sol+'_scl1 = velodat[i*14+11]\ns'+sol+'_scl2 = velodat[i*14+12]\ns'+sol+'_scl3 = velodat[i*14+13]\n')

    # NEAREST ND INTERPOLATION

    iz = np.size(sol_d)
    jz = np.size(tsr_d)

    # Transferring parameter arrays to data matrix
    dataset = np.zeros((14*iz,jz))

    for i in range(iz):
        sol = str(i+1)
        exec('dataset[i*14] = s'+sol+'_men')
        exec('dataset[i*14+1] = s'+sol+'_sdv1')
        exec('dataset[i*14+2] = s'+sol+'_sdv2')
        exec('dataset[i*14+3] = s'+sol+'_sdv3')
        exec('dataset[i*14+4] = s'+sol+'_sdv4')
        exec('dataset[i*14+5] = s'+sol+'_rat')
        exec('dataset[i*14+6] = s'+sol+'_wdt')
        exec('dataset[i*14+7] = s'+sol+'_spr1')
        exec('dataset[i*14+8] = s'+sol+'_spr2')
        exec('dataset[i*14+9] = s'+sol+'_spr3')
        exec('dataset[i*14+10] = s'+sol+'_spr4')
        exec('dataset[i*14+11] = s'+sol+'_scl1')
        exec('dataset[i*14+12] = s'+sol+'_scl2')
        exec('dataset[i*14+13] = s'+sol+'_scl3')

    # Transferring data matrix to TSR, soldity, and parameter arrays
    points1 = np.array([])
    points2 = np.array([])
    points3 = np.array([])
    values = np.array([])
    for i in range(jz):
        for j in range(iz):
            for k in range(14):
                points1 = np.append(points1,tsr_d[i])
                points2 = np.append(points2,sol_d[j])
                points3 = np.append(points3,k+1)
                values = np.append(values,dataset[j*14+k,i])

    # Identifying nearest ND point to CFD data set
    met = 'nearest'
    men = griddata((points1,points2,points3),values,(tsr,solidity,1),method=met)
    sdv1 = griddata((points1,points2,points3),values,(tsr,solidity,2),method=met)
    sdv2 = griddata((points1,points2,points3),values,(tsr,solidity,3),method=met)
    sdv3 = griddata((points1,points2,points3),values,(tsr,solidity,4),method=met)
    sdv4 = griddata((points1,points2,points3),values,(tsr,solidity,5),method=met)
    rat = griddata((points1,points2,points3),values,(tsr,solidity,6),method=met)
    wdt = griddata((points1,points2,points3),values,(tsr,solidity,7),method=met)
    spr1 = griddata((points1,points2,points3),values,(tsr,solidity,8),method=met)
    spr2 = griddata((points1,points2,points3),values,(tsr,solidity,9),method=met)
    spr3 = griddata((points1,points2,points3),values,(tsr,solidity,10),method=met)
    spr4 = griddata((points1,points2,points3),values,(tsr,solidity,11),method=met)
    scl1 = griddata((points1,points2,points3),values,(tsr,solidity,12),method=met)
    scl2 = griddata((points1,points2,points3),values,(tsr,solidity,13),method=met)
    scl3 = griddata((points1,points2,points3),values,(tsr,solidity,14),method=met)

    # Identifying nearest TSR and solidity to given variables
    pointst = np.copy(points1)/tsr - 1.
    pointss = np.copy(points2)/solidity - 1.
    ti = np.argmin(fabs(pointst))
    si = np.argmin(fabs(pointss))
    tsrn = points1[ti]
    soln = points2[si]

    # Creating arrays of the parameters
    men = np.array([men])
    sdv = np.array([sdv1, sdv2, sdv3, sdv4])
    rat = np.array([rat])
    wdt = np.array([wdt])
    spr = np.array([spr1, spr2, spr3, spr4])
    scl = np.array([scl1, scl2, scl3])

    return men,sdv,rat,wdt,spr,scl,tsrn,soln


def overlay(xt,ys,tsr,sol):
    ntsr = np.size(tsr)
    nsol = np.size(sol)

    dtsr = 'f1 = 0.'
    for i in range(ntsr):
        tpow = '+'+str(tsr[i])+'*xt**'+str(ntsr-1-i)
        dtsr = dtsr+tpow

    dsol = 'f2 = 0.'
    for i in range(nsol):
        spow = '+'+str(sol[i])+'*ys**'+str(nsol-1-i)
        dsol = dsol+spow

    exec(dtsr)
    exec(dsol)
    return f1*f2


def velocity2(tsr,solidity):
    ment = np.array([0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00])
    sdv1t = np.array([-5.49968724e-03,9.61937568e-02,-5.54670889e-01,2.49196469e+00])
    sdv2t = np.array([1.19177362e-03,-2.97620621e-03,3.24915036e-04,8.88376997e-03])
    sdv3t = np.array([-5.32550486e-01,1.36796928e+01,-1.36290158e+02,4.88975685e+02])
    sdv4t = np.array([-3.00669046e-01,3.62985562e+00,-1.14149652e+01,1.42138555e+01])
    ratt = np.array([7.79857758e-01,-9.72827823e+00,3.57499750e+01,-1.77061536e+01])
    tnst = np.array([2.82323961e-02,-4.01600512e-01,1.49054030e+00,9.59240910e-01])
    spr1t = np.array([4.75980490e-03,-8.13347861e-02,5.00350642e-01,2.56944711e+00])
    spr2t = np.array([1.03772552e-03,-2.43729195e-02,-9.16216810e-02,2.39188975e-02])
    spr3t = np.array([2.00064577e-01,-1.05319231e+01,9.91461423e+01,-2.93841893e+02])
    spr4t = np.array([2.79059938e-01,-3.23201793e+00,1.15980962e+01,-5.87869620e+00])
    scl1t = np.array([1.55734596e-03,-7.23612213e-03,-9.40544723e-03,1.06726267e+00])
    scl2t = np.array([-5.20639533e-04,4.97916878e-04,1.17604650e-01,-2.55380782e-02])
    scl3t = np.array([-2.88231112e-02,2.72948674e-01,-4.01043675e-01,6.63626086e-01])
    mens = np.array([0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00])
    sdv1s = np.array([-3.72447988e-01,-1.98990833e-01,5.93913665e-02,5.14618501e-01])
    sdv2s = np.array([2.39556464e-01,8.88174262e-02,-3.95793113e-01,1.45938185e-01])
    sdv3s = np.array([1.86558487e+03,2.56113269e+02,-2.92022347e+02,2.59231995e+01])
    sdv4s = np.array([-3.70159535e+00,5.82354220e+00,-2.91402824e+00,4.59207298e-01])
    rats = np.array([-9.71306609e+02,1.34584010e+03,-3.72179434e+02,3.57344765e+01])
    tnss = np.array([-9.57522091e+00,1.76555533e+01,-1.51151040e+01,9.93503682e+00])
    spr1s = np.array([-3.10057176e+00,4.95103006e+00,-1.18363880e+00,4.31493696e-01])
    spr2s = np.array([-3.84261851e-01,3.23932705e+00,-7.84005884e-01,8.41289090e-02])
    spr3s = np.array([-2.40700774e+02,5.62918841e+02,-4.37946480e+02,1.16580190e+02])
    spr4s = np.array([-1.89233223e+02,2.95952219e+02,-1.04029722e+02,5.66542903e+00])
    scl1s = np.array([-7.00286183e-01,1.27314224e+00,-7.64015056e-01,5.71113239e-01])
    scl2s = np.array([-2.94973127e-02,4.51269816e-01,6.77950402e-01,7.15663584e-02])
    scl3s = np.array([2.91410061e+00,1.10869146e+01,-1.76449554e+01,6.10637931e+00])


    ment = np.array([0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00])
    sdv1t = np.array([4.80702257e-04,8.03681634e-02,-2.36330244e-01,7.29863231e-01])
    sdv2t = np.array([2.11285630e-04,-3.16101679e-04,-1.26270146e-02,6.57959715e-02])
    sdv3t = np.array([2.22716598e-01,-9.04256829e-01,-2.49105240e-01,1.19056348e+01])
    sdv4t = np.array([-8.04309278e-03,1.32885513e-01,-8.19131933e-01,2.11400704e+00])
    ratt = np.array([2.43433069e-01,-2.72968444e+00,8.62008306e+00,1.85250887e+00])
    tnst = np.array([1.59041935e-02,-2.32827659e-01,8.94192129e-01,8.86760524e-01])
    spr1t = np.array([1.17590965e-03,-1.34828411e-02,5.15300934e-02,-4.35065337e-02])
    spr2t = np.array([-1.11474542e-02,1.34769685e-01,-2.04538671e-01,5.45610682e-01])
    spr3t = np.array([-5.24259374e-02,4.79679128e-01,-7.05366229e-01,-3.36233596e+00])
    spr4t = np.array([7.87168914e-02,-8.17323003e-01,2.50119222e+00,2.96791153e-01])
    scl1t = np.array([2.11976958e-03,-2.87168152e-02,1.43251744e-01,2.10307240e-01])
    scl2t = np.array([9.03104823e-04,-1.34569375e-02,1.42934606e-01,-5.21819331e-02])
    scl3t = np.array([-4.50618047e-02,5.51811159e-01,-1.13434268e+00,1.21880973e+00])
    mens = np.array([0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00])
    sdv1s = np.array([-4.27178664e-01,-6.05098542e-01,1.20187349e+00,1.22003706e+00])
    sdv2s = np.array([1.35941967e-01,-1.09329420e-01,-2.70402400e-01,2.83947389e-01])
    sdv3s = np.array([3.33777885e+02,-5.53838971e+01,1.38126011e+02,9.11675842e+00])
    sdv4s = np.array([-4.63841806e-04,3.59121619e+00,-5.96370397e+00,-1.36286088e-01])
    rats = np.array([-1.61498707e+02,2.43675212e+02,-3.05919484e+01,5.38578311e+00])
    tnss = np.array([1.44203370e+00,-5.77706129e+00,1.64940309e+00,6.33296935e+00])
    spr1s = np.array([8.58623428e-01,-5.35149908e-01,9.44068982e-02,1.20078844e-03])
    spr2s = np.array([2.00123699e+00,-4.29269027e+00,3.67931690e+00,-1.27519122e-01])
    spr3s = np.array([-1.14940716e+01,2.40033024e+01,1.47769155e+00,2.79175293e+00])
    spr4s = np.array([-2.20808946e+01,4.85630658e+01,-4.50017925e+00,9.87259794e-01])
    scl1s = np.array([4.72026239e-01,-7.05406384e-01,1.40982351e-01,9.75327984e-01])
    scl2s = np.array([1.90852734e-01,5.25471845e-02,1.08647122e+00,6.03490045e-02])
    scl3s = np.array([-2.98976936e+00,1.26754049e+01,-1.27888937e+01,4.75385199e+00])


    ment = np.array([0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00])
    sdv1t = np.array([0.00000000e+00,-2.14161553e-03,-4.21556328e-03, 1.48152720e-01])
    sdv2t = np.array([0.00000000e+00, 4.57006047e-03, 2.24110496e-02, 4.92413935e-01])
    sdv3t = np.array([-5.28369268e-03,-9.89053847e-02,-4.21273359e-01, 1.45167187e+01])
    sdv4t = np.array([0.00000000e+00, 1.79452090e-02,-8.88797924e-02, 2.03963797e-03])
    ratt = np.array([7.57992485e-01,-1.07514453e+01, 2.84371459e+01, 7.21470835e+01])
    tnst = np.array([1.36090201e-01,-1.41933971e+00, 4.95510763e+00,-4.10649797e+00])
    spr1t = np.array([0.00000000e+00, 4.15053874e-02,-1.21436212e-01, 9.01176495e-02])
    spr2t = np.array([-9.51987791e-05,-2.42813363e-03,-3.14921447e-02, 4.27419977e-01])
    spr3t = np.array([-8.57669739e-01, 8.15743937e-02, 9.68978378e+00,-1.67991922e+02])
    spr4t = np.array([-2.31984687e-02,-1.98556659e-02, 8.14305825e-01, 6.04599309e+00])
    scl1t = np.array([0.00000000e+00, 1.09103873e-02,-1.57044023e-01, 5.81535656e-01])
    scl2t = np.array([0.00000000e+00, 1.85914844e-02,-1.75816877e-02, 5.47783767e-01])
    scl3t = np.array([-2.40071938e-01, 1.51012911e+00, 5.34643308e-01, 7.09568655e+00])
    mens = np.array([0.00000000e+00,0.00000000e+00,0.00000000e+00,0.00000000e+00])
    sdv1s = np.array([-8.72642010e-02,-2.04665371e-02, 7.25709019e-02, 1.94550404e-01])
    sdv2s = np.array([-6.24110794e-01,-2.49915441e-01, 7.80670041e-03, 4.60330311e-01])
    sdv3s = np.array([-4.15481030e+00, 3.02940787e+00, 7.80647011e+00, 2.12571927e+01])
    sdv4s = np.array([8.28483625e-01, 3.19263507e-01,-8.78338367e-01,-1.76082841e-01])
    rats = np.array([-1.24839749e+02, 1.70780462e+02,-4.54202150e+01, 7.23787704e+00])
    tnss = np.array([-3.39515268e+01, 4.51152171e+01,-1.23404907e+01, 1.56008702e+00])
    spr1s = np.array([-1.76549905e+00, 1.18934484e+00, 1.00749787e+00, 1.43673970e-01])
    spr2s = np.array([1.03616521e+00, 9.63793185e-02, 9.14732634e-02, 2.60096317e-01])
    spr3s = np.array([-2.83189876e+00,-5.95341735e+00,-2.37347088e+00, 2.43989236e+01])
    spr4s = np.array([7.25152289e+00, 1.33620627e+01, 2.60022973e+01, 1.16939523e+01])
    scl1s = np.array([2.81442143e+00,-3.79403885e+00, 1.46089819e+00, 2.39420036e-02])
    scl2s = np.array([-1.12001008e+00, 6.18471752e-01,-3.16135317e-01, 4.16299140e-01])
    scl3s = np.array([6.32022229e+00,-9.48437842e+01, 8.86383530e+01, 1.56930523e+01])

    men = overlay(tsr,solidity,ment,mens)
    sdv1 = overlay(tsr,solidity,sdv1t,sdv1s)
    sdv2 = overlay(tsr,solidity,sdv2t,sdv2s)
    sdv3 = overlay(tsr,solidity,sdv3t,sdv3s)
    sdv4 = overlay(tsr,solidity,sdv4t,sdv4s)
    rat = overlay(tsr,solidity,ratt,rats)
    tns = overlay(tsr,solidity,tnst,tnss)
    spr1 = overlay(tsr,solidity,spr1t,spr1s)
    spr2 = overlay(tsr,solidity,spr2t,spr2s)
    spr3 = overlay(tsr,solidity,spr3t,spr3s)
    spr4 = overlay(tsr,solidity,spr4t,spr4s)
    scl1 = overlay(tsr,solidity,scl1t,scl1s)
    scl2 = overlay(tsr,solidity,scl2t,scl2s)
    scl3 = overlay(tsr,solidity,scl3t,scl3s)

    # Creating arrays of the parameters
    men = np.array([men])
    sdv = np.array([sdv1, sdv2, sdv3, sdv4])
    rat = np.array([rat])
    wdt = np.array([tns])
    spr = np.array([spr1, spr2, spr3, spr4])
    scl = np.array([scl1, scl2, scl3])

    return men,sdv,rat,wdt,spr,scl