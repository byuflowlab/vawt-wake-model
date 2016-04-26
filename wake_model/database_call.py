import csv
from os import path
import numpy as np
from scipy.interpolate import RectBivariateSpline,SmoothBivariateSpline,griddata

# from matplotlib import rcParams
# rcParams['font.family'] = 'Times New Roman'


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
    
    # BIVARIATE SPLINE FITTING
    
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
        array of the mean parameter (3 values)
    spr : array
        array of the spread parameter (3 values)
    scl : array
        array of the scale parameter (3 values)
    rat : array
        array of the rate parameter (2 values)
    tns : array
        array of the translation parameter (2 values)
    """
    # Reading in csv file (vorticity database)
    basepath = path.join(path.dirname(path.realpath(__file__)), 'data')
    # fdata = basepath + path.sep + 'velodatabase_SMG.csv'
    fdata = basepath + path.sep + 'velodatabase_SMG4.csv'
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
        
        exec('s'+sol+'_men1 = velodat[i*14]\ns'+sol+'_men2 = velodat[i*14+1]\ns'+sol+'_men3 = velodat[i*14+2]\ns'+sol+'_spr1 = velodat[i*14+3]\ns'+sol+'_spr2 = velodat[i*14+4]\ns'+sol+'_spr3 = velodat[i*14+5]\ns'+sol+'_spr4 = velodat[i*14+6]\ns'+sol+'_scl1 = velodat[i*14+7]\ns'+sol+'_scl2 = velodat[i*14+8]\ns'+sol+'_scl3 = velodat[i*14+9]\ns'+sol+'_rat1 = velodat[i*14+10]\ns'+sol+'_rat2 = velodat[i*14+11]\ns'+sol+'_tns1 = velodat[i*14+12]\ns'+sol+'_tns2 = velodat[i*14+13]\n')

    spline = True
    # spline = False
    grid = True
    grid = False

    if spline == True:
        # BIVARIATE SPLINE FITTING

        iz = np.size(sol_d)
        jz = np.size(tsr_d)

        # Initializing rectangular matrices
        z_men1 = np.zeros((iz,  jz))
        z_men2 = np.zeros((iz, jz))
        z_men3 = np.zeros((iz, jz))
        z_spr1 = np.zeros((iz, jz))
        z_spr2 = np.zeros((iz, jz))
        z_spr3 = np.zeros((iz, jz))
        z_spr4 = np.zeros((iz, jz))
        z_scl1 = np.zeros((iz, jz))
        z_scl2 = np.zeros((iz, jz))
        z_scl3 = np.zeros((iz, jz))
        z_rat1 = np.zeros((iz, jz))
        z_rat2 = np.zeros((iz, jz))
        z_tns1 = np.zeros((iz, jz))
        z_tns2 = np.zeros((iz, jz))

        # Transferring raw data into rectangular matrices
        for i in range(iz):
            for j in range(jz):
                sol = str(i+1)
                exec('z_men1[i,j] = s'+sol+'_men1[j]')
                exec('z_men2[i,j] = s'+sol+'_men2[j]')
                exec('z_men3[i,j] = s'+sol+'_men3[j]')
                exec('z_spr1[i,j] = s'+sol+'_spr1[j]')
                exec('z_spr2[i,j] = s'+sol+'_spr2[j]')
                exec('z_spr3[i,j] = s'+sol+'_spr3[j]')
                exec('z_spr4[i,j] = s'+sol+'_spr4[j]')
                exec('z_scl1[i,j] = s'+sol+'_scl1[j]')
                exec('z_scl2[i,j] = s'+sol+'_scl2[j]')
                exec('z_scl3[i,j] = s'+sol+'_scl3[j]')
                exec('z_rat1[i,j] = s'+sol+'_rat1[j]')
                exec('z_rat2[i,j] = s'+sol+'_rat2[j]')
                exec('z_tns1[i,j] = s'+sol+'_tns1[j]')
                exec('z_tns2[i,j] = s'+sol+'_tns2[j]')

        # Creating a rectangular bivariate spline of the parameter data
        # if tsr > 2.:
        #     xs = 3
        #     ys = 3
        # else:
        #     xs = 1
        #     ys = 1
        xs = 3
        ys = 3
        sm = 0
        s_men1 = RectBivariateSpline(sol_d, tsr_d, z_men1, kx=xs, ky=ys, s=sm)
        s_men2 = RectBivariateSpline(sol_d, tsr_d, z_men2, kx=xs, ky=ys, s=sm)
        s_men3 = RectBivariateSpline(sol_d, tsr_d, z_men3, kx=xs, ky=ys, s=sm)
        s_spr1 = RectBivariateSpline(sol_d, tsr_d, z_spr1, kx=xs, ky=ys, s=sm)
        s_spr2 = RectBivariateSpline(sol_d, tsr_d, z_spr2, kx=xs, ky=ys, s=sm)
        s_spr3 = RectBivariateSpline(sol_d, tsr_d, z_spr3, kx=xs, ky=ys, s=sm)
        s_spr4 = RectBivariateSpline(sol_d, tsr_d, z_spr4, kx=xs, ky=ys, s=sm)
        s_scl1 = RectBivariateSpline(sol_d, tsr_d, z_scl1, kx=xs, ky=ys, s=sm)
        s_scl2 = RectBivariateSpline(sol_d, tsr_d, z_scl2, kx=xs, ky=ys, s=sm)
        s_scl3 = RectBivariateSpline(sol_d, tsr_d, z_scl3, kx=xs, ky=ys, s=sm)
        s_rat1 = RectBivariateSpline(sol_d, tsr_d, z_rat1, kx=xs, ky=ys, s=sm)
        s_rat2 = RectBivariateSpline(sol_d, tsr_d, z_rat2, kx=xs, ky=ys, s=sm)
        s_tns1 = RectBivariateSpline(sol_d, tsr_d, z_tns1, kx=xs, ky=ys, s=sm)
        s_tns2 = RectBivariateSpline(sol_d, tsr_d, z_tns2, kx=xs, ky=ys, s=sm)

        # Selecting the specific parameters to use for TSR and solidity
        men1 = s_men1(solidity, tsr)
        men2 = s_men2(solidity, tsr)
        men3 = s_men3(solidity, tsr)
        spr1 = s_spr1(solidity, tsr)
        spr2 = s_spr2(solidity, tsr)
        spr3 = s_spr3(solidity, tsr)
        spr4 = s_spr4(solidity, tsr)
        scl1 = s_scl1(solidity, tsr)
        scl2 = s_scl2(solidity, tsr)
        scl3 = s_scl3(solidity, tsr)
        rat1 = s_rat1(solidity, tsr)
        rat2 = s_rat2(solidity, tsr)
        tns1 = s_tns1(solidity, tsr)
        tns2 = s_tns2(solidity, tsr)

        # Creating arrays of the parameters
        men = np.array([men1[0, 0], men2[0, 0], men3[0, 0]])
        spr = np.array([spr1[0, 0], spr2[0, 0], spr3[0, 0], spr4[0, 0]])
        scl = np.array([scl1[0, 0], scl2[0, 0], scl3[0, 0]])
        rat = np.array([rat1[0, 0], rat2[0, 0]])
        tns = np.array([tns1[0, 0], tns2[0, 0]])

    elif grid == True:
        iz = np.size(sol_d)
        jz = np.size(tsr_d)

        dataset = np.zeros((14*iz,jz))

        for i in range(iz):
            sol = str(i+1)
            exec('dataset[i*14] = s'+sol+'_men1')
            exec('dataset[i*14+1] = s'+sol+'_men2')
            exec('dataset[i*14+2] = s'+sol+'_men3')
            exec('dataset[i*14+3] = s'+sol+'_spr1')
            exec('dataset[i*14+4] = s'+sol+'_spr2')
            exec('dataset[i*14+5] = s'+sol+'_spr3')
            exec('dataset[i*14+6] = s'+sol+'_spr4')
            exec('dataset[i*14+7] = s'+sol+'_scl1')
            exec('dataset[i*14+8] = s'+sol+'_scl2')
            exec('dataset[i*14+9] = s'+sol+'_scl3')
            exec('dataset[i*14+10] = s'+sol+'_rat1')
            exec('dataset[i*14+11] = s'+sol+'_rat2')
            exec('dataset[i*14+12] = s'+sol+'_tns1')
            exec('dataset[i*14+13] = s'+sol+'_tns2')

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

        met = 'nearest'
        # met = 'linear'
        men1 = griddata((points1,points2,points3),values,(tsr,solidity,1),method=met)
        men2 = griddata((points1,points2,points3),values,(tsr,solidity,2),method=met)
        men3 = griddata((points1,points2,points3),values,(tsr,solidity,3),method=met)
        spr1 = griddata((points1,points2,points3),values,(tsr,solidity,4),method=met)
        spr2 = griddata((points1,points2,points3),values,(tsr,solidity,5),method=met)
        spr3 = griddata((points1,points2,points3),values,(tsr,solidity,6),method=met)
        spr4 = griddata((points1,points2,points3),values,(tsr,solidity,7),method=met)
        scl1 = griddata((points1,points2,points3),values,(tsr,solidity,8),method=met)
        scl2 = griddata((points1,points2,points3),values,(tsr,solidity,9),method=met)
        scl3 = griddata((points1,points2,points3),values,(tsr,solidity,10),method=met)
        rat1 = griddata((points1,points2,points3),values,(tsr,solidity,11),method=met)
        rat2 = griddata((points1,points2,points3),values,(tsr,solidity,12),method=met)
        tns1 = griddata((points1,points2,points3),values,(tsr,solidity,13),method=met)
        tns2 = griddata((points1,points2,points3),values,(tsr,solidity,14),method=met)

        # Creating arrays of the parameters
        men = np.array([men1, men2, men3])
        spr = np.array([spr1, spr2, spr3, spr4])
        scl = np.array([scl1, scl2, scl3])
        rat = np.array([rat1, rat2])
        tns = np.array([tns1, tns2])




    return men, spr, scl, rat, tns


def quad(tsr, solidity):
    """
    Using quadratic distribution parameters to define the velocity strength and shape
    
    Parameters
    ----------
    tsr : float
        tip-speed ratio
    solidity : float
        turbine solidity
    
    Returns
    ----------
    scl : array
        array of the scale parameter (3 values)
    trn : array
        array of the translation parameter (3 values)
    """
    # Reading in csv file (vorticity database)
    # basepath = path.join(path.dirname(path.realpath(__file__)), 'data')
    # fdata = basepath + path.sep + 'velodatabase_2.csv'
    fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/data/velodatabase3_4_edit.csv'
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
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_scl1 = velodat[i*6]\ns'+sol+'_scl2 = velodat[i*6+1]\ns'+sol+'_scl3 = velodat[i*6+2]\ns'+sol+'_trn1 = velodat[i*6+3]\ns'+sol+'_trn2 = velodat[i*6+4]\ns'+sol+'_trn3 = velodat[i*6+5]\n')
    
    # BIVARIATE SPLINE FITTING
    
    iz = np.size(sol_d)
    jz = np.size(tsr_d)
    
    # Initializing rectangular matrices
    z_scl1 = np.zeros((iz, jz))
    z_scl2 = np.zeros((iz, jz))
    z_scl3 = np.zeros((iz, jz))
    z_trn1 = np.zeros((iz, jz))
    z_trn2 = np.zeros((iz, jz))
    z_trn3 = np.zeros((iz, jz))
    
    # Transferring raw data into rectangular matrices
    for i in range(iz):
        for j in range(jz):
            sol = str(i+1)
            exec('z_scl1[i, j] = s'+sol+'_scl1[j]')
            exec('z_scl2[i, j] = s'+sol+'_scl2[j]')
            exec('z_scl3[i, j] = s'+sol+'_scl3[j]')
            exec('z_trn1[i, j] = s'+sol+'_trn1[j]')
            exec('z_trn2[i, j] = s'+sol+'_trn2[j]')
            exec('z_trn3[i, j] = s'+sol+'_trn3[j]')
    
    # Creating a rectangular bivariate spline of the parameter data
    s_scl1 = RectBivariateSpline(sol_d, tsr_d, z_scl1)
    s_scl2 = RectBivariateSpline(sol_d, tsr_d, z_scl2)
    s_scl3 = RectBivariateSpline(sol_d, tsr_d, z_scl3)
    s_trn1 = RectBivariateSpline(sol_d, tsr_d, z_trn1)
    s_trn2 = RectBivariateSpline(sol_d, tsr_d, z_trn2)
    s_trn3 = RectBivariateSpline(sol_d, tsr_d, z_trn3)
    
    # Selecting the specific parameters to use for TSR and solidity
    scl1 = s_scl1(solidity, tsr)
    scl2 = s_scl2(solidity, tsr)
    scl3 = s_scl3(solidity, tsr)
    trn1 = s_trn1(solidity, tsr)
    trn2 = s_trn2(solidity, tsr)
    trn3 = s_trn3(solidity, tsr)
    
    # Creating arrays of the parameters
    scl = np.array([scl1[0, 0], scl2[0, 0], scl3[0, 0]])
    trn = np.array([trn1[0, 0], trn2[0, 0], trn3[0, 0]])

    return scl, trn
