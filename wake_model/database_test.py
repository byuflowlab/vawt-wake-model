import numpy as np
import matplotlib.pyplot as plt
import csv
from os import path
from scipy.interpolate import RectBivariateSpline


def vorticity(tsr,solidity):
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
    basepath = path.join(path.dirname(path.realpath(__file__)),'data')
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
            sol_d = np.append(sol_d,float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            vortdat = np.vstack([vortdat,raw]) # adding entry to vorticity database array
        i += 1
    f.close()
    
    vortdat = np.delete(vortdat,(0),axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    vortdat = vortdat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_loc1 = vortdat[i*10]\ns'+sol+'_loc2 = vortdat[i*10+1]\ns'+sol+'_loc3 = vortdat[i*10+2]\ns'+sol+'_spr1 = vortdat[i*10+3]\ns'+sol+'_spr2 = vortdat[i*10+4]\ns'+sol+'_skw1 = vortdat[i*10+5]\ns'+sol+'_skw2 = vortdat[i*10+6]\ns'+sol+'_scl1 = vortdat[i*10+7]\ns'+sol+'_scl2 = vortdat[i*10+8]\ns'+sol+'_scl3 = vortdat[i*10+9]\n')
    
    # BIVARIATE SPLINE FITTING
    
    iz = np.size(sol_d)
    jz = np.size(tsr_d)
    
    # Initializing rectangular matrices
    Z_loc1 = np.zeros((iz,jz))
    Z_loc2 = np.zeros((iz,jz))
    Z_loc3 = np.zeros((iz,jz))
    Z_spr1 = np.zeros((iz,jz))
    Z_spr2 = np.zeros((iz,jz))
    Z_skw1 = np.zeros((iz,jz))
    Z_skw2 = np.zeros((iz,jz))
    Z_scl1 = np.zeros((iz,jz))
    Z_scl2 = np.zeros((iz,jz))
    Z_scl3 = np.zeros((iz,jz))
    
    # Transferring raw data into rectangular matrices
    for i in range(iz):
        for j in range(jz):
            sol = str(i+1)
            exec('Z_loc1[i,j] = s'+sol+'_loc1[j]')
            exec('Z_loc2[i,j] = s'+sol+'_loc2[j]')
            exec('Z_loc3[i,j] = s'+sol+'_loc3[j]')
            exec('Z_spr1[i,j] = s'+sol+'_spr1[j]')
            exec('Z_spr2[i,j] = s'+sol+'_spr2[j]')
            exec('Z_skw1[i,j] = s'+sol+'_skw1[j]')
            exec('Z_skw2[i,j] = s'+sol+'_skw2[j]')
            exec('Z_scl1[i,j] = s'+sol+'_scl1[j]')
            exec('Z_scl2[i,j] = s'+sol+'_scl2[j]')
            exec('Z_scl3[i,j] = s'+sol+'_scl3[j]')
    
    # Creating a rectangular bivariate spline of the parameter data
    s_loc1 = RectBivariateSpline(sol_d,tsr_d,Z_loc1)
    s_loc2 = RectBivariateSpline(sol_d,tsr_d,Z_loc2)
    s_loc3 = RectBivariateSpline(sol_d,tsr_d,Z_loc3)
    s_spr1 = RectBivariateSpline(sol_d,tsr_d,Z_spr1)
    s_spr2 = RectBivariateSpline(sol_d,tsr_d,Z_spr2)
    s_skw1 = RectBivariateSpline(sol_d,tsr_d,Z_skw1)
    s_skw2 = RectBivariateSpline(sol_d,tsr_d,Z_skw2)
    s_scl1 = RectBivariateSpline(sol_d,tsr_d,Z_scl1)
    s_scl2 = RectBivariateSpline(sol_d,tsr_d,Z_scl2)
    s_scl3 = RectBivariateSpline(sol_d,tsr_d,Z_scl3)
    
    # Selecting the specific parameters to use for TSR and solidity
    loc1 = s_loc1(solidity,tsr)
    loc2 = s_loc2(solidity,tsr)
    loc3 = s_loc3(solidity,tsr)
    spr1 = s_spr1(solidity,tsr)
    spr2 = s_spr2(solidity,tsr)
    skw1 = s_skw1(solidity,tsr)
    skw2 = s_skw2(solidity,tsr)
    scl1 = s_scl1(solidity,tsr)
    scl2 = s_scl2(solidity,tsr)
    scl3 = s_scl3(solidity,tsr)
    
    # Creating arrays of the parameters
    loc = np.array([loc1[0,0],loc2[0,0],loc3[0,0]])
    spr = np.array([spr1[0,0],spr2[0,0]])
    skw = np.array([skw1[0,0],skw2[0,0]])
    scl = np.array([scl1[0,0],scl2[0,0],scl3[0,0]])
    
    return loc,spr,skw,scl


def velocity(tsr,solidity):
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
    # basepath = path.join(path.dirname(path.realpath(__file__)),'data')
    # fdata = basepath + path.sep + 'velodatabase_2.csv'
    fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/data/velodatabase.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw,0)
            velodat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d,float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            velodat = np.vstack([velodat,raw]) # adding entry to vorticity database array
        i += 1
    f.close()
    
    velodat = np.delete(velodat,(0),axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    velodat = velodat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('s'+sol+'_men1 = velodat[i*13]\ns'+sol+'_men2 = velodat[i*13+1]\ns'+sol+'_men3 = velodat[i*13+2]\ns'+sol+'_spr1 = velodat[i*13+3]\ns'+sol+'_spr2 = velodat[i*13+4]\ns'+sol+'_spr3 = velodat[i*13+5]\ns'+sol+'_scl1 = velodat[i*13+6]\ns'+sol+'_scl2 = velodat[i*13+7]\ns'+sol+'_scl3 = velodat[i*13+8]\ns'+sol+'_rat1 = velodat[i*13+9]\ns'+sol+'_rat2 = velodat[i*13+10]\ns'+sol+'_tns1 = velodat[i*13+11]\ns'+sol+'_tns2 = velodat[i*13+12]\n')
    
    # BIVARIATE SPLINE FITTING
    
    iz = np.size(sol_d)
    jz = np.size(tsr_d)
    
    # Initializing rectangular matrices
    Z_men1 = np.zeros((iz,jz))
    Z_men2 = np.zeros((iz,jz))
    Z_men3 = np.zeros((iz,jz))
    Z_spr1 = np.zeros((iz,jz))
    Z_spr2 = np.zeros((iz,jz))
    Z_spr3 = np.zeros((iz,jz))
    Z_scl1 = np.zeros((iz,jz))
    Z_scl2 = np.zeros((iz,jz))
    Z_scl3 = np.zeros((iz,jz))
    Z_rat1 = np.zeros((iz,jz))
    Z_rat2 = np.zeros((iz,jz))
    Z_tns1 = np.zeros((iz,jz))
    Z_tns2 = np.zeros((iz,jz))
    
    # Transferring raw data into rectangular matrices
    for i in range(iz):
        for j in range(jz):
            sol = str(i+1)
            exec('Z_men1[i,j] = s'+sol+'_men1[j]')
            exec('Z_men2[i,j] = s'+sol+'_men2[j]')
            exec('Z_men3[i,j] = s'+sol+'_men3[j]')
            exec('Z_spr1[i,j] = s'+sol+'_spr1[j]')
            exec('Z_spr2[i,j] = s'+sol+'_spr2[j]')
            exec('Z_spr3[i,j] = s'+sol+'_spr3[j]')
            exec('Z_scl1[i,j] = s'+sol+'_scl1[j]')
            exec('Z_scl2[i,j] = s'+sol+'_scl2[j]')
            exec('Z_scl3[i,j] = s'+sol+'_scl3[j]')
            exec('Z_rat1[i,j] = s'+sol+'_rat1[j]')
            exec('Z_rat2[i,j] = s'+sol+'_rat2[j]')
            exec('Z_tns1[i,j] = s'+sol+'_tns1[j]')
            exec('Z_tns2[i,j] = s'+sol+'_tns2[j]')
    
    # Creating a rectangular bivariate spline of the parameter data
    s_men1 = RectBivariateSpline(sol_d,tsr_d,Z_men1)
    s_men2 = RectBivariateSpline(sol_d,tsr_d,Z_men2)
    s_men3 = RectBivariateSpline(sol_d,tsr_d,Z_men3)
    s_spr1 = RectBivariateSpline(sol_d,tsr_d,Z_spr1)
    s_spr2 = RectBivariateSpline(sol_d,tsr_d,Z_spr2)
    s_spr3 = RectBivariateSpline(sol_d,tsr_d,Z_spr3)
    s_scl1 = RectBivariateSpline(sol_d,tsr_d,Z_scl1)
    s_scl2 = RectBivariateSpline(sol_d,tsr_d,Z_scl2)
    s_scl3 = RectBivariateSpline(sol_d,tsr_d,Z_scl3)
    s_rat1 = RectBivariateSpline(sol_d,tsr_d,Z_rat1)
    s_rat2 = RectBivariateSpline(sol_d,tsr_d,Z_rat2)
    s_tns1 = RectBivariateSpline(sol_d,tsr_d,Z_tns1)
    s_tns2 = RectBivariateSpline(sol_d,tsr_d,Z_tns2)
    
    # Selecting the specific parameters to use for TSR and solidity
    men1 = s_men1(solidity,tsr)
    men2 = s_men2(solidity,tsr)
    men3 = s_men3(solidity,tsr)
    spr1 = s_spr1(solidity,tsr)
    spr2 = s_spr2(solidity,tsr)
    spr3 = s_spr3(solidity,tsr)
    scl1 = s_scl1(solidity,tsr)
    scl2 = s_scl2(solidity,tsr)
    scl3 = s_scl3(solidity,tsr)
    rat1 = s_rat1(solidity,tsr)
    rat2 = s_rat2(solidity,tsr)
    tns1 = s_tns1(solidity,tsr)
    tns2 = s_tns2(solidity,tsr)
    
    # Creating arrays of the parameters
    men = np.array([men1[0,0],men2[0,0],men3[0,0]])
    spr = np.array([spr1[0,0],spr2[0,0],spr3[0,0]])
    scl = np.array([scl1[0,0],scl2[0,0],scl3[0,0]])
    rat = np.array([rat1[0,0],rat2[0,0]])
    tns = np.array([tns1[0,0],tns2[0,0]])

    return men,spr,scl,rat,tns

##Main
# print velocity(4.,0.25)

tsr = np.linspace(1.5,7,23)
solidity = 0.25

men1 = np.zeros_like(tsr)
men2 = np.zeros_like(tsr)
men3 = np.zeros_like(tsr)
spr1 = np.zeros_like(tsr)
spr2 = np.zeros_like(tsr)
spr3 = np.zeros_like(tsr)
scl1 = np.zeros_like(tsr)
scl2 = np.zeros_like(tsr)
scl3 = np.zeros_like(tsr)
rat1 = np.zeros_like(tsr)
rat2 = np.zeros_like(tsr)
tns1 = np.zeros_like(tsr)
tns2 = np.zeros_like(tsr)

for i in range(np.size(tsr)):
    men,spr,scl,rat,tns = velocity(tsr[i],solidity)
    men1[i] = men[0]
    men2[i] = men[1]
    men3[i] = men[2]
    spr1[i] = spr[0]
    spr2[i] = spr[1]
    spr3[i] = spr[2]
    scl1[i] = scl[0]
    scl2[i] = scl[1]
    scl3[i] = scl[2]
    rat1[i] = rat[0]
    rat2[i] = rat[1]
    tns1[i] = tns[0]
    tns2[i] = tns[1]
    

plt.figure()
plt.subplot(3,5,1)
plt.plot(tsr,men1)
plt.subplot(3,5,6)
plt.plot(tsr,men2)
plt.subplot(3,5,11)
plt.plot(tsr,men3)
plt.subplot(3,5,2)
plt.plot(tsr,spr1)
plt.subplot(3,5,7)
plt.plot(tsr,spr2)
plt.subplot(3,5,12)
plt.plot(tsr,spr3)
plt.subplot(3,5,3)
plt.plot(tsr,scl1)
plt.subplot(3,5,8)
plt.plot(tsr,scl2)
plt.subplot(3,5,13)
plt.plot(tsr,scl3)
plt.subplot(3,5,4)
plt.plot(tsr,rat1)
plt.subplot(3,5,9)
plt.plot(tsr,rat2)
plt.subplot(3,5,5)
plt.plot(tsr,tns1)
plt.subplot(3,5,10)
plt.plot(tsr,tns2)

plt.show()