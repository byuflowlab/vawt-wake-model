import csv
from os import path
import numpy as np
from numpy import pi,exp,fabs,sqrt
import matplotlib.pyplot as plt
from scipy.integrate import dblquad
from scipy.interpolate import RectBivariateSpline


from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'


def vorticity(show):
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
    # basepath = path.join(path.dirname(path.realpath(__file__)),'data')
    # fdata = basepath + path.sep + 'errordatabase.csv'
    fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/validation/error_cfd_vort.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw,0)
            errordat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d,float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            errordat = np.vstack([errordat,raw]) # adding entry to vorticity database array
        i += 1
    f.close()
    
    errordat = np.delete(errordat,(0),axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    errordat = errordat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('error2D'+sol+' = errordat[i*12]\nerror4D'+sol+' = errordat[i*12+1]\nerror6D'+sol+' = errordat[i*12+2]\nerror8D'+sol+' = errordat[i*12+3]\nerror10D'+sol+' = errordat[i*12+4]\nerror15D'+sol+' = errordat[i*12+5]\nstd2D'+sol+' = errordat[i*12+6]\nstd4D'+sol+' = errordat[i*12+7]\nstd6D'+sol+' = errordat[i*12+8]\nstd8D'+sol+' = errordat[i*12+9]\nstd10D'+sol+' = errordat[i*12+10]\nstd15D'+sol+' = errordat[i*12+11]\n')
    
    errorplot = np.zeros((np.size(sol_d),np.size(tsr_d)))
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        for j in range(np.size(tsr_d)):
            exec('errorplot[i,j] = fabs(np.average([error2D'+sol+'[j],error4D'+sol+'[j],error6D'+sol+'[j],error8D'+sol+'[j],error10D'+sol+'[j],error15D'+sol+'[j]]))')
    
    fs = 15
    fig = plt.figure(figsize=(8,5))
    # fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
    lb = 0. # lower bound on velocity to display
    ub = 1.0 # upper bound on velocity to display
    ran = 50 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
    CS = plt.contourf(tsr_d,sol_d,errorplot,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.parula) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel('% Error',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    CB.ax.set_aspect(20)
    plt.xlabel('TSR',fontsize=fs)
    plt.ylabel('Solidity',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    if show == True:
        plt.show()
    
    return 0


def velocity(show):
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
    # fdata = basepath + path.sep + 'velodatabase.csv'
    fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/validation/error_cfd_SMG2.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw,0)
            errordat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d,float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            errordat = np.vstack([errordat,raw]) # adding entry to vorticity database array
        i += 1
    f.close()
    
    errordat = np.delete(errordat,(0),axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    errordat = errordat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('error2D'+sol+' = errordat[i*12]\nerror4D'+sol+' = errordat[i*12+1]\nerror6D'+sol+' = errordat[i*12+2]\nerror8D'+sol+' = errordat[i*12+3]\nerror10D'+sol+' = errordat[i*12+4]\nerror15D'+sol+' = errordat[i*12+5]\nstd2D'+sol+' = errordat[i*12+6]\nstd4D'+sol+' = errordat[i*12+7]\nstd6D'+sol+' = errordat[i*12+8]\nstd8D'+sol+' = errordat[i*12+9]\nstd10D'+sol+' = errordat[i*12+10]\nstd15D'+sol+' = errordat[i*12+11]\n')
    
    errorplot = np.zeros((np.size(sol_d),np.size(tsr_d)))
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        for j in range(np.size(tsr_d)):
            exec('errorplot[i,j] = fabs(np.average([error2D'+sol+'[j],error4D'+sol+'[j],error6D'+sol+'[j],error8D'+sol+'[j],error10D'+sol+'[j],error15D'+sol+'[j]]))')
    
    fs = 15
    fig = plt.figure(figsize=(8,5))
    # fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
    lb = 0. # lower bound on velocity to display
    ub = 1.0 # upper bound on velocity to display
    ran = 50 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
    CS = plt.contourf(tsr_d,sol_d,errorplot,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.parula) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel('% Error',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    CB.ax.set_aspect(20)
    plt.xlabel('TSR',fontsize=fs)
    plt.ylabel('Solidity',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    if show == True:
        plt.show()
    
    return 0


def quad(show):
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
    # basepath = path.join(path.dirname(path.realpath(__file__)),'data')
    # fdata = basepath + path.sep + 'velodatabase_2.csv'
    fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/validation/error_cfd_quad.csv'
    f = open(fdata)
    csv_f = csv.reader(f)
    
    i = 0
    sol_d = np.array([])
    for row in csv_f:
        if i == 0:
            raw = row
            raw = np.delete(raw,0)
            errordat = raw
            tsr_d = raw # range of tip-speed ratios included
        if row[0] == 'solidity':
            sol_d = np.append(sol_d,float(row[1])) # range of solidities included
        elif row[0] != 'TSR' and row[0] != 'solidity':
            raw = row
            raw = np.delete(raw,0)
            errordat = np.vstack([errordat,raw]) # adding entry to vorticity database array
        i += 1
    f.close()
    
    errordat = np.delete(errordat,(0),axis=0) # eliminating first row used as a placeholder
    tsr_d = tsr_d.astype(np.float) # converting tip-speed ratio entries into floats
    errordat = errordat.astype(np.float) # converting vorticity database entries into floats
    
    # Creating arrays for each EMG parameter
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        
        exec('error2D'+sol+' = errordat[i*12]\nerror4D'+sol+' = errordat[i*12+1]\nerror6D'+sol+' = errordat[i*12+2]\nerror8D'+sol+' = errordat[i*12+3]\nerror10D'+sol+' = errordat[i*12+4]\nerror15D'+sol+' = errordat[i*12+5]\nstd2D'+sol+' = errordat[i*12+6]\nstd4D'+sol+' = errordat[i*12+7]\nstd6D'+sol+' = errordat[i*12+8]\nstd8D'+sol+' = errordat[i*12+9]\nstd10D'+sol+' = errordat[i*12+10]\nstd15D'+sol+' = errordat[i*12+11]\n')
    
    errorplot = np.zeros((np.size(sol_d),np.size(tsr_d)))
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        for j in range(np.size(tsr_d)):
            exec('errorplot[i,j] = fabs(np.average([error2D'+sol+'[j],error4D'+sol+'[j],error6D'+sol+'[j],error8D'+sol+'[j],error10D'+sol+'[j],error15D'+sol+'[j]]))')
    
    fs = 15
    fig = plt.figure(figsize=(8,5))
    # fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
    lb = 0. # lower bound on velocity to display
    ub = 1.0 # upper bound on velocity to display
    ran = 50 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
    CS = plt.contourf(tsr_d,sol_d,errorplot,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.parula) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel('% Error',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    CB.ax.set_aspect(20)
    plt.xlabel('TSR',fontsize=fs)
    plt.ylabel('Solidity',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    if show == True:
        plt.show()
    
    return 0

## Main
# vorticity(False)
velocity(False)
# quad(False)

plt.show()