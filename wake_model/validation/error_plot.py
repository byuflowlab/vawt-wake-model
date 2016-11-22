import csv
from os import path
import numpy as np
from numpy import fabs
import matplotlib.pyplot as plt
import scipy.ndimage


from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'


def vorticity(show):
    # Reading in csv file (vorticity database)
    basepath = path.join(path.dirname(path.realpath(__file__)))
    # fdata = basepath + path.sep + 'error_cfd_vort_EMG.csv'
    # fdata = basepath + path.sep + 'error_cfd_vort_EMG_deficit_abs.csv'
    # fdata = basepath + path.sep + 'error_cfd_vort_EMG_deficit_rel.csv'
    fdata = basepath + path.sep + 'error_cfd_vort_EMG_deficit_rms.csv'
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
        
        exec('error2D'+sol+' = errordat[i*14]\nerror4D'+sol+' = errordat[i*14+1]\nerror6D'+sol+' = errordat[i*14+2]\nerror8D'+sol+' = errordat[i*14+3]\nerror10D'+sol+' = errordat[i*14+4]\nerror15D'+sol+' = errordat[i*14+5]\nerror20D'+sol+' = errordat[i*14+6]\nstd2D'+sol+' = errordat[i*14+7]\nstd4D'+sol+' = errordat[i*14+8]\nstd6D'+sol+' = errordat[i*14+9]\nstd8D'+sol+' = errordat[i*14+10]\nstd10D'+sol+' = errordat[i*14+11]\nstd15D'+sol+' = errordat[i*14+12]\nstd20D'+sol+' = errordat[i*14+13]')
    
    errorplot = np.zeros((np.size(sol_d),np.size(tsr_d)))
    for i in range(np.size(sol_d)):
        sol = str(i+1)
        for j in range(np.size(tsr_d)):
            exec('errorplot[i,j] = fabs(np.average([error2D'+sol+'[j],error4D'+sol+'[j],error6D'+sol+'[j],error8D'+sol+'[j],error10D'+sol+'[j],error15D'+sol+'[j],error20D'+sol+'[j]]))')
    
    fs = 15
    fig = plt.figure(figsize=(8,5))
    fig.subplots_adjust(bottom=.12)#,left=.05,right=1.0)
    lb = 0. # lower bound on velocity to display
    ub = 0.15 # upper bound on velocity to display
    ran = 50 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar

    zoomval = 3.
    data1 = scipy.ndimage.interpolation.zoom(tsr_d, zoom=zoomval)
    data2 = scipy.ndimage.interpolation.zoom(sol_d, zoom=zoomval)
    data3 = scipy.ndimage.interpolation.zoom(errorplot, zoom=zoomval)

    CS = plt.contourf(data1,data2,data3,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.parula) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel(r'RMS Error Normalized by $U_\infty$',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    CB.ax.set_aspect(20)
    plt.xlabel('TSR',fontsize=fs)
    plt.ylabel('Solidity',fontsize=fs)
    # plt.title('Average RMS Error')
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    
    if show == True:
        plt.savefig('/Users/ning1/Documents/FLOW Lab/tingey-2016-vawt-wake-model/journal_version/images/error_plot.pdf')
        # plt.savefig('/Users/ning1/Documents/FLOW Lab/error_plot_rms.png')
        plt.show()
    
    return 0

## Main
vorticity(True)
