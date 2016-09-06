import numpy as np
import matplotlib.pyplot as plt
import csv
from VAWT_Wake_Model import velocity_field
from database_call import vorticity,velocity,velocity2
from numpy import fabs

from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

r = 0.5
v = 1.

flip = 1.


x15 = flip*np.linspace(-1.7,1.7,21)
x15r = flip*np.linspace(-1.7,1.7,100)

y15c = np.array([1.06,1.06,1.07,1.09,1.1,1.11,1.09,0.91,0.65,0.56,0.66,0.62,0.8,1.05,1.12,1.115,1.115,1.11,1.105,1.1,1.095])
y15o = np.array([0.98,0.98,0.985,0.99,0.995,1.005,0.96,0.73,0.5,0.44,0.54,0.51,0.66,0.9,1.01,1.0,0.99,0.985,0.98,0.98,0.97])

rad = 0.515
dia = 2*rad
velf = 16.0
sol = 0.5
tsr = 1.6
rom15 = np.zeros_like(x15r)
rom15t = np.zeros_like(x15)
error_test = np.zeros_like(x15)

for k in range(1):
    # Choose whether CFD vorticity or velocity data will be used as the basis
    if k == 1:
        cfd_data = 'velo'
        sol = 0.5
        tsr = 1.6
    elif k == 0:
        cfd_data = 'velo2'
        sol = 0.5
        tsr = 1.6
    
    if cfd_data == 'vort':
        loc,spr,skw,scl = vorticity(tsr,sol)
        param = np.array([loc,spr,skw,scl])
        
    elif cfd_data == 'velo':
        men1,sdv1,rat1,wdt1,spr1,scl1,tsrn1,_ = velocity(tsr-0.1249,sol)
        men2,sdv2,rat2,wdt2,spr2,scl2,tsrn2,_ = velocity(tsr+0.1249,sol)
        if sol >= 0.35:
            men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,sol-0.1249)
            men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,sol+0.1249)
        elif sol >=0.25:
            men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,sol-0.049)
            men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,sol+0.1249)
        else:
            men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,sol-0.049)
            men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,sol+0.049)
        if tsrn1 == tsrn2:
            p = 0.
        else:
            p = (tsr-tsrn1)/(tsrn2-tsrn1)
        if soln1 == soln2:
            q = 0.
        else:
            q = (sol-soln1)/(soln2-soln1)

        # import time
        # print tsrn1,tsrn2
        # print soln1,soln2
        # print p,q
        # time.sleep(10)
        param = np.array([men1,sdv1,rat1,wdt1,spr1,scl1,men2,sdv2,rat2,wdt2,spr2,scl2,men3,sdv3,rat3,wdt3,spr3,scl3,men4,sdv4,rat4,wdt4,spr4,scl4,p,q])

    elif cfd_data == 'velo2':
        spr1,pow1,pow2,spr2,skw,scl1,scl2,scl3 = velocity2(tsr,sol)
        param = np.array([spr1,pow1,pow2,spr2,skw,scl1,scl2,scl3])
        # print param
        # import time
        #
        # time.sleep(10)

    for i in range(np.size(x15)):
        rom15t[i] = velocity_field(0.,0.,1.5*dia,x15[i]*dia,velf,dia,tsr,sol,cfd_data,param)
        error_test[i] = (rom15t[i]-y15o[i])/y15o[i]
    error = np.average(fabs(error_test))
    errorstd = np.std(fabs(error_test))
    for i in range(np.size(rom15)):
        rom15[i] = velocity_field(0.,0.,1.5*dia,x15r[i]*dia,velf,dia,tsr,sol,cfd_data,param)
        print i
    
    fs = 15
    
    plt.figure(1)
    # plt.subplot(2,3,1)
    # plt.plot(x15,y15c,'r.',label='Experimental (closed)')
    # plt.plot(x15,y15o,'b.',label='Experimental (open)')
    if k == 0:
        plt.plot(x15,y15o,'b.',label='Experimental')
        plt.plot(x15r,rom15,'m-',label='Vorticity')
        plt.xlim(-1.75,1.75)
        # plt.ylim(0.3,1.4)
        plt.xlabel('$y/D$',fontsize=fs)
        plt.ylabel(r'$u/U_\infty$',fontsize=fs)
        plt.xticks(fontsize=fs)
        plt.yticks(fontsize=fs)
    elif k == 1:
        plt.plot(x15r,rom15,'g-',label='Velocity')
        plt.legend(loc=1,fontsize=fs)
    # plt.text(-0.5,0.9,'X/D = 1.5')
    # print '1.5 modc',(min(rom15)-min(y15c))/min(y15c)
    print '1.5 modo',(min(rom15)-min(y15o))/min(y15o),error,errorstd
plt.show()

# Vort
# 1.5 modo -0.151387906156 0.0974744457882 0.0791398227522
# Velo- val
# 1.5 modo 0.156755962966 0.0837934819892 0.142086247374
# Velo- flip
# 1.5 modo 0.156755962966 0.0531856178905 0.0699217546663


