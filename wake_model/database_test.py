import numpy as np
import matplotlib.pyplot as plt
from database_call import vorticity,velocity,quad
import csv
from os import path
from scipy.interpolate import RectBivariateSpline



##Main
# print velocity(4.,0.25)
type = 'vort'
type = 'velo'
# type = 'quad'


tsr = np.linspace(1.5,7,230)
# solidity = np.array([.15])
solidity = np.array([0.15,0.25,0.5,0.75,1.0])
col = np.array(['b','g','y','r','m'])


if type == 'vort':
    loc1 = np.zeros_like(tsr)
    loc2 = np.zeros_like(tsr)
    loc3 = np.zeros_like(tsr)
    spr1 = np.zeros_like(tsr)
    spr2 = np.zeros_like(tsr)
    skw1 = np.zeros_like(tsr)
    skw2 = np.zeros_like(tsr)
    scl1 = np.zeros_like(tsr)
    scl2 = np.zeros_like(tsr)
    scl3 = np.zeros_like(tsr)
    
    for i in range(np.size(tsr)):
        loc,spr,skw,scl = vorticity(tsr[i],solidity)
        loc1[i] = loc[0]
        loc2[i] = loc[1]
        loc3[i] = loc[2]
        spr1[i] = spr[0]
        spr2[i] = spr[1]
        skw1[i] = skw[0]
        skw2[i] = skw[1]
        scl1[i] = scl[0]
        scl2[i] = scl[1]
        scl3[i] = scl[2]
        
    
    plt.figure()
    plt.subplot(3,4,1)
    plt.plot(tsr,loc1)
    plt.subplot(3,4,5)
    plt.plot(tsr,loc2)
    plt.subplot(3,4,9)
    plt.plot(tsr,loc3)
    plt.subplot(3,4,2)
    plt.plot(tsr,spr1)
    plt.subplot(3,4,6)
    plt.plot(tsr,spr2)
    plt.subplot(3,4,3)
    plt.plot(tsr,skw1)
    plt.subplot(3,4,7)
    plt.plot(tsr,skw2)
    plt.subplot(3,4,4)
    plt.plot(tsr,scl1)
    plt.subplot(3,4,8)
    plt.plot(tsr,scl2)
    plt.subplot(3,4,12)
    plt.plot(tsr,scl3)


elif type == 'velo':
    # for k in range(np.size(solidity)):
    #     men1 = np.zeros_like(tsr)
    #     men2 = np.zeros_like(tsr)
    #     men3 = np.zeros_like(tsr)
    #     spr1 = np.zeros_like(tsr)
    #     spr2 = np.zeros_like(tsr)
    #     spr3 = np.zeros_like(tsr)
    #     spr4 = np.zeros_like(tsr)
    #     scl1 = np.zeros_like(tsr)
    #     scl2 = np.zeros_like(tsr)
    #     scl3 = np.zeros_like(tsr)
    #     rat1 = np.zeros_like(tsr)
    #     rat2 = np.zeros_like(tsr)
    #     tns1 = np.zeros_like(tsr)
    #     tns2 = np.zeros_like(tsr)
    #
    #     for i in range(np.size(tsr)):
    #         men,spr,scl,rat,tns = velocity(tsr[i],solidity[k])
    #         men1[i] = men[0]
    #         men2[i] = men[1]
    #         men3[i] = men[2]
    #         spr1[i] = spr[0]
    #         spr2[i] = spr[1]
    #         spr3[i] = spr[2]
    #         spr4[i] = spr[3]
    #         scl1[i] = scl[0]
    #         scl2[i] = scl[1]
    #         scl3[i] = scl[2]
    #         rat1[i] = rat[0]
    #         rat2[i] = rat[1]
    #         tns1[i] = tns[0]
    #         tns2[i] = tns[1]
    #
    #     plt.figure()
    #     plt.subplot(4,5,1)
    #     plt.plot(tsr,men1,col[k])
    #     plt.subplot(4,5,6)
    #     plt.plot(tsr,men2,col[k])
    #     plt.subplot(4,5,11)
    #     plt.plot(tsr,men3,col[k])
    #     plt.subplot(4,5,2)
    #     plt.plot(tsr,spr1,col[k])
    #     plt.subplot(4,5,7)
    #     plt.plot(tsr,spr2,col[k])
    #     plt.subplot(4,5,12)
    #     plt.plot(tsr,spr3,col[k])
    #     plt.subplot(4,5,17)
    #     plt.plot(tsr,spr4,col[k])
    #     plt.subplot(4,5,3)
    #     plt.plot(tsr,scl1,col[k])
    #     plt.subplot(4,5,8)
    #     plt.plot(tsr,scl2,col[k])
    #     plt.subplot(4,5,13)
    #     plt.plot(tsr,scl3,col[k])
    #     plt.subplot(4,5,4)
    #     plt.plot(tsr,rat1,col[k])
    #     plt.subplot(4,5,9)
    #     plt.plot(tsr,rat2,col[k])
    #     plt.subplot(4,5,5)
    #     plt.plot(tsr,tns1,col[k])
    #     plt.subplot(4,5,10)
    #     plt.plot(tsr,tns2,col[k])



    for k in range(15):
        men11 = np.zeros_like(solidity)
        men22 = np.zeros_like(solidity)
        men33 = np.zeros_like(solidity)
        spr11 = np.zeros_like(solidity)
        spr22 = np.zeros_like(solidity)
        spr33 = np.zeros_like(solidity)
        spr44 = np.zeros_like(solidity)
        scl11 = np.zeros_like(solidity)
        scl22 = np.zeros_like(solidity)
        scl33 = np.zeros_like(solidity)
        rat11 = np.zeros_like(solidity)
        rat22 = np.zeros_like(solidity)
        tns11 = np.zeros_like(solidity)
        tns22 = np.zeros_like(solidity)

        for i in range(np.size(solidity)):
            men,spr,scl,rat,tns = velocity(tsr[k],solidity[i])
            men11[i] = men[0]
            men22[i] = men[1]
            men33[i] = men[2]
            spr11[i] = spr[0]
            spr22[i] = spr[1]
            spr33[i] = spr[2]
            spr44[i] = spr[3]
            scl11[i] = scl[0]
            scl22[i] = scl[1]
            scl33[i] = scl[2]
            rat11[i] = rat[0]
            rat22[i] = rat[1]
            tns11[i] = tns[0]
            tns22[i] = tns[1]

        plt.figure()
        plt.subplot(4,5,1)
        plt.plot(solidity,men11)
        plt.subplot(4,5,6)
        plt.plot(solidity,men22)
        plt.subplot(4,5,11)
        plt.plot(solidity,men33)
        plt.subplot(4,5,2)
        plt.plot(solidity,spr11)
        plt.subplot(4,5,7)
        plt.plot(solidity,spr22)
        plt.subplot(4,5,12)
        plt.plot(solidity,spr33)
        plt.subplot(4,5,17)
        plt.plot(solidity,spr44)
        plt.subplot(4,5,3)
        plt.plot(solidity,scl11)
        plt.subplot(4,5,8)
        plt.plot(solidity,scl22)
        plt.subplot(4,5,13)
        plt.plot(solidity,scl33)
        plt.subplot(4,5,4)
        plt.plot(solidity,rat11)
        plt.subplot(4,5,9)
        plt.plot(solidity,rat22)
        plt.subplot(4,5,5)
        plt.plot(solidity,tns11)
        plt.subplot(4,5,10)
        plt.plot(solidity,tns22)

    
elif type == 'quad':
    scl1 = np.zeros_like(tsr)
    scl2 = np.zeros_like(tsr)
    scl3 = np.zeros_like(tsr)
    trn1 = np.zeros_like(tsr)
    trn2 = np.zeros_like(tsr)
    trn3 = np.zeros_like(tsr)
    
    for i in range(np.size(tsr)):
        scl,trn = quad(tsr[i],solidity)
        scl1[i] = scl[0]
        scl2[i] = scl[1]
        scl3[i] = scl[2]
        trn1[i] = trn[0]
        trn2[i] = trn[1]
        trn3[i] = trn[2]
        
    plt.figure()
    plt.subplot(3,2,1)
    plt.plot(tsr,scl1)
    plt.subplot(3,2,3)
    plt.plot(tsr,scl2)
    plt.subplot(3,2,5)
    plt.plot(tsr,scl3)
    plt.subplot(3,2,2)
    plt.plot(tsr,trn1)
    plt.subplot(3,2,4)
    plt.plot(tsr,trn2)
    plt.subplot(3,2,6)
    plt.plot(tsr,trn3)
    
plt.show()