import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from database_call import vorticity,velocity,velocity2
import csv
from os import path
from scipy.interpolate import RectBivariateSpline



##Main
# print velocity(4.,0.25)
type = 'vort'
type = 'velo'
# type = 'velo2'
# type = 'quad'


tsr = np.linspace(1.5,7,23)
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
    X,Y = np.meshgrid(tsr,solidity)
    n = 5
    m = 23
    men1 = np.zeros((n,m))
    sdv1 = np.zeros((n,m))
    sdv2 = np.zeros((n,m))
    sdv3 = np.zeros((n,m))
    sdv4 = np.zeros((n,m))
    rat1 = np.zeros((n,m))
    wdt1 = np.zeros((n,m))
    spr1 = np.zeros((n,m))
    spr2 = np.zeros((n,m))
    spr3 = np.zeros((n,m))
    spr4 = np.zeros((n,m))
    scl1 = np.zeros((n,m))
    scl2 = np.zeros((n,m))
    scl3 = np.zeros((n,m))


    for i in range(np.size(solidity)):
        for j in range(np.size(tsr)):
            men,sdv,rat,wdt,spr,scl,_,_ = velocity(X[i,j],Y[i,j])
            men1[i,j] = men[0]
            sdv1[i,j] = sdv[0]
            sdv2[i,j] = sdv[1]
            sdv3[i,j] = sdv[2]
            sdv4[i,j] = sdv[3]
            rat1[i,j] = rat[0]
            wdt1[i,j] = wdt[0]
            spr1[i,j] = spr[0]
            spr2[i,j] = spr[1]
            spr3[i,j] = spr[2]
            spr4[i,j] = spr[3]
            scl1[i,j] = scl[0]
            scl2[i,j] = scl[1]
            scl3[i,j] = scl[2]

        fig = plt.figure(1)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,men1,'b.')
        ax.set_xlabel('men1')
        fig = plt.figure(2)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,sdv1,'b.')
        ax.set_xlabel('sdv1')
        fig = plt.figure(3)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,sdv2,'b.')
        ax.set_xlabel('sdv2')
        fig = plt.figure(4)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,sdv3,'b.')
        ax.set_xlabel('sdv3')
        fig = plt.figure(5)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,sdv4,'b.')
        ax.set_xlabel('sdv4')
        fig = plt.figure(6)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,rat1,'b.')
        ax.set_xlabel('rat1')
        fig = plt.figure(7)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,wdt1,'b.')
        ax.set_xlabel('wdt1')
        fig = plt.figure(8)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,spr1,'b.')
        ax.set_xlabel('spr1')
        fig = plt.figure(9)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,spr2,'b.')
        ax.set_xlabel('spr2')
        fig = plt.figure(10)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,spr3,'b.')
        ax.set_xlabel('spr3')
        fig = plt.figure(11)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,spr4,'b.')
        ax.set_xlabel('spr4')
        fig = plt.figure(12)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,scl1,'b.')
        ax.set_xlabel('scl1')
        fig = plt.figure(13)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,scl2,'b.')
        ax.set_xlabel('scl2')
        fig = plt.figure(14)
        ax = fig.gca(projection='3d')
        pnt = ax.scatter(X,Y,scl3,'b.')
        ax.set_xlabel('scl3')



    # for k in range(15):
    #     men11 = np.zeros_like(solidity)
    #     sdv11 = np.zeros_like(solidity)
    #     sdv22 = np.zeros_like(solidity)
    #     sdv33 = np.zeros_like(solidity)
    #     sdv44 = np.zeros_like(solidity)
    #     rat11 = np.zeros_like(solidity)
    #     wdt11 = np.zeros_like(solidity)
    #     spr11 = np.zeros_like(solidity)
    #     spr22 = np.zeros_like(solidity)
    #     spr33 = np.zeros_like(solidity)
    #     spr44 = np.zeros_like(solidity)
    #     scl11 = np.zeros_like(solidity)
    #     scl22 = np.zeros_like(solidity)
    #     scl33 = np.zeros_like(solidity)
    #
    #     for i in range(np.size(solidity)):
    #         men,sdv,rat,wdt,spr,scl,_,_ = velocity(tsr[k],solidity[i])
    #         men11[i] = men[0]
    #         sdv11[i] = sdv[0]
    #         sdv22[i] = sdv[1]
    #         sdv33[i] = sdv[2]
    #         sdv44[i] = sdv[3]
    #         rat11[i] = rat[0]
    #         wdt11[i] = wdt[0]
    #         spr11[i] = spr[0]
    #         spr22[i] = spr[1]
    #         spr33[i] = spr[2]
    #         spr44[i] = spr[3]
    #         scl11[i] = scl[0]
    #         scl22[i] = scl[1]
    #         scl33[i] = scl[2]
    #
    #     plt.figure()
    #     plt.subplot(4,6,1)
    #     plt.plot(solidity,men11,col[0])
    #     plt.subplot(4,6,2)
    #     plt.plot(solidity,sdv11,col[0])
    #     plt.subplot(4,6,8)
    #     plt.plot(solidity,sdv22,col[0])
    #     plt.subplot(4,6,14)
    #     plt.plot(solidity,sdv33,col[0])
    #     plt.subplot(4,6,20)
    #     plt.plot(solidity,sdv44,col[0])
    #     plt.subplot(4,6,3)
    #     plt.plot(solidity,rat11,col[0])
    #     plt.subplot(4,6,4)
    #     plt.plot(solidity,wdt11,col[0])
    #     plt.subplot(4,6,5)
    #     plt.plot(solidity,spr11,col[0])
    #     plt.subplot(4,6,11)
    #     plt.plot(solidity,spr22,col[0])
    #     plt.subplot(4,6,17)
    #     plt.plot(solidity,spr33,col[0])
    #     plt.subplot(4,6,23)
    #     plt.plot(solidity,spr44,col[0])
    #     plt.subplot(4,6,6)
    #     plt.plot(solidity,scl11,col[0])
    #     plt.subplot(4,6,12)
    #     plt.plot(solidity,scl22,col[0])
    #     plt.subplot(4,6,18)
    #     plt.plot(solidity,scl33,col[0])


elif type == 'velo2':
    for k in range(np.size(solidity)):
        men1 = np.zeros_like(tsr)
        sdv1 = np.zeros_like(tsr)
        sdv2 = np.zeros_like(tsr)
        sdv3 = np.zeros_like(tsr)
        sdv4 = np.zeros_like(tsr)
        rat1 = np.zeros_like(tsr)
        wdt1 = np.zeros_like(tsr)
        spr1 = np.zeros_like(tsr)
        spr2 = np.zeros_like(tsr)
        spr3 = np.zeros_like(tsr)
        spr4 = np.zeros_like(tsr)
        scl1 = np.zeros_like(tsr)
        scl2 = np.zeros_like(tsr)
        scl3 = np.zeros_like(tsr)

        for i in range(np.size(tsr)):
            men,sdv,rat,wdt,spr,scl = velocity2(tsr[i],solidity[k])
            men1[i] = men[0]
            sdv1[i] = sdv[0]
            sdv2[i] = sdv[1]
            sdv3[i] = sdv[2]
            sdv4[i] = sdv[3]
            rat1[i] = rat[0]
            wdt1[i] = wdt[0]
            spr1[i] = spr[0]
            spr2[i] = spr[1]
            spr3[i] = spr[2]
            spr4[i] = spr[3]
            scl1[i] = scl[0]
            scl2[i] = scl[1]
            scl3[i] = scl[2]

        plt.figure()
        plt.subplot(4,6,1)
        plt.plot(tsr,men1,col[k])
        plt.subplot(4,6,2)
        plt.plot(tsr,sdv1,col[k])
        plt.subplot(4,6,8)
        plt.plot(tsr,sdv2,col[k])
        plt.subplot(4,6,14)
        plt.plot(tsr,sdv3,col[k])
        plt.subplot(4,6,20)
        plt.plot(tsr,sdv4,col[k])
        plt.subplot(4,6,3)
        plt.plot(tsr,rat1,col[k])
        plt.subplot(4,6,4)
        plt.plot(tsr,wdt1,col[k])
        plt.subplot(4,6,5)
        plt.plot(tsr,spr1,col[k])
        plt.subplot(4,6,11)
        plt.plot(tsr,spr2,col[k])
        plt.subplot(4,6,17)
        plt.plot(tsr,spr3,col[k])
        plt.subplot(4,6,23)
        plt.plot(tsr,spr4,col[k])
        plt.subplot(4,6,6)
        plt.plot(tsr,scl1,col[k])
        plt.subplot(4,6,12)
        plt.plot(tsr,scl2,col[k])
        plt.subplot(4,6,18)
        plt.plot(tsr,scl3,col[k])



    # for k in range(15):
    #     men11 = np.zeros_like(solidity)
    #     sdv11 = np.zeros_like(solidity)
    #     sdv22 = np.zeros_like(solidity)
    #     sdv33 = np.zeros_like(solidity)
    #     sdv44 = np.zeros_like(solidity)
    #     rat11 = np.zeros_like(solidity)
    #     wdt11 = np.zeros_like(solidity)
    #     spr11 = np.zeros_like(solidity)
    #     spr22 = np.zeros_like(solidity)
    #     spr33 = np.zeros_like(solidity)
    #     spr44 = np.zeros_like(solidity)
    #     scl11 = np.zeros_like(solidity)
    #     scl22 = np.zeros_like(solidity)
    #     scl33 = np.zeros_like(solidity)
    #
    #     for i in range(np.size(solidity)):
    #         men,sdv,rat,wdt,spr,scl = velocity2(tsr[k],solidity[i])
    #         men11[i] = men[0]
    #         sdv11[i] = sdv[0]
    #         sdv22[i] = sdv[1]
    #         sdv33[i] = sdv[2]
    #         sdv44[i] = sdv[3]
    #         rat11[i] = rat[0]
    #         wdt11[i] = wdt[0]
    #         spr11[i] = spr[0]
    #         spr22[i] = spr[1]
    #         spr33[i] = spr[2]
    #         spr44[i] = spr[3]
    #         scl11[i] = scl[0]
    #         scl22[i] = scl[1]
    #         scl33[i] = scl[2]
    #
    #     plt.figure()
    #     plt.subplot(4,6,1)
    #     plt.plot(solidity,men11,col[0])
    #     plt.subplot(4,6,2)
    #     plt.plot(solidity,sdv11,col[0])
    #     plt.subplot(4,6,8)
    #     plt.plot(solidity,sdv22,col[0])
    #     plt.subplot(4,6,14)
    #     plt.plot(solidity,sdv33,col[0])
    #     plt.subplot(4,6,20)
    #     plt.plot(solidity,sdv44,col[0])
    #     plt.subplot(4,6,3)
    #     plt.plot(solidity,rat11,col[0])
    #     plt.subplot(4,6,4)
    #     plt.plot(solidity,wdt11,col[0])
    #     plt.subplot(4,6,5)
    #     plt.plot(solidity,spr11,col[0])
    #     plt.subplot(4,6,11)
    #     plt.plot(solidity,spr22,col[0])
    #     plt.subplot(4,6,17)
    #     plt.plot(solidity,spr33,col[0])
    #     plt.subplot(4,6,23)
    #     plt.plot(solidity,spr44,col[0])
    #     plt.subplot(4,6,6)
    #     plt.plot(solidity,scl11,col[0])
    #     plt.subplot(4,6,12)
    #     plt.plot(solidity,scl22,col[0])
    #     plt.subplot(4,6,18)
    #     plt.plot(solidity,scl33,col[0])

    
plt.show()