from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp,sqrt,pi,fabs

n = 100
x = np.linspace(-2,2,n)
y = np.linspace(0,150./6.,n)

men = 0.0706029475998
sdv = 0.361485432916
rat = 29.2307972266
spr = 20.8552963734
scl1 = 131.542734985
scl2 = 0.000245565522327
scl3 = 15.2204942377

men = -0.0275297605532
sdv = 0.687598282579
rat = 17.0583226917
spr = 12.153515003
scl1 = 0.378233548012
scl2 = 0.135586038933
scl3 = 19.9479281363

# men = 0.242939177452
# sdv = 0.511597681187
# rat = 32.0060071153
# spr = 16.5348856064
# scl1 = 3.09220361845
# scl2 = 0.00199400483863
# scl3 = 20.2838888522

# men = 0.0740949336081
# sdv = 0.369311970256
# rat = 29.6854166068
# spr = 19.7709398653
# scl1 = 266.998161146
# scl2 = 0.00011972640774
# scl3 = 15.5600231975

men = 0.0766567990473
sdv = 0.370036128166
rat = 10.8615667761
spr = 8.39081017485
scl1 = 1.0
scl2 = 0.0393076036
scl3 = 10.8396908876

men = -0.0275296976446
sdv1 = 0.261654851058
sdv2 = -9.63893914504e-08
sdv3 = 19.9449868059
sdv4 = 0.687598665934
rat = 17.0583314564
spr = 12.1535216043
scl1 = 0.378233474572
scl2 = 0.135586064238
scl3 = 19.9479225419

men = -0.0264524653241
sdv1 = 0.0645418855488
sdv2 = 0.205589658547
sdv3 = 29.6181325216
sdv4 = 0.00746603948577
rat1 = -1.32348928602
rat2 = 24.9137627006
spr1 = -1.1720762744
spr2 = 19.0426787529
scl1 = 0.0655078542251
scl2 = 0.217055323455
scl3 = 50.0

men = -0.0234253600051
sdv1 = 6.54745750152
sdv2 = 0.00609915951455
sdv3 = 17.2706153977
sdv4 = 0.0
rat = 30.5434995612
tns = 15.6651700448
spr1 = 0.21939104869
spr2 = 0.216838482049
spr3 = 3.79650578382
spr4 = 1.06611237486
scl1 = 0.547915533139
scl2 = 0.113773089751
scl3 = 15.6659553933

# sdv_v = sdv3*sdv2*sdv1*np.exp(sdv2*y)*np.exp(sdv1)*np.exp(-sdv1*np.exp(sdv2*y)) + sdv4
# plt.figure()
# plt.plot(y,sdv_v)
#
# rat_v = rat1*y + rat2
# plt.figure()
# plt.plot(y,rat_v)
#
# spr_v = spr1*y + spr2
# plt.figure()
# plt.plot(y,spr_v)
#
# scl_v = scl3*scl2*scl1*np.exp(scl2*y)*np.exp(scl1)*np.exp(-scl1*np.exp(scl2*y))
# plt.figure()
# plt.plot(y,scl_v)
#
#
# plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
G = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        sdv_v = sdv3*sdv2*sdv1*exp(sdv2*Y[i,j])*exp(sdv1)*exp(-sdv1*exp(sdv2*Y[i,j]))+sdv4

        spr_v = spr3*spr2*spr1*exp(spr2*Y[i,j])*exp(spr1)*exp(-spr1*exp(spr2*Y[i,j]))+spr4
        # spr_v = 1.

        f1 = -1./(sdv_v*sqrt(2.*pi))*exp(-((X[i,j]/spr_v)-men)**2/(2.*sdv_v**2))*(1./(1.+exp(rat*fabs((X[i,j]/spr_v))-tns)))
        f2 = scl3*scl2*scl1*exp(scl2*Y[i,j])*exp(scl1)*exp(-scl1*exp(scl2*Y[i,j]))
        G[i,j] = f1*f2 + 1.


surf = ax.plot_surface(X, Y, G, rstride=1, cstride=1, cmap=cm.parula, linewidth=0, antialiased=False)

plt.show()