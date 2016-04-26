from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
import numpy as np
from numpy import exp,sqrt,pi,fabs

n = 100
x = np.linspace(-2,2,n)
y = np.linspace(0,100.,n)

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

# gom = scl3*scl2*scl1*np.exp(scl2*y)*np.exp(scl1)*np.exp(-scl1*np.exp(scl2*y))
#
# plt.plot(y,gom)
# plt.show()

fig = plt.figure()
ax = fig.gca(projection='3d')
X, Y = np.meshgrid(x, y)
G = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        g1 = -1./(sdv*np.sqrt(2.*np.pi))*np.exp(-(X[i,j]-men)**2/(2.*sdv**2))*(1./(1.+exp(rat*fabs(X[i,j])-spr)))
        g2 = scl3*scl2*scl1*np.exp(scl2*Y[i,j])*np.exp(scl1)*np.exp(-scl1*np.exp(scl2*Y[i,j]))
        G[i,j] = g1*g2 + 1.


surf = ax.plot_surface(X, Y, G, rstride=1, cstride=1, cmap=cm.parula, linewidth=0, antialiased=False)

plt.show()