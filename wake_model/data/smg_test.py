import numpy as np
from numpy import sqrt, pi, exp, fabs
import matplotlib.pyplot as plt

x = np.linspace(-10,10,1000)

color = np.array(['r','orange','y','g','b'])


m = np.ones(5)*-5.
sd = np.ones(5)*3.
r = np.ones(5)*1.
t = np.ones(5)*10.
sc = np.ones(5)*1.

mp = m[0]
sdp = sd[0]
rp = r[0]
tp = t[0]
scp = sc[0]

# m = np.array([-5.,-1.,0.,1.,5.])
sd = np.array([0.5,1.,2.,3.,4.])
# r = np.array([0.5,1.,6.,8.,10.])
# t = np.array([1.,5.,10.,20.,50.])
# sc = np.array([0.1,0.5,0.7,1.,10.])

def smg(x,m,sd,r,t,sc):
    return (-1./(sd*sqrt(2.*pi))*exp(-(x-m)**2/(2.*sd**2))*(1./(1.+exp(r*fabs(x)-t))))*sc + 1.


plt.figure()
for i in range(5):
    plt.plot(smg(x,m[i],sd[i],r[i],t[i],sc[i]),x,color[i])
    
p = (-1./(sdp*sqrt(2.*pi))*exp(-(x-mp)**2/(2.*sdp**2)))*scp+1.
plt.plot(p,x,'k--')

s = (1./(1.+exp(rp*fabs(x)-tp)))
plt.plot(s,x,'k:')

plt.show()
