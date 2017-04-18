import numpy as np
import matplotlib.pyplot as plt
from VAWT_Wake_Model import velocity_field

from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

velf = 0.15 # free stream wind speed (m/s)
dia = .12 # turbine diameter (m)
tsr = 2.14 # tip speed ratio
# tsr = 3.85 # tip speed ratio
rot = tsr*velf/(dia/2.)
B = 2. # number of blades
chord = 0.02 # chord lenth (m)
solidity = (chord*B)/(dia/2.)
r = dia*2.

xt = 0.
yt = 0.

veltype = 'x'

# Data from Brochier, et al. (provided by Dr. Shamsoddin and Dr. Porte-Agel)
yr_a = np.array([-1.2991,-1.1776,-0.8318,-0.6262,-0.4393,-0.3178,-0.1589,0.,0.0841,0.1776,0.3551,0.5140,0.6636,0.8598,1.0280,1.1589,1.3458])
uU_a = np.array([1.6837,1.6721,1.2375,0.9849,0.7508,0.6245,0.4757,0.3937,0.3785,0.3893,0.3664,0.3399,0.4208,0.5905,0.9565,1.3709,1.6072])
yr_b = np.array([-1.2500,-1.1667,-0.9815,-0.8056,-0.6389,-0.5000,-0.2778,-0.1574,-0.0741,0.0185,0.2130,0.3426,0.5185,0.6759,0.8519,1.0093,1.1852,1.3333])
uU_b = np.array([1.7222,1.6222,1.4407,1.2037,1.0074,0.8556,0.6481,0.5370,0.4852,0.4630,0.4185,0.3963,0.4259,0.4704,0.6815,1.0185,1.3667,1.5593])
yr_c = np.array([-1.2123,-1.1296,-0.9657,-0.7947,-0.6588,-0.4522,-0.2979,-0.1623,0.0134,0.1534,0.3228,0.4864,0.6797,0.8496,1.0222,1.1853,1.3615])
uU_c = np.array([1.5897,1.5414,1.4103,1.2448,1.0517,0.9276,0.7690,0.5690,0.5172,0.4207,0.4172,0.4793,0.6414,0.8517,1.1241,1.3724,1.5310])
yr_d = np.array([-1.3035,-1.2146,-1.0159,-0.8538,-0.6720,-0.5194,-0.3390,-0.1682,-0.0190,0.1472,0.3087,0.4695,0.6482,0.7693,0.9650,1.1016,1.2899])
uU_d = np.array([1.5000,1.4680,1.3760,1.2360,1.0920,0.9480,0.8320,0.7120,0.6400,0.6160,0.6880,0.7760,0.8960,1.0000,1.1680,1.3560,1.4760])
yr_e = np.array([-1.2427,-1.1262,-0.9515,-0.7961,-0.6311,-0.4660,-0.3204,-0.1359, 0.,0.1553,0.3398,0.5437,0.7087,0.8738,1.0388,1.1845,1.3786,1.5340])
uU_e = np.array([1.3739,1.3742,1.3014,1.1940,1.0867,1.0177,0.9488,0.8722,0.8417,0.8343,0.8577,0.9197,1.0008,1.0895,1.2091,1.3209,1.3597,1.3062])

xr_a = 2.5
xr_b = 3.33
xr_c = 5.
xr_d = 8.33
xr_e = 11.67

N = 50
rom_y = np.linspace(-1.5,1.5,N)
vel_a = np.zeros_like(rom_y)
vel_b = np.zeros_like(rom_y)
vel_c = np.zeros_like(rom_y)
vel_d = np.zeros_like(rom_y)
vel_e = np.zeros_like(rom_y)

for i in range(N):
    vel_a[i] = velocity_field(xt,yt,xr_a*r,rom_y[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    vel_b[i] = velocity_field(xt,yt,xr_b*r,rom_y[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    vel_c[i] = velocity_field(xt,yt,xr_c*r,rom_y[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    vel_d[i] = velocity_field(xt,yt,xr_d*r,rom_y[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    vel_e[i] = velocity_field(xt,yt,xr_e*r,rom_y[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    print i + 1

fig = plt.figure(2,figsize=(12.5,6))
fig.subplots_adjust(left=.05,right=.86,wspace=.36,hspace=.35)
plt.subplot(2,3,1)
plt.plot(yr_a/2.,uU_a,'b.')
plt.plot(rom_y,vel_a,'g-')
plt.xlim(-1.5,1.5)
plt.ylim(0.15,2.)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.45,1.75,r'$x/D$ = 1.25')

plt.subplot(2,3,2)
plt.plot(yr_b/2.,uU_b,'b.')
plt.plot(rom_y,vel_b,'g-')
plt.xlim(-1.5,1.5)
plt.ylim(0.15,2.)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.45,1.75,r'$x/D$ = 1.665')

plt.subplot(2,3,3)
plt.plot(yr_c/2.,uU_c,'b.',label='Experiment')
plt.plot(rom_y,vel_c,'g-',label='Model')
plt.xlim(-1.5,1.5)
plt.ylim(0.15,2.)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.45,1.75,r'$x/D$ = 2.5')
plt.legend(loc="upper left", bbox_to_anchor=(1,1))

plt.subplot(2,3,4)
plt.plot(yr_d/2.,uU_d,'b.')
plt.plot(rom_y,vel_d,'g-')
plt.xlim(-1.5,1.5)
plt.ylim(0.15,2.)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.45,1.75,r'$x/D$ = 4.165')

plt.subplot(2,3,5)
plt.plot(yr_e/2.,uU_e,'b.')
plt.plot(rom_y,vel_e,'g-')
plt.xlim(-1.5,1.5)
plt.ylim(0.15,2.)
plt.xlabel('$y/D$')
plt.ylabel(r'$u/U_\infty$')
plt.text(-0.45,1.75,r'$x/D$ = 5.835')

plt.show()
