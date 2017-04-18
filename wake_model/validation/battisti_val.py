import numpy as np
import matplotlib.pyplot as plt
from VAWT_Wake_Model import velocity_field
from numpy import fabs

from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

r = 0.5 # radius
v = 1.0 # velocity

veltype = 'x'
errortype = 'abs'
errortype = 'rel'
epsilon = 1e-3
errortype = 'rms'

direction = 'horz'
direction = 'vert'

x15 = np.linspace(-1.7,1.7,21)
x15r = np.linspace(-1.75,1.75,100)

y15c = np.array([1.06,1.06,1.07,1.09,1.1,1.11,1.09,0.91,0.65,0.56,0.66,0.62,0.8,1.05,1.12,1.115,1.115,1.11,1.105,1.1,1.095])
y15o = np.array([0.98,0.98,0.985,0.99,0.995,1.005,0.96,0.73,0.5,0.44,0.54,0.51,0.66,0.9,1.01,1.0,0.99,0.985,0.98,0.98,0.97])

rad = 0.515
dia = 2*rad
velf = 16.0
sol = 0.5
tsr = 1.6
rot = tsr*velf/rad
B = 3
chord = sol*rad/B

rom15 = np.zeros_like(x15r)
rom15t = np.zeros_like(x15)
error_test = np.zeros_like(x15)

for i in range(np.size(x15)):
    rom15t[i] = velocity_field(0.,0.,1.5*dia,x15[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)

    if errortype == 'abs':
        error_test[i] = fabs((1.-rom15t[i])-(1.-y15o[i]))
    elif errortype == 'rel':
        if fabs(1.-y15o[i]) >= epsilon:
            error_test[i] = fabs(((1.-rom15t[i])-(1.-y15o[i]))/(1.-y15o[i]))
        else:
            error_test[i] = fabs((1.-rom15t[i])-(1.-y15o[i]))
    elif errortype == 'rms':
        error_test[i] = ((1.-rom15t[i])-(1.-y15o[i]))**2.
    print 'error calc',i+1,'of',np.size(x15)

if errortype == 'abs' or errortype == 'rel':
    error = np.average(fabs(error_test))
    errorstd = np.std(fabs(error_test))
elif errortype == 'rms':
    error = np.sqrt(np.average(error_test))
    errorstd = 1.
for i in range(np.size(rom15)):
    rom15[i] = velocity_field(0.,0.,1.5*dia,x15r[i]*dia,velf,dia,rot,chord,B,param=None,veltype=veltype)
    print 'plot point',i+1,'of',np.size(x15r)

fs = 19 # journal
#fs = 20 # thesis

fig = plt.figure(1)
fig.subplots_adjust(bottom=.12)
if direction == 'horz':
    # plt.plot(y15c,x15,'r.',label='Experimental (closed)')
    # plt.plot(y15o,x15,'g.',label='Experimental (open)')
    plt.plot(y15o,x15,'.',color='k',label='Experimental')
    plt.plot(rom15,x15r,'r-',label='Model')
    plt.ylim(-1.75,1.75)
    plt.xlim(0.3,1.2)
    plt.ylabel('$y/D$',fontsize=fs)
    plt.xlabel(r'$u/U_\infty$',fontsize=fs)
elif direction == 'vert':
    # plt.plot(x15,y15c,'r.',label='Experimental (closed)')
    # plt.plot(x15,y15o,'g.',label='Experimental (open)')
    plt.plot(x15,y15o,'.',color='k',label='Experimental')
    plt.plot(x15r,rom15,'r-',label='Model')
    plt.xlim(-1.75,1.75)
    plt.ylim(0.3,1.4)
    plt.xlabel('$y/D$',fontsize=fs)
    plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
if direction == 'horz':
    plt.legend(loc=2,fontsize=fs)
elif direction == 'vert':
    plt.legend(loc=1,fontsize=fs)
# plt.text(-0.5,0.9,'X/D = 1.5')

print '----1.5 modo----'
if errortype == 'abs':
    print '\n Average Absolute Error'
elif errortype == 'rel':
    print '\n Average Relative Error'
elif errortype == 'rms':
    print '\n Root Mean Squared Error'
print '\nAverage Error:',error
print 'Average Stand Dev:',errorstd
low = y15o.min()
lowrom = rom15.min()
if errortype == 'abs' or errortype == 'rms':
     mdferror = fabs((1.-lowrom)-(1.-low))
elif errortype == 'rel':
    mdferror = fabs(((1.-lowrom)-(1.-low))/(1.-low))
print 'Maximum Deficit Error:',mdferror


plt.show()
