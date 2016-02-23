import numpy as np
import matplotlib.pyplot as plt
import csv
from VAWT_Wake_Model import velocity_field

from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

r = 0.5
v = 1.


x15 = np.linspace(-1.7,1.7,21)
x15r = np.linspace(-1.7,1.7,100)

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

for i in range(np.size(x15)):
    rom15t[i] = velocity_field(0.,0.,1.5*dia,x15[i]*dia,velf,dia,tsr,sol)
    error_test[i] = (rom15t[i]-y15o[i])/y15o[i]
error = np.average(error_test)
errorstd = np.std(error_test)
for i in range(np.size(rom15)):
    rom15[i] = velocity_field(0.,0.,1.5*dia,x15r[i]*dia,velf,dia,tsr,sol)
    print i

fs = 15

plt.figure()
# plt.subplot(2,3,1)
# plt.plot(x15,y15c,'r.',label='Experimental (closed)')
plt.plot(x15,y15o,'b.',label='Experimental (open)')
plt.plot(x15r,rom15,'g-',label='Model')
plt.xlim(-1.75,1.75)
plt.ylim(0.3,1.4)
plt.xlabel('$y/D$',fontsize=fs)
plt.ylabel(r'$u/U_\infty$',fontsize=fs)
plt.legend(loc=1,fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
# plt.text(-0.5,0.9,'X/D = 1.5')
# print '1.5 modc',(min(rom15)-min(y15c))/min(y15c)
print '1.5 modo',(min(rom15)-min(y15o))/min(y15o),error,errorstd
plt.show()


# 1.5 modo -0.0511120658493 0.0466994707498 0.117194753988

