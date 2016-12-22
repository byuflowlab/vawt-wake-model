import numpy as np
from numpy import log,fabs
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

fs = 18

## Grid Convergence
#extabs: cell:1777931, cell_turb:36152

cell = np.array([387391,477697,629345,1418137,2372024,4768290])#,4750946])
cell_turb = np.array([24385,27460,31954,43681,62579,110335])#,92991])
cp = np.array([0.321406528,0.42425241,0.431833702,0.438355123,0.441759825,0.43911712])#,0.436701608]) #b4001:b5001
# cp = np.array([0.325914367,0.427952694,0.435097835,0.441709186,0.445323805,0.440465415]) #b3501:b5001


# cell = np.array([481295,631815,1419944,2375520,4757450])
# cell_turb = np.array([27470,31957,44795,65691,99495])
# cp = np.array([0.42899673,0.43194673,0.439188886,0.442006833,0.437654271])

cell = np.array([387391,477697,629345,1427721,2384838,4802978])#,4768290])#,4750946])
cell_turb = np.array([24385,27460,31954,52572,75009,145023])#,110335])#,92991])
cp = np.array([0.321406528,0.42425241,0.431833702,0.443527399,0.446379056,0.444569346])#,0.43911712])#

cell = np.array([387391,477697,629345,1418137,2372024,4802978])#,4750946])
cell_turb = np.array([24385,27460,31954,43681,62579,145023])#,92991])
cp = np.array([0.321406528,0.42425241,0.431833702,0.438355123,0.441759825,0.444569346])#,0.436701608]) #b4001:b5001

# 0.42730301, 0.422036135

p = log((0.438355123-0.441759825)/(0.441759825-0.444569346))/log(1.4)
print 'p',p
for i in range(5):
    print cp[i+1] + (cp[i+1]-cp[i])/(1.4**p - 1)
    # print 0.444569346 + (0.444569346-0.441759825)/(1.4**2 - 1)
    
facts = 1.25
gc1 = facts*fabs((0.444569346-0.441759825)/0.444569346)/(1.4**p -1)
gc2 = facts*fabs((0.441759825-0.438355123)/0.441759825)/(1.4**p -1)
print 'error',gc1
print gc2/(1.4**p*gc1)
xerr = np.linspace(300000,10000000,100)
yerr = xerr*0.0 + (0.444569346 + (0.444569346-0.441759825)/(1.4**2 - 1))
xerrb = 629345.
yerrb = 0.444569346 + (0.444569346-0.441759825)/(1.4**2 - 1)

## Plotting
plt.figure(1)
# plt.subplot(1,2,1)
plt.semilogx(cell,cp,'bo-')
plt.plot(629345,0.431833702,'ro',markersize=8)
plt.plot(xerr,yerr,'r-')
plt.errorbar(xerrb,yerrb,yerr=gc1,color='r')
plt.xlabel('Cell Count',fontsize=fs)
plt.ylabel('Power Coefficient ($C_P$)',fontsize=fs)
plt.text(4000000,0.45,'Converged Value',fontsize=fs-4)
plt.text(550000,0.475,'Error Band',rotation=90,fontsize=fs-4)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
# plt.gca().invert_xaxis()
plt.xlim(300000,10000000)

# plt.subplot(1,2,2)
# plt.semilogx(cell_turb,cp,'ro-')
# plt.xlabel('Cell Count (Turbine)',fontsize=fs)
# plt.ylabel('Power Coefficient ($C_P$)',fontsize=fs)
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)

# plt.figure(2)
# plt.plot((1./cell)*1e6,cp,'bo-')
# plt.plot((1./629345.)*1e6,0.431833702,'ro',markersize=8)
# plt.xlabel('(Cell Count)$^{-1}$x$10^6$',fontsize=fs)
# plt.ylabel('Power Coefficient ($C_P$)',fontsize=fs)
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)

## Time Convergence
delta = np.array([1.5,1.0])
cp_vel = np.array([])

# vel1 = np.array([5.174577,5.178434])
# vel2 = np.array([6.204054,6.206861])

# plt.figure()
# plt.subplot(1,2,1)
# plt.plot(delta,vel1,'b')
# plt.xlabel('Time Step (s)',fontsize=fs)
# plt.ylabel('Sample Velocity 1 (m/s)',fontsize=fs)
# plt.gca().invert_xaxis()
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)
# plt.subplot(1,2,2)
# plt.plot(delta,vel2,'r')
# plt.xlabel('Time Step (s)',fontsize=fs)
# plt.ylabel('Sample Velocity 2 (m/s)',fontsize=fs)
# plt.gca().invert_xaxis()
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)

plt.show()
