import numpy as np
from numpy import log,fabs
import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

fs = 18

# Grid Convergence

#grid1,grid2,grid3,grid4,grid5,grid6,grid7
cell = np.array([387391,477697,629345,1418137,2372024,4802978])
cell_turb = np.array([24385,27460,31954,43681,62579,145023])
cp = np.array([0.321406528,0.42425241,0.431833702,0.438355123,0.441759825,0.444569346]) #b4001:b5001

p = log((0.438355123-0.441759825)/(0.441759825-0.444569346))/log(1.4)
print 'p',p
for i in range(5):
    print cp[i+1] + (cp[i+1]-cp[i])/(1.4**p - 1)
    
facts = 1.25
gc1 = facts*fabs((0.444569346-0.441759825)/0.444569346)/(1.4**p -1)
gc2 = facts*fabs((0.441759825-0.438355123)/0.441759825)/(1.4**p -1)
print 'error',gc1
print gc2/(1.4**p*gc1)
xerr = np.linspace(300000,10000000,100)
yerr = xerr*0.0 + (0.444569346 + (0.444569346-0.441759825)/(1.4**2 - 1))
xerrb = 629345.
yerrb = 0.444569346 + (0.444569346-0.441759825)/(1.4**2 - 1)

# Plotting
plt.figure(1)
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
plt.xlim(300000,10000000)

plt.show()