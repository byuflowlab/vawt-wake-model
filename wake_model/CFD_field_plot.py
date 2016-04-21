import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import re
import csv
from scipy.interpolate import SmoothBivariateSpline, griddata
import matplotlib.mlab as ml
from matplotlib.tri import Triangulation, UniformTriRefiner,CubicTriInterpolator
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'

fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/velocity.csv'
dia = 6.
wind = 15.

for i in range(101):
    name = str(i+1)
    exec('pos'+name+' = np.array([])')
    exec('velo'+name+' = np.array([])')

f = open(fdata)

csv_f = csv.reader(f)

i = 0
for row in csv_f:
    for j in range(101):
        npos = str(j+1)
        nrowe = str(2*j)
        nrowo = str(2*j+1)
        nulls = 'null'

        exec('if row['+nrowe+'] != nulls and i != 0:\n\tpos'+npos+' = np.append(pos'+npos+',float(row['+nrowe+']))')
        exec('if row['+nrowo+'] != nulls and i != 0:\n\tvelo'+npos+' = np.append(velo'+npos+',float(row['+nrowo+']))')
    i += 1

f.close()
print 'DATA IMPORTED'

region = 0.#0.3
for i in range(101):
    name = str(i+1)

    #Ordering the data numerically by the position
    exec('pos'+name+', velo'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', velo'+name+'))))')
    #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
    exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvelo'+name+' = np.copy(velo'+name+'_0)')
    
    # exec('delvec = np.array([])\nfor j in range(np.size(pos'+name+')):\n\tif pos'+name+'[j] < dia/2.+region and pos'+name+'[j] > -dia/2.-region:\n\t\tdelvec = np.append(delvec,j)\ndelvec = delvec[::-1]\nfor j in range(np.size(delvec)):\n\tvelo'+name+' = np.delete(velo'+name+',delvec[j])\n\tpos'+name+' = np.delete(pos'+name+',delvec[j])')

max = 0
for i in range(101):
    name = str(i+1)
    exec('pos'+name+' = pos'+name+'/dia')
    exec('velo'+name+' = velo'+name+'/wind')
    
    exec('if np.size(pos'+name+') > max:\n\tmax = np.size(pos'+name+')')

x = np.linspace(-18,102,101)/dia

xs = np.array([])
ys = np.array([])
vels = np.array([])

for i in range(101):
    name = str(i+1)
    
    exec('xs = np.append(xs,np.ones_like(pos'+name+')*x[i])')
    exec('ys = np.append(ys,pos'+name+')')
    exec('vels = np.append(vels,velo'+name+')')

points = np.zeros((np.size(xs),2))
for i in range(np.size(xs)):
    points[i,0] = xs[i]
    points[i,1] = ys[i]
    


N = 100 #number of points for plotting/interpolation

xi = np.linspace(xs.min(), xs.max(), N)
yi = np.linspace(ys.min(), ys.max(), 7*N)
evalp = np.zeros((np.size(xi),2))
for i in range(np.size(xi)):
    points[i,0] = xi[i]
    points[i,1] = yi[i]

# [X,Y] = np.meshgrid(xi,yi)
zi = ml.griddata(xs, ys, vels, xi, yi, 'linear')#, method='cubic')

# triang = Triangulation(xs, ys)
# refiner = UniformTriRefiner(triang)
# tri_refi, z_test_refi = refiner.refine_field(vels, subdiv=3)

# print 'xi = np.array(',xi.tolist(),')'
# print 'yi = np.array(',yi.tolist(),')'
# print 'zi = np.array(',zi.tolist(),')'

fs = 25
fig = plt.figure(figsize=(19,5))
fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
lb = 0.15 # lower bound on velocity to display
ub = 1.15 # upper bound on velocity to display
ran = 32 # number of contours between the velocity bounds
bounds = np.linspace(lb,ub,ran)
v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
CS = plt.contourf(xi,yi,zi,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
# CS = plt.tricontourf(tri_refi, z_test_refi,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
CB = plt.colorbar(CS, ticks=v) # creating colorbar
CB.ax.set_ylabel(r'$u/U_\infty$',fontsize=fs)
CB.ax.tick_params(labelsize=fs)
# CB.ax.set_aspect(20)
plt.xlabel('$x/D$',fontsize=fs)
plt.ylabel('$y/D$',fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
circ = plt.Circle((0,0),0.5,edgecolor='k',fill=False)
plt.gca().add_patch(circ)
plt.xlim(-3,17)
plt.ylim(-2.5,2.5)
plt.show()


# surface = SmoothBivariateSpline(xs,ys,vels)
# xp = np.copy(x)
# yp = np.linspace(min(ys),max(ys),101)








# [X,Y] = np.meshgrid(xs,ys)
# 
# Z = np.zeros_like(X)
# for i in range(np.size(xp)):
#     for j in range(np.size(yp)):
#         Z[i,j] = surface(X[i,j],Y[i,j])
# fs = 18
# plt.figure()
# lb = 0.15 # lower bound on velocity to display
# ub = 1.15 # upper bound on velocity to display
# ran = 32 # number of contours between the velocity bounds
# bounds = np.linspace(lb,ub,ran)
# v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
# CS = plt.contourf(X/dia,Y,Z,ran,cmap=plt.cm.coolwarm)#,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
# CB = plt.colorbar(CS)#, ticks=v) # creating colorbar
# CB.ax.set_ylabel(r'$u/U_\infty$',fontsize=fs)
# CB.ax.tick_params(labelsize=fs)
# CB.ax.set_aspect(40)
# plt.xlabel('$x/D$',fontsize=fs)
# plt.ylabel('$y/D$',fontsize=fs)
# plt.xticks(fontsize=fs)
# plt.yticks(fontsize=fs)
# plt.xlim(-3,17)
# plt.ylim(-2.5,2.5)
# circ = plt.Circle((0./dia,0./dia),0.5,edgecolor='k',fill=False)
# plt.gca().add_patch(circ)
# plt.show()












# #Ordering the data numerically by the position
# exec('pos'+name+', velo'+name+' = (list(t) for t in zip(*sorted(zip(pos'+name+', velo'+name+'))))')
# #STAR-CCM+ data contained repeated values; this creates new sets of data with repeats eliminated
# exec('pos'+name+'_0 = np.array([])\nvelo'+name+'_0 = np.array([])\nfor i in range(np.size(pos'+name+')):\n\tif pos'+name+'[i] not in pos'+name+'_0:\n\t\tpos'+name+'_0 = np.append(pos'+name+'_0,pos'+name+'[i])\n\t\tvelo'+name+'_0 = np.append(velo'+name+'_0,velo'+name+'[i])\npos'+name+' = np.copy(pos'+name+'_0)\nvelo'+name+' = np.copy(velo'+name+'_0)')
# 
# 
# 
# #3D Plotting
# x,y,z = zip(*holder)
# # fig = plt.figure()
# # ax = Axes3D(fig)
# # ax.plot(x,y,z)
# 
# [X,Y] = np.meshgrid(x,y)
# 
# Z = np.zeros_like(X)
# for i in range(np.size(x)):
#     for j in range(np.size(y)):
#         Z[i,j] = z[i]
# 
# plt.figure()
# print x,y,z
# plt.contourf(X,Y,Z)
# plt.show()