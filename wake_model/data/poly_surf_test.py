import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def surf(x,y):
    #cub
    a = 0.
    b = 0.
    c = 0.
    d = 0.
    #quad
    e = .1
    f = .1
    g = .1
    #lin
    h = 10.
    i = 1.
    j = 1.
    
    return a*x**3 + b*y**3 + c*x**2*y + d*x*y**2 + e*x**2 + f*y**2 + g*x*y + h*x + i*y + j

n = 100
x = np.linspace(1.5,7,n)
y = np.linspace(.15,1.,n)


X,Y = np.meshgrid(x,y)
Z = np.zeros((n,n))

for i in range(n):
    for j in range(n):
        Z[i,j] = surf(X[i,j],Y[i,j])
        
fig = plt.figure(1)
ax = fig.gca(projection='3d')
# pnt = ax.scatter(X,Y,Z)
surf = ax.plot_surface(X,Y,Z,rstride=1,cstride=1,color='g',linewidth=0,antialiased=True,alpha=0.5)#cmap=cm.parula


plt.show()