import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def surf(x,y):
    a = -1.
    b = 1.
    c = 1.
    d = 1.
    e = 1.
    f = 1.
    g = 1.
    h = 1.
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
pnt = ax.scatter(X,Y,Z)


plt.show()