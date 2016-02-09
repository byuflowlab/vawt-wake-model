import numpy as np
import matplotlib.pyplot as plt
from numpy import pi
from VAWT_Wake_Model import velocity_field

# Enter the values desired

velf = 15.0 # free stream wind speed (m/s)
dia = 6.0  # turbine diameter (m)
tsr = 4.0  # tip speed ratio
B = 3. # number of blades
chord = 0.25 # chord lenth (m)
solidity = (chord*B)/(dia/2.)

# Enter the positions of the turbine and velocity calculation
xt = 0. # downstream position of turbine (m)
yt = 0. # later position of turbine (m)
x0 = 5. # downstream position of velocity calculation from turbine (m)
y0 = 0. # lateral position of velocity calculation from turbine (m)

vel = velocity_field(xt,yt,xt + x0,yt + y0,velf,dia,tsr,solidity)

print '\nNormalized velocity at (',x0,',',y0,') from the turbine =',vel,'\n' # output velocity (normalized by free stream wind speed)

## Plotting
fs = 15 # font size for plots

# Option to plot velocity profiles
vel_slice = True
vel_slice = False # comment this out if desired on

# Option to plot a full velocity domain
plot_dist = True
plot_dist = False # comment this out if desired on

# Plotting velocity profiles
if vel_slice == True:
    leng = 100 # data points in the velocity profile
    wide = 2.0*dia # width of the profile
    
    d_lab1 = str(wide/dia) # y-axis label
    d_lab2 = str(wide/(2*dia)) # y-axis label
    
    x = np.array([2*dia,4*dia,6*dia,8*dia,10*dia,15*dia]) # plotting at 2D, 4D, 6D, 8D, 10D, and 15D (changeable)
    y = np.linspace(-wide,wide,leng)
    
    color = np.array(['b','c','g','y','r','m']) # identifying six colors to use for differentiation
    
    iterp = 0
    for i in range(int(np.size(x))):
        vel = np.array([])
        val = str(x[i]/dia)
        lab = '$x/D$ = '+val
        for j in range(int(np.size(y))):
            velp = velocity_field(xt,yt,x[i],y[j],velf,dia,tsr,solidity)
            vel = np.append(vel,velp)
            iterp += 1
            print 'Vel Slice ('+str(iterp)+' of '+str(leng*np.size(x))+')'
        plt.figure(1)
        plt.plot(vel,y,color[i],label=lab)
    
    tix = np.array([-wide,-wide/2.,0.,wide/2.,wide])
    tixl = np.array([d_lab1,d_lab2,'0.0',d_lab2,d_lab1])
    # plt.legend(loc="upper left",bbox_to_anchor=(1, 1),fontsize=fs) # legend off to the side
    plt.legend(loc=2,fontsize=fs) # legend in the plot
    plt.xlim(0.,1.2)
    plt.ylim(-wide,wide)
    plt.xticks(fontsize=fs)
    plt.yticks(tix,tixl,fontsize=fs)
    plt.xlabel(r'$u/U_\infty$', fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)

    # Plotting full velocity domain
if plot_dist == True:
    xi = -1.0*dia # starting point in downstream direction
    xf = 10.0*dia # ending point in downstream direction
    yd = -2.0*dia # lateral extent on down side
    yu = 2.0*dia # lateral extent on up side
    
    N = 100 # N**2 = number of data points in domain
    
    xp = np.linspace(xi,xf,N)
    yp = np.linspace(yd,yu,N)
    [X, Y] = np.meshgrid(xp,yp)
    VEL = np.zeros((N, N)) # initiallizing velocity data point array
    
    iter = 0
    for i in range(N):
        for j in range(N):
            VEL[i,j] = velocity_field(xt,yt,X[i,j],Y[i,j],velf,dia,tsr,solidity)
            iter = iter +1
            print 'Plot ('+str(iter)+' of '+str(N*N)+')'
    
    plt.figure(2)
    lb = 0.15 # lower bound on velocity to display
    ub = 1.15 # upper bound on velocity to display
    ran = 200 # number of contours between the velocity bounds
    bounds = np.linspace(lb,ub,ran)
    v = np.linspace(lb,ub,5) # setting the number of tick marks on colorbar
    CS = plt.contourf(X/dia,Y/dia,VEL,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.jet) # plotting the contour plot
    CB = plt.colorbar(CS, ticks=v) # creating colorbar
    CB.ax.set_ylabel(r'$u/U_\infty$',fontsize=fs)
    CB.ax.tick_params(labelsize=fs)
    plt.xlabel('$x/D$',fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    circ = plt.Circle((xt/dia,yt/dia),0.5,edgecolor='k',fill=False)
    plt.gca().add_patch(circ)

plt.show()
