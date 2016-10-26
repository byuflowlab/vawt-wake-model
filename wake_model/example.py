import numpy as np
import matplotlib.pyplot as plt
from VAWT_Wake_Model import velocity_field
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import time

start = time.time()

# Option to plot velocity profiles
vel_slice = True
vel_slice = False # comment this out if desired on

# Option to plot a full velocity domain
plot_dist = True
# plot_dist = False # comment this out if desired on

# Enter the values desired

velf = 15.0 # free stream wind speed (m/s)
dia = 6.  # turbine diameter (m)
tsr = 4.  # tip speed ratio
B = 3. # number of blades
chord = 0.25 # chord lenth (m)
rot = tsr*velf/(dia/2.)

# Enter the positions of the turbine and velocity calculation
xt = 0. # downstream position of turbine in flow domain (m)
yt = 0. # later position of turbine in flow domain (m)
x0 = 24. # downstream distance from turbine for velocity calculation (m)
y0 = 0. # lateral distance from turbine for velocity calculation (m)

veltype = 'vort'

veltype = 'all'
# veltype = 'x'
# veltype = 'y'
# veltype = 'ind'

integration = 'simp'
m = 220
n = 200
# integration = 'gskr'

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

vel = velocity_field(xt,yt,x0,y0,velf,dia,rot,chord,B,param=None,veltype=veltype,integration=integration,m=m,n=n)
print '\nNormalized velocity at (',x0,',',y0,') =',vel,'\n' # output velocity (normalized by free stream wind speed)
pointval1 = 1.
pointval2 = 0.
pointval3 = 0.

## Plotting
fs = 25 # 18# font size for plots

# Plotting velocity profiles
if vel_slice == True:
    leng = 100. # data points in the velocity profile
    wide = 2.0*dia # width of the profile

    d_lab1 = str(wide/dia) # y-axis label
    d_lab2 = str(wide/(2*dia)) # y-axis label

    x = np.array([2*dia,4*dia,6*dia,8*dia,10*dia,15*dia]) # plotting at 2D, 4D, 6D, 8D, 10D, and 15D (changeable)
    y = np.linspace(-wide,wide,leng)

    pointval2 = leng*np.size(x)

    color = np.array(['b','c','g','y','r','m']) # identifying six colors to use for differentiation

    iterp = 0
    for i in range(int(np.size(x))):
        vel = np.array([])
        val = str(x[i]/dia)
        lab = '$x/D$ = '+val
        for j in range(int(np.size(y))):
            velp = velocity_field(xt,yt,x[i],y[j],velf,dia,rot,chord,B,param=None,veltype=veltype,integration=integration,m=m,n=n)
            vel = np.append(vel,velp)
            iterp += 1
            print 'Vel Slice ('+str(iterp)+' of '+str(pointval2)+')'
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
    xi = -3.*dia # starting point in downstream direction
    xf = 17.0*dia # ending point in downstream direction
    yd = -2.5*dia # lateral extent on down side
    yu = 2.5*dia # lateral extent on up side

    N = 100 # N**2 = number of data points in domain
    pointval3 = N*N

    xp = np.linspace(xi,xf,N)
    yp = np.linspace(yd,yu,N)
    [X, Y] = np.meshgrid(xp,yp)
    VEL = np.zeros((N, N)) # initiallizing velocity data point array
    VELy = np.zeros((N, N))

    iter = 0
    for i in range(N):
        for j in range(N):
            if veltype == 'all' or veltype == 'x' or veltype == 'y' or veltype == 'velfort':
                VEL[i,j] = velocity_field(xt,yt,X[i,j],Y[i,j],velf,dia,rot,chord,B,param=None,veltype=veltype,integration=integration,m=m,n=n)
            elif veltype == 'ind':
                velfd = velocity_field(xt,yt,X[i,j],Y[i,j],velf,dia,rot,chord,B,param=None,veltype=veltype,integration=integration,m=m,n=n)
                VEL[i,j] = velfd[0]
                VELy[i,j] = velfd[1]
            elif veltype == 'vort':
                VEL[i,j] = velocity_field(xt,yt,X[i,j],Y[i,j],velf,dia,rot,chord,B,param=None,veltype=veltype,integration=integration,m=m,n=n)
            iter = iter +1
            print 'Plot ('+str(iter)+' of '+str(N*N)+')'

    if veltype == 'all' or veltype == 'x' or veltype == 'y' or veltype == 'velfort':
        fig = plt.figure(2,figsize=(19,5))
        fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
        if veltype == 'all' or veltype == 'x' or veltype == 'velfort':
            lb = 0.15 # lower bound on velocity to display
            ub = 1.15 # upper bound on velocity to display

            # lb = -0.75 # lower bound on velocity to display ind x-vel
            # ub = 0.4 # upper bound on velocity to display ind x-vel
            # lb = 0.0 # lower bound on velocity to display ind x-vel
            # ub = 25./15. # upper bound on velocity to display ind x-vel

        elif veltype == 'y':
            lb = -0.35 # lower bound on velocity to display
            ub = 0.35 # upper bound on velocity to display
        ran = 32 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,VEL,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        CB.ax.set_ylabel(r'$u/U_\infty$',fontsize=fs)
        CB.ax.tick_params(labelsize=fs)
        CB.ax.set_aspect(20)
    elif veltype == 'ind':
        fig = plt.figure(2,figsize=(19,5))
        # fig.subplots_adjust(bottom=.16,left=.05,right=0.81)
        # CS = plt.quiver(X/dia,Y/dia, VEL, VELy)
        # plt.quiverkey(Q,1.1,0.8,1,'1 m/s',fontproperties={'size':fs})
        fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
        speed = np.sqrt(VEL*VEL + VELy*VELy)
        CS = plt.streamplot(X/dia, Y/dia, VEL, VELy, density=2, color=speed, cmap=plt.cm.coolwarm)
        CB = fig.colorbar(CS.lines,ticks=np.linspace(0.0,0.75,6))
        CB.ax.set_ylabel(r'$velocity mag/U_\infty$',fontsize=fs)
        CB.ax.tick_params(labelsize=fs)
        CB.ax.set_aspect(20)
    elif veltype == 'vort':
        fig = plt.figure(2,figsize=(19,5))
        fig.subplots_adjust(bottom=.16,left=.05,right=1.0)
        lb = -1.0/rot # lower bound on velocity to display
        ub = 1.0/rot # upper bound on velocity to display
        ran = 32 # number of contours between the velocity bounds
        bounds = np.linspace(lb,ub,ran)
        v = np.linspace(lb,ub,6) # setting the number of tick marks on colorbar
        CS = plt.contourf(X/dia,Y/dia,VEL,ran,vmax=ub,vmin=lb,levels=bounds,cmap=plt.cm.coolwarm) # plotting the contour plot
        CB = plt.colorbar(CS, ticks=v) # creating colorbar
        CB.ax.set_ylabel(r'$\gamma/\Omega$',fontsize=fs)
        CB.ax.tick_params(labelsize=fs)
        CB.ax.set_aspect(20)
    plt.xlabel('$x/D$',fontsize=fs)
    plt.ylabel('$y/D$',fontsize=fs)
    plt.xticks(fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.xlim(xi/dia,xf/dia)
    plt.ylim(yd/dia,yu/dia)
    circ = plt.Circle((xt/dia,yt/dia),0.5,edgecolor='k',fill=False)
    plt.gca().add_patch(circ)


    # if veltype == 'all':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/mag-vel.png')
    # elif veltype == 'x':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/x-vel_newplot.png')
    # elif veltype == 'y':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/y-vel.png')
    # elif veltype == 'ind':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/ind-vel.png')
    # elif veltype == 'vort':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/vort-plot.png')
    # elif veltype == 'velfort':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/mag-velfort_xind.png')

sec = time.time()-start
pointval = pointval1+pointval2+pointval3
print '\nRun Statistics:\n\tCalculated',int(pointval),'points'
print '\tRan for:\n\t\t',sec,'seconds\n\t\t',sec/60.,'minutes\n\t',sec/pointval,'seconds per point (average)'

plt.show()