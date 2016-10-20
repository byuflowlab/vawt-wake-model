import numpy as np
import matplotlib.pyplot as plt
from database_call import vorticity,velocity,vorticity2,velocity2
from VAWT_Wake_Model import velocity_field
from matplotlib import rcParams
rcParams['font.family'] = 'Times New Roman'
import time
import csv

def parameterval(tsr,sol,coef):
    a = coef[0]
    b = coef[1]
    c = coef[2]
    d = coef[3]
    e = coef[4]
    f = coef[5]
    g = coef[6]
    h = coef[7]
    i = coef[8]
    j = coef[9]

    return a + b*tsr + c*sol + d*tsr**2 + e*tsr*sol + f*sol**2 + g*tsr**3 + h*tsr**2*sol + i*tsr*sol**2 + j*sol**3

start = time.time()

# Enter the values desired

velf = 15.0 # free stream wind speed (m/s)
dia = 6.  # turbine diameter (m)
tsr = 4.  # tip speed ratio
B = 3. # number of blades
chord = 0.25 # chord lenth (m)
solidity = (chord*B)/(dia/2.)

# Enter the positions of the turbine and velocity calculation
xt = 0. # downstream position of turbine in flow domain (m)
yt = 0. # later position of turbine in flow domain (m)
x0 = 24. # downstream distance from turbine for velocity calculation (m)
y0 = 0. # lateral distance from turbine for velocity calculation (m)

# Choose whether CFD vorticity or velocity data will be used as the basis
cfd_data = 'vort'
cfd_data = 'vort2'
# cfd_data = 'velo'
# cfd_data = 'velo2'

veltype = 'all'
veltype = 'x'
veltype = 'y'
veltype = 'ind'
veltype = 'vort'
veltype = 'velfort'

if cfd_data == 'vort':
    loc,spr,skw,scl = vorticity(tsr,solidity)
    param = np.array([loc,spr,skw,scl])

    # print param
    # import time
    #
    # time.sleep(10)
    
elif cfd_data == 'velo':
    men1,sdv1,rat1,wdt1,spr1,scl1,tsrn1,_ = velocity(tsr-0.1249,solidity)
    men2,sdv2,rat2,wdt2,spr2,scl2,tsrn2,_ = velocity(tsr+0.1249,solidity)
    if solidity >= 0.35:
        men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,solidity-0.1249)
        men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,solidity+0.1249)
    elif solidity >=0.25:
        men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,solidity-0.049)
        men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,solidity+0.1249)
    else:
        men3,sdv3,rat3,wdt3,spr3,scl3,_,soln1 = velocity(tsr,solidity-0.049)
        men4,sdv4,rat4,wdt4,spr4,scl4,_,soln2 = velocity(tsr,solidity+0.049)
    if tsrn1 == tsrn2:
        p = 0.
    else:
        p = (tsr-tsrn1)/(tsrn2-tsrn1)
    if soln1 == soln2:
        q = 0.
    else:
        q = (solidity-soln1)/(soln2-soln1)

    param = np.array([men1,sdv1,rat1,wdt1,spr1,scl1,men2,sdv2,rat2,wdt2,spr2,scl2,men3,sdv3,rat3,wdt3,spr3,scl3,men4,sdv4,rat4,wdt4,spr4,scl4,p,q])

elif cfd_data == 'vort2':
    loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3 = vorticity2(tsr,solidity)
    param = np.array([loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3])

    # print param
    # import time
    #
    # time.sleep(10)

elif cfd_data == 'velo2':
    spr1,pow1,pow2,spr2,skw,scl1,scl2,scl3 = velocity2(tsr,solidity)
    param = np.array([spr1,pow1,pow2,spr2,skw,scl1,scl2,scl3])



vel = velocity_field(xt,yt,xt + x0,yt + y0,velf,dia,tsr,solidity,cfd_data,param,veltype)

print '\nNormalized velocity at (',x0,',',y0,') from the turbine =',vel,'\n' # output velocity (normalized by free stream wind speed)

## Plotting
fs = 25 # 18# font size for plots

# Option to plot velocity profiles
vel_slice = True
vel_slice = False # comment this out if desired on

# Option to plot a full velocity domain
plot_dist = True
# plot_dist = False # comment this out if desired on

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
            velp = velocity_field(xt,yt,x[i],y[j],velf,dia,tsr,solidity,cfd_data,param,veltype)
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
    xi = -3.*dia # starting point in downstream direction
    xf = 17.0*dia # ending point in downstream direction
    yd = -2.5*dia # lateral extent on down side
    yu = 2.5*dia # lateral extent on up side

    N = 100 # N**2 = number of data points in domain

    xp = np.linspace(xi,xf,N)
    yp = np.linspace(yd,yu,N)
    [X, Y] = np.meshgrid(xp,yp)
    VEL = np.zeros((N, N)) # initiallizing velocity data point array
    VELy = np.zeros((N, N))

    # vel = np.zeros((10000,3))
    # iter = 0
    # for i in range(N):
    #     for j in range(N):
    #         vel[iter,2],m,n = velocity_field(xt,yt,xp[i],yp[j],velf,dia,tsr,solidity,cfd_data,param,veltype)
    #         vel[iter,0] = xp[i]
    #         vel[iter,1] = yp[j]
    #         iter += 1
    #         print 'CSV ('+str(iter)+' of '+str(N*N)+')'+str(m)+'x'+str(n)
    #
    # if veltype == 'all':
    #     fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/data/vel_all-plot.csv'
    # elif veltype == 'velfort':
    #     fdata = '/Users/ning1/Documents/Flow Lab/VAWTWakeModel/wake_model/data/vel_fort-plot_'+str(m)+'x'+str(n)+'.csv'
    #
    # with open(fdata,'w') as fp:
    #     a = csv.writer(fp)
    #
    #     data = np.array(['x','y','vel'])
    #
    #
    #     vals = np.zeros((10000,3)).astype(np.str)
    #     for i in range(N*N):
    #         vals[i,0] = vel[i,0]
    #         vals[i,1] = vel[i,1]
    #         vals[i,2] = vel[i,2]
    #
    #     data = np.vstack([data,vals])
    #
    #     a.writerows(data)



    iter = 0
    for i in range(N):
        for j in range(N):
            if veltype == 'all' or veltype == 'x' or veltype == 'y' or veltype == 'velfort':
                VEL[i,j],m,n = velocity_field(xt,yt,X[i,j],Y[i,j],velf,dia,tsr,solidity,cfd_data,param,veltype)
            elif veltype == 'ind':
                velfd,m,n = velocity_field(xt,yt,X[i,j],Y[i,j],velf,dia,tsr,solidity,cfd_data,param,veltype)
                VEL[i,j] = velfd[0]
                VELy[i,j] = velfd[1]
            elif veltype == 'vort':
                VEL[i,j],m,n = velocity_field(xt,yt,X[i,j],Y[i,j],velf,dia,tsr,solidity,cfd_data,param,veltype)
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
        rot = tsr*velf/(dia/2.)
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


    coef0 = np.array( [0.0025703809856661534, -0.0007386258659065129, 0.004595508188667984, 0.000380123563204793, -0.0005090098755683027, 0.005744581813281894, -4.103393770815313e-05, -0.0014146918534486358, -0.013975958482495927, 0.0] )
    coef1 = np.array( [-0.5047504670963536, 0.23477391362058556, 0.8414256436198028, -0.04252528528617351, -0.06962875967504166, -0.6566907653208429, 0.002839318332370807, 0.00571803958194812, 0.0070744372783060295, 0.22805286438890995] )
    coef2 = np.array( [0.2878345841026334, 0.11512552658662782, 0.7303949879914625, -0.007035517839387948, -0.18284850673545897, -0.5241921153256568, -0.0003704899921255296, 0.010972527139685873, 0.04380801537377295, 0.1724129349605399] )
    coef3 = np.array( [0.08234816067475287, -0.03530687906626052, -0.3662863944976986, 0.003240141344532779, 0.12172015102204112, 0.2993048183466721, 0.0, -0.009253185586804007, -0.057469126406649716, -0.07257633583877886] )
    coef4 = np.array( [-0.07083579909945328, 0.016182024377569406, 0.1985436342461859, 0.0017738254727425816, -0.09111094817943823, -0.06561408122153217, -0.0005115133402638633, 0.009434288536679505, 0.022392136905926813, 0.0] )
    coef5 = np.array( [-1.6712830849073221, 1.5625053380692426, -6.180392756736983, -0.20407668040293722, -4.6476103643607685, 29.380064536220306, 0.0, 0.7502978877582536, -0.16358232641365608, -19.937609244085568] )
    coef6 = np.array( [-3.423561091777921, -9.228795430171687, 86.95722105482042, 2.772872601988039, -11.968168333741515, -150.61261090270446, -0.24715316589674527, 0.5283723108899993, 4.537286811245538, 82.50581844010263] )
    coef7 = np.array( [-0.19815381951708524, 0.08438758133540872, 1.2650146439483734, -0.007606115512168328, -0.2747023984740461, -0.8844640101378567, 0.0, 0.01870057580949183, 0.0699898278743648, 0.2794360008051127] )
    coef8 = np.array( [2.3932787625531815, -2.020874419612962, -8.938221963838357, 0.576323845480877, 2.8782448498416944, 16.598492450314534, -0.04746016700352029, -0.197101203594028, -1.3860007472886064, -8.289767128060362] )
    coef9 = np.array( [104.40501489600803, -29.942999569370276, -174.42008279158216, 3.708514822202037, 25.14336546356742, 132.35546551746415, -0.16479555172343271, -1.351556690339512, -6.721810844025761, -40.39565289044579] )

    x0t = np.linspace(0,xf/dia,100)
    loc1 = parameterval(tsr,solidity,coef0)
    loc2 = parameterval(tsr,solidity,coef1)
    loc3 = parameterval(tsr,solidity,coef2)
    spr1 = parameterval(tsr,solidity,coef3)
    spr2 = parameterval(tsr,solidity,coef4)
    skw1 = parameterval(tsr,solidity,coef5)
    skw2 = parameterval(tsr,solidity,coef6)
    scl1 = parameterval(tsr,solidity,coef7)
    scl2 = parameterval(tsr,solidity,coef8)
    scl3 = parameterval(tsr,solidity,coef9)
    plt.figure(20)
    plt.subplot(2,2,1)
    plt.plot(x0t,loc1*x0t**2+loc2*x0t+loc3)
    plt.title('LOC')
    plt.text(2,0,str(loc1)+'x^2 + '+str(loc2)+'x + '+str(loc3))
    plt.subplot(2,2,2)
    plt.plot(x0t,spr1*x0t+spr2)
    plt.title('SPR')
    plt.text(2,0,str(spr1)+'x + '+str(spr2))
    plt.subplot(2,2,3)
    plt.plot(x0t,skw1*x0t+skw2)
    plt.title('SKW')
    plt.text(2,0,str(skw1)+'x + '+str(skw2))
    plt.subplot(2,2,4)
    plt.plot(x0t,scl1/(1.0 + np.exp(scl2*(x0t - scl3))))
    plt.title('SCL')
    plt.text(2,0,str(scl1)+'/(1 + exp('+str(scl2)+'(x - '+str(scl3)+'))')

    # if veltype == 'all':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/mag-vel.png')
    # elif veltype == 'x':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/x-vel.png')
    # elif veltype == 'y':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/y-vel.png')
    # elif veltype == 'ind':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/ind-vel.png')
    # elif veltype == 'vort':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/vort-plot.png')
    # elif veltype == 'velfort':
    #     plt.savefig('/Users/ning1/Documents/FLOW Lab/mag-velfort_xind.png')

print time.time()-start

plt.show()


