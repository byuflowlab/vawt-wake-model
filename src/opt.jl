#Pkg.build("Snopt")
using Snopt
using vawtwake
using BPM
using PyCall
using vawtac
@pyimport matplotlib.pyplot as plt


#Replacing Eric's Globals with a struct
struct globals
    turb_dia
    dia
    rot
    chord
    twist
    delta
    B
    H
    rho
    mu
    Vinf
    windroseDirections
    windFrequencies
    obs
    Hub
    af
    power_iso_tot
    Vpn
    Vnn
    Vtp
    Cpp
    Cpn
    ntheta
    thetavec
    useAC::Bool
    SPLlim::Float64
end #globals


function full_obj(xopt,gl)
    #Get Number of turbines
    nturb=length(gl.dia)
    #Global Variables?
    x=xopt[1:nturb]
    y=xopt[nturb+1:end]

    nwind = length(gl.windroseDirections)  #numberof wind directions

    power_turb= zeros(nturb)
    power_dir = zeros(nwind)

    #Reordering rotation directions to a matrix of wind directions
    rotw = zeros((nwind,nturb))
    k=1
    for i=1:nwind
        for j=1:nturb
            rotw[i,j]=gl.rot[k]
            k+=1
        end
    end

    winddir_turb=zeros(nwind)
    for d=1:nwind
        println("Wind Direction #",d)
        #adjusting coordinate system for wind direction
        winddir_turb[d] = 270-gl.windroseDirections[d]
        if winddir_turb[d] < 0
            winddir_turb[d] +=360
        end
        winddir_turb_rad=pi*winddir_turb[d]/180
        xw=x*cos(-winddir_turb_rad) - y*sin(-winddir_turb_rad)
        yw = x*sin(-winddir_turb_rad) + y*cos(-winddir_turb_rad)

        tic()
        #calculate wake velocity components
        wakex,wakey = vawt_wake(xw,yw,gl.dia,rotw[d,:],gl.ntheta,gl.chord,gl.B,gl.Vinf)
        toc()
        #Calculate Power
        for i=1:nturb
            power_turb[i]=vawt_power(i,gl.dia,rotw[d,:],gl.ntheta,gl.chord,gl.H,gl.B,gl.Vinf,gl.twist,gl.delta,gl.rho,wakex,wakey,gl)
        end
        power_dir[d]=sum(power_turb)*gl.windFrequencies[d]
        #Calculating Noise (dB)
        SPL_d=bpm_noise(x,y,gl.windroseDirections[d],rotw[d],wakex,wakey,gl)
        if d==1
            SPL=SPL_d
        else
            append!(SPL,SPL_d)
        end
    end
    power=sum(power_dir)
    percent=power/gl.power_iso_tot
    mSPL=maximum(SPL)
    pit=gl.power_iso_tot
    println("Power: $power W (Isolated: $pit W; $percent)   Max SPL: $mSPL dB")

    #Scale Values for Optimizer
    obj=-1*power/1e3
    SPL=(SPL-gl.SPLlim)/10
    #Calculating separation between turbines
    sep = sep_func(x,y,gl.turb_dia)
    funcs=obj
    fail=false
    #Combine Constants into 1 group
    c=[SPL,sep]
    return funcs,c,fail

end #obj_func

function vawt_wake(xw,yw,dia,rotw,ntheta,chord,B,Vinf)
    t=length(xw) #number of turbines
    wakex=Array{Float64}(1)
    wakey=Array{Float64}(1)
    for i=1:t
        #Copy Variables (to remove element)
        xt=copy(xw)
        yt=copy(yw)
        diat=copy(dia)
        rott=copy(rotw)
        #Remove Stacked Turbines
        deleteat!(xt,i)
        deleteat!(yt,i)
        deleteat!(diat,i)
        deleteat!(rott,i)

        wakexd,wakexy = vawtwake.overlap(ntheta,xt,yt,diat,rott,chord,B,xw[i],yw[i],diat,Vinf,false)
        append!(wakex,wakexd)
        append!(wakey,wakexy)
        if i==1
            deleteat!(wakex,1)
            deleteat!(wakey,1)
        end
    end
    return wakex,wakey
end #vawt_wake

function vawt_power(i,dia,rotw,ntheta,chord,H,B,Vinf,twist,delta,rho,wakext,wakeyt,gl)
    wakex=zeros(gl.ntheta)
    wakey=zeros(gl.ntheta)
    for j=1:ntheta
        wakex[j]=wakext[j+ntheta*(i-1)]
        wakey[j]=wakeyt[j+ntheta*(i-1)]
    end
    if gl.useAC==true
        Omega=rotw[i]
        env=vawtac.Environment(gl.Vinf,gl.rho,gl.mu)
        turbines=Array{vawtac.Turbine}(1)
        turbines[1]=vawtac.Turbine(gl.turb_dia/2,gl.chord,gl.twist,gl.delta,gl.B,gl.af,Omega,0.0,0.0)
        _,Cp,_,_,_,_= vawtac.actuatorcylinder(turbines,env,gl.ntheta)

        power_turb=(0.5*rho*Vinf^3)*(dia[i]*H)*Cp

    else
        power_turb,Cp=vawtwake.powercalc(gl.thetavec,gl.Vinf,wakex,wakey,gl.Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,rotw[i],dia[i]/2.,H,af_data,cl_data,cd_data,twist,rho,interp)
    end
    return power_turb[1]
end #vawt_power

function bpm_noise(turbineX,turbineY,winddir,rot,wakex,wakey,gl)

    nobs,_=size(gl.obs)
    SPL=zeros(nobs)

    noise_corr=1
    nu = 1.78e-5
    c0 = 343.2
    psi = 14.0
    AR = 5.
    rad = gl.turb_dia/2.

    div = 5

    c=ones(div)*gl.chord
    c1=c*0.5
    alpha=ones(div)*0.0
    high=linspace(0,gl.H,div+1)

    for i=1:nobs
        SPL[i]=BPM.turbinepos_VAWT(gl.ntheta,turbineX,turbineY,gl.obs[i,:],winddir,gl.B,gl.Hub,high,rad,c,c1,alpha,nu,c0,psi,AR,noise_corr,gl.rot,gl.Vinf,wakex,wakey)
    end
    return SPL
end #bpm_noise

function sep_func(x,y,turb_dia)

    space=1.65 # turbine diameters apart

    n=length(x)

    sep=zeros((n-1)*n/2)

    k=1
    for i=1:n
        for j=i+1:n
            sep[k]=(x[j]-x[i])^2+(y[j]-y[i])^2
            k+=1
        end
    end
    return sep-(space*turb_dia)^2
end #sep_func

function optimize()
    #run optimization
    optimize=false

    #plot results

    #save Results

    useAC=false
    #Model Setup Parameters
    SPLlim = 100.           # sound pressure level limit of observers
    rotdir_spec = "cn"      # rotation direction (cn- counter-rotating, co- co-rotating)
    ntheta = 72             # number of points around blade flight path
    nRows = 2               # number of paired group rows
    nCols = 2               # number of paired group columns

    windroseDirections=[205.,225.,245.]
    windFrequencies=[0.25,0.5,0.25]
    #One Directional Test
    #windroseDirections=[225]
    #windFrequencies=[1]

    nwind=length(windFrequencies)

    #Define Turbine Specifications
    Vinf = 8.                           # free stream velocity (m/s)
    turb_dia = 1.2                      # turbine diameter (m)
    tsrd = 2.625                        # tip-speed ratio
    turb_rot = tsrd*Vinf/(turb_dia/2.)  # turbine rotation rate (rad/sec)
    twist = 0.0                         # blade twist angle (rad)
    delta = 0.0                         # blade curvature angle (rad)
    B = 3                               # number of turbine blades
    chord = 0.128                       # chord length (m)
    H = 6.1                             # turbine blade height (m)
    Hub = 2.                          # turbine hub height (m)

    #Constants
    rho = 1.225                         # air density (kg/m^3)
    mu = 1.7894e-5                      # fluid viscosity (kg/ms)

    thetavec = zeros(ntheta)
    for i=1:ntheta
        thetavec[i]=(2*pi/ntheta)*(i+1)-(2*pi/ntheta)/2
    end

    #setup initial turbine positions
    # setup initial turbine positions
    grid_start = 2.             # location of starting corner (m)
    pair_sep = 3.               # separation distance between pairs (m)
    group_sep = pair_sep + 10.  # separation distance between paired groups (m)

    #Limits
    xlim=[0,20]
    ylim=[0,20]


    nRC=nRows*nCols
    x01=zeros(nRC)
    y01=zeros(nRC)
    x02=zeros(nRC)
    y02=zeros(nRC)

    x1=linspace(grid_start,grid_start+group_sep*(nCols-1),nCols)
    y1=linspace(grid_start,grid_start+group_sep*(nRows-1),nRows)
    x2=linspace(grid_start+pair_sep,grid_start+pair_sep+group_sep*(nCols-1),nCols)

    k=1
    for i=1:nRows
        for j=1:nCols
            x01[k]=x1[j]
            x01[k] = x1[j]
            y01[k] = y1[i]
            x02[k] = x2[j]
            y02[k] = y1[i]
            # y02[k] = y2[i]
            k += 1
        end
    end
    x0=x01
    y0=y01
    append!(x0,x02)
    append!(y0,y02)
    nturb=length(x0)

    println("x0:",x0)
    println("y0:",y0)

    dia=ones(nturb)*turb_dia

    if rotdir_spec == "co"
        rot1=ones(length(x01))*turb_rot
        rot2=ones(length(x02))*-turb_rot
        rot=rot1
        append!(rot,rot2)
    else
        rot=ones(nturb)*turb_rot
    end
    println("rot:",rot)

    #Duplicate rot for each wind direction
    for i=1:nwind
        append!(rot,rot)
    end

    # specifying boundary locations
    spaceval = 2.
    xlow = 0.
    xupp = max(x2[end],y1[end]) + 2.
    # xupp = max(x2[-1],y2[-1]) + 2.
    ylow = 0.
    yupp = max(x2[end],y1[end]) + 2.
    # yupp = max(x2[-1],y2[-1]) + 2.

    # specifying observer locations
    grid_x = xupp/2.
    grid_y = yupp/2.
    grid_radius = round(Int64,(sqrt((grid_x)^2 + (grid_y)^2))) + 4.

    nobs = 8
    obs_theta = linspace(-pi,pi,nobs+1)
    obs = zeros(nobs,3)
    for i=1:nobs
        obs[i,1] = grid_x + grid_radius*cos(obs_theta[i])
        obs[i,2] = grid_y + grid_radius*sin(obs_theta[i])
        obs[i,3] = 2.
    end
    println(nobs," observers around a radius of ",grid_radius)

    path,_ = splitdir(@__FILE__)
    foildata="$path/../data/airfoils/du06w200.dat" #saved in data/airfoils/
    af=vawtac.readaerodyn(foildata)

    #AC Calculation
    Omega=0.0
    turbines=Array{vawtac.Turbine}(1)
    turbines[1] = vawtac.Turbine(turb_dia/2.,chord,twist,delta,B,af,Omega,0.0,0.0)
    env=vawtac.Environment(Vinf,rho,mu)

    Vpn=0.0
    Vnn=0.0
    Vtp=0.0
    Cpp=0.0
    Cpn=0.0

    _,Cp_iso,_,_,_,_ = vawtac.actuatorcylinder(turbines,env,ntheta)

    power_iso= (0.5*rho*Vinf^3)*(dia[1]*H)*Cp_iso
    power_iso_tot=power_iso*nturb


    useAC=true

    #Setup Global Struct
    globes=globals(turb_dia,dia,rot,chord,twist,delta,B,H,rho,mu,Vinf,windroseDirections,
    windFrequencies,obs,Hub,af,power_iso_tot[1],Vpn,Vnn,Vtp,Cpp,Cpn,ntheta,
    thetavec,useAC,SPLlim)

    #Test OBJ Function
    println("Testing OBJ Function")
    xf=x0
    append!(xf,y0)
    tic()
    a=full_obj(xf,globes)
    toc()

    #Optimize
    #Limits
    #lb=
    #ub=



end #optimize

optimize()

function powercalc(n,wake_x,wake_y,Omega,r,gl)
    n=gl.ntheta
    #local Variables
    W20=zeros(n)
    phi0=zeros(n)
    alpha0=zeros(n)
    cl0=zeros(n)
    cd0=zeros(n)
    ct0=zeros(n)
    Vn=zeros(n)
    Vt=zeros(n)
    W2=zeros(n)
    phi=zeros(n)
    alpha=zeros(n)
    cl=zeros(n)
    cd=zeros(n)
    ct=zeros(n)
    Cpl=zeros(n)

    if (Omega>=0.0)
        for i=1:n
            W20[i]=gl.Vnp[i]^2+gl.Vtp[i]^2
            phi0[i]=atan2(gl.Vnp[i],gl.Vtp[i])
            alpha0[i]=phi0[i]-gl.twist

            ct0[i]=cl0[i]*sin(alpha0[i])-cd0[i]*cos(alpha0[i])

            #Correct Normal/tangential velocities with wake velocities
            Vn[i] = gl.Vnp[i] + wake_x[i]*sin(gl.thetavec[i]) - wake_y[i]*cos(gl.thetavec[i])
            Vt[i] = gl.Vtp[i] + wake_x[i]*cos(gl.thetavec[i]) + wake_y[i]*sin(gl.thetavec[i])
            W2[i] = Vn[i]^2 + Vt[i]^2

            #Compute new inflow angle
            phi[i]=atan2(Vn[i],Vt[i])
            alpha[i]=phi[i]-gl.twist

            #Airfoil
            cl[i], cd[i] = gl.af(alpha[i])
            #Compute new tangential force coefficient
            ct[i]=cl[i]*sin(alpha[i])-cd[i]/ct0[i]

            #Provide relative correction to power coefficient
            Cpl[i]=gl.Cpp[i]*W2[i]/W20[i]*ct[i]/ct0[i]
        end
    else
        for i=1:n
            #calculate baseline values
            W20[i] = Vnn[i]^2 + Vtn[i]^2
            phi0[i] = atan2(Vnn[i],Vtn[i])
            alpha0[i] = phi0[i] - twist
            cl0[i], cd0[i] = gl.af(alpha[i])

            ct0[i] = cl0[i]*sin(alpha0[i]) - cd0[i]*cos(alpha0[i])

            # correct normal/tangential velocities with wake velocities
            Vn[i] = Vnn[i] + wake_x[i]*sin(thetavec[i]) - wake_y[i]*cos(thetavec[i])
            Vt[i] = Vtn[i] - wake_x[i]*cos(thetavec[i]) - wake_y[i]*sin(thetavec[i])
            W2[i] = Vn[i]^2 + Vt[i]^2

            # compute new inflow angle
            phi[i] = atan2(Vn[i],Vt[i])
            alpha[i] = phi[i] - twist

            # airfoil
            cl[i], cd[i] = gl.af(alpha[i])

            # compute new tangential force coefficient
            ct[i] = cl[i]*sin(phi[i]) - cd[i]*cos(phi[i])

            # provide relative correction to power coefficient
            Cpl[i] = Cpn[i]*W2[i]/W20[i]*ct[i]/ct0[i]
        end
    end
    return P,Cp

end #powercalc
