#Pkg.build("Snopt")
using Snopt
using vawtwake
#@everywhere include("vawtwake.jl")
using BPM
using PyCall
using vawtac
@pyimport matplotlib.pyplot as plt
using JLD


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
        tic()
        println("Wind Direction #",d)
        #adjusting coordinate system for wind direction
        winddir_turb[d] = 270-gl.windroseDirections[d]
        if winddir_turb[d] < 0
            winddir_turb[d] +=360
        end
        winddir_turb_rad=pi*winddir_turb[d]/180
        xw=x*cos(-winddir_turb_rad) - y*sin(-winddir_turb_rad)
        yw = x*sin(-winddir_turb_rad) + y*cos(-winddir_turb_rad)

        #calculate wake velocity components
        tic()
        wakex,wakey = vawt_wake(xw,yw,gl.dia,rotw[d,:],gl.ntheta,gl.chord,gl.B,gl.Vinf)
        toc()
        tic()
        #Calculate Power
        power_turb=vawt_power(0.0,gl.dia,rotw[d,:],gl.ntheta,gl.chord,gl.H,gl.B,gl.Vinf,gl.twist,gl.delta,gl.rho,wakex,wakey,gl,xw,yw)
        toc()
        #for i=1:nturb
        #    power_turb[i]=vawt_power(i,gl.dia,rotw[d,:],gl.ntheta,gl.chord,gl.H,gl.B,gl.Vinf,gl.twist,gl.delta,gl.rho,wakex,wakey,gl)
        #end
        power_dir[d]=sum(power_turb)*gl.windFrequencies[d]
        #Calculating Noise (dB)
        SPL_d=bpm_noise(x,y,gl.windroseDirections[d],rotw[d],wakex,wakey,gl)
        if d==1
            SPL=SPL_d
        else
            append!(SPL,SPL_d)
        end
        toc()
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
    #c=[SPL,sep]
    c=SPL
    append!(c,sep)
    cbounds=circularBoundary(x,y,5)
    append!(c,cbounds)
    #println(funcs,c,fail)
    return funcs,c,fail

end #obj_func

function circularBoundary(turbineX,turbineY,circle_radius)
    nTurb=length(turbineX)
    cbounds=zeros(nTurb)
    for i=1:nTurb
        R=sqrt((turbineX[i]-0)^2+(turbineY[i]-0)^2)
        cbounds[i]=circle_radius-R
    end
    return cbounds
end #circularBoundary

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

function vawt_power(k,dia,rotw,ntheta,chord,H,B,Vinf,twist,delta,rho,wakext,wakeyt,gl,xw,yw)
    wakex=zeros(gl.ntheta)
    wakey=zeros(gl.ntheta)
    nturb=length(rotw)
    power_turb=zeros(nturb)
    #for j=1:ntheta
    #    wakex[j]=wakext[j+ntheta*(i-1)]
    #    wakey[j]=wakeyt[j+ntheta*(i-1)]
    #end
    if gl.useAC==true
        env=vawtac.Environment(Vinf,gl.rho,gl.mu)
        turbines=Array{vawtac.Turbine}(nturb)
        for i=1:nturb
            turbines[i]=vawtac.Turbine(gl.dia[i]/2,gl.chord,gl.twist,gl.delta,gl.B,gl.af,rotw[i],xw[i],yw[i])
        end
        _,Cp,_,_,_,_= vawtac.actuatorcylinder(turbines,env,gl.ntheta)
        for i=1:nturb
            power_turb[i]=(0.5*rho*Vinf^3)*(gl.dia[i]*H)*Cp[i]
        end
        #power_turb=
    else
        power_turb,Cp=vawtwake.powercalc(gl.thetavec,gl.Vinf,wakex,wakey,gl.Vnp,Vnn,Vtp,Vtn,Cpp,Cpn,rotw[i],dia[i]/2.,H,af_data,cl_data,cd_data,twist,rho,interp)
    end
    return power_turb
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
    return (space*turb_dia)^2-sep
end #sep_func

function dglobes()
    useAC=false
    #Model Setup Parameters
    SPLlim = 100.           # sound pressure level limit of observers
    rotdir_spec = "cn"      # rotation direction (cn- counter-rotating, co- co-rotating)
    ntheta = 72             # number of points around blade flight path
    nRows = 2               # number of paired group rows
    nCols = 2               # number of paired group columns

    #windroseDirections=[205.,225.,245.]
    #windFrequencies=[0.25,0.5,0.25]
    #One Directional Test
    windroseDirections=[225]
    windFrequencies=[1]

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
    pair_sep = 2.               # separation distance between pairs (m)
    group_sep = pair_sep + 4.  # separation distance between paired groups (m)

    #Limits
    xlim=[1,10]
    ylim=[1,10]


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
    Omega=turb_rot
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


    xf=x0
    append!(xf,y0)
    return globes,xf
end #dglobes

function optimize(opt)
    #run optimization
    optimize=opt

    #plot results
    xlim=[-5,5]
    ylim=[-5,5]
    #save Results

    #Load Model
    globes,xf=dglobes()

    if optimize==false
        #Test OBJ Function
        println("Testing OBJ Function")
        tic()
        xopt,c,info=full_obj(xf,globes)
        toc()
        println(c)
        xopt=xf
        fopt=0.0
        info="Skipped Function"
    else
    #Warm Start
    xf=[1.0, -3.6962, 1.0, 3.09605, 4.5, 2.91555, 1.09723, -3.41997, 2.67624, 3.41189, -4.48557, 3.43793, 1.5, 1.02541,-2.7729, 4.00314]
    nturb=8
    #Optimize
    sz=ones(nturb)
    #Convert Limits to Snopt format
    lb=sz*xlim[1]
    append!(lb,sz*ylim[1])
    ub=sz*xlim[2]
    append!(ub,sz*ylim[2])

    options = Dict{String, Any}()
    options["Derivative option"] = 0

    #Define Objective Function (no parameters)
    sobj(xx)=full_obj(xx,globes)

    xopt,fopt,info=snopt(sobj,xf,lb,ub,options)

    println(xopt)
    println(fopt)
    println(info)
    end
    return xopt,globes
end #optimize

function contour(globes,x)
    #Grid Resolution
    gs=20

    nturb=floor(Int,length(x)/2)
    #Plot Boundary
    bxlim=[1,10]
    bylim=[1,10]

    xls=linspace(bxlim[1],bxlim[2],gs)
    yls=linspace(bylim[1],bylim[2],gs)

    grid=zeros(nturb,gs,gs,5)
    xss=x

    for t=1:nturb
        for i=1:gs
            for j=1:gs
                println([t,i,j])
                xxs=copy(x)
                #Modify turbine x and y
                xxs[t]=xls[i]
                xxs[t+nturb]=yls[j]

                #Test Seperation
                sep=sep_func(xxs[1:nturb],xxs[nturb+1:end],globes.turb_dia)
                #println(maximum(sep))
                #println(xxs)
                if maximum(sep)>(0.65*1.2)
                    tf=0
                    tc=zeros(nturb+1)
                    append!(tc,sep)
                else
                    tf,tc,~=full_obj(xxs,globes)
                end
                grid[t,i,j,1]=tf
                grid[t,i,j,2]=maximum(tc[1:nturb])
                grid[t,i,j,3]=maximum(tc[nturb+1:end])
                grid[t,i,j,4]=xls[i]
                grid[t,i,j,5]=yls[j]
            end
        end
    end
    println(grid)
    save("grid2.jld", "data", grid)


end

x,globes=optimize(false)
#x=[1.0, 8.6962, 1.0, 3.09605, 3.09744, 9.91555, 1.09723, 5.41997, 2.67624, 3.41189, 8.48557, 3.43793, 1.0, 1.02541, 5.7729, 4.00314]

println(x)
#contour(globes,x)
