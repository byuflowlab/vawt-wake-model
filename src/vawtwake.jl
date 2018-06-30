module vawtwake
using SpecialFunctions
using QuadGK
using CSV
using Dierckx
export vfieldx,vfieldy,vorticitystrength,model_gen,velocity_field,overlap,airfoil_data
path,_ = splitdir(@__FILE__)

struct Arguments
        x0::Float64
        y0::Float64
        dia::Float64
        loc1::Float64
        loc2::Float64
        loc3::Float64
        spr1::Float64
        spr2::Float64
        skw1::Float64
        skw2::Float64
        scl1::Float64
        scl2::Float64
        scl3::Float64
        ybound1::Float64
        ybound2::Float64
        rot::Float64
end #Arguments

#Calculating EMG parameter values based on given polynomial surface
function parametereval(tsr,sol,coef)
    a = coef[1]
    b = coef[2]
    c = coef[3]
    d = coef[4]
    e = coef[5]
    f = coef[6]
    g = coef[7]
    h = coef[8]
    i = coef[9]
    j = coef[10]

    val = a + b*tsr + c*sol + d*tsr^2 + e*tsr*sol + f*sol^2 + g*tsr^3 +
    h*tsr^2*sol + i*tsr*sol^2 + j*sol^3
    return val
end #parametereval

function EMGdist(y,loc,spr,skw,scl)
    gam_skew = scl*skw/2.0*exp(skw/2.0*(2.0*loc+skw*spr^2-2.0*y))*(1.0-erf((loc + skw*spr^2.0 - y)/(sqrt(2.0)*spr)))
    return gam_skew
end #EMGdist

function modifyparams(x,y,Args)
        xd=x/Args.dia
        #limiting parameter components to create expected behavior
        if Args.loc1>-0.001 #ensure concave down
                loc1d=-0.001
        else
                loc1d=Args.loc1
        end
        if Args.loc2<0.01       #ensure slight incrase moving downstream
                loc2d=0.01
        else
                loc2d=Args.loc2
        end
        if Args.loc3<0.48       #ensure wake originating from edge of turbine
                loc3d=0.48
        else
                loc3d=Args.loc3
        end
        #EMG Location
        loc=loc1d*xd*xd+loc2d*xd+loc3d

        if Args.spr1 >-0.001
                spr1d=-0.001
        else
                spr1d=Args.spr1
        end
        if Args.spr2>0.0
                spr2d=0.0
        else
                spr2d=Args.spr2
        end

        #EMG spread
        spr=spr1d*xd+spr2d

        skw1d = Args.skw1       #no limitataions necessary
        if Args.skw2>0.0        #Ensure value does not begin positive
                skw2d=0.0
        else
                skw2d=Args.skw2
        end
        #EMG Skew
        skw=skw1d*xd+skw2d

        if Args.scl1<0.0        #ensure positive maximum vorticity strength
                scl1d=0.0
        else
                scl1d=Args.scl1
        end
        if Args.scl2<0.05       #ensure decay moving downstream
                scl2d=0.05
        else
                scl2d=Args.scl2
        end
        if Args.scl3<0.0        #ensure decay occurs downstream
                scl3d=0.0
        else
                scl3d = Args.scl3
        end
        #EMG Scale
        scl=scl1d/(1.0+exp(scl2d*(xd-scl3d)))

        #Limiting parameters to the maximum values the EMG distributino can handle
        if loc<0.2
                loc=0.2
        end
        if spr<-0.5
                spr=-0.5
        elseif (spr>-0.001)
                spr=-0.001
        end
        if skw>0.0
                skw=0.0
        end
        return loc,spr,skw,scl
end #modifyparams

function vorticitystrength(x,y,Args)
        yd=y/Args.dia
        loc,spr,skw,scl=modifyparams(x,y,Args)
        g1=EMGdist(yd,loc,spr,skw,scl)
        g2=EMGdist(yd,-loc,-spr,-skw,-scl)
        gam_lat=(g1-g2)
        return gam_lat
end #vorticitystrength

function vorticitystrengthx(x,y,Args)
        yd=y/Args.dia
        loc,spr,skw,scl=modifyparams(x,y,Args)
        kpp = skw*spr*spr
        #Exponentially Modified Gaussian Distribution
        gam_lat =  (1.0/(2.0*spr))*scl*skw*(exp(-(loc-yd)*(loc-yd)/(2.0*spr*spr))*
        sqrt(2.0/pi) + exp(-(loc+yd)*(loc+yd)/(2.0*spr*spr))*sqrt(2.0/pi) +
        exp(0.5*skw*(2.0*loc + kpp - 2.0*y)*skw*spr*(-1.0 +
        erf((loc + kpp - yd)/(sqrt(2.0)*spr)))) + exp(0.5*skw*(2.0*
        loc + kpp + 2.0*y)*skw*spr*(-1.0 + erf((loc + kpp +
        yd)/(sqrt(2.0)*spr)))))
        return gam_lat
end #vorticitystrengthx

function vorticitystrengthy(x,y,Args)
    #Normalize x and y by diameter
    xd=x/Args.dia
    yd=y/Args.dia
    loc,spr,skw,scl=modifyparams(x,y,Args)
    g1=EMGdist(yd,loc,spr,skw,scl)
    g2=EMGdist(yd,-loc,-spr,-skw,-scl)
    gam_lat=(g1-g2)
    return gam_lat
end #vorticitystrengthy

function integrandx(x,y,Args)
        gammav=vorticitystrength(x,y,Args)
        num=(y-Args.y0)
        den=((x-Args.x0)*(x-Args.x0))+((y-Args.y0)*(y-Args.y0))
        inte=gammav*num/den
        return inte
end #integrandx

function integrandy(x,y,Args)
        gammav=vorticitystrength(x,y,Args)
        num=Args.x0-x
        den=(x-Args.x0)*(x-Args.x0)+(y-Args.y0)*(y-Args.y0)
        inte=gammav*num/den
        return inte
end #integrandy

function fxy(x,y,Args)
        results=integrandx(x,y,Args)
        return results
end #fxy

function fyx(x,y,Args)
        results=integrandy(x,y,Args)
        return results
end #fyx

function fx(x,Args)
        #Integration Tolerances
        epsabs=1.49e-8
        epsrel=1.49e-8
        #Create Inner Function
        funcypass(y)=fxy(x,y,Args)
        #Inner Integral
        result,tol=quadgk(funcypass,Args.ybound1,Args.ybound2;reltol=epsrel,abstol=epsabs)
        return result
end #fx

function fy(x,Args)
        #Integration Tolerances
        epsabs=1.49e-8
        epsrel=1.49e-8
        #Create Inner Function
        funcxpass(y)=fyx(x,y,Args)
        #Inner Integral
        result,tol=quadgk(funcxpass,Args.ybound1,Args.ybound2;reltol=epsrel,abstol=epsabs)
        return result

end #fy

function vfieldx(Args)
        #Integration Properties
        epsabs=1.49e-8
        epsrel=1.49e-8
        #limits of intigration
        xbounds=(0,(Args.scl3+5.0)*Args.dia)
        #Integration
        funcxpass(x)=fx(x,Args)
        #Outer Integral
        result,tol=quadgk(funcxpass,xbounds[1],xbounds[2];reltol=epsrel,abstol=epsabs)
        return result
end #vfieldx

function vfieldy(Args)
        #Integration Properties
        epsabs=1.49e-8
        epsrel=1.49e-8
        #limits of intigration
        xbounds=(0,(Args.scl1+5.0)*Args.dia)
        #Integration
        funcypass(x)=fy(x,Args)
        #Outer Integral
        result,tol=quadgk(funcypass,xbounds[1],xbounds[2];reltol=epsrel,abstol=epsabs)
        return result
end #vfieldy

function model_gen(x0,y0,tsr,solidity,dia,rot)
    #Get file locatoin (absolute)
    #path,_ = splitdir(@__FILE__)
    fileloc="$path/../data/VAWTPolySurfaceCoef_pub.csv" #Not sure if this is the best way
    csvdata=CSV.read(fileloc,delim=',')
    #Project Point Onto Surface
    loc1=_parameterval(tsr,solidity,csvdata[:,1])
    loc2=_parameterval(tsr,solidity,csvdata[:,2])
    loc3=_parameterval(tsr,solidity,csvdata[:,3])
    spr1=_parameterval(tsr,solidity,csvdata[:,4])
    spr2=_parameterval(tsr,solidity,csvdata[:,5])
    skw1=_parameterval(tsr,solidity,csvdata[:,6])
    skw2=_parameterval(tsr,solidity,csvdata[:,7])
    scl1=_parameterval(tsr,solidity,csvdata[:,8])
    scl2=_parameterval(tsr,solidity,csvdata[:,9])
    scl3=_parameterval(tsr,solidity,csvdata[:,10])
    #Prepare Model Arguments Struct
    xbound = (scl3+5.0)*dia
    Args=Arguments(x0,y0,dia,loc1,loc2,loc3,spr1,spr2,skw1,skw2,scl1,scl2,scl3,-dia,dia,rot)
    return Args
end #model_gen

function velocity_field(xt,yt,x0,y0,Vinf,dia,rot,chord,B,veltype,param="None")
        rad=dia/2.0
        tsr=rad*abs(rot)/Vinf
        solidity = (chord*B)/rad
        x0t=x0-xt
        x0y=y0-yt

    if param=="None"
            Args=model_gen(x0,y0,tsr,solidity,dia,rot)
    else
            Args=param
    end


    if veltype=="vort"
        #Vorticity Caclculation (No Integration)
        if x0t < 0
            vel=0.0
        else
            vel=vorticitystrength(x0t,y0t,Args)
        end
    else
        if veltype== "all" || veltype=="x" || veltype=="ind"
            vel_x = vfieldx(Args)
            vel_xs=(vel_x[1]*abs(Args.rot))/(2*pi)
        end
        if veltype== "all" || veltype=="y" || veltype=="ind"
            vel_y = vfieldy(Args)
            vel_ys=(vel_y[1]*abs(Args.rot))/(2*pi)
        end
        if veltype=="all"
            vel=sqrt((vel_xs+Vinf)^2+vel_ys^2)/Vinf
        elseif veltype=="x"
            vel=(vel_xs+Vinf)/Vinf
        elseif veltype=="y"
            vel=vel_ys/Vinf
        elseif veltype=="ind"
            vel=[vel_xs,vel_ys]/Vinf
        else
            throw(DomainError)
            println("Error in veltype Option")
        end
    end
    return vel
end #velocity_field

function _parameterval(tsr,sol,coef)
    #calculate coefficent based on given tsr and sol on surface
    a = coef[1]
    b = coef[2]
    c = coef[3]
    d = coef[4]
    e = coef[5]
    f = coef[6]
    g = coef[7]
    h = coef[8]
    i = coef[9]
    j = coef[10]

    surf = a + b*tsr + c*sol + d*tsr^2 + e*tsr*sol + f*sol^2 + g*tsr^3 + h*tsr^2*sol + i*tsr*sol^2 + j*sol^3
    return surf
end #_parameterval

function pvawt(p,xt,yt,xd,yd,Vinf,diat,rott,chord,B,veltype,param)
        a=zeros(p,2)
        for k=1:p
                wake=vawtwake.velocity_field(xt,yt,xd[k],yd[k],Vinf,diat,rott,chord,B,veltype,param)
                a[k,1]=wake[1]
                a[k,2]=wake[2]
        end
        return a
end

function overlap(p,xt,yt,diat,rott,chord,B,x0,y0,dia,Vinf,pointcalc,param="None",veltype="ind")
        np=nprocs()
        t=length(xt)
        xd=zeros(p)
        yd=zeros(p)
        velx=zeros(p)
        vely=zeros(p)
        velx_int=zeros(p)
        vely_int=zeros(p)
        #Find points around the flight path of the blades

        for i=1:p
                if pointcalc==false
                        theta=(2*pi/p)*i-(2*pi/p)/2
                        xd[i]=x0-sin(theta)*(dia[1]/2.0)
                        yd[i]=y0+cos(theta)*(dia[1]/2.0)
                elseif pointcalc==true
                        xd[1]=0
                        yd[1]=0
                end
        end
        intex=zeros(p)
        intey=zeros(p)

        if t==1 #Coupled Confiugration
                if pointcalc==false
                        for j=1:p
                                wake=velocity_field(xt[1],yt[1],xd[j],yd[j],Vinf,diat[1],chord,B,veltype,param)
                                velx[j]=wake[1]*Vinf
                                vely[j]=wake[2]*Vinf
                        end
                elseif pointcalc==true
                        wake = velocity_field(xt[1],yt[1],xd[1],yd[1],Vinf,diat[1],rott[1],chord,B,veltype,param)
                        velx[1] = wake[1]*Vinf
                        vely[1] = wake[2]*Vinf
                end
        else
                if pointcalc==false
                        if np==1
                                for j=1:t
                                        for k=1:p
                                                wake=velocity_field(xt[j],yt[j],xd[k],yd[k],Vinf,diat[j],rott[j],chord,B,veltype,param)
                                                velx_int[k] = -wake[1]
                                                vely_int[k] = wake[2]
                                                if (velx_int[k] >= 0.0)
                                                        intex[k] = intex[k] + (velx_int[k])^2
                                                else
                                                        intex[k] = intex[k] - (velx_int[k])^2
                                                end
                                                if (vely_int[k] >= 0.0)
                                                        intey[k] = intey[k] + (vely_int[k])^2
                                                else
                                                        intey[k] = intey[k] - (vely_int[k])^2
                                                end
                                        end
                                end
                        elseif np>1
                                wake=zeros(t,p,2)
                                #Do the Calculation First
                                wake=@parallel vcat for j=1:t
                                    pvawt(p,xt[j],yt[j],xd,yd,Vinf,diat[j],rott[j],chord,B,veltype,param)
                                end
                                #Then Perform numbers
                                for j=1:t
                                        for k=1:p
                                                velx_int[1] = -wake[(j-1)*p+k,1]
                                                vely_int[1] = wake[(j-1)*p+k,2]

                                                # sum of squares of velocity deficits
                                                if (velx_int[1] >= 0.0)
                                                        intex[1] = intex[1] + (velx_int[1])^2
                                                else
                                                        intex[1] = intex[1] - (velx_int[1])^2
                                                end
                                                if (vely_int[1] >= 0.0)
                                                        intey[1] = intey[1] + (vely_int[1])^2
                                                else
                                                        intey[1] = intey[1] - (vely_int[1])^2
                                                end
                                        end
                                end

                        end
        elseif pointcalc==true
                        for j=1:t
                                wake = velocity_field(xt[j],yt[j],xd[1],yd[1],Vinf,diat[j],rott[j],chord,B,veltype,param)
                                velx_int[1] = -wake[1]
                                vely_int[1] = wake[2]

                                # sum of squares of velocity deficits
                                if (velx_int[1] >= 0.0)
                                        intex[1] = intex[1] + (velx_int[1])^2
                                else
                                        intex[1] = intex[1] - (velx_int[1])^2
                                end
                                if (vely_int[1] >= 0.0)
                                        intey[1] = intey[1] + (vely_int[1])^2
                                else
                                        intey[1] = intey[1] - (vely_int[1])^2
                                end
                        end
                end
                #square root of sum squares
                for l=1:p
                        if (intex[l] >= 0.0)
                                velx[l] = -Vinf*(sqrt(intex[l]))
                        else
                                velx[l] = Vinf*(sqrt(abs(intex[l])))
                        end
                        if (intey[l] >= 0.0)
                                vely[l] = Vinf*(sqrt(intey[l]))
                        else
                                vely[l] = -Vinf*(sqrt(abs(intey[l])))
                        end
                end
        end
        return velx,vely
end #overlap

function airfoil_data(file)
        fileloc="$path/../data/airfoils/$file"
        af_data=Array{Float64}(0)
        cl_data=Array{Float64}(0)
        cd_data=Array{Float64}(0)
        open(fileloc) do f
                for i in enumerate(eachline(f))
                        count=i[1]
                        str=i[2]
                        if count>13
                                line=split(str)
                                if line[1]=="EOT"

                                else
                                        af_data=append!(af_data,float(line[1]))
                                        cl_data=append!(cl_data,float(line[2]))
                                        cd_data=append!(cd_data,float(line[3]))
                                end
                        end
                end
        end
        return af_data,cl_data,cd_data
        #Smoothing data with a spline


        clsmooth=Dierckx.Spline1D(af_data,cl_data,s=0.1)
        cdsmooth=Dierckx.Spline1D(af_data,cd_data,s=0.001)

        #Reassigning imported airfoild data with the smoothed airfoil data
        af_data_smooth = linspace(af_data[0],af_data[end],4000)
        cl_data_smooth = evaluate(clsmooth,af_data_smooth)
        cd_data_smooth = evaluate(cdsmooth,af_data_smooth)

        #Return Smoothed data
        return af_data_smooth,cl_data_smooth,cd_data_smooth
end #airfoil data


end #module
