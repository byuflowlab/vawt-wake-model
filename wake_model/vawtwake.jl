module vawtwake
using SpecialFunctions
using QuadGK
using CSV
export vfieldx,vfieldy,vorticitystrength,model_gen,velocity_field

struct Arguments
        xt
        yt
        dia
        loc1
        loc2
        loc3
        spr1
        spr2
        skw1
        skw2
        scl1
        scl2
        scl3
        ybound1
        ybound2
        rot
end

function interpolate(n,x,y,xval)
    for i=1:n
        if (xval <x[i])
            yval = y[i-1]+(xval - x[i-1])*((y[i]-y[i-1])/(x[i]-x[i-1]))
            break
        elseif (xval==x[i])
            yval = y[i]
        end
    end
    return yval
end

#Cubic Splin interpolation setup
function splineint(n,x,y,xval)
    for i=1:n
        if (xval<x[i])
            if (i==2)
                x1=x[1]
                x2=x[2]
                x3=x[3]
                y1=y[1]
                y2=y[2]
                y3=y[3]
                yval = cubspline(x1,x2,x3,y1,y2,y3,xval)
            elseif (i==n)
                x1=x[n-2]
                x2=x[n-1]
                x3=x[n]
                y1=y[n-2]
                y2=y[n-1]
                y3=y[n]
                yval=cubspline(x1,x2,x3,y1,y2,y3,xval)
            else
                if (xval <= (x[i]+x[i-1])/2.0)
                    x1=x[i-2]
                    x2=x[i-1]
                    x3=x[i]
                    y1=y[i-2]
                    y2=y[i-1]
                    y3=y[i]
                    yval=cubspline(x1,x2,x3,y1,y2,y3,xval)
                else
                    x1=x[i-1]
                    x2=x[i]
                    x3=x[i+1]
                    y1=y[i-1]
                    y2=y[i]
                    y3=y[i+1]
                    yval=cubsline(x1,x2,x3,y1,y2,y3,xval)
                end
            end
            return yval
        elseif (xval==x[i])
            yval=y[i]
            return yval
        end
    end
end #splineint

#cubic spline interpolation (for airfoil data)

function cubspline(x1,x2,x3,y1,y2,y3,xval)
    #solve tridiagonal linear equation system
    a11=2.0/(x2-x1)
    a12=1.0/(x2-x1)
    a13=0.0
    a21=1.0/(x2-x1)
    a22=2.0*((1.0/(x2-x1))+(1.0/(x3-x2)))
    a23=1.0/(x3-x2)
    a31=0.0
    a32=1.0/(x3-x2)
    a33=2.0/(x3-x2)
    b1=3.0*(y2-y1)/(x20x1)^2
    b2=3.0*(((y2-y1)/(x2-x1)^2)+((y3-y2)/(x3-x2)^2))
    b3=3.0*(y3-y2)/(x3-x2)^2

    #solve using inverse matrix method
    bot = a11*a22*a33 + a12*a23*a31 + a13*a21*a32 - a13*a22*a31 - a12*a21*a33 - a11*a23*a32
    if (xval < x2)
      xtop = b1*a22*a33 + a12*a23*b3 + a13*b2*a32 - a13*a22*b3 - a12*b2*a33 - b1*a23*a32
      ytop = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - a13*b2*a31 - b1*a21*a33 - a11*a23*b3

      k1 = xtop/bot
      k2 = ytop/bot

      a = k1*(x2-x1) - (y2-y1)
      b = -k2*(x2-x1) + (y2-y1)
      t = (xval-x1)/(x2-x1)

      yval = (1.0 - t)*y1 + t*y2 + t*(1.0 - t)*(a*(1.0 - t) + b*t)
    else
      ytop = a11*b2*a33 + b1*a23*a31 + a13*a21*b3 - a13*b2*a31 - b1*a21*a33 - a11*a23*b3
      ztop = a11*a22*b3 + a12*b2*a31 + b1*a21*a32 - b1*a22*a31 - a12*a21*b3 - a11*b2*a32

      k2 = ytop/bot
      k3 = ztop/bot

      a = k2*(x3-x2) - (y3-y2)
      b = -k3*(x3-x2) + (y3-y2)
      t = (xval-x2)/(x3-x2)

      yval = (1.0 - t)*y2 + t*y3 + t*(1.0 - t)*(a*(1.0 - t) + b*t)
    end

    return yval
end #cubspline

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
    gam_skew = scl*skw/2.0*exp(skw/2.0*(2.0*loc+skw*spr^2-2.0*y))*(1.0-erf((loc + skw*spr^2.0 - y)/(sqrt(2.0*spr))))
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
        num=(y-Args.yt)
        den=((x-Args.xt)*(x-Args.xt))+((y-Args.yt)*(y-Args.yt))
        inte=gammav*num/den
        return inte
end #integrandx

function integrandy(x,y,Args)
        gammav=vorticitystrength(x,y,Args)
        num=x0-Args.xt
        den=(x-Args.xt)*(x-Args.xt)+(y-Args.yt)*(y-Args.yt)
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
        xbounds=(0,(Args.scl1+5.0)*Args.dia)
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
    fileloc=pwd()*"/data/VAWTPolySurfaceCoef_pub.csv" #Not sure if this is the best way
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

function velocity_field(x0,y0,Args,veltype)

    x0t=Args.xt-x0
    y0t=Args.yt-y0

    #Eric Includes 2 options here
    #Caclaulate Surface coefficents
    #Load Surface coefficents

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
            println("Error in veltype Option")
        end
        return vel
    end

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

end #module
