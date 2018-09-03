using PyCall
using JLD
@pyimport matplotlib.pyplot as plt
@pyimport matplotlib as mpl
@pyimport numpy as np

nt=360
g=load("grid2.jld")
grid=g["data"]
println(size(grid))

#Grid Save Test
nturb=8
turb_dia=1.2
gs=40
#grid=zeros(nturb,gs,gs,3)
#tic()
#for t=1:nturb
#    for i=1:gs
#        for j=1:gs
#            grid[t,i,j,1]=i*t+j
#            grid[t,i,j,2]=i-j
#            grid[t,i,j,3]=i-t*j^2
#        end
#    end
#end
#toc()

#println(grid)
t=1

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

#x=[2.0, 8.0, 2.0, 8.0, 4.0, 10.0, 4.0, 10.0, 2.0, 2.0, 8.0, 8.0, 2.0, 2.0, 8.0, 8.0]
x=[1.0, 8.6962, 1.0, 3.09605, 3.09744, 9.91555, 1.09723, 5.41997, 2.67624, 3.41189, 8.48557, 3.43793, 1.0, 1.02541, 5.7729, 4.00314]


nr=360
#gs2=20
xls=np.linspace(1,10,gs)
yls=np.linspace(1,10,gs)
X,Y=np.meshgrid(xls,yls)
#xls2=np.linspace(1,10,gs2)
#yls2=np.linspace(1,10,gs2)
#test=zeros(gs2,gs2)

#X2,Y2=np.meshgrid(xls2,yls2)
for t=1:nturb
    plt.figure()
    #fig, axs = plt.subplots(1,2)
    X,Y=np.meshgrid(xls,yls)
    xf=x[1:nturb]
    yf=x[nturb+1:end]


    #for i=1:gs2
    #    for j=1:gs2
    #        xxs=copy(x)
    #        #Modify turbine x and y
    #        xxs[t]=xls2[i]
    #        xxs[t+nturb]=yls2[j]
    #
    #        #Test Seperation
    #        sep=sep_func(xxs[1:nturb],xxs[nturb+1:end],1.2/1.65)
    #        test[i,j]=maximum(sep)
    #    end
    #end
    lev=np.linspace(3.5,5.1,100)
    CS=plt.contourf(X,Y,-grid[t,:,:,1],levels=lev)
    #pcm=mpl.pcolormesh(x,y,Z,vmin=-6,vmax=-2)
    #plt.contour(X,Y,grid[t,:,:,2],levels=[100],colors="k")
    plt.contour(X,Y,grid[t,:,:,3],colors="r",levels=[0])#,linetype="--")
    #plt.contour(X2,Y2,test,colors="r",levels=[0])
    plt.colorbar(CS)
    #Draw Turbines
    tt=linspace(0,360,nr+1)
    yz=cosd.(tt)*turb_dia/2
    xz=sind.(tt)*turb_dia/2
    for i=1:nturb
        if i==t
            plt.plot(xz+yf[i],yz+xf[i],"r")
        else
            plt.plot(xz+yf[i],yz+xf[i],"k")
        end
    end
    #plt.colorbar()
end
plt.show()
