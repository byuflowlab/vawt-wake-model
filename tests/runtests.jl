using vawtwake
using Base.Test

@testset "Generic VAWT Nearfield" begin
##################################################################################
##################################################################################
##################################################################################

#Test Farm Validation
Vinf=15.0               #free stream wind speed (m/s)
dia=6.0                 #turbine diameter (m)
tsr=4.                  #tip speed ratio ((dia/2)*rot/Vinf)
B=3.                    #number of blades
chord=0.25              #chord length (m)
rot=tsr*Vinf/(dia/2.)   #rotation rate (rad/s)

rad = dia/2.
tsr=rad*abs(rot)/Vinf
solidity = (chord*B)/rad

#Turbine and velocity calulation locations
xt = 0.         # downstream position of turbine in flow domain (m)
yt = 0.         # lateral position of turbine in flow domain (m)
x0 = 12.        # downstream distance for velocity calculation (m)
y0 = 0.         # lateral distance for velocity calculation (m)

veltype = "all"

#Create Model
res=vawtwake.velocity_field(xt,yt,x0,y0,Vinf,dia,rot,chord,B,veltype)
@test isapprox(res,0.284053,atol=1e-6)
end #Nearfield test

@testset "Generic VAWT Farfield" begin
##################################################################################
##################################################################################
##################################################################################

#Test Farm Validation
Vinf=15.0               #free stream wind speed (m/s)
dia=6.0                 #turbine diameter (m)
tsr=4.                  #tip speed ratio ((dia/2)*rot/Vinf)
B=3.                    #number of blades
chord=0.25              #chord length (m)
rot=tsr*Vinf/(dia/2.)   #rotation rate (rad/s)

rad = dia/2.
tsr=rad*abs(rot)/Vinf
solidity = (chord*B)/rad

#Turbine and velocity calulation locations
xt = 0.         # downstream position of turbine in flow domain (m)
yt = 0.         # lateral position of turbine in flow domain (m)
x0 = 50.        # downstream distance for velocity calculation (m)
y0 = 20.         # lateral distance for velocity calculation (m)

veltype = "all"

#Create Model
res=vawtwake.velocity_field(xt,yt,x0,y0,Vinf,dia,rot,chord,B,veltype)
@test isapprox(res,1.04282,atol=1e-6)
end #farfield test
