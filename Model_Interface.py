import numpy as np
import matplotlib.pyplot as plt

from VAWT_Wake_Model import velocity_field

# Inputting parameters
x0 = 5.0 # downstream position
y0 = 0.0 # lateral position

velf = 15.0 # free stream wind speed (m/s)
dia = 6.0 # turbine diameter (m)
tsr = 4.0 # turbine tip-speed ratio
B = 3 # number of turbine blades
chord = 0.25 # turbine blade chord length
solidity = (chord*B)/(dia/2.)

vel = velocity_field(x0,y0,velf,dia,tsr,solidity)

print vel