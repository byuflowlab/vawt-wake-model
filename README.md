# VAWT Wake Model

## Installation instructions

- system requirements: gfortran, python 2.7, numpy, scipy
- navigate to the directory and run the following command in the terminal: 
```
$ f2py -c  --opt=-O2 -m _vortrun vorticity.f90
```   

## Running the Python code

This python code can be run from another file using:
```python
from VAWT_Wake_Model import velocity_field
velocity_field(x0,y0,velf,dia,tsr,solidity)  # velocity calculation at any point (x0,y0) for a given free stream wind speed, turbine diameter, tip-speed ratio, and solidity
``` 

An example code is available to see how to call the wake model code and calculate a noramlized velocity at a given location. Plotting of a velocity profile at a specific downstream distance as well as plotting the entire flow domain is also demonstrated in the example.

The complete data set of the wake vorticity calculations used to produce this model is available to access at:
https://figshare.com/articles/Parameterized_Wake_Model/2059947
