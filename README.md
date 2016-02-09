# VAWT Wake Model

## Installation instructions

- system requirements: gfortran (using MinGW for Windows in order to use the commands here), python 2.7, numpy, scipy
- navigate to the directory and run the following command in the terminal to build the Fortran code:

Mac
```
$ cd wake_model
$ f2py -c  --opt=-O2 -m _vortmodel vorticity.f90
```

Windows
```
cd wake_model
python <\your\path\to\f2py.py> -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _vortmodel vorticity.f90
```
(most likely C:\Python27\Scripts\f2py.py)

- a unit test can then be run using the command:
```
$ cd tests
$ python test.py
```

## Running the Python code

This python code can be run from another file using:
```python
from VAWT_Wake_Model import velocity_field
velocity_field(x0,y0,velf,dia,tsr,solidity)  # velocity calculation at any point (x0,y0) for a given free stream wind speed, turbine diameter, tip-speed ratio, and solidity
```

An example code is available to see how to call the wake model code and calculate a normalized velocity at a given location. Plotting of a velocity profile at a specific downstream distance as well as plotting the entire flow domain is also demonstrated in the example.

The complete data set of the wake vorticity calculations used to produce this model is available to access at:
https://figshare.com/articles/Parameterized_Wake_Model/2059947

## Referencing the Model

To reference any part of this repository, use the reference of the paper entitled "Parameterized Vertical-Axis Wind Turbine Wake Model Using CFD Vorticity Data" found at http://arc.aiaa.org/doi/10.2514/6.2016-1730.
