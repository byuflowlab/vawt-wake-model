# VAWT Wake Model

Parameterized VAWT Wake Model using CFD vorticity data

Developed by Eric Tingey at Brigham Young University, 2015-2017

This code models the wake behind a vertical-axis wind turbine based on parameters like tip-speed ratio, solidity and wind speed by using databases of curve fit coefficients to produce velocity information. The model uses CFD data obtained from STAR-CCM+ of simulated turbines to make the wake model as accurate as possible.

Only valid for tip-speed ratios between 1.5 and 7.0 and solidities between 0.15 and 1.0. Reynolds numbers should also be around the range of 600,000 to 6,000,000.


## Installation instructions

- system requirements: gfortran (using MinGW for Windows in order to use the commands here), python 2.7, numpy, scipy, matplotlib, joblib
- run:
```
python setup.py install
```
- or navigate to the directory and run the following command in the terminal to build the Fortran code:

Mac
```
$ cd wake_model
$ f2py -c  --opt=-O2 -m _vawtwake VAWT_Wake_Model.f90
```

Windows
```
cd wake_model
python <\your\path\to\f2py.py> -c --opt=-O2 --compiler=mingw32 --fcompiler=gfortran -m _vawtwake VAWT_Wake_Model.f90
```
(<\your\path\to\f2py.py>: most likely C:\Python27\Scripts\f2py.py)
MinGW can be installed using the download installer from http://mingw.org/.
Using the installer, use the basic setup and install msys-base, mingw32-base, mingw32-gcc-g++, mingw32-gcc-objc, and mingw-developer-tools. Ensure that gcc.exe is installed in C:\MinGW\bin\. Set the path variable of 'C:\MinGW\bin\' so it can be recognized by the computer. An example of this process can be viewed at https://www.youtube.com/watch?v=DHekr3EtDOA.

- a unit test can then be run using the command:
```
$ cd tests
$ python test.py
```

## Running the Python code

This python code can be run from another file using:
```python
from VAWT_Wake_Model import velocity_field
velocity_field(x0,y0,velf,dia,rot,chord,B,param=None,veltype='all',integration='simp',m=220,n=200)
```

An example code is available to see how to call the wake model code and calculate a normalized velocity at a given location. Plotting of a velocity profile at a specific downstream distance as well as plotting the entire flow domain is also demonstrated in the example.

The complete data set of the wake vorticity calculations used to produce this model is available to access at:
https://figshare.com/articles/Parameterized_Wake_Model/2059947

## Referencing the Model

To reference any part of this repository, use the reference of the paper entitled "Parameterized Vertical-Axis Wind Turbine Wake Model Using CFD Vorticity Data" found at http://arc.aiaa.org/doi/10.2514/6.2016-1730.
