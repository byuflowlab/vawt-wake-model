VAWT Wake Model
-----------------------------
Installation instructions for MAC
-----------------------------
- system requirements: gfortran, python 2.7, numpy, scipy
- put all files in desired directory
- navigate to the directory and run the following command in the terminal: $ f2py -c  --opt=-O2 -m _vortrun vorticity.f90
    
Running the Python code
-----------------------------
This python code can be run from another file using:
  - from VAWT_Wake_Model import velocity_field
  - velocity_field(x0,y0,velf,dia,tsr,solidity) --> velocity calculation at any point (x0,y0) for a given free stream wind speed, turbine diameter, tip-speed ratio, and solidity

It can also be run from the file itself under the main file section. Specify the parameters of the VAWT according to a specific case desired. Plotting of a velocity profile at a specific downstream distance as well as plotting the entire flow domain is available to use.
