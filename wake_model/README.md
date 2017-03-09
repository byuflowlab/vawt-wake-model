This folder contains the following files:

VAWT_Wake_Model.py: the VAWT wake model that uses the trend fitting from the CFD wake data

Integrate.py: file that uses the Scipy.integrate package and calls only necessary subroutines

Velocity_Example.py: an example code that calls VAWT_Wake_Model.py and returns velocity plots and results

ACsingle.py: the actuator cylinder code that converges the self-induced velocities in the turbine region for power calculations

Wind_Farm_Layout.py: a code that accounts for wake overlap of multiple turbines and calculates power output of the wind farm

Optimizer: a code that uses Wind_Farm_Layout.py to optimize a simple layout of turbines using SNOPT (from pyoptsparse); note that SNOPT and pyoptsparse must be functional to use this code

VAWT_Wake_Model.f90: the Fortran version of the wake model code; build this using the instructions found in the main directory

Cross_Validate_EMG: the polynomial surface fitting of the EMG coefficients used in the VAWT wake model
