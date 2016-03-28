This folder contains the following files:

VAWT_Wake_Model.py: the VAWT wake model that uses the trend fitting from the CFD vorticity data

database_call.py: a code that identifies the EMG and SMG parameters for a given tip-speed ratio and solidity

example.py: a code that can be run that calls VAWT_Wake_Model.py and returns and plots the results

Wind_Farm_Layout.py: a code that accounts for wake overlap of multiple turbines and calculates power output of the wind farm

Optimizer: a code that uses Wind_Farm_Layout.py to optimize a simple layout of turbines using SNOPT (from pyoptsparse); note that SNOPT and pyoptsparse must be functional to use this code

vorticity.f90: the Fortran code that runs the vorticity calculations for VAWT_Wake_Model.py; build this using the instructions found in the main directory
