This folder contains the following files:

VAWT_Wake_Model.py: the VAWT wake model based originally on cross-validated fitting of CFD wake data

VAWT_Wake_Model.f90: the Fortran version of the wake model code including a simplified integration method for converting vorticity into velocity, blade aerodynamic force calculations, wake overlap calculations, and power calculations

Velocity_Example.py: an example code that calls VAWT_Wake_Model.py and returns velocity plots and results

Cross_Validate_EMG: the polynomial surface fitting of the EMG coefficients used in the VAWT wake model

ACsingle.py: the actuator cylinder code that converges the self-induced velocities in the turbine region for power calculations

Power_Example.py: an example code that calculates power based on multiple wake interaction and correction factors of the actuator cylinder model to create a Cp/TSR curve, a contour plot and polar plot of the normalized power of coupled VAWTs, and a contour plot of the velocity overlap of multiple wakes

Optimizer: a code that uses Wind_Farm_Layout.py to optimize a simple layout of turbines using SNOPT (from pyoptsparse); note that SNOPT and pyoptsparse must be functional to use this code

Poisson_Solver.py: an example of using the Poisson equation and a discrete sine transform to calculate the velocity field from the vorticity data (still under development for more accurate results)
