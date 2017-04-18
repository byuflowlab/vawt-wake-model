# VAWT Wake Model Testing

This folder contains various files used for validation of the VAWT wake model. Cross Validation Results.xlsx contains the results for the final cross validated polynomial surfaces and the chosen surface highlighted in yellow. The files error_cfd_vort...csv contain comparisons between the CFD and the wake model using different methods of error calculation. Python files are also used to calculate the error between the CFD and the wake model (errorcalc_full.py) and plot the results (error_plot.py). The folders Grid Convergence and MATLAB contain files relating to aspects of the grid convergence study conducted of the CFD model and the files used as a basis of the cross-validation.

The wake data used in this validation is from:

Tescione, G., Ragni, D., He, C., Ferreira, C. S., and van Bussel, G., “Near wake flow analysis of a vertical axis wind turbine by stereoscopic particle image velocimetry,” Renewable Energy, Vol. 70, 2014, pp. 47-61. (tescione_val.py with test_cfd.csv)

Battisti, L., Zanne, L., Dell’Anna, S., Dossena, V., Persico, G., and Paradiso, B., “Aerodynamic Measurements on a Vertical Axis Wind Turbine in a Large Scale Wind Tunnel,” Journal of Energy Resources Technology, Vol. 133, September 2011. (battisti_val.py)

Shamsoddin, S., and Porte-Agel, F., “Large eddy simulation of vertical axis wind turbine wakes,” Energies, Vol. 7, 2011, pp. 8990-8912. (brochier_val.py)
