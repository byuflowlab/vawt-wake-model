
#include <gsl/gsl_integration.h>

struct Arguments{
  // Defines all parameters and integration methods for inner and outer intergation
  double x0;
  double y0;
  double dia;
  double loc1;
  double loc2;
  double loc3;
  double spr1;
  double spr2;
  double skw1;
  double skw2;
  double scl1;
  double scl2;
  double scl3;
  double ybound1;
  double ybound2;
  double xvalue;
  double workspacesize;
  double imethod;
  gsl_integration_workspace *giw;
  int allocationsize;
  double* workspace1;
  double* workspace2;
};


// Function definitions that are visible in Python
double velocity_fieldx_c(double * in_array,int size);
double EMGdists(double x,double mu,double sigma,double lamda,double scale);
double velocity_fieldy_c(double * in_array,int size);
