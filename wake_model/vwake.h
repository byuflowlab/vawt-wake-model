
#include <gsl/gsl_integration.h>

struct Arguments{
  // Amplitude
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
};

double tester(double n);
double velocity_fieldx_c(double * in_array,int size);
double EMGdists(double x,double mu,double sigma,double lamda,double scale);

double integrandxext(double x,double y,double x0,double y0,double dia,double loc1,double loc2,double loc3,double spr1,double spr2,double skw1,double skw2,double scl1,double scl2,double scl3);
