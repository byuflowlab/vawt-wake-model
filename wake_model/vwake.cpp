#import <iostream>
#include <math.h>
#include "vwake.h"

// Function Declarations for Functions that cannot be called from Python
double integrandx(double y,double x,Arguments Args);
double integrandy(double y,double x,Arguments Args);
double fx(double x,void *p);
double fxy(double x,void *p);
double fy(double x,void *p);
double fyx(double x,void *p);
double argtest(double x,double y,Arguments Args);
double modifyParams(double x, double y, Arguments Args);
double vorticitystrength(double x,double y,Arguments Args);
Arguments unpack(double *in_array);
double vorticitystrengthx(double x,double y,Arguments Args);
double vorticitystrengthy(double x,double y,Arguments Args);

double EMGdists(double x,double mu,double sigma,double lamda,double scale){
  double lss = lamda*sigma*sigma;
  double p1=exp((lamda/2.0)*(2.0*mu+lss-2*x));
  double p2=1.0-erf((mu+lss-x)/(sqrt(2.0)*sigma));
  double EMG=(lamda/2.0)*p1*p2*scale;
  return EMG;
}
void modifyParams(double x,double y,Arguments Args,double params[]){
  double xd  = x/Args.dia; //Normalizing x by the diameter
  //double yd = y/Args.dia; // Normalizing y by the diameter - unused

  // Limiting parameter components to create expected behavior
  double loc1d;
  if (Args.loc1 > -0.001){ // ensure concave down
    loc1d=-0.001;
  }else{
    loc1d=Args.loc1;
  }
  double loc2d;
  if (Args.loc2<0.01){   // ensure slight increase mocing downstream
    loc2d = 0.01;
  }else{
    loc2d=Args.loc2;
  }
  double loc3d;
  if (Args.loc3 < 0.48){ // ensure wake originating from edge of turbine
    loc3d=0.48;
  }else{
    loc3d=Args.loc3;
  }
  double loc = loc1d*xd*xd + loc2d*xd + loc3d; // EMG Location

  double spr1d;
  if (Args.spr1 > -0.001){   // ensure drecrease in value (more spread downstream)
    spr1d= -0.001;
  }else{
    spr1d=Args.spr1;
  }
  double spr2d;
  if (Args.spr2>0.0){      // Ensure calue does not begin positive
    spr2d=0.0;
  }else{
    spr2d=Args.spr2;
  }

  double spr=spr1d*xd+spr2d; // EMG Spread

  double skw1d = Args.skw1;//no limitations necssary
  double skw2d;
  if (Args.skw2>0.0){ // Ensure value does not begin positive
    skw2d=0.0;
  }else{
    skw2d=Args.skw2;
  }

  double skw=skw1d*xd+skw2d; // EMG Skew

  double scl1d;
  if (Args.scl1 <0.0){   // ensure positive maximum vorticity strength
    scl1d=0.0;
  }else{
    scl1d=Args.scl1;
  }
  double scl2d;
  if (Args.scl2<0.05){   // ensure decay moving downstream
    scl2d=0.05;
  }else{
    scl2d=Args.scl2;
  }
  double scl3d;
  if (Args.scl3<0.0){    //ensure decay occurs downstream
    scl3d=0.0;
  } else{
    scl3d=Args.scl3;
  }

  double scl= scl1d/(1.0+exp(scl2d*(xd-scl3d))); // EMG Scale

  // Limiting parameters to the maximum values the EMG distributino can handle
  if (loc < 0.2){
    loc=0.2;
  }
  if (spr<-0.5){
    spr=-0.5;
  }else if (spr >-0.001){
    spr=-0.001;
  }
  if (skw>0.0){
    skw=0.0;
  }
  params[0]=loc;
  params[1]=spr;
  params[2]=skw;
  params[3]=scl;
  //params={loc,spr,skw,scl};
  //return params;
}
double vorticitystrength(double x,double y,Arguments Args){
  //double xd  = x/Args.dia; //Normalizing x by the diameter
  double yd = y/Args.dia; // Normalizing y by the diameter

  double params[4];
  modifyParams(x,y,Args,params);
  double loc = params[0];
  double spr = params[1];
  double skw = params[2];
  double scl = params[3];
  double g1=EMGdists(yd,loc,spr,skw,scl);
  double g2=EMGdists(yd,-loc,-spr,-skw,-scl);
  double gam_lat = (g1-g2);
  return gam_lat;
}
double vorticitystrengthx(double x,double y,Arguments Args){
  //double xd  = x/Args.dia; //Normalizing x by the diameter
  double yd = y/Args.dia; // Normalizing y by the diameter
  double pi = 3.1415926535897932;
  double params[4];
  modifyParams(x,y,Args,params);
  double loc = params[0];
  double spr = params[1];
  double skw = params[2];
  double scl = params[3];
  double kpp = skw*spr*spr;

  //double gam_lat=EMGdists(yd,loc,spr,skw,scl);

  //double gam_lat = (1.0/(2.0*spr))*scl*skw*(exp(-(loc-yd)*(loc-yd)/(2.0*spr*spr))*sqrt(2.0/pi) + exp(-(loc+yd)*(loc+yd)/(2.*spr*spr))*sqrt(2.0/pi) +exp(0.5*skw*(2.0*loc + skw*spr*spr - 2.0*y)*skw*spr*(-1.0 + erf((loc + skw*spr*spr - yd)/(sqrt(2.0)*spr)))) + exp(0.5*skw*(2.0* loc + skw*spr*spr + 2.0*y)*skw*spr*(-1.0 + erf((loc + skw*spr*spr +yd)/(sqrt(2.0)*spr)))));
  double gam_lat =  (1.0/(2.0*spr))*scl*skw*(exp(-(loc-yd)*(loc-yd)/(2.0*spr*spr))*\
  sqrt(2.0/pi) + exp(-(loc+yd)*(loc+yd)/(2.0*spr*spr))*sqrt(2.0/pi) + \
  exp(0.5*skw*(2.0*loc + kpp - 2.0*y)*skw*spr*(-1.0 + \
  erf((loc + kpp - yd)/(sqrt(2.0)*spr)))) + exp(0.5*skw*(2.0*
  loc + kpp + 2.0*y)*skw*spr*(-1.0 + erf((loc + kpp + \
  yd)/(sqrt(2.0)*spr)))));
  //return yd;

  return gam_lat;
}
double vorticitystrengthy(double x,double y,Arguments Args){
  //double xd  = x/Args.dia; //Normalizing x by the diameter
  double yd = y/Args.dia; // Normalizing y by the diameter
  double params[4];
  modifyParams(x,y,Args,params);
  double loc = params[0];
  double spr = params[1];
  double skw = params[2];
  double scl = params[3];

  double g1=EMGdists(yd,loc,spr,skw,scl);
  double g2=EMGdists(yd,-loc,-spr,-skw,-scl);
  double gam_lat = (g1-g2);

  return gam_lat;
}
double integrandx(double x,double y,Arguments Args){
  double gammav=vorticitystrength(x,y,Args);
  double num = (y-Args.y0);
  double den = ((x-Args.x0)*(x-Args.x0))+((y-Args.y0)*(y-Args.y0));
  double inte = gammav*num/den;
  //double inte = gammav*((y-Args.y0)/(((x-Args.x0)*(x-Args.x0))+((y-Args.y0)*(y-Args.y0))));
  //inte = Args.dia;
  return inte;
}
double integrandy(double x,double y,Arguments Args){
  double gammav=vorticitystrength(x,y,Args);
  double inte = gammav*((Args.x0-x)/((x-Args.x0)*(x-Args.x0)+(y-Args.y0)*(y-Args.y0)));
  return inte;
}
Arguments ptoArgs(Arguments Args,void *p){
  // converts GSL Void pointer into Arguments Struct
    double *in_array=(double *) p;
    Args.x0=in_array[0];
    Args.y0=in_array[1];
    Args.dia=in_array[2];
    Args.loc1=in_array[3];
    Args.loc2=in_array[4];
    Args.loc3=in_array[5];
    Args.spr1=in_array[6];
    Args.spr2=in_array[7];
    Args.skw1=in_array[8];
    Args.skw2=in_array[9];
    Args.scl1=in_array[10];
    Args.scl2=in_array[11];
    Args.scl3=in_array[12];
    Args.ybound1=in_array[13];
    Args.ybound2=in_array[14];
    Args.xvalue=in_array[15];
    Args.workspacesize=in_array[16];
    return Args;
}
double fx(double x,void *p){
  Arguments Args = *(Arguments *)p;
  Args.xvalue = x;
  gsl_function F;
  F.function = &fxy;
  F.params = &Args;

  // Initialise values to put the result in
  double result;
  double abserror;
  double epsabs=1.49e-8;
  double epsrel=1.49e-8;

  gsl_integration_qag(&F, Args.ybound1, Args.ybound2, epsabs, epsrel, Args.workspacesize, Args.imethod, Args.giw, &result, &abserror);
  //result=integrandx(5.0,5.0,Args);
  return result;
}
double fxy(double y,void *p){
  Arguments Args = *(Arguments *)p;
  double results;
  results = integrandx(Args.xvalue,y,Args);
  return results;
}
double fy(double x,void *p){
  Arguments Args = *(Arguments*)p;
  Args.xvalue = x;
  gsl_function F;
  F.function = &fyx;
  F.params = &Args;
  // Initialise values to put the result in
  double result;
  double abserror;
  double epsabs=1.49e-8;
  double epsrel=1.49e-8;

  gsl_integration_qag(&F, Args.ybound1, Args.ybound2, epsabs, epsrel, Args.workspacesize, Args.imethod, Args.giw, &result, &abserror);
  return result;
}
double fyx(double y,void *p){
  Arguments Args = *(Arguments*)p;
  double results;
  results = integrandy(Args.xvalue,y,Args);
  return results;
}
double functest(double x,double y,Arguments Args){
  // This is a function tester
  //double value=integrandy(x,y,Args);
  double value = vorticitystrength(x,y,Args);
  return value;
}
Arguments unpack(double *in_array,double allocationsize){
  // Unpack Numpy Array into Arguments Struct
  Arguments Args;
  Args.x0=in_array[2];
  Args.y0=in_array[3];
  Args.dia=in_array[4];
  Args.loc1=in_array[5];
  Args.loc2=in_array[6];
  Args.loc3=in_array[7];
  Args.spr1=in_array[8];
  Args.spr2=in_array[9];
  Args.skw1=in_array[10];
  Args.skw2=in_array[11];
  Args.scl1=in_array[12];
  Args.scl2=in_array[13];
  Args.scl3=in_array[14];
  Args.workspacesize=allocationsize;
  // Set Ybounds in Struct for passing into GSL
  Args.ybound1=-1.0*Args.dia;
  Args.ybound2=Args.dia;
  Args.giw= gsl_integration_workspace_alloc(allocationsize);
  // Return Struct
  return Args;
}
double velocity_fieldx_c(double * in_array,int size){
    // Unpack Pyton (numpy) Values
    // These Two Arguments are only used for debugging
    double x=in_array[0];
    double y=in_array[1];

    int allocationsize=10000; // maximum only affects memory useage, but to little will fail
    //Arguments Args;
    Arguments  Args;
    Args = unpack(in_array,allocationsize);
    // Set bounds of integration
    double xbounds[2]={0,(Args.scl3+5.0)*Args.dia};

    // Set G-K Points
    int imethod = 1;
    /*
    int imethod = 1; // 15pt
    int imethod = 2; // 21 pt
    int imethod = 3; // 31 pt
    int imethod = 4; // 41 pt
    int imethod = 5; // 51 pt
    int imethod = 6; // 61 pt
    */


    Args.imethod=imethod;

    // Allocate integration workspace
    gsl_integration_workspace *giw = gsl_integration_workspace_alloc(allocationsize);

    // Create GSL function
    gsl_function F;
    F.function = &fx;
    F.params = &Args;

    // Initialise values to put the result in
    double result;
    double abserror;
    double epsabs=1.49e-8;
    double epsrel=1.49e-8;
    // Perform integration
    gsl_integration_qag(&F, xbounds[0], xbounds[1], epsabs, epsrel, allocationsize, Args.imethod, giw, &result, &abserror);

    // Free the integration workspace
    gsl_integration_workspace_free(giw);
    gsl_integration_workspace_free(Args.giw);
    //result=0;
    //result = functest(x,y,Args);
   return result;
}
double velocity_fieldy_c(double * in_array,int size){
    // Unpack Pyton (numpy) Values

    double x=in_array[0];
    double y=in_array[1];
    int allocationsize=10000; // maximum only affects memory useage, but to little will fail
    //Arguments Args;
    Arguments Args;
    Args = unpack(in_array,allocationsize);
    // Set bounds of integration
    double xbounds[2]={0,(Args.scl3+5.0)*Args.dia};

    // Set G-K Points
    int imethod = 1;
    /*
    int imethod = 1; // 15pt
    int imethod = 2; // 21 pt
    int imethod = 3; // 31 pt
    int imethod = 4; // 41 pt
    int imethod = 5; // 51 pt
    int imethod = 6; // 61 pt
    */

    Args.imethod=imethod;

    // Allocate integration workspace

    gsl_integration_workspace *giw = gsl_integration_workspace_alloc(allocationsize);

    // Create GSL function
    gsl_function F;
    F.function = &fy;
    F.params = &Args;

    // Initialise values to put the result in
    double result;
    double abserror;
    double epsabs=1.49e-8;
    double epsrel=1.49e-8;
    // Perform integration
    gsl_integration_qag(&F, xbounds[0], xbounds[1], epsabs, epsrel, allocationsize, Args.imethod, giw, &result, &abserror);

    // Free the integration workspace
    gsl_integration_workspace_free(giw);
    gsl_integration_workspace_free(Args.giw);
    //result=0;
    //result = functest(x,y,Args);

   return result;
}
double integrandxext(double y,double x,double x0,double y0,double dia,double loc1,double loc2,double loc3,double spr1,double spr2,double skw1,double skw2,double scl1,double scl2,double scl3){
  Arguments Args;
  Args.x0=x0;
  Args.y0=y0;
  Args.dia=dia;
  Args.loc1=loc1;
  Args.loc2=loc2;
  Args.loc3=loc3;
  Args.spr1=spr1;
  Args.spr2=spr2;
  Args.skw1=skw1;
  Args.skw2=skw2;
  Args.scl1=scl1;
  Args.scl2=scl2;
  Args.scl3=scl3;
  double inte = integrandx(x,y,Args);
  return inte;
}
