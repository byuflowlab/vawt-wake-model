#import <iostream>
#include <math.h>
#include "vwake.h"
//
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
double rombergintegration(double (*F)(double,void*),double a,double b,Arguments Args,double error,double *workspace1,double *workspace2);


double EMGdists(double x,double mu,double sigma,double lamda,double scale){
  // Exponentially modified Gaussian distribution Implementation
  double lss = lamda*sigma*sigma; // Stored Valued to Save on computation time
  double p1=exp((lamda/2.0)*(2.0*mu+lss-2*x)); // part 1
  double p2=1.0-erf((mu+lss-x)/(sqrt(2.0)*sigma)); // part 2
  double EMG=(lamda/2.0)*p1*p2*scale; // combine together and scale
  return EMG;
}
void modifyParams(double x,double y,Arguments Args,double params[]){
  double xd  = x/Args.dia; //Normalizing x by the diameter

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
}
double vorticitystrength(double x,double y,Arguments Args){
  double yd = y/Args.dia; // Normalizing y by the diameter

  double params[4]; // Initalize Array that is returned by reference
  modifyParams(x,y,Args,params); // Ensure correct behavior
  // Unpack Results
  double loc = params[0];
  double spr = params[1];
  double skw = params[2];
  double scl = params[3];
  // Calculate Vorticity Strength
  double g1=EMGdists(yd,loc,spr,skw,scl);
  double g2=EMGdists(yd,-loc,-spr,-skw,-scl);
  double gam_lat = (g1-g2);
  return gam_lat;
}
double vorticitystrengthx(double x,double y,Arguments Args){
  double yd = y/Args.dia; // Normalizing y by the diameter
  double pi = 3.141592653589793238462643383;

  double params[4];// Initalize Array that is returned by reference
  modifyParams(x,y,Args,params); // Ensure correct behavior
  // Unpack Results
  double loc = params[0];
  double spr = params[1];
  double skw = params[2];
  double scl = params[3];
  double kpp = skw*spr*spr;

  // Calculate Voriticty Strength X
  double gam_lat =  (1.0/(2.0*spr))*scl*skw*(exp(-(loc-yd)*(loc-yd)/(2.0*spr*spr))*\
  sqrt(2.0/pi) + exp(-(loc+yd)*(loc+yd)/(2.0*spr*spr))*sqrt(2.0/pi) + \
  exp(0.5*skw*(2.0*loc + kpp - 2.0*y)*skw*spr*(-1.0 + \
  erf((loc + kpp - yd)/(sqrt(2.0)*spr)))) + exp(0.5*skw*(2.0*
  loc + kpp + 2.0*y)*skw*spr*(-1.0 + erf((loc + kpp + \
  yd)/(sqrt(2.0)*spr)))));
  return gam_lat;
}
double vorticitystrengthy(double x,double y,Arguments Args){
  double yd = y/Args.dia; // Normalizing y by the diameter
  double params[4];  // Initalize Array that is returned by reference
  modifyParams(x,y,Args,params);// Ensure correct behavior
  // Unpack Results
  double loc = params[0];
  double spr = params[1];
  double skw = params[2];
  double scl = params[3];
  // Calculate Vorticity Strength
  double g1=EMGdists(yd,loc,spr,skw,scl);
  double g2=EMGdists(yd,-loc,-spr,-skw,-scl);
  double gam_lat = (g1-g2);

  return gam_lat;
}
double integrandx(double x,double y,Arguments Args){
  //Calculates the X integrand
  double gammav=vorticitystrength(x,y,Args);
  double num = (y-Args.y0);
  double den = ((x-Args.x0)*(x-Args.x0))+((y-Args.y0)*(y-Args.y0));
  double inte = gammav*num/den;
  return inte;
}
double integrandy(double x,double y,Arguments Args){
  //Calculates the Y integrand
  double gammav=vorticitystrength(x,y,Args);
  double num = (Args.x0-x);
  double den = (x-Args.x0)*(x-Args.x0)+(y-Args.y0)*(y-Args.y0);
  double inte=gammav*num/den;
  return inte;
}

double fx(double x,void *p){
  // Outside X Integral Function
  Arguments Args = *(Arguments *)p; // Convert void pointer to Arguments struct
  Args.xvalue = x; // Store x in Argument Struct
  // Define Inner Integrand Function
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
  // X Inner Integrand function
  Arguments Args = *(Arguments *)p; // Convert void pointer to Arguments struct
  double results;
  results = integrandx(Args.xvalue,y,Args);
  return results; // Return Integrand x
}
double fxr(double x,void *p){
  // Outside X Integral Function
  Arguments Args = *(Arguments *)p; // Convert void pointer to Arguments struct
  Args.xvalue = x; // Store x in Argument Struct
  // Define Inner Integrand Function

  // Initialise values to put the result in
  double result=0.0;
  //double abserror;
  double epsabs=1.49e-8;
  //double epsrel=1.49e-8;

  // Perform Romberg Integration
  result=rombergintegration(*fxr,Args.ybound1, Args.ybound2,Args,epsabs,Args.workspace1,Args.workspace2);
  return result;
}
double fxyr(double y,void *p){
  // X Inner Integrand function
  Arguments Args = *(Arguments *)p; // Convert void pointer to Arguments struct
  double results;
  results = integrandx(Args.xvalue,y,Args);
  return results; // Return Integrand x
}
double fy(double x,void *p){
  // Outside of Y double integrand
  Arguments Args = *(Arguments*)p; // Convert void pointer to Arguments struct
  Args.xvalue = x; // Store x in Argument Struct

  // Define Inner Integrand Function
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
  // Y Inner Integrand function
  Arguments Args = *(Arguments*)p; // Convert void pointer to Arguments struct
  double results;
  results = integrandy(Args.xvalue,y,Args);
  return results; // Return Integrand y
}
double fyr(double x,void *p){
  // Outside Y Integral Function
  Arguments Args = *(Arguments *)p; // Convert void pointer to Arguments struct
  Args.xvalue = x; // Store x in Argument Struct
  // Define Inner Integrand Function

  // Initialise values to put the result in
  double result=0.0;
  //double abserror;
  double epsabs=1.49e-8;
  //double epsrel=1.49e-8;

  // Perform Romberg Integration
  result=rombergintegration(*fyr,Args.ybound1, Args.ybound2,Args,epsabs,Args.workspace1,Args.workspace2);

  return result;
}
double fyxr(double y,void *p){
  // Y Inner Integrand function
  Arguments Args = *(Arguments *)p; // Convert void pointer to Arguments struct
  double results;
  results = integrandy(Args.xvalue,y,Args);
  return results; // Return Integrand x
}
double functest(double x,double y,Arguments Args){
  // This is a function tester
  //double value=integrandy(x,y,Args);
  double value = vorticitystrength(x,y,Args);
  return value;
}
Arguments unpack(double *in_array,double allocationsize,int imethod){
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
  if (imethod<7){
    // Allocation Inner Integral workspace to minimize memory functions (speed)
    Args.giw= gsl_integration_workspace_alloc(allocationsize);
  }
  // Return Struct
  return Args;
}

double rombergintegration(double (*F)(double,void*),double a,double b,Arguments Args,double error,double* workspace1,double* workspace2){
      // Use Allocated Memory
      //double R1[sizework] = workspace1;
      //double R2[sizework]= workspace2;
      //double R1[max_steps], R2[max_steps]; //buffers
    int max_steps = Args.allocationsize;
    double acc = 1.49e-8;
    double *Rp = &workspace1[0], *Rc = &workspace2[0]; //Rp is previous row, Rc is current row
    double h = (b-a); //step size
    Rp[0] = (F(a,&Args) + F(b,&Args))*h*.5; //first trapezoidal step
    /*
     for(size_t i = 1; i < max_steps; ++i){
        h /= 2.;
        double c = 0;
        size_t ep = 1 << (i-1); //2^(n-1)
        for(size_t j = 1; j <= ep; ++j){
           c += F(a+(2*j-1)*h,&Args);
        }
        Rc[0] = h*c + .5*Rp[0]; //R(i,0)

        for(size_t j = 1; j <= i; ++j){
           double n_k = pow(4, j);
           Rc[j] = (n_k*Rc[j-1] - Rp[j-1])/(n_k-1); //compute R(i,j)
        }

        if(i > 1 && fabs(Rp[i-1]-Rc[i]) < acc){
          return Rc[i-1];
        }

        //swap Rn and Rc as we only need the last row
        double *rt = Rp;
        Rp = Rc;
        Rc = rt;
     }
     */
     return Rp[max_steps-1]; //return our best guess
}

double velocity_fieldx_c(double * in_array,int size){
    // Unpack Pyton (numpy) Values

    // Set G-K Points
    int imethod = 2;
    /*
    GSL Guass-Konrod
    int imethod = 1; // 15pt
    int imethod = 2; // 21 pt
    int imethod = 3; // 31 pt
    int imethod = 4; // 41 pt
    int imethod = 5; // 51 pt
    int imethod = 6; // 61 pt
    7 for Romberg
    */

    int allocationsize=10000; // maximum only affects memory useage, but too little will fail


    // These Two Arguments are only used for debugging
    //double x=in_array[0]; // To be removed
    //double y=in_array[1]; // To be removed

    //Convert Input to Parameter Struct
    Arguments  Args;
    Args = unpack(in_array,allocationsize,imethod);

    Args.imethod=imethod;// Set to Arguments Struct

    // Initialise values to put the result in
    double result=0.0;
    //double abserror;
    double epsabs=1.49e-8;
    double epsrel=1.49e-8;

    if(Args.imethod<7){

      // Set bounds of integration
      double xbounds[2]={0,(Args.scl3+5.0)*Args.dia};

      // Allocate integration workspace
      gsl_integration_workspace *giw = gsl_integration_workspace_alloc(allocationsize);

      // Create GSL function
      gsl_function F;
      F.function = &fx;
      F.params = &Args;

      // Initialise values to put the result in
      double abserror;

      // Perform integration
      gsl_integration_qag(&F, xbounds[0], xbounds[1], epsabs, epsrel, allocationsize, Args.imethod, giw, &result, &abserror);

      // Free the integration workspace
      gsl_integration_workspace_free(giw);
      gsl_integration_workspace_free(Args.giw); // Free inner integration workspace
      //result = functest(x,y,Args);
    }
    if(Args.imethod==7){
      // Set bounds of integration
      double xbounds[2]={0,(Args.scl3+5.0)*Args.dia};

      // Allocate integration workspace
      double workspace1[allocationsize];
      double workspace2[allocationsize];



      result=rombergintegration(*fxr,xbounds[0],xbounds[1],Args,epsabs,workspace2,workspace1);

    }
   return result;
}
double velocity_fieldy_c(double * in_array,int size){
  // Unpack Pyton (numpy) Values

  // Set G-K Points
  int imethod = 2;
  /*
  GSL Guass-Konrod
  int imethod = 1; // 15pt
  int imethod = 2; // 21 pt
  int imethod = 3; // 31 pt
  int imethod = 4; // 41 pt
  int imethod = 5; // 51 pt
  int imethod = 6; // 61 pt
  int imethod = 7; // for Romberg
  */
  int allocationsize=10000; // maximum only affects memory useage, but to little will fail

  // These Two Arguments are only used for debugging
  //double x=in_array[0]; // To be removed
  //double y=in_array[1]; // To be removed

  Arguments  Args;
  Args = unpack(in_array,allocationsize,imethod);

  Args.imethod=imethod;// Set to Arguments Struct

  // Initialise values to put the result in
  double result;
  double abserror;
  double epsabs=1.49e-6;
  double epsrel=1.49e-8;


  if(Args.imethod<7){

    //Convert Input to Parameter Struct

    // Set bounds of integration
    double xbounds[2]={0,(Args.scl3+5.0)*Args.dia};

    // Allocate integration workspace
    gsl_integration_workspace *giw = gsl_integration_workspace_alloc(allocationsize);

    // Create GSL function
    gsl_function F;
    F.function = &fy;
    F.params = &Args;


    // Perform integration
    gsl_integration_qag(&F, xbounds[0], xbounds[1], epsabs, epsrel, allocationsize, Args.imethod, giw, &result, &abserror);

    // Free the integration workspace
    gsl_integration_workspace_free(giw);
    gsl_integration_workspace_free(Args.giw); // Free inner integration workspace
    //result = functest(x,y,Args);
  }
  if(Args.imethod==7){

    // Set bounds of integration
    double xbounds[2]={0,(Args.scl3+5.0)*Args.dia};


    // Allocate integration workspace
    double workspace1[allocationsize];
    double workspace2[allocationsize];
    double workspace3[allocationsize];
    double workspace4[allocationsize];
    Args.workspace1=workspace3;
    Args.workspace2=workspace4;

    result=rombergintegration(*fyr,xbounds[0],xbounds[1],Args,epsabs,workspace2,workspace1);

  }
 return result;
}
