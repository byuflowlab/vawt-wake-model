%module vwake
%{
#define SWIG_FILE_WITH_INIT
#include "vwake.h"
%}

%include "vwake.h"

/*  include the numpy typemaps */
%include "numpy.i"
/*  need this for correct module initialization */
%init %{
    import_array();
%}

/*  typemaps for the two arrays, the second will be modified in-place */
%apply (double* IN_ARRAY1, int DIM1) {(double * in_array, int size_in)}

/*  Wrapper for velocity_fieldx that massages the types */
%inline %{
    /*  takes as input two numpy arrays */
    double velocity_fieldx(double * in_array, int size_in) {
        /*  calls the original funcion, providing only the size of the first */
        return velocity_fieldx_c(in_array, size_in);
    }
%}
%inline %{
    /*  takes as input two numpy arrays */
    double velocity_fieldy(double * in_array, int size_in) {
        /*  calls the original funcion, providing only the size of the first */
        return velocity_fieldy_c(in_array, size_in);
    }
%}
