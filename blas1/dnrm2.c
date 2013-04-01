
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dnrm2_(int *n, LONG DOUBLE *x, int *incx)
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qnrm2(int *n, LONG DOUBLE *x, int *incx)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qnrm2_(int *n, LONG DOUBLE *x, int *incx)
#endif

{
  LONG DOUBLE norm = 0.0;
  int    k,i;

  if (*incx != 1) {
    k = 0;
    for ( i=0; i<*n; i++ ) {
      norm += x[k]*x[k];
      k    += *incx;
    }
  } else {
    for ( i=0; i<*n; i++ ) {
      norm += x[i]*x[i];
    }
  }
  norm = sqrt(norm);
  return norm;
}

