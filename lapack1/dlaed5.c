#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaed5_(int *i, LONG DOUBLE *d, LONG DOUBLE *z, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaed5(int *i, LONG DOUBLE *d, LONG DOUBLE *z, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaed5_(int *i, LONG DOUBLE *d, LONG DOUBLE *z, 
#endif

	LONG DOUBLE *delta, LONG DOUBLE *rho, LONG DOUBLE *dlam)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Oak Ridge National Lab, Argonne National Lab, 
  
       Courant Institute, NAG Ltd., and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    This subroutine computes the I-th eigenvalue of a symmetric rank-one 
  
    modification of a 2-by-2 diagonal matrix   

               diag( D )  +  RHO *  Z * transpose(Z) .   

    The diagonal elements in the array D are assumed to satisfy   

               D(i) < D(j)  for  i < j .   

    We also assume RHO > 0 and that the Euclidean norm of the vector   
    Z is one.   

    Arguments   
    =========   

    I      (input) INTEGER   
           The index of the eigenvalue to be computed.  I = 1 or I = 2.   

    D      (input) LONG DOUBLE PRECISION array, dimension (2)   
           The original eigenvalues.  We assume D(1) < D(2).   

    Z      (input) LONG DOUBLE PRECISION array, dimension (2)   
           The components of the updating vector.   

    DELTA  (output) LONG DOUBLE PRECISION array, dimension (2)   
           The vector DELTA contains the information necessary   
           to construct the eigenvectors.   

    RHO    (input) LONG DOUBLE PRECISION   
           The scalar in the symmetric updating formula.   

    DLAM   (output) LONG DOUBLE PRECISION   
           The computed lambda_I, the I-th updated eigenvalue.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE temp, b, c, w, del, tau;


#define DELTA(I) delta[(I)-1]
#define Z(I) z[(I)-1]
#define D(I) d[(I)-1]


    del = D(2) - D(1);
    if (*i == 1) {
	w = *rho * 2. * (Z(2) * Z(2) - Z(1) * Z(1)) / del + 1.;
	if (w > 0.) {
	    b = del + *rho * (Z(1) * Z(1) + Z(2) * Z(2));
	    c = *rho * Z(1) * Z(1) * del;

/*           B > ZERO, always */

	    d__1 = b * b - c * 4.;
            tau = c * 2. / (b + sqrt(( ABS(d__1))));
	    *dlam = D(1) + tau;
	    DELTA(1) = -Z(1) / tau;
	    DELTA(2) = Z(2) / (del - tau);
	} else {
	    b = -del + *rho * (Z(1) * Z(1) + Z(2) * Z(2));
	    c = *rho * Z(2) * Z(2) * del;
	    if (b > 0.) {
		tau = c * -2. / (b + sqrt(b * b + c * 4.));
	    } else {
		tau = (b - sqrt(b * b + c * 4.)) / 2.;
	    }
	    *dlam = D(2) + tau;
	    DELTA(1) = -Z(1) / (del + tau);
	    DELTA(2) = -Z(2) / tau;
	}
	temp = sqrt(DELTA(1) * DELTA(1) + DELTA(2) * DELTA(2));
	DELTA(1) /= temp;
	DELTA(2) /= temp;
    } else {

/*     Now I=2 */

	b = -del + *rho * (Z(1) * Z(1) + Z(2) * Z(2));
	c = *rho * Z(2) * Z(2) * del;
	if (b > 0.) {
	    tau = (b + sqrt(b * b + c * 4.)) / 2.;
	} else {
	    tau = c * 2. / (-b + sqrt(b * b + c * 4.));
	}
	*dlam = D(2) + tau;
	DELTA(1) = -Z(1) / (del + tau);
	DELTA(2) = -Z(2) / tau;
	temp = sqrt(DELTA(1) * DELTA(1) + DELTA(2) * DELTA(2));
	DELTA(1) /= temp;
	DELTA(2) /= temp;
    }
    return;

/*     End OF DLAED5 */

} /* dlaed5_ */

