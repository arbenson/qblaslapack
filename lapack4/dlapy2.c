
#include <math.h>
#define MIN(a,b)           ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)           ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)             ( ((a)<0.0)   ? -(a) : (a) )



#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlapy2_(LONG DOUBLE *x, LONG DOUBLE *y)
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlapy2(LONG DOUBLE *x, LONG DOUBLE *y)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlapy2_(LONG DOUBLE *x, LONG DOUBLE *y)
#endif

{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary 
  
    overflow.   

    Arguments   
    =========   

    X       (input) LONG DOUBLE PRECISION   
    Y       (input) LONG DOUBLE PRECISION   
            X and Y specify the values x and y.   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    LONG DOUBLE ret_val, d__1;
    /* Local variables */
    static LONG DOUBLE xabs, yabs, w, z;



    xabs = ABS(*x);
    yabs = ABS(*y);
    w = MAX(xabs,yabs);
    z = MIN(xabs,yabs);
    if (z == 0.) {
	ret_val = w;
    } else {
/* Computing 2nd power */
	d__1 = z / w;
	ret_val = w * sqrt(d__1 * d__1 + 1.);
    }
    return ret_val;

/*     End of DLAPY2 */

} /* dlapy2_ */

