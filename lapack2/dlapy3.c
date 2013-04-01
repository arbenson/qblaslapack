#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlapy3_(LONG DOUBLE *x, LONG DOUBLE *y, LONG DOUBLE *z)
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlapy3(LONG DOUBLE *x, LONG DOUBLE *y, LONG DOUBLE *z)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlapy3_(LONG DOUBLE *x, LONG DOUBLE *y, LONG DOUBLE *z)
#endif

{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause   
    unnecessary overflow.   

    Arguments   
    =========   

    X       (input) LONG DOUBLE PRECISION   
    Y       (input) LONG DOUBLE PRECISION   
    Z       (input) LONG DOUBLE PRECISION   
            X, Y and Z specify the values x, y and z.   

    ===================================================================== 
*/
/* >>Start of File<<   
       System generated locals */
    LONG DOUBLE ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE xabs, yabs, zabs, w;



    xabs = ABS(*x);
    yabs = ABS(*y);
    zabs = ABS(*z);
/* Computing MAX */
    d__1 = MAX(xabs,yabs);
    w = MAX(d__1,zabs);
    if (w == 0.) {
	ret_val = 0.;
    } else {
/* Computing 2nd power */
	d__1 = xabs / w;
/* Computing 2nd power */
	d__2 = yabs / w;
/* Computing 2nd power */
	d__3 = zabs / w;
	ret_val = w * sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
    }
    return ret_val;

/*     End of DLAPY3 */

} /* dlapy3_ */

