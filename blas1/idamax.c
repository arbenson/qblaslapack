
/*  -- translated by f2c (version 19940927).
*/

#include <math.h>
#define ABS(a) ( ((a)<0.0)   ? -(a) : (a) )


#ifdef PETSC_PREFIX_SUFFIX
int idamax_(int *n, LONG DOUBLE *dx, int *incx)
#endif
#ifdef Q_C_PREFIX_SUFFIX
int iqamax(int *n, LONG DOUBLE *dx, int *incx)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
int iqamax_(int *n, LONG DOUBLE *dx, int *incx)
#endif

{


    /* System generated locals */
    int ret_val, i__1;
    LONG DOUBLE d__1;

    /* Local variables */
    static LONG DOUBLE dmax__;
    static int i, ix;


/*     finds the index of element having max. absolute value.   
       jack dongarra, linpack, 3/11/78.   
       modified 3/93 to return if incx .le. 0.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DX(I) dx[(I)-1]


    ret_val = 0;
    if (*n < 1 || *incx <= 0) {
	return ret_val;
    }
    ret_val = 1;
    if (*n == 1) {
	return ret_val;
    }
    if (*incx == 1) {
	goto L20;
    }

/*        code for increment not equal to 1 */

    ix = 1;
    dmax__ = ABS(DX(1));
    ix += *incx;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if ((d__1 = DX(ix), ABS(d__1)) <= dmax__) {
	    goto L5;
	}
	ret_val = i;
	dmax__ = (d__1 = DX(ix), ABS(d__1));
L5:
	ix += *incx;
/* L10: */
    }
    return ret_val;

/*        code for increment equal to 1 */

L20:
    dmax__ = ABS(DX(1));
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	if ((d__1 = DX(i), ABS(d__1)) <= dmax__) {
	    goto L30;
	}
	ret_val = i;
	dmax__ = (d__1 = DX(i), ABS(d__1));
L30:
	;
    }
    return ret_val;
} /* idamax_ */

