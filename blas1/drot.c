#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void drot_(int *n, LONG DOUBLE *dx, int *incx, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qrot(int *n, LONG DOUBLE *dx, int *incx, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qrot_(int *n, LONG DOUBLE *dx, int *incx, 
#endif

	LONG DOUBLE *dy, int *incy, LONG DOUBLE *c, LONG DOUBLE *s)
{


    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i;
    static LONG DOUBLE dtemp;
    static int ix, iy;


/*     applies a plane rotation.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*       code for unequal increments or equal increments not equal   
           to 1 */

    ix = 1;
    iy = 1;
    if (*incx < 0) {
	ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0) {
	iy = (-(*n) + 1) * *incy + 1;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp = *c * DX(ix) + *s * DY(iy);
	DY(iy) = *c * DY(iy) - *s * DX(ix);
	DX(ix) = dtemp;
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return;

/*       code for both increments equal to 1 */

L20:
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	dtemp = *c * DX(i) + *s * DY(i);
	DY(i) = *c * DY(i) - *s * DX(i);
	DX(i) = dtemp;
/* L30: */
    }
    return;
} /* drot_ */

