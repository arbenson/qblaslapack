
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void daxpy_(int *n, LONG DOUBLE *da, LONG DOUBLE *dx, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qaxpy(int *n, LONG DOUBLE *da, LONG DOUBLE *dx, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qaxpy_(int *n, LONG DOUBLE *da, LONG DOUBLE *dx, 
#endif

	int *incx, LONG DOUBLE *dy, int *incy)
{


    /* System generated locals */
    int i__1;

    /* Local variables */
    static int i, m, ix, iy, mp1;


/*     constant times a vector plus a vector.   
       uses unrolled loops for increments equal to one.   
       jack dongarra, linpack, 3/11/78.   
       modified 12/3/93, array(1) declarations changed to array(*)   


    
   Parameter adjustments   
       Function Body */
#define DY(I) dy[(I)-1]
#define DX(I) dx[(I)-1]


    if (*n <= 0) {
	return;
    }
    if (*da == 0.) {
	return;
    }
    if (*incx == 1 && *incy == 1) {
	goto L20;
    }

/*        code for unequal increments or equal increments   
            not equal to 1 */

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
	DY(iy) += *da * DX(ix);
	ix += *incx;
	iy += *incy;
/* L10: */
    }
    return;

/*        code for both increments equal to 1   


          clean-up loop */

L20:
    m = *n % 4;
    if (m == 0) {
	goto L40;
    }
    i__1 = m;
    for (i = 1; i <= m; ++i) {
	DY(i) += *da * DX(i);
/* L30: */
    }
    if (*n < 4) {
	return;
    }
L40:
    mp1 = m + 1;
    i__1 = *n;
    for (i = mp1; i <= *n; i += 4) {
	DY(i) += *da * DX(i);
	DY(i + 1) += *da * DX(i + 1);
	DY(i + 2) += *da * DX(i + 2);
	DY(i + 3) += *da * DX(i + 3);
/* L50: */
    }
    return;
} /* daxpy_ */

