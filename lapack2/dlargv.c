#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlargv_(int *n, LONG DOUBLE *x, int *incx, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlargv(int *n, LONG DOUBLE *x, int *incx, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlargv_(int *n, LONG DOUBLE *x, int *incx, 
#endif

	LONG DOUBLE *y, int *incy, LONG DOUBLE *c, int *incc)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARGV generates a vector of real plane rotations, determined by   
    elements of the real vectors x and y. For i = 1,2,...,n   

       (  c(i)  s(i) ) ( x(i) ) = ( a(i) )   
       ( -s(i)  c(i) ) ( y(i) ) = (   0  )   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of plane rotations to be generated.   

    X       (input/output) LONG DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCX)   
            On entry, the vector x.   
            On exit, x(i) is overwritten by a(i), for i = 1,...,n.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    Y       (input/output) LONG DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCY)   
            On entry, the vector y.   
            On exit, the sines of the plane rotations.   

    INCY    (input) INTEGER   
            The increment between elements of Y. INCY > 0.   

    C       (output) LONG DOUBLE PRECISION array, dimension (1+(N-1)*INCC)   
            The cosines of the plane rotations.   

    INCC    (input) INTEGER   
            The increment between elements of C. INCC > 0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static int i;
    static LONG DOUBLE w;
    static int ic, ix, iy;
    static LONG DOUBLE xi, yi, tt;


#define C(I) c[(I)-1]
#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]


    ix = 1;
    iy = 1;
    ic = 1;
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	xi = X(ix);
	yi = Y(iy);
	if (xi == 0.) {
	    C(ic) = 0.;
	    Y(iy) = 1.;
	    X(ix) = yi;
	} else {
/* Computing MAX */
	    d__1 = ABS(xi), d__2 = ABS(yi);
	    w = MAX(d__1,d__2);
	    xi /= w;
	    yi /= w;
	    tt = sqrt(xi * xi + yi * yi);
	    C(ic) = xi / tt;
	    Y(iy) = yi / tt;
	    X(ix) = w * tt;
	}
	ix += *incx;
	iy += *incy;
	ic += *incc;
/* L10: */
    }
    return;

/*     End of DLARGV */

} /* dlargv_ */

