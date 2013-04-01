#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlartv_(int *n, LONG DOUBLE *x, int *incx, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlartv(int *n, LONG DOUBLE *x, int *incx, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlartv_(int *n, LONG DOUBLE *x, int *incx, 
#endif

	LONG DOUBLE *y, int *incy, LONG DOUBLE *c, LONG DOUBLE *s, int *
	incc)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLARTV applies a vector of real plane rotations to elements of the   
    real vectors x and y. For i = 1,2,...,n   

       ( x(i) ) := (  c(i)  s(i) ) ( x(i) )   
       ( y(i) )    ( -s(i)  c(i) ) ( y(i) )   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of plane rotations to be applied.   

    X       (input/output) LONG DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCX)   
            The vector x.   

    INCX    (input) INTEGER   
            The increment between elements of X. INCX > 0.   

    Y       (input/output) LONG DOUBLE PRECISION array,   
                           dimension (1+(N-1)*INCY)   
            The vector y.   

    INCY    (input) INTEGER   
            The increment between elements of Y. INCY > 0.   

    C       (input) LONG DOUBLE PRECISION array, dimension (1+(N-1)*INCC)   
            The cosines of the plane rotations.   

    S       (input) LONG DOUBLE PRECISION array, dimension (1+(N-1)*INCC)   
            The sines of the plane rotations.   

    INCC    (input) INTEGER   
            The increment between elements of C and S. INCC > 0.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    /* Local variables */
    static int i, ic, ix, iy;
    static LONG DOUBLE xi, yi;


#define S(I) s[(I)-1]
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
	X(ix) = C(ic) * xi + S(ic) * yi;
	Y(iy) = C(ic) * yi - S(ic) * xi;
	ix += *incx;
	iy += *incy;
	ic += *incc;
/* L10: */
    }
    return;

/*     End of DLARTV */

} /* dlartv_ */

