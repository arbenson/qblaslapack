#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspr_(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspr(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspr_(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif

	LONG DOUBLE *x, int *incx, LONG DOUBLE *ap)
{


    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int info;
    static LONG DOUBLE temp;
    static int i, j, k;
    extern long int lsame_(char *, char *);
    static int kk, ix, jx, kx;
    extern /* Subroutine */ void xerbla_(char *, int *);


/*  Purpose   
    =======   

    DSPR    performs the symmetric rank 1 operation   

       A := alpha*x*x' + A,   

    where alpha is a real scalar, x is an n element vector and A is an   
    n by n symmetric matrix, supplied in packed form.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the matrix A is supplied in the packed   
             array AP as follows:   

                UPLO = 'U' or 'u'   The upper triangular part of A is   
                                    supplied in AP.   

                UPLO = 'L' or 'l'   The lower triangular part of A is   
                                    supplied in AP.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    ALPHA  - LONG DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    X      - LONG DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ABS( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    AP     - LONG DOUBLE PRECISION array of DIMENSION at least   
             ( ( n*( n + 1 ) )/2 ).   
             Before entry with  UPLO = 'U' or 'u', the array AP must   
             contain the upper triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )   
             and a( 2, 2 ) respectively, and so on. On exit, the array   
             AP is overwritten by the upper triangular part of the   
             updated matrix.   
             Before entry with UPLO = 'L' or 'l', the array AP must   
             contain the lower triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )   
             and a( 3, 1 ) respectively, and so on. On exit, the array   
             AP is overwritten by the lower triangular part of the   
             updated matrix.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define AP(I) ap[(I)-1]
#define X(I) x[(I)-1]


    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    }
    if (info != 0) {
	xerbla_("DSPR  ", &info);
	return;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return;
    }

/*     Set the start point in X if the increment is not unity. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of the array AP 
  
       are accessed sequentially with one pass through AP. */

    kk = 1;
    if (lsame_(uplo, "U")) {

/*        Form  A  when upper triangle is stored in AP. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0.) {
		    temp = *alpha * X(j);
		    k = kk;
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			AP(k) += X(i) * temp;
			++k;
/* L10: */
		    }
		}
		kk += j;
/* L20: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    ix = kx;
		    i__2 = kk + j - 1;
		    for (k = kk; k <= kk+j-1; ++k) {
			AP(k) += X(ix) * temp;
			ix += *incx;
/* L30: */
		    }
		}
		jx += *incx;
		kk += j;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when lower triangle is stored in AP. */

	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0.) {
		    temp = *alpha * X(j);
		    k = kk;
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			AP(k) += X(i) * temp;
			++k;
/* L50: */
		    }
		}
		kk = kk + *n - j + 1;
/* L60: */
	    }
	} else {
	    jx = kx;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    ix = jx;
		    i__2 = kk + *n - j;
		    for (k = kk; k <= kk+*n-j; ++k) {
			AP(k) += X(ix) * temp;
			ix += *incx;
/* L70: */
		    }
		}
		jx += *incx;
		kk = kk + *n - j + 1;
/* L80: */
	    }
	}
    }

    return;

/*     End of DSPR  . */

} /* dspr_ */

