#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsyr2_(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsyr2(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsyr2_(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif

	LONG DOUBLE *x, int *incx, LONG DOUBLE *y, int *incy, 
	LONG DOUBLE *a, int *lda)
{


    /* System generated locals */
    int  i__1, i__2;

    /* Local variables */
    static int info;
    static LONG DOUBLE temp1, temp2;
    static int i, j;
    extern long int lsame_(char *, char *);
    static int ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ void xerbla_(char *, int *);


/*  Purpose   
    =======   

    DSYR2  performs the symmetric rank 2 operation   

       A := alpha*x*y' + alpha*y*x' + A,   

    where alpha is a scalar, x and y are n element vectors and A is an n 
  
    by n symmetric matrix.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the array A is to be referenced as   
             follows:   

                UPLO = 'U' or 'u'   Only the upper triangular part of A   
                                    is to be referenced.   

                UPLO = 'L' or 'l'   Only the lower triangular part of A   
                                    is to be referenced.   

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

    Y      - LONG DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ABS( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y.   
             Unchanged on exit.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
             Unchanged on exit.   

    A      - LONG DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular part of the symmetric matrix and the strictly   
             lower triangular part of A is not referenced. On exit, the   
             upper triangular part of the array A is overwritten by the   
             upper triangular part of the updated matrix.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular part of the symmetric matrix and the strictly   
             upper triangular part of A is not referenced. On exit, the   
             lower triangular part of the array A is overwritten by the   
             lower triangular part of the updated matrix.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             MAX( 1, n ).   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 5;
    } else if (*incy == 0) {
	info = 7;
    } else if (*lda < MAX(1,*n)) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DSYR2 ", &info);
	return;
    }

/*     Quick return if possible. */

    if (*n == 0 || *alpha == 0.) {
	return;
    }

/*     Set up the start points in X and Y if the increments are not both 
  
       unity. */

    if (*incx != 1 || *incy != 1) {
	if (*incx > 0) {
	    kx = 1;
	} else {
	    kx = 1 - (*n - 1) * *incx;
	}
	if (*incy > 0) {
	    ky = 1;
	} else {
	    ky = 1 - (*n - 1) * *incy;
	}
	jx = kx;
	jy = ky;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the triangular part   
       of A. */

    if (lsame_(uplo, "U")) {

/*        Form  A  when A is stored in the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0. || Y(j) != 0.) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L10: */
		    }
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0. || Y(jy) != 0.) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = kx;
		    iy = ky;
		    i__2 = j;
		    for (i = 1; i <= j; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L30: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L40: */
	    }
	}
    } else {

/*        Form  A  when A is stored in the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(j) != 0. || Y(j) != 0.) {
		    temp1 = *alpha * Y(j);
		    temp2 = *alpha * X(j);
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(i) * temp1 
				+ Y(i) * temp2;
/* L50: */
		    }
		}
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0. || Y(jy) != 0.) {
		    temp1 = *alpha * Y(jy);
		    temp2 = *alpha * X(jx);
		    ix = jx;
		    iy = jy;
		    i__2 = *n;
		    for (i = j; i <= *n; ++i) {
			A(i,j) = A(i,j) + X(ix) * temp1 
				+ Y(iy) * temp2;
			ix += *incx;
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		jy += *incy;
/* L80: */
	    }
	}
    }

    return;

/*     End of DSYR2 . */

} /* dsyr2_ */

