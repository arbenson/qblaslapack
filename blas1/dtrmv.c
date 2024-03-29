#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtrmv_(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtrmv(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtrmv_(char *uplo, char *trans, char *diag, int *n, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *x, int *incx)
{


    /* System generated locals */
    int  i__1, i__2;

    /* Local variables */
    static int info;
    static LONG DOUBLE temp;
    static int i, j;
    extern long int lsame_(char *, char *);
    static int ix, jx, kx;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int nounit;


/*  Purpose   
    =======   

    DTRMV  performs one of the matrix-vector operations   

       x := A*x,   or   x := A'*x,   

    where x is an n element vector and  A is an n by n unit, or non-unit, 
  
    upper or lower triangular matrix.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   x := A*x.   

                TRANS = 'T' or 't'   x := A'*x.   

                TRANS = 'C' or 'c'   x := A'*x.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - LONG DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             MAX( 1, n ).   
             Unchanged on exit.   

    X      - LONG DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ABS( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element vector x. On exit, X is overwritten with the   
             tranformed vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
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

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < MAX(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("DTRMV ", &info);
	return;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return;
    }

    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (lsame_(trans, "N")) {

/*        Form  x := A*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(j) != 0.) {
			temp = X(j);
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    X(i) += temp * A(i,j);
/* L10: */
			}
			if (nounit) {
			    X(j) *= A(j,j);
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(jx) != 0.) {
			temp = X(jx);
			ix = kx;
			i__2 = j - 1;
			for (i = 1; i <= j-1; ++i) {
			    X(ix) += temp * A(i,j);
			    ix += *incx;
/* L30: */
			}
			if (nounit) {
			    X(jx) *= A(j,j);
			}
		    }
		    jx += *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (X(j) != 0.) {
			temp = X(j);
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    X(i) += temp * A(i,j);
/* L50: */
			}
			if (nounit) {
			    X(j) *= A(j,j);
			}
		    }
/* L60: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    if (X(jx) != 0.) {
			temp = X(jx);
			ix = kx;
			i__1 = j + 1;
			for (i = *n; i >= j+1; --i) {
			    X(ix) += temp * A(i,j);
			    ix -= *incx;
/* L70: */
			}
			if (nounit) {
			    X(jx) *= A(j,j);
			}
		    }
		    jx -= *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := A'*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = X(j);
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = j - 1; i >= 1; --i) {
			temp += A(i,j) * X(i);
/* L90: */
		    }
		    X(j) = temp;
/* L100: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    temp = X(jx);
		    ix = jx;
		    if (nounit) {
			temp *= A(j,j);
		    }
		    for (i = j - 1; i >= 1; --i) {
			ix -= *incx;
			temp += A(i,j) * X(ix);
/* L110: */
		    }
		    X(jx) = temp;
		    jx -= *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(j);
		    if (nounit) {
			temp *= A(j,j);
		    }
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			temp += A(i,j) * X(i);
/* L130: */
		    }
		    X(j) = temp;
/* L140: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(jx);
		    ix = jx;
		    if (nounit) {
			temp *= A(j,j);
		    }
		    i__2 = *n;
		    for (i = j + 1; i <= *n; ++i) {
			ix += *incx;
			temp += A(i,j) * X(ix);
/* L150: */
		    }
		    X(jx) = temp;
		    jx += *incx;
/* L160: */
		}
	    }
	}
    }

    return;

/*     End of DTRMV . */

} /* dtrmv_ */

