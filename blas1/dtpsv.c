#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtpsv_(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtpsv(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtpsv_(char *uplo, char *trans, char *diag, int *n, 
#endif

	LONG DOUBLE *ap, LONG DOUBLE *x, int *incx)
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
    static long int nounit;


/*  Purpose   
    =======   

    DTPSV  solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix, supplied in packed form. 
  

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   A'*x = b.   

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

    AP     - LONG DOUBLE PRECISION array of DIMENSION at least   
             ( ( n*( n + 1 ) )/2 ).   
             Before entry with  UPLO = 'U' or 'u', the array AP must   
             contain the upper triangular matrix packed sequentially,   
             column by column, so that AP( 1 ) contains a( 1, 1 ),   
             AP( 2 ) and AP( 3 ) contain a( 1, 2 ) and a( 2, 2 )   
             respectively, and so on.   
             Before entry with UPLO = 'L' or 'l', the array AP must   
             contain the lower triangular matrix packed sequentially,   
             column by column, so that AP( 1 ) contains a( 1, 1 ),   
             AP( 2 ) and AP( 3 ) contain a( 2, 1 ) and a( 3, 1 )   
             respectively, and so on.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced, but are assumed to be unity.   
             Unchanged on exit.   

    X      - LONG DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ABS( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

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
#define AP(I) ap[(I)-1]


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
    } else if (*incx == 0) {
	info = 7;
    }
    if (info != 0) {
	xerbla_("DTPSV ", &info);
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

/*     Start the operations. In this version the elements of AP are   
       accessed sequentially with one pass through AP. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    if (X(j) != 0.) {
			if (nounit) {
			    X(j) /= AP(kk);
			}
			temp = X(j);
			k = kk - 1;
			for (i = j - 1; i >= 1; --i) {
			    X(i) -= temp * AP(k);
			    --k;
/* L10: */
			}
		    }
		    kk -= j;
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		    if (X(jx) != 0.) {
			if (nounit) {
			    X(jx) /= AP(kk);
			}
			temp = X(jx);
			ix = jx;
			i__1 = kk - j + 1;
			for (k = kk - 1; k >= kk-j+1; --k) {
			    ix -= *incx;
			    X(ix) -= temp * AP(k);
/* L30: */
			}
		    }
		    jx -= *incx;
		    kk -= j;
/* L40: */
		}
	    }
	} else {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(j) != 0.) {
			if (nounit) {
			    X(j) /= AP(kk);
			}
			temp = X(j);
			k = kk + 1;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    X(i) -= temp * AP(k);
			    ++k;
/* L50: */
			}
		    }
		    kk += *n - j + 1;
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    if (X(jx) != 0.) {
			if (nounit) {
			    X(jx) /= AP(kk);
			}
			temp = X(jx);
			ix = jx;
			i__2 = kk + *n - j;
			for (k = kk + 1; k <= kk+*n-j; ++k) {
			    ix += *incx;
			    X(ix) -= temp * AP(k);
/* L70: */
			}
		    }
		    jx += *incx;
		    kk += *n - j + 1;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame_(uplo, "U")) {
	    kk = 1;
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(j);
		    k = kk;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
			temp -= AP(k) * X(i);
			++k;
/* L90: */
		    }
		    if (nounit) {
			temp /= AP(kk + j - 1);
		    }
		    X(j) = temp;
		    kk += j;
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		    temp = X(jx);
		    ix = kx;
		    i__2 = kk + j - 2;
		    for (k = kk; k <= kk+j-2; ++k) {
			temp -= AP(k) * X(ix);
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
			temp /= AP(kk + j - 1);
		    }
		    X(jx) = temp;
		    jx += *incx;
		    kk += j;
/* L120: */
		}
	    }
	} else {
	    kk = *n * (*n + 1) / 2;
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		    temp = X(j);
		    k = kk;
		    i__1 = j + 1;
		    for (i = *n; i >= j+1; --i) {
			temp -= AP(k) * X(i);
			--k;
/* L130: */
		    }
		    if (nounit) {
			temp /= AP(kk - *n + j);
		    }
		    X(j) = temp;
		    kk -= *n - j + 1;
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		    temp = X(jx);
		    ix = kx;
		    i__1 = kk - (*n - (j + 1));
		    for (k = kk; k >= kk-(*n-(j+1)); --k) {
			temp -= AP(k) * X(ix);
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
			temp /= AP(kk - *n + j);
		    }
		    X(jx) = temp;
		    jx -= *incx;
		    kk -= *n - j + 1;
/* L160: */
		}
	    }
	}
    }

    return;

/*     End of DTPSV . */

} /* dtpsv_ */

