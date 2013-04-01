#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspmv_(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspmv(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspmv_(char *uplo, int *n, LONG DOUBLE *alpha, 
#endif

	LONG DOUBLE *ap, LONG DOUBLE *x, int *incx, LONG DOUBLE *beta, 
	LONG DOUBLE *y, int *incy)
{


    /* System generated locals */
    int i__1, i__2;

    /* Local variables */
    static int info;
    static LONG DOUBLE temp1, temp2;
    static int i, j, k;
    extern long int lsame_(char *, char *);
    static int kk, ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ void xerbla_(char *, int *);


/*  Purpose   
    =======   

    DSPMV  performs the matrix-vector operation   

       y := alpha*A*x + beta*y,   

    where alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric matrix, supplied in packed form.   

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

    AP     - LONG DOUBLE PRECISION array of DIMENSION at least   
             ( ( n*( n + 1 ) )/2 ).   
             Before entry with UPLO = 'U' or 'u', the array AP must   
             contain the upper triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 1, 2 )   
             and a( 2, 2 ) respectively, and so on.   
             Before entry with UPLO = 'L' or 'l', the array AP must   
             contain the lower triangular part of the symmetric matrix   
             packed sequentially, column by column, so that AP( 1 )   
             contains a( 1, 1 ), AP( 2 ) and AP( 3 ) contain a( 2, 1 )   
             and a( 3, 1 ) respectively, and so on.   
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

    BETA   - LONG DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - LONG DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*ABS( INCY ) ).   
             Before entry, the incremented array Y must contain the n   
             element vector y. On exit, Y is overwritten by the updated   
             vector y.   

    INCY   - INTEGER.   
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   
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
#define Y(I) y[(I)-1]
#define X(I) x[(I)-1]
#define AP(I) ap[(I)-1]


    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*incx == 0) {
	info = 6;
    } else if (*incy == 0) {
	info = 9;
    }
    if (info != 0) {
	xerbla_("DSPMV ", &info);
	return;
    }

/*     Quick return if possible. */

    if (*n == 0 || (*alpha == 0. && *beta == 1.)) {
	return;
    }

/*     Set up the start points in  X  and  Y. */

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

/*     Start the operations. In this version the elements of the array AP 
  
       are accessed sequentially with one pass through AP.   

       First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = 0.;
/* L10: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = *n;
		for (i = 1; i <= *n; ++i) {
		    Y(iy) = *beta * Y(iy);
		    iy += *incy;
/* L40: */
		}
	    }
	}
    }
    if (*alpha == 0.) {
	return;
    }
    kk = 1;
    if (lsame_(uplo, "U")) {

/*        Form  y  when AP contains the upper triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		k = kk;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    Y(i) += temp1 * AP(k);
		    temp2 += AP(k) * X(i);
		    ++k;
/* L50: */
		}
		Y(j) = Y(j) + temp1 * AP(kk + j - 1) + *alpha * temp2;
		kk += j;
/* L60: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.;
		ix = kx;
		iy = ky;
		i__2 = kk + j - 2;
		for (k = kk; k <= kk+j-2; ++k) {
		    Y(iy) += temp1 * AP(k);
		    temp2 += AP(k) * X(ix);
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		Y(jy) = Y(jy) + temp1 * AP(kk + j - 1) + *alpha * temp2;
		jx += *incx;
		jy += *incy;
		kk += j;
/* L80: */
	    }
	}
    } else {

/*        Form  y  when AP contains the lower triangle. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		Y(j) += temp1 * AP(kk);
		k = kk + 1;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    Y(i) += temp1 * AP(k);
		    temp2 += AP(k) * X(i);
		    ++k;
/* L90: */
		}
		Y(j) += *alpha * temp2;
		kk += *n - j + 1;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.;
		Y(jy) += temp1 * AP(kk);
		ix = jx;
		iy = jy;
		i__2 = kk + *n - j;
		for (k = kk + 1; k <= kk+*n-j; ++k) {
		    ix += *incx;
		    iy += *incy;
		    Y(iy) += temp1 * AP(k);
		    temp2 += AP(k) * X(ix);
/* L110: */
		}
		Y(jy) += *alpha * temp2;
		jx += *incx;
		jy += *incy;
		kk += *n - j + 1;
/* L120: */
	    }
	}
    }

    return;

/*     End of DSPMV . */

} /* dspmv_ */

