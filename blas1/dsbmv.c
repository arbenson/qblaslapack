#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsbmv_(char *uplo, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsbmv(char *uplo, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsbmv_(char *uplo, int *n, int *k, LONG DOUBLE *
#endif

	alpha, LONG DOUBLE *a, int *lda, LONG DOUBLE *x, int *incx, 
	LONG DOUBLE *beta, LONG DOUBLE *y, int *incy)
{


    /* System generated locals */
    int i__1, i__2, i__3, i__4;

    /* Local variables */
    static int info;
    static LONG DOUBLE temp1, temp2;
    static int i, j, l;
    extern long int lsame_(char *, char *);
    static int kplus1, ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ void xerbla_(char *, int *);


/*  Purpose   
    =======   

    DSBMV  performs the matrix-vector  operation   

       y := alpha*A*x + beta*y,   

    where alpha and beta are scalars, x and y are n element vectors and   
    A is an n by n symmetric band matrix, with k super-diagonals.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the upper or lower   
             triangular part of the band matrix A is being supplied as   
             follows:   

                UPLO = 'U' or 'u'   The upper triangular part of A is   
                                    being supplied.   

                UPLO = 'L' or 'l'   The lower triangular part of A is   
                                    being supplied.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    K      - INTEGER.   
             On entry, K specifies the number of super-diagonals of the   
             matrix A. K must satisfy  0 .le. K.   
             Unchanged on exit.   

    ALPHA  - LONG DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - LONG DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )   
             by n part of the array A must contain the upper triangular   
             band part of the symmetric matrix, supplied column by   
             column, with the leading diagonal of the matrix in row   
             ( k + 1 ) of the array, the first super-diagonal starting at 
  
             position 2 in row k, and so on. The top left k by k triangle 
  
             of the array A is not referenced.   
             The following program segment will transfer the upper   
             triangular part of a symmetric band matrix from conventional 
  
             full matrix storage to band storage:   

                   DO 20, J = 1, N   
                      M = K + 1 - J   
                      DO 10, I = MAX( 1, J - K ), J   
                         A( M + I, J ) = matrix( I, J )   
                10    CONTINUE   
                20 CONTINUE   

             Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )   
             by n part of the array A must contain the lower triangular   
             band part of the symmetric matrix, supplied column by   
             column, with the leading diagonal of the matrix in row 1 of 
  
             the array, the first sub-diagonal starting at position 1 in 
  
             row 2, and so on. The bottom right k by k triangle of the   
             array A is not referenced.   
             The following program segment will transfer the lower   
             triangular part of a symmetric band matrix from conventional 
  
             full matrix storage to band storage:   

                   DO 20, J = 1, N   
                      M = 1 - J   
                      DO 10, I = J, MIN( N, J + K )   
                         A( M + I, J ) = matrix( I, J )   
                10    CONTINUE   
                20 CONTINUE   

             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             ( k + 1 ).   
             Unchanged on exit.   

    X      - LONG DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*ABS( INCX ) ).   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - LONG DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta.   
             Unchanged on exit.   

    Y      - LONG DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*ABS( INCY ) ).   
             Before entry, the incremented array Y must contain the   
             vector y. On exit, Y is overwritten by the updated vector y. 
  

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
#define X(I) x[(I)-1]
#define Y(I) y[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (*n < 0) {
	info = 2;
    } else if (*k < 0) {
	info = 3;
    } else if (*lda < *k + 1) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    } else if (*incy == 0) {
	info = 11;
    }
    if (info != 0) {
	xerbla_("DSBMV ", &info);
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

/*     Start the operations. In this version the elements of the array A 
  
       are accessed sequentially with one pass through A.   

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
    if (lsame_(uplo, "U")) {

/*        Form  y  when upper triangle of A is stored. */

	kplus1 = *k + 1;
	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		l = kplus1 - j;
/* Computing MAX */
		i__2 = 1, i__3 = j - *k;
		i__4 = j - 1;
		for (i = MAX(1,j-*k); i <= j-1; ++i) {
		    Y(i) += temp1 * A(l+i,j);
		    temp2 += A(l+i,j) * X(i);
/* L50: */
		}
		Y(j) = Y(j) + temp1 * A(kplus1,j) + *alpha * temp2;
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
		l = kplus1 - j;
/* Computing MAX */
		i__4 = 1, i__2 = j - *k;
		i__3 = j - 1;
		for (i = MAX(1,j-*k); i <= j-1; ++i) {
		    Y(iy) += temp1 * A(l+i,j);
		    temp2 += A(l+i,j) * X(ix);
		    ix += *incx;
		    iy += *incy;
/* L70: */
		}
		Y(jy) = Y(jy) + temp1 * A(kplus1,j) + *alpha * 
			temp2;
		jx += *incx;
		jy += *incy;
		if (j > *k) {
		    kx += *incx;
		    ky += *incy;
		}
/* L80: */
	    }
	}
    } else {

/*        Form  y  when lower triangle of A is stored. */

	if (*incx == 1 && *incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(j);
		temp2 = 0.;
		Y(j) += temp1 * A(1,j);
		l = 1 - j;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = MIN(i__4,i__2);
		for (i = j + 1; i <= MIN(*n,j+*k); ++i) {
		    Y(i) += temp1 * A(l+i,j);
		    temp2 += A(l+i,j) * X(i);
/* L90: */
		}
		Y(j) += *alpha * temp2;
/* L100: */
	    }
	} else {
	    jx = kx;
	    jy = ky;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp1 = *alpha * X(jx);
		temp2 = 0.;
		Y(jy) += temp1 * A(1,j);
		l = 1 - j;
		ix = jx;
		iy = jy;
/* Computing MIN */
		i__4 = *n, i__2 = j + *k;
		i__3 = MIN(i__4,i__2);
		for (i = j + 1; i <= MIN(*n,j+*k); ++i) {
		    ix += *incx;
		    iy += *incy;
		    Y(iy) += temp1 * A(l+i,j);
		    temp2 += A(l+i,j) * X(ix);
/* L110: */
		}
		Y(jy) += *alpha * temp2;
		jx += *incx;
		jy += *incy;
/* L120: */
	    }
	}
    }

    return;

/*     End of DSBMV . */

} /* dsbmv_ */

