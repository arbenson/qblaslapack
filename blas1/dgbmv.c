#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))

/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbmv_(char *trans, int *m, int *n, int *kl, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbmv(char *trans, int *m, int *n, int *kl, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbmv_(char *trans, int *m, int *n, int *kl, 
#endif

	int *ku, LONG DOUBLE *alpha, LONG DOUBLE *a, int *lda, 
	LONG DOUBLE *x, int *incx, LONG DOUBLE *beta, LONG DOUBLE *y, 
	int *incy)
{


    /* System generated locals */
    int  i__1, i__2, i__3, i__4, i__5, i__6;

    /* Local variables */
    static int info;
    static LONG DOUBLE temp;
    static int lenx, leny, i, j, k;
    extern long int lsame_(char *, char *);
    static int ix, iy, jx, jy, kx, ky;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static int kup1;


/*  Purpose   
    =======   

    DGBMV  performs one of the matrix-vector operations   

       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   

    where alpha and beta are scalars, x and y are vectors and A is an   
    m by n band matrix, with kl sub-diagonals and ku super-diagonals.   

    Parameters   
    ==========   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the operation to be performed as   
             follows:   

                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   

                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   

                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

             Unchanged on exit.   

    M      - INTEGER.   
             On entry, M specifies the number of rows of the matrix A.   
             M must be at least zero.   
             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the number of columns of the matrix A. 
  
             N must be at least zero.   
             Unchanged on exit.   

    KL     - INTEGER.   
             On entry, KL specifies the number of sub-diagonals of the   
             matrix A. KL must satisfy  0 .le. KL.   
             Unchanged on exit.   

    KU     - INTEGER.   
             On entry, KU specifies the number of super-diagonals of the 
  
             matrix A. KU must satisfy  0 .le. KU.   
             Unchanged on exit.   

    ALPHA  - LONG DOUBLE PRECISION.   
             On entry, ALPHA specifies the scalar alpha.   
             Unchanged on exit.   

    A      - LONG DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry, the leading ( kl + ku + 1 ) by n part of the   
             array A must contain the matrix of coefficients, supplied   
             column by column, with the leading diagonal of the matrix in 
  
             row ( ku + 1 ) of the array, the first super-diagonal   
             starting at position 2 in row ku, the first sub-diagonal   
             starting at position 1 in row ( ku + 2 ), and so on.   
             Elements in the array A that do not correspond to elements   
             in the band matrix (such as the top left ku by ku triangle) 
  
             are not referenced.   
             The following program segment will transfer a band matrix   
             from conventional full matrix storage to band storage:   

                   DO 20, J = 1, N   
                      K = KU + 1 - J   
                      DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )   
                         A( K + I, J ) = matrix( I, J )   
                10    CONTINUE   
                20 CONTINUE   

             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             ( kl + ku + 1 ).   
             Unchanged on exit.   

    X      - LONG DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( n - 1 )*ABS( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*ABS( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   
             Unchanged on exit.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   

    BETA   - LONG DOUBLE PRECISION.   
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   
             Unchanged on exit.   

    Y      - LONG DOUBLE PRECISION array of DIMENSION at least   
             ( 1 + ( m - 1 )*ABS( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*ABS( INCY ) ) otherwise.   
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
    if (! lsame_(trans, "N") && ! lsame_(trans, "T") && ! 
	    lsame_(trans, "C")) {
	info = 1;
    } else if (*m < 0) {
	info = 2;
    } else if (*n < 0) {
	info = 3;
    } else if (*kl < 0) {
	info = 4;
    } else if (*ku < 0) {
	info = 5;
    } else if (*lda < *kl + *ku + 1) {
	info = 8;
    } else if (*incx == 0) {
	info = 10;
    } else if (*incy == 0) {
	info = 13;
    }
    if (info != 0) {
	xerbla_("DGBMV ", &info);
	return;
    }

/*     Quick return if possible. */

    if (*m == 0 || *n == 0 || (*alpha == 0. && *beta == 1.)) {
	return;
    }

/*     Set  LENX  and  LENY, the lengths of the vectors x and y, and set 
  
       up the start points in  X  and  Y. */

    if (lsame_(trans, "N")) {
	lenx = *n;
	leny = *m;
    } else {
	lenx = *m;
	leny = *n;
    }
    if (*incx > 0) {
	kx = 1;
    } else {
	kx = 1 - (lenx - 1) * *incx;
    }
    if (*incy > 0) {
	ky = 1;
    } else {
	ky = 1 - (leny - 1) * *incy;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through the band part of A.   

       First form  y := beta*y. */

    if (*beta != 1.) {
	if (*incy == 1) {
	    if (*beta == 0.) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(i) = 0.;
/* L10: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(i) = *beta * Y(i);
/* L20: */
		}
	    }
	} else {
	    iy = ky;
	    if (*beta == 0.) {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
		    Y(iy) = 0.;
		    iy += *incy;
/* L30: */
		}
	    } else {
		i__1 = leny;
		for (i = 1; i <= leny; ++i) {
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
    kup1 = *ku + 1;
    if (lsame_(trans, "N")) {

/*        Form  y := alpha*A*x + y. */

	jx = kx;
	if (*incy == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    k = kup1 - j;
/* Computing MAX */
		    i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__4 = MIN(i__5,i__6);
		    for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i) {
			Y(i) += temp * A(k+i,j);
/* L50: */
		    }
		}
		jx += *incx;
/* L60: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		if (X(jx) != 0.) {
		    temp = *alpha * X(jx);
		    iy = ky;
		    k = kup1 - j;
/* Computing MAX */
		    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
		    i__5 = *m, i__6 = j + *kl;
		    i__3 = MIN(i__5,i__6);
		    for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i) {
			Y(iy) += temp * A(k+i,j);
			iy += *incy;
/* L70: */
		    }
		}
		jx += *incx;
		if (j > *ku) {
		    ky += *incy;
		}
/* L80: */
	    }
	}
    } else {

/*        Form  y := alpha*A'*x + y. */

	jy = ky;
	if (*incx == 1) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp = 0.;
		k = kup1 - j;
/* Computing MAX */
		i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__2 = MIN(i__5,i__6);
		for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i) {
		    temp += A(k+i,j) * X(i);
/* L90: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
/* L100: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		temp = 0.;
		ix = kx;
		k = kup1 - j;
/* Computing MAX */
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__4 = MIN(i__5,i__6);
		for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i) {
		    temp += A(k+i,j) * X(ix);
		    ix += *incx;
/* L110: */
		}
		Y(jy) += *alpha * temp;
		jy += *incy;
		if (j > *ku) {
		    kx += *incx;
		}
/* L120: */
	    }
	}
    }

    return;

/*     End of DGBMV . */

} /* dgbmv_ */

