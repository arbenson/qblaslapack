#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlagtf_(int *n, LONG DOUBLE *a, LONG DOUBLE *lambda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlagtf(int *n, LONG DOUBLE *a, LONG DOUBLE *lambda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlagtf_(int *n, LONG DOUBLE *a, LONG DOUBLE *lambda, 
#endif

	LONG DOUBLE *b, LONG DOUBLE *c, LONG DOUBLE *tol, LONG DOUBLE *d, int 
	*in, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLAGTF factorizes the matrix (T - lambda*I), where T is an n by n   
    tridiagonal matrix and lambda is a scalar, as   

       T - lambda*I = PLU,   

    where P is a permutation matrix, L is a unit lower tridiagonal matrix 
  
    with at most one non-zero sub-diagonal elements per column and U is   
    an upper triangular matrix with at most two non-zero super-diagonal   
    elements per column.   

    The factorization is obtained by Gaussian elimination with partial   
    pivoting and implicit row scaling.   

    The parameter LAMBDA is included in the routine so that DLAGTF may   
    be used, in conjunction with DLAGTS, to obtain eigenvectors of T by   
    inverse iteration.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix T.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, A must contain the diagonal elements of T.   

            On exit, A is overwritten by the n diagonal elements of the   
            upper triangular matrix U of the factorization of T.   

    LAMBDA  (input) LONG DOUBLE PRECISION   
            On entry, the scalar lambda.   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, B must contain the (n-1) super-diagonal elements of 
  
            T.   

            On exit, B is overwritten by the (n-1) super-diagonal   
            elements of the matrix U of the factorization of T.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, C must contain the (n-1) sub-diagonal elements of   
            T.   

            On exit, C is overwritten by the (n-1) sub-diagonal elements 
  
            of the matrix L of the factorization of T.   

    TOL     (input) LONG DOUBLE PRECISION   
            On entry, a relative tolerance used to indicate whether or   
            not the matrix (T - lambda*I) is nearly singular. TOL should 
  
            normally be chose as approximately the largest relative error 
  
            in the elements of T. For example, if the elements of T are   
            correct to about 4 significant figures, then TOL should be   
            set to about 5*10**(-4). If TOL is supplied as less than eps, 
  
            where eps is the relative machine precision, then the value   
            eps is used in place of TOL.   

    D       (output) LONG DOUBLE PRECISION array, dimension (N-2)   
            On exit, D is overwritten by the (n-2) second super-diagonal 
  
            elements of the matrix U of the factorization of T.   

    IN      (output) INTEGER array, dimension (N)   
            On exit, IN contains details of the permutation matrix P. If 
  
            an interchange occurred at the kth step of the elimination,   
            then IN(k) = 1, otherwise IN(k) = 0. The element IN(n)   
            returns the smallest positive int j such that   

               ABS( u(j,j) ).le. norm( (T - lambda*I)(j) )*TOL,   

            where norm( A(j) ) denotes the sum of the absolute values of 
  
            the jth row of the matrix A. If no such j exists then IN(n)   
            is returned as zero. If IN(n) is returned as positive, then a 
  
            diagonal element of U is small, indicating that   
            (T - lambda*I) is singular or nearly singular,   

    INFO    (output)   
            = 0   : successful exit   
            .lt. 0: if INFO = -k, the kth argument had an illegal value   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1, d__2;
    /* Local variables */
    static LONG DOUBLE temp, mult;
    static int k;
    static LONG DOUBLE scale1, scale2;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE tl;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE eps, piv1, piv2;


#define IN(I) in[(I)-1]
#define D(I) d[(I)-1]
#define C(I) c[(I)-1]
#define B(I) b[(I)-1]
#define A(I) a[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -1;
	i__1 = -(*info);
	xerbla_("DLAGTF", &i__1);
	return;
    }

    if (*n == 0) {
	return;
    }

    A(1) -= *lambda;
    IN(*n) = 0;
    if (*n == 1) {
	if (A(1) == 0.) {
	    IN(1) = 1;
	}
	return;
    }


#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Epsilon");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Epsilon");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Epsilon");
#endif


    tl = MAX(*tol,eps);
    scale1 = ABS(A(1)) + ABS(B(1));
    i__1 = *n - 1;
    for (k = 1; k <= *n-1; ++k) {
	A(k + 1) -= *lambda;
	scale2 = (d__1 = C(k), ABS(d__1)) + (d__2 = A(k + 1), ABS(d__2));
	if (k < *n - 1) {
	    scale2 += (d__1 = B(k + 1), ABS(d__1));
	}
	if (A(k) == 0.) {
	    piv1 = 0.;
	} else {
	    piv1 = (d__1 = A(k), ABS(d__1)) / scale1;
	}
	if (C(k) == 0.) {
	    IN(k) = 0;
	    piv2 = 0.;
	    scale1 = scale2;
	    if (k < *n - 1) {
		D(k) = 0.;
	    }
	} else {
	    piv2 = (d__1 = C(k), ABS(d__1)) / scale2;
	    if (piv2 <= piv1) {
		IN(k) = 0;
		scale1 = scale2;
		C(k) /= A(k);
		A(k + 1) -= C(k) * B(k);
		if (k < *n - 1) {
		    D(k) = 0.;
		}
	    } else {
		IN(k) = 1;
		mult = A(k) / C(k);
		A(k) = C(k);
		temp = A(k + 1);
		A(k + 1) = B(k) - mult * temp;
		if (k < *n - 1) {
		    D(k) = B(k + 1);
		    B(k + 1) = -mult * D(k);
		}
		B(k) = temp;
		C(k) = mult;
	    }
	}
	if (MAX(piv1,piv2) <= tl && IN(*n) == 0) {
	    IN(*n) = k;
	}
/* L10: */
    }
    if ((d__1 = A(*n), ABS(d__1)) <= scale1 * tl && IN(*n) == 0) {
	IN(*n) = *n;
    }

    return;

/*     End of DLAGTF */

} /* dlagtf_ */

