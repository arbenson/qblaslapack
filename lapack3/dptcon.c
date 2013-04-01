#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dptcon_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qptcon(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qptcon_(int *n, LONG DOUBLE *d, LONG DOUBLE *e, 
#endif

	LONG DOUBLE *anorm, LONG DOUBLE *rcond, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPTCON computes the reciprocal of the condition number (in the   
    1-norm) of a real symmetric positive definite tridiagonal matrix   
    using the factorization A = L*D*L**T or A = U**T*D*U computed by   
    DPTTRF.   

    Norm(inv(A)) is computed by a direct method, and the reciprocal of   
    the condition number is computed as   
                 RCOND = 1 / (ANORM * norm(inv(A))).   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from the   
            factorization of A, as computed by DPTTRF.   

    E       (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) off-diagonal elements of the unit bidiagonal factor 
  
            U or L from the factorization of A,  as computed by DPTTRF.   

    ANORM   (input) LONG DOUBLE PRECISION   
            The 1-norm of the original matrix A.   

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is the   
            1-norm of inv(A) computed in this routine.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The method used is described in Nicholas J. Higham, "Efficient   
    Algorithms for Computing the Condition Number of a Tridiagonal   
    Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1;
    /* Local variables */
    static int i, ix;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE ainvnm;



#define WORK(I) work[(I)-1]
#define E(I) e[(I)-1]
#define D(I) d[(I)-1]


    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*anorm < 0.) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPTCON", &i__1);
	return;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return;
    } else if (*anorm == 0.) {
	return;
    }

/*     Check that D(1:N) is positive. */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (D(i) <= 0.) {
	    return;
	}
/* L10: */
    }

/*     Solve M(A) * x = e, where M(A) = (m(i,j)) is given by   

          m(i,j) =  ABS(A(i,j)), i = j,   
          m(i,j) = -ABS(A(i,j)), i .ne. j,   

       and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.   

       Solve M(L) * x = e. */

    WORK(1) = 1.;
    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	WORK(i) = WORK(i - 1) * (d__1 = E(i - 1), ABS(d__1)) + 1.;
/* L20: */
    }

/*     Solve D * M(L)' * x = b. */

    WORK(*n) /= D(*n);
    for (i = *n - 1; i >= 1; --i) {
	WORK(i) = WORK(i) / D(i) + WORK(i + 1) * (d__1 = E(i), ABS(d__1));
/* L30: */
    }

/*     Compute AINVNM = MAX(x(i)), 1<=i<=n. */


#ifdef PETSC_PREFIX_SUFFIX
    ix = idamax_(n, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    ix = iqamax(n, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    ix = iqamax_(n, &WORK(1), &c__1);
#endif

    ainvnm = (d__1 = WORK(ix), ABS(d__1));

/*     Compute the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

    return;

/*     End of DPTCON */

} /* dptcon_ */

