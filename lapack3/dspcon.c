#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspcon_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspcon(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspcon_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif

	ipiv, LONG DOUBLE *anorm, LONG DOUBLE *rcond, LONG DOUBLE *work, int 
	*iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSPCON estimates the reciprocal of the condition number (in the   
    1-norm) of a real symmetric packed matrix A using the factorization   
    A = U*D*U**T or A = L*D*L**T computed by DSPTRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the details of the factorization are stored 
  
            as an upper or lower triangular matrix.   
            = 'U':  Upper triangular, form is A = U*D*U**T;   
            = 'L':  Lower triangular, form is A = L*D*L**T.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The block diagonal matrix D and the multipliers used to   
            obtain the factor U or L as computed by DSPTRF, stored as a   
            packed triangular matrix.   

    IPIV    (input) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D   
            as determined by DSPTRF.   

    ANORM   (input) LONG DOUBLE PRECISION   
            The 1-norm of the original matrix A.   

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an   
            estimate of the 1-norm of inv(A) computed in this routine.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (2*N)   

    IWORK    (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1;
    /* Local variables */
    static int kase, i;
    extern long int lsame_(char *, char *);
    static long int upper;
    static int ip;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *, int *), xerbla_(char *, int *);
    static LONG DOUBLE ainvnm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsptrs_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsptrs(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsptrs_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *);



#define IWORK(I) iwork[(I)-1]
#define WORK(I) work[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*anorm < 0.) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPCON", &i__1);
	return;
    }

/*     Quick return if possible */

    *rcond = 0.;
    if (*n == 0) {
	*rcond = 1.;
	return;
    } else if (*anorm <= 0.) {
	return;
    }

/*     Check that the diagonal matrix D is nonsingular. */

    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

	ip = *n * (*n + 1) / 2;
	for (i = *n; i >= 1; --i) {
	    if (IPIV(i) > 0 && AP(ip) == 0.) {
		return;
	    }
	    ip -= i;
/* L10: */
	}
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

	ip = 1;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (IPIV(i) > 0 && AP(ip) == 0.) {
		return;
	    }
	    ip = ip + *n - i + 1;
/* L20: */
	}
    }

/*     Estimate the 1-norm of the inverse. */

    kase = 0;
L30:

#ifdef PETSC_PREFIX_SUFFIX
    dlacon_(n, &WORK(*n + 1), &WORK(1), &IWORK(1), &ainvnm, &kase);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlacon(n, &WORK(*n + 1), &WORK(1), &IWORK(1), &ainvnm, &kase);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlacon_(n, &WORK(*n + 1), &WORK(1), &IWORK(1), &ainvnm, &kase);
#endif

    if (kase != 0) {

/*        Multiply by inv(L*D*L') or inv(U*D*U'). */


#ifdef PETSC_PREFIX_SUFFIX
	dsptrs_(uplo, n, &c__1, &AP(1), &IPIV(1), &WORK(1), n, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsptrs(uplo, n, &c__1, &AP(1), &IPIV(1), &WORK(1), n, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsptrs_(uplo, n, &c__1, &AP(1), &IPIV(1), &WORK(1), n, info);
#endif

	goto L30;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

    return;

/*     End of DSPCON */

} /* dspcon_ */

