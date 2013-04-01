#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgtcon_(char *norm, int *n, LONG DOUBLE *dl, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgtcon(char *norm, int *n, LONG DOUBLE *dl, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgtcon_(char *norm, int *n, LONG DOUBLE *dl, 
#endif

	LONG DOUBLE *d, LONG DOUBLE *du, LONG DOUBLE *du2, int *ipiv, 
	LONG DOUBLE *anorm, LONG DOUBLE *rcond, LONG DOUBLE *work, int *
	iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGTCON estimates the reciprocal of the condition number of a real   
    tridiagonal matrix A using the LU factorization as computed by   
    DGTTRF.   

    An estimate is obtained for norm(inv(A)), and the reciprocal of the   
    condition number is computed as RCOND = 1 / (ANORM * norm(inv(A))).   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies whether the 1-norm condition number or the   
            infinity-norm condition number is required:   
            = '1' or 'O':  1-norm;   
            = 'I':         Infinity-norm.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    DL      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) multipliers that define the matrix L from the   
            LU factorization of A as computed by DGTTRF.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the upper triangular matrix U from 
  
            the LU factorization of A.   

    DU      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) elements of the first superdiagonal of U.   

    DU2     (input) LONG DOUBLE PRECISION array, dimension (N-2)   
            The (n-2) elements of the second superdiagonal of U.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    ANORM   (input) LONG DOUBLE PRECISION   
            If NORM = '1' or 'O', the 1-norm of the original matrix A.   
            If NORM = 'I', the infinity-norm of the original matrix A.   

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal of the condition number of the matrix A,   
            computed as RCOND = 1/(ANORM * AINVNM), where AINVNM is an   
            estimate of the 1-norm of inv(A) computed in this routine.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (2*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1;
    /* Local variables */
    static int kase, kase1, i;
    extern long int lsame_(char *, char *);

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
    static long int onenrm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgttrs_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgttrs(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgttrs_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *,
	     LONG DOUBLE *, int *, int *);



#define IWORK(I) iwork[(I)-1]
#define WORK(I) work[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define DU2(I) du2[(I)-1]
#define DU(I) du[(I)-1]
#define D(I) d[(I)-1]
#define DL(I) dl[(I)-1]


    *info = 0;
    onenrm = *(unsigned char *)norm == '1' || lsame_(norm, "O");
    if (! onenrm && ! lsame_(norm, "I")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*anorm < 0.) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGTCON", &i__1);
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

/*     Check that D(1:N) is non-zero. */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	if (D(i) == 0.) {
	    return;
	}
/* L10: */
    }

    ainvnm = 0.;
    if (onenrm) {
	kase1 = 1;
    } else {
	kase1 = 2;
    }
    kase = 0;
L20:

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
	if (kase == kase1) {

/*           Multiply by inv(U)*inv(L). */


#ifdef PETSC_PREFIX_SUFFIX
	    dgttrs_("No transpose", n, &c__1, &DL(1), &D(1), &DU(1), &DU2(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgttrs("No transpose", n, &c__1, &DL(1), &D(1), &DU(1), &DU2(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgttrs_("No transpose", n, &c__1, &DL(1), &D(1), &DU(1), &DU2(1), 
#endif

		    &IPIV(1), &WORK(1), n, info);
	} else {

/*           Multiply by inv(L')*inv(U'). */


#ifdef PETSC_PREFIX_SUFFIX
	    dgttrs_("Transpose", n, &c__1, &DL(1), &D(1), &DU(1), &DU2(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgttrs("Transpose", n, &c__1, &DL(1), &D(1), &DU(1), &DU2(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgttrs_("Transpose", n, &c__1, &DL(1), &D(1), &DU(1), &DU2(1), &
#endif

		    IPIV(1), &WORK(1), n, info);
	}
	goto L20;
    }

/*     Compute the estimate of the reciprocal condition number. */

    if (ainvnm != 0.) {
	*rcond = 1. / ainvnm / *anorm;
    }

    return;

/*     End of DGTCON */

} /* dgtcon_ */

