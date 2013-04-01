#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dptsvx_(char *fact, int *n, int *nrhs, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qptsvx(char *fact, int *n, int *nrhs, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qptsvx_(char *fact, int *n, int *nrhs, 
#endif

	LONG DOUBLE *d, LONG DOUBLE *e, LONG DOUBLE *df, LONG DOUBLE *ef, 
	LONG DOUBLE *b, int *ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *
	rcond, LONG DOUBLE *ferr, LONG DOUBLE *berr, LONG DOUBLE *work, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DPTSVX uses the factorization A = L*D*L**T to compute the solution   
    to a real system of linear equations A*X = B, where A is an N-by-N   
    symmetric positive definite tridiagonal matrix and X and B are   
    N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'N', the matrix A is factored as A = L*D*L**T, where L   
       is a unit lower bidiagonal matrix and D is diagonal.  The   
       factorization can also be regarded as having the form   
       A = U**T*D*U.   

    2. The factored form of A is used to compute the condition number   
       of the matrix A.  If the reciprocal of the condition number is   
       less than machine precision, steps 3 and 4 are skipped.   

    3. The system of equations is solved for X using the factored form   
       of A.   

    4. Iterative refinement is applied to improve the computed solution   
       matrix and calculate error bounds and backward error estimates   
       for it.   

    Arguments   
    =========   

    FACT    (input) CHARACTER*1   
            Specifies whether or not the factored form of A has been   
            supplied on entry.   
            = 'F':  On entry, DF and EF contain the factored form of A.   
                    D, E, DF, and EF will not be modified.   
            = 'N':  The matrix A will be copied to DF and EF and   
                    factored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the tridiagonal matrix A.   

    E       (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) subdiagonal elements of the tridiagonal matrix A.   

    DF      (input or output) LONG DOUBLE PRECISION array, dimension (N)   
            If FACT = 'F', then DF is an input argument and on entry   
            contains the n diagonal elements of the diagonal matrix D   
            from the L*D*L**T factorization of A.   
            If FACT = 'N', then DF is an output argument and on exit   
            contains the n diagonal elements of the diagonal matrix D   
            from the L*D*L**T factorization of A.   

    EF      (input or output) LONG DOUBLE PRECISION array, dimension (N-1)   
            If FACT = 'F', then EF is an input argument and on entry   
            contains the (n-1) subdiagonal elements of the unit   
            bidiagonal factor L from the L*D*L**T factorization of A.   
            If FACT = 'N', then EF is an output argument and on exit   
            contains the (n-1) subdiagonal elements of the unit   
            bidiagonal factor L from the L*D*L**T factorization of A.   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The N-by-NRHS right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            If INFO = 0, the N-by-NRHS solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= MAX(1,N).   

    RCOND   (output) LONG DOUBLE PRECISION   
            The reciprocal condition number of the matrix A.  If RCOND   
            is less than the machine precision (in particular, if   
            RCOND = 0), the matrix is singular to working precision.   
            This condition is indicated by a return code of INFO > 0,   
            and the solution and error bounds are not computed.   

    FERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).   

    BERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in any   
            element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, and i is   
                  <= N  the leading minor of order i of A is not   
                  positive definite, so the factorization could not be   
                  completed unless i = N, and the solution and error   
                  bounds could not be computed.   
                  = N+1 RCOND is less than machine precision.  The   
                  factorization has been completed, but the matrix is   
                  singular to working precision, and the solution and   
                  error bounds have not been computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    /* Local variables */
    extern long int lsame_(char *, char *);
    static LONG DOUBLE anorm;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static long int nofact;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacpy_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacpy_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlanst_(char *, int *, LONG DOUBLE *, LONG DOUBLE *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dptcon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qptcon(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qptcon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif


#ifdef PETSC_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), dptrfs_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qptrfs(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qptrfs_(
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), dpttrf_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qpttrf(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qpttrf_(
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *), dpttrs_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *), qpttrs(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *), qpttrs_(
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, int *);



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define DF(I) df[(I)-1]
#define EF(I) ef[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    if (! nofact && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < MAX(1,*n)) {
	*info = -9;
    } else if (*ldx < MAX(1,*n)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPTSVX", &i__1);
	return;
    }

    if (nofact) {

/*        Compute the L*D*L' (or U'*D*U) factorization of A. */


#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &D(1), &c__1, &DF(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &D(1), &c__1, &DF(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &D(1), &c__1, &DF(1), &c__1);
#endif

	if (*n > 1) {
	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__1, &E(1), &c__1, &EF(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__1, &E(1), &c__1, &EF(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__1, &E(1), &c__1, &EF(1), &c__1);
#endif

	}

#ifdef PETSC_PREFIX_SUFFIX
	dpttrf_(n, &DF(1), &EF(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qpttrf(n, &DF(1), &EF(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qpttrf_(n, &DF(1), &EF(1), info);
#endif


/*        Return if INFO is non-zero. */

	if (*info != 0) {
	    if (*info > 0) {
		*rcond = 0.;
	    }
	    return;
	}
    }

/*     Compute the norm of the matrix A. */


#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlanst_("1", n, &D(1), &E(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlanst("1", n, &D(1), &E(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlanst_("1", n, &D(1), &E(1));
#endif


/*     Compute the reciprocal of the condition number of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dptcon_(n, &DF(1), &EF(1), &anorm, rcond, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qptcon(n, &DF(1), &EF(1), &anorm, rcond, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qptcon_(n, &DF(1), &EF(1), &anorm, rcond, &WORK(1), info);
#endif


/*     Return if the matrix is singular to working precision. */


#ifdef PETSC_PREFIX_SUFFIX
    if (*rcond < dlamch_("Epsilon")) {
#endif
#ifdef Q_C_PREFIX_SUFFIX
    if (*rcond < qlamch("Epsilon")) {
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    if (*rcond < qlamch_("Epsilon")) {
#endif

	*info = *n + 1;
	return;
    }

/*     Compute the solution vectors X. */


#ifdef PETSC_PREFIX_SUFFIX
    dlacpy_("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlacpy("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlacpy_("Full", n, nrhs, &B(1,1), ldb, &X(1,1), ldx);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    dpttrs_(n, nrhs, &DF(1), &EF(1), &X(1,1), ldx, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpttrs(n, nrhs, &DF(1), &EF(1), &X(1,1), ldx, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpttrs_(n, nrhs, &DF(1), &EF(1), &X(1,1), ldx, info);
#endif


/*     Use iterative refinement to improve the computed solutions and   
       compute error bounds and backward error estimates for them. */


#ifdef PETSC_PREFIX_SUFFIX
    dptrfs_(n, nrhs, &D(1), &E(1), &DF(1), &EF(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qptrfs(n, nrhs, &D(1), &E(1), &DF(1), &EF(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qptrfs_(n, nrhs, &D(1), &E(1), &DF(1), &EF(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), info);
#endif


    return;

/*     End of DPTSVX */

} /* dptsvx_ */

