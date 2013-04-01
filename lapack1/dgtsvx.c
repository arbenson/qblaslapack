#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgtsvx_(char *fact, char *trans, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgtsvx(char *fact, char *trans, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgtsvx_(char *fact, char *trans, int *n, int *
#endif

	nrhs, LONG DOUBLE *dl, LONG DOUBLE *d, LONG DOUBLE *du, LONG DOUBLE *dlf, 
	LONG DOUBLE *df, LONG DOUBLE *duf, LONG DOUBLE *du2, int *ipiv, 
	LONG DOUBLE *b, int *ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *
	rcond, LONG DOUBLE *ferr, LONG DOUBLE *berr, LONG DOUBLE *work, int *
	iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGTSVX uses the LU factorization to compute the solution to a real   
    system of linear equations A * X = B or A**T * X = B,   
    where A is a tridiagonal matrix of order N and X and B are N-by-NRHS 
  
    matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'N', the LU decomposition is used to factor the matrix A 
  
       as A = L * U, where L is a product of permutation and unit lower   
       bidiagonal matrices and U is upper triangular with nonzeros in   
       only the main diagonal and first two superdiagonals.   

    2. The factored form of A is used to estimate the condition number   
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
            = 'F':  DLF, DF, DUF, DU2, and IPIV contain the factored   
                    form of A; DL, D, DU, DLF, DF, DUF, DU2 and IPIV   
                    will not be modified.   
            = 'N':  The matrix will be copied to DLF, DF, and DUF   
                    and factored.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) subdiagonal elements of A.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of A.   

    DU      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) superdiagonal elements of A.   

    DLF     (input or output) LONG DOUBLE PRECISION array, dimension (N-1)   
            If FACT = 'F', then DLF is an input argument and on entry   
            contains the (n-1) multipliers that define the matrix L from 
  
            the LU factorization of A as computed by DGTTRF.   

            If FACT = 'N', then DLF is an output argument and on exit   
            contains the (n-1) multipliers that define the matrix L from 
  
            the LU factorization of A.   

    DF      (input or output) LONG DOUBLE PRECISION array, dimension (N)   
            If FACT = 'F', then DF is an input argument and on entry   
            contains the n diagonal elements of the upper triangular   
            matrix U from the LU factorization of A.   

            If FACT = 'N', then DF is an output argument and on exit   
            contains the n diagonal elements of the upper triangular   
            matrix U from the LU factorization of A.   

    DUF     (input or output) LONG DOUBLE PRECISION array, dimension (N-1)   
            If FACT = 'F', then DUF is an input argument and on entry   
            contains the (n-1) elements of the first superdiagonal of U. 
  

            If FACT = 'N', then DUF is an output argument and on exit   
            contains the (n-1) elements of the first superdiagonal of U. 
  

    DU2     (input or output) LONG DOUBLE PRECISION array, dimension (N-2)   
            If FACT = 'F', then DU2 is an input argument and on entry   
            contains the (n-2) elements of the second superdiagonal of   
            U.   

            If FACT = 'N', then DU2 is an output argument and on exit   
            contains the (n-2) elements of the second superdiagonal of   
            U.   

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains the pivot indices from the LU factorization of A as 
  
            computed by DGTTRF.   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains the pivot indices from the LU factorization of A;   
            row i of the matrix was interchanged with row IPIV(i).   
            IPIV(i) will always be either i or i+1; IPIV(i) = i indicates 
  
            a row interchange was not required.   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The N-by-NRHS right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            If INFO = 0, the N-by-NRHS solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= MAX(1,N).   

    RCOND   (output) LONG DOUBLE PRECISION   
            The estimate of the reciprocal condition number of the matrix 
  
            A.  If RCOND is less than the machine precision (in   
            particular, if RCOND = 0), the matrix is singular to working 
  
            precision.  This condition is indicated by a return code of   
            INFO > 0, and the solution and error bounds are not computed. 
  

    FERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (3*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, and i is   
                  <= N:  U(i,i) is exactly zero.  The factorization   
                         has not been completed unless i = N, but the   
                         factor U is exactly singular, so the solution   
                         and error bounds could not be computed.   
                 = N+1:  RCOND is less than machine precision.  The   
                         factorization has been completed, but the   
                         matrix is singular to working precision, and   
                         the solution and error bounds have not been   
                         computed.   

    ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static char norm[1];
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
    extern LONG DOUBLE dlamch_(char *), dlangt_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *), dlangt_(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *), dlangt_(char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *);
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

#ifdef PETSC_PREFIX_SUFFIX
	    xerbla_(char *, int *), dgtcon_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    xerbla_(char *, int *), qgtcon(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    xerbla_(char *, int *), qgtcon_(char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *,

#ifdef PETSC_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), dgtrfs_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), qgtrfs(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), qgtrfs_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     int *, int *), dgttrf_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *, int *), qgttrf(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *, int *), qgttrf_(int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, int *);
    static long int notran;

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



#define DL(I) dl[(I)-1]
#define D(I) d[(I)-1]
#define DU(I) du[(I)-1]
#define DLF(I) dlf[(I)-1]
#define DF(I) df[(I)-1]
#define DUF(I) duf[(I)-1]
#define DU2(I) du2[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    notran = lsame_(trans, "N");
    if (! nofact && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, 
	    "C")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*ldb < MAX(1,*n)) {
	*info = -14;
    } else if (*ldx < MAX(1,*n)) {
	*info = -16;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGTSVX", &i__1);
	return;
    }

    if (nofact) {

/*        Compute the LU factorization of A. */


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
	    dcopy_(&i__1, &DL(1), &c__1, &DLF(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__1, &DL(1), &c__1, &DLF(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__1, &DL(1), &c__1, &DLF(1), &c__1);
#endif

	    i__1 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__1, &DU(1), &c__1, &DUF(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__1, &DU(1), &c__1, &DUF(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__1, &DU(1), &c__1, &DUF(1), &c__1);
#endif

	}

#ifdef PETSC_PREFIX_SUFFIX
	dgttrf_(n, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgttrf(n, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgttrf_(n, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), info);
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

    if (notran) {
	*(unsigned char *)norm = '1';
    } else {
	*(unsigned char *)norm = 'I';
    }

#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlangt_(norm, n, &DL(1), &D(1), &DU(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlangt(norm, n, &DL(1), &D(1), &DU(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlangt_(norm, n, &DL(1), &D(1), &DU(1));
#endif


/*     Compute the reciprocal of the condition number of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dgtcon_(norm, n, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), &anorm, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgtcon(norm, n, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), &anorm, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgtcon_(norm, n, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), &anorm, 
#endif

	    rcond, &WORK(1), &IWORK(1), info);

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
    dgttrs_(trans, n, nrhs, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), &X(1,1), ldx, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgttrs(trans, n, nrhs, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), &X(1,1), ldx, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgttrs_(trans, n, nrhs, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(1), &X(1,1), ldx, info);
#endif


/*     Use iterative refinement to improve the computed solutions and   
       compute error bounds and backward error estimates for them. */


#ifdef PETSC_PREFIX_SUFFIX
    dgtrfs_(trans, n, nrhs, &DL(1), &D(1), &DU(1), &DLF(1), &DF(1), &DUF(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgtrfs(trans, n, nrhs, &DL(1), &D(1), &DU(1), &DLF(1), &DF(1), &DUF(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgtrfs_(trans, n, nrhs, &DL(1), &D(1), &DU(1), &DLF(1), &DF(1), &DUF(1), &
#endif

	    DU2(1), &IPIV(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), 
	    &BERR(1), &WORK(1), &IWORK(1), info);

    return;

/*     End of DGTSVX */

} /* dgtsvx_ */

