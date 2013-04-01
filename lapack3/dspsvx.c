#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspsvx_(char *fact, char *uplo, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspsvx(char *fact, char *uplo, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspsvx_(char *fact, char *uplo, int *n, int *
#endif

	nrhs, LONG DOUBLE *ap, LONG DOUBLE *afp, int *ipiv, LONG DOUBLE *b, 
	int *ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *rcond, 
	LONG DOUBLE *ferr, LONG DOUBLE *berr, LONG DOUBLE *work, int *iwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSPSVX uses the diagonal pivoting factorization A = U*D*U**T or   
    A = L*D*L**T to compute the solution to a real system of linear   
    equations A * X = B, where A is an N-by-N symmetric matrix stored   
    in packed format and X and B are N-by-NRHS matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'N', the diagonal pivoting method is used to factor A as 
  
          A = U * D * U**T,  if UPLO = 'U', or   
          A = L * D * L**T,  if UPLO = 'L',   
       where U (or L) is a product of permutation and unit upper (lower) 
  
       triangular matrices and D is symmetric and block diagonal with   
       1-by-1 and 2-by-2 diagonal blocks.   

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
            = 'F':  On entry, AFP and IPIV contain the factored form of   
                    A.  AP, AFP and IPIV will not be modified.   
            = 'N':  The matrix A will be copied to AFP and factored.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The upper or lower triangle of the symmetric matrix A, packed 
  
            columnwise in a linear array.  The j-th column of A is stored 
  
            in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. 
  
            See below for further details.   

    AFP     (input or output) LONG DOUBLE PRECISION array, dimension   
                              (N*(N+1)/2)   
            If FACT = 'F', then AFP is an input argument and on entry   
            contains the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T as computed by DSPTRF, stored as 
  
            a packed triangular matrix in the same storage format as A.   

            If FACT = 'N', then AFP is an output argument and on exit   
            contains the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T as computed by DSPTRF, stored as 
  
            a packed triangular matrix in the same storage format as A.   

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains details of the interchanges and the block structure 
  
            of D, as determined by DSPTRF.   
            If IPIV(k) > 0, then rows and columns k and IPIV(k) were   
            interchanged and D(k,k) is a 1-by-1 diagonal block.   
            If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and   
            columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) 
  
            is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =   
            IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were   
            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains details of the interchanges and the block structure 
  
            of D, as determined by DSPTRF.   

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
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0 and <= N: if INFO = i, D(i,i) is exactly zero.  The   
                 factorization has been completed, but the block diagonal 
  
                 matrix D is exactly singular, so the solution and error 
  
                 bounds could not be computed.   
            = N+1: the block diagonal matrix D is nonsingular, but RCOND 
  
                 is less than machine precision.  The factorization has   
                 been completed, but the matrix is singular to working   
                 precision, so the solution and error bounds have not   
                 been computed.   

    Further Details   
    ===============   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the symmetric matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = aji)   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

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
    extern LONG DOUBLE dlansp_(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlansp(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlansp_(char *, char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dspcon_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qspcon(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qspcon_(char *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), dsprfs_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qsprfs(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qsprfs_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,

#ifdef PETSC_PREFIX_SUFFIX
	     int *, int *), dsptrf_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	     int *, int *), qsptrf(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	     int *, int *), qsptrf_(char *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), dsptrs_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), qsptrs(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), qsptrs_(char *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, int *);



#define AP(I) ap[(I)-1]
#define AFP(I) afp[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    nofact = lsame_(fact, "N");
    if (! nofact && ! lsame_(fact, "F")) {
	*info = -1;
    } else if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*ldb < MAX(1,*n)) {
	*info = -9;
    } else if (*ldx < MAX(1,*n)) {
	*info = -11;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPSVX", &i__1);
	return;
    }

    if (nofact) {

/*        Compute the factorization A = U*D*U' or A = L*D*L'. */

	i__1 = *n * (*n + 1) / 2;

#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(&i__1, &AP(1), &c__1, &AFP(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(&i__1, &AP(1), &c__1, &AFP(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(&i__1, &AP(1), &c__1, &AFP(1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dsptrf_(uplo, n, &AFP(1), &IPIV(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsptrf(uplo, n, &AFP(1), &IPIV(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsptrf_(uplo, n, &AFP(1), &IPIV(1), info);
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
    anorm = dlansp_("I", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlansp("I", uplo, n, &AP(1), &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlansp_("I", uplo, n, &AP(1), &WORK(1));
#endif


/*     Compute the reciprocal of the condition number of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dspcon_(uplo, n, &AFP(1), &IPIV(1), &anorm, rcond, &WORK(1), &IWORK(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qspcon(uplo, n, &AFP(1), &IPIV(1), &anorm, rcond, &WORK(1), &IWORK(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qspcon_(uplo, n, &AFP(1), &IPIV(1), &anorm, rcond, &WORK(1), &IWORK(1), 
#endif

	    info);

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
    dsptrs_(uplo, n, nrhs, &AFP(1), &IPIV(1), &X(1,1), ldx, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsptrs(uplo, n, nrhs, &AFP(1), &IPIV(1), &X(1,1), ldx, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsptrs_(uplo, n, nrhs, &AFP(1), &IPIV(1), &X(1,1), ldx, info);
#endif


/*     Use iterative refinement to improve the computed solutions and   
       compute error bounds and backward error estimates for them. */


#ifdef PETSC_PREFIX_SUFFIX
    dsprfs_(uplo, n, nrhs, &AP(1), &AFP(1), &IPIV(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), &IWORK(1), info)
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsprfs(uplo, n, nrhs, &AP(1), &AFP(1), &IPIV(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), &IWORK(1), info)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsprfs_(uplo, n, nrhs, &AP(1), &AFP(1), &IPIV(1), &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1), &IWORK(1), info)
#endif

	    ;

    return;

/*     End of DSPSVX */

} /* dspsvx_ */

