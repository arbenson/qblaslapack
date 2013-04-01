#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsysvx_(char *fact, char *uplo, int *n, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsysvx(char *fact, char *uplo, int *n, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsysvx_(char *fact, char *uplo, int *n, int *
#endif

	nrhs, LONG DOUBLE *a, int *lda, LONG DOUBLE *af, int *ldaf, 
	int *ipiv, LONG DOUBLE *b, int *ldb, LONG DOUBLE *x, int *
	ldx, LONG DOUBLE *rcond, LONG DOUBLE *ferr, LONG DOUBLE *berr, 
	LONG DOUBLE *work, int *lwork, int *iwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYSVX uses the diagonal pivoting factorization to compute the   
    solution to a real system of linear equations A * X = B,   
    where A is an N-by-N symmetric matrix and X and B are N-by-NRHS   
    matrices.   

    Error bounds on the solution and a condition estimate are also   
    provided.   

    Description   
    ===========   

    The following steps are performed:   

    1. If FACT = 'N', the diagonal pivoting method is used to factor A.   
       The form of the factorization is   
          A = U * D * U**T,  if UPLO = 'U', or   
          A = L * D * L**T,  if UPLO = 'L',   
       where U (or L) is a product of permutation and unit upper (lower) 
  
       triangular matrices, and D is symmetric and block diagonal with   
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
            = 'F':  On entry, AF and IPIV contain the factored form of   
                    A.  AF and IPIV will not be modified.   
            = 'N':  The matrix A will be copied to AF and factored.   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            The symmetric matrix A.  If UPLO = 'U', the leading N-by-N   
            upper triangular part of A contains the upper triangular part 
  
            of the matrix A, and the strictly lower triangular part of A 
  
            is not referenced.  If UPLO = 'L', the leading N-by-N lower   
            triangular part of A contains the lower triangular part of   
            the matrix A, and the strictly upper triangular part of A is 
  
            not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    AF      (input or output) LONG DOUBLE PRECISION array, dimension (LDAF,N) 
  
            If FACT = 'F', then AF is an input argument and on entry   
            contains the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T as computed by DSYTRF.   

            If FACT = 'N', then AF is an output argument and on exit   
            returns the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L from the factorization   
            A = U*D*U**T or A = L*D*L**T.   

    LDAF    (input) INTEGER   
            The leading dimension of the array AF.  LDAF >= MAX(1,N).   

    IPIV    (input or output) INTEGER array, dimension (N)   
            If FACT = 'F', then IPIV is an input argument and on entry   
            contains details of the interchanges and the block structure 
  
            of D, as determined by DSYTRF.   
            If IPIV(k) > 0, then rows and columns k and IPIV(k) were   
            interchanged and D(k,k) is a 1-by-1 diagonal block.   
            If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and   
            columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) 
  
            is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =   
            IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were   
            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.   

            If FACT = 'N', then IPIV is an output argument and on exit   
            contains details of the interchanges and the block structure 
  
            of D, as determined by DSYTRF.   

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

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of WORK.  LWORK >= 3*N, and for best performance   
            LWORK >= N*NB, where NB is the optimal blocksize for   
            DSYTRF.   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, and i is   
                 <= N: D(i,i) is exactly zero.  The factorization has   
                       has been completed, but the block diagonal   
                       matrix D is exactly singular, so the solution and 
  
                       error bounds could not be computed.   
                 = N+1: the block diagonal matrix D is nonsingular, but   
                       RCOND is less than machine precision.  The   
                       factorization has been completed, but the matrix   
                       is singular to working precision, so the solution 
  
                       and error bounds have not been computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    extern long int lsame_(char *, char *);
    static LONG DOUBLE anorm;

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
    extern LONG DOUBLE dlansy_(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlansy(char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlansy_(char *, char *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsycon_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsycon(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsycon_(char *, int *, LONG DOUBLE *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *), dsyrfs_(char *, int *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *), qsyrfs(char *, int *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *), qsyrfs_(char *, int *, int 
#endif

	    *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dsytrf_(char *, int *, LONG DOUBLE *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsytrf(char *, int *, LONG DOUBLE *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsytrf_(char *, int *, LONG DOUBLE *, int *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), dsytrs_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), qsytrs(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, int *), qsytrs_(char *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, int *, 
	    LONG DOUBLE *, int *, int *);


#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define AF(I,J) af[(I)-1 + ((J)-1)* ( *ldaf)]
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
    } else if (*lda < MAX(1,*n)) {
	*info = -6;
    } else if (*ldaf < MAX(1,*n)) {
	*info = -8;
    } else if (*ldb < MAX(1,*n)) {
	*info = -11;
    } else if (*ldx < MAX(1,*n)) {
	*info = -13;
    } else if (*lwork < *n << 1) {
	*info = -18;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYSVX", &i__1);
	return;
    }

    if (nofact) {

/*        Compute the factorization A = U*D*U' or A = L*D*L'. */


#ifdef PETSC_PREFIX_SUFFIX
	dlacpy_(uplo, n, n, &A(1,1), lda, &AF(1,1), ldaf);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacpy(uplo, n, n, &A(1,1), lda, &AF(1,1), ldaf);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacpy_(uplo, n, n, &A(1,1), lda, &AF(1,1), ldaf);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dsytrf_(uplo, n, &AF(1,1), ldaf, &IPIV(1), &WORK(1), lwork, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsytrf(uplo, n, &AF(1,1), ldaf, &IPIV(1), &WORK(1), lwork, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsytrf_(uplo, n, &AF(1,1), ldaf, &IPIV(1), &WORK(1), lwork, 
#endif

		info);

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
    anorm = dlansy_("I", uplo, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlansy("I", uplo, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlansy_("I", uplo, n, &A(1,1), lda, &WORK(1));
#endif


/*     Compute the reciprocal of the condition number of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dsycon_(uplo, n, &AF(1,1), ldaf, &IPIV(1), &anorm, rcond, &WORK(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsycon(uplo, n, &AF(1,1), ldaf, &IPIV(1), &anorm, rcond, &WORK(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsycon_(uplo, n, &AF(1,1), ldaf, &IPIV(1), &anorm, rcond, &WORK(1), 
#endif

	    &IWORK(1), info);

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
    dsytrs_(uplo, n, nrhs, &AF(1,1), ldaf, &IPIV(1), &X(1,1), ldx, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsytrs(uplo, n, nrhs, &AF(1,1), ldaf, &IPIV(1), &X(1,1), ldx, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsytrs_(uplo, n, nrhs, &AF(1,1), ldaf, &IPIV(1), &X(1,1), ldx, 
#endif

	    info);

/*     Use iterative refinement to improve the computed solutions and   
       compute error bounds and backward error estimates for them. */


#ifdef PETSC_PREFIX_SUFFIX
    dsyrfs_(uplo, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &IPIV(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsyrfs(uplo, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &IPIV(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsyrfs_(uplo, n, nrhs, &A(1,1), lda, &AF(1,1), ldaf, &IPIV(1), 
#endif

	    &B(1,1), ldb, &X(1,1), ldx, &FERR(1), &BERR(1), &WORK(1)
	    , &IWORK(1), info);

    return;

/*     End of DSYSVX */

} /* dsysvx_ */

