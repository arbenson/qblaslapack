#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dposv_(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qposv(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qposv_(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif

	*a, int *lda, LONG DOUBLE *b, int *ldb, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPOSV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N symmetric positive definite matrix and X and B   
    are N-by-NRHS matrices.   

    The Cholesky decomposition is used to factor A as   
       A = U**T* U,  if UPLO = 'U', or   
       A = L * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is a lower triangular   
    matrix.  The factored form of A is then used to solve the system of   
    equations A * X = B.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i of A is not   
                  positive definite, so the factorization could not be   
                  completed, and the solution has not been computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dpotrf_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpotrf(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpotrf_(
#endif

	    char *, int *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dpotrs_(char *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qpotrs(char *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qpotrs_(char *, int *, int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *, int *);



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPOSV ", &i__1);
	return;
    }

/*     Compute the Cholesky factorization A = U'*U or A = L*L'. */


#ifdef PETSC_PREFIX_SUFFIX
    dpotrf_(uplo, n, &A(1,1), lda, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpotrf(uplo, n, &A(1,1), lda, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpotrf_(uplo, n, &A(1,1), lda, info);
#endif

    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dpotrs_(uplo, n, nrhs, &A(1,1), lda, &B(1,1), ldb, info)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qpotrs(uplo, n, nrhs, &A(1,1), lda, &B(1,1), ldb, info)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qpotrs_(uplo, n, nrhs, &A(1,1), lda, &B(1,1), ldb, info)
#endif

		;

    }
    return;

/*     End of DPOSV */

} /* dposv_ */

