#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsysv_(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsysv(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsysv_(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif

	*a, int *lda, int *ipiv, LONG DOUBLE *b, int *ldb, 
	LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYSV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N symmetric matrix and X and B are N-by-NRHS   
    matrices.   

    The diagonal pivoting method is used to factor A as   
       A = U * D * U**T,  if UPLO = 'U', or   
       A = L * D * L**T,  if UPLO = 'L',   
    where U (or L) is a product of permutation and unit upper (lower)   
    triangular matrices, and D is symmetric and block diagonal with   
    1-by-1 and 2-by-2 diagonal blocks.  The factored form of A is then   
    used to solve the system of equations A * X = B.   

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

            On exit, if INFO = 0, the block diagonal matrix D and the   
            multipliers used to obtain the factor U or L from the   
            factorization A = U*D*U**T or A = L*D*L**T as computed by   
            DSYTRF.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D, as 
  
            determined by DSYTRF.  If IPIV(k) > 0, then rows and columns 
  
            k and IPIV(k) were interchanged, and D(k,k) is a 1-by-1   
            diagonal block.  If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0,   
            then rows and columns k-1 and -IPIV(k) were interchanged and 
  
            D(k-1:k,k-1:k) is a 2-by-2 diagonal block.  If UPLO = 'L' and 
  
            IPIV(k) = IPIV(k+1) < 0, then rows and columns k+1 and   
            -IPIV(k) were interchanged and D(k:k+1,k:k+1) is a 2-by-2   
            diagonal block.   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of WORK.  LWORK >= 1, and for best performance   
            LWORK >= N*NB, where NB is the optimal blocksize for   
            DSYTRF.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, D(i,i) is exactly zero.  The factorization 
  
                 has been completed, but the block diagonal matrix D is   
                 exactly singular, so the solution could not be computed. 
  

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dsytrf_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qsytrf(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qsytrf_(
#endif

	    char *, int *, LONG DOUBLE *, int *, int *, LONG DOUBLE 

#ifdef PETSC_PREFIX_SUFFIX
	    *, int *, int *), dsytrs_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    *, int *, int *), qsytrs(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    *, int *, int *), qsytrs_(char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, 
	    int *, int *);


#define IPIV(I) ipiv[(I)-1]
#define WORK(I) work[(I)-1]

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
	*info = -8;
    } else if (*lwork < 1) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYSV ", &i__1);
	return;
    }

/*     Compute the factorization A = U*D*U' or A = L*D*L'. */


#ifdef PETSC_PREFIX_SUFFIX
    dsytrf_(uplo, n, &A(1,1), lda, &IPIV(1), &WORK(1), lwork, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qsytrf(uplo, n, &A(1,1), lda, &IPIV(1), &WORK(1), lwork, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qsytrf_(uplo, n, &A(1,1), lda, &IPIV(1), &WORK(1), lwork, info);
#endif

    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dsytrs_(uplo, n, nrhs, &A(1,1), lda, &IPIV(1), &B(1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsytrs(uplo, n, nrhs, &A(1,1), lda, &IPIV(1), &B(1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsytrs_(uplo, n, nrhs, &A(1,1), lda, &IPIV(1), &B(1,1), ldb,
#endif

		 info);

    }
    return;

/*     End of DSYSV */

} /* dsysv_ */

