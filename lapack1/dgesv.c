#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgesv_(int *n, int *nrhs, LONG DOUBLE *a, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgesv(int *n, int *nrhs, LONG DOUBLE *a, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgesv_(int *n, int *nrhs, LONG DOUBLE *a, int 
#endif

	*lda, int *ipiv, LONG DOUBLE *b, int *ldb, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGESV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N matrix and X and B are N-by-NRHS matrices.   

    The LU decomposition with partial pivoting and row interchanges is   
    used to factor A as   
       A = P * L * U,   
    where P is a permutation matrix, L is unit lower triangular, and U is 
  
    upper triangular.  The factored form of A is then used to solve the   
    system of equations A * X = B.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of linear equations, i.e., the order of the   
            matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the N-by-N coefficient matrix A.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            The pivot indices that define the permutation matrix P;   
            row i of the matrix was interchanged with row IPIV(i).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the N-by-NRHS matrix of right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization 
  
                  has been completed, but the factor U is exactly   
                  singular, so the solution could not be computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgetrf_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgetrf(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgetrf_(int *, int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *, int *), xerbla_(char *, int *), dgetrs_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *, int *), xerbla_(char *, int *), qgetrs(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *, int *), xerbla_(char *, int *), qgetrs_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, int *);


#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    } else if (*ldb < MAX(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGESV ", &i__1);
	return;
    }

/*     Compute the LU factorization of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dgetrf_(n, n, &A(1,1), lda, &IPIV(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgetrf(n, n, &A(1,1), lda, &IPIV(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgetrf_(n, n, &A(1,1), lda, &IPIV(1), info);
#endif

    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dgetrs_("No transpose", n, nrhs, &A(1,1), lda, &IPIV(1), &B(1,1), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgetrs("No transpose", n, nrhs, &A(1,1), lda, &IPIV(1), &B(1,1), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgetrs_("No transpose", n, nrhs, &A(1,1), lda, &IPIV(1), &B(1,1), ldb, info);
#endif

    }
    return;

/*     End of DGESV */

} /* dgesv_ */

