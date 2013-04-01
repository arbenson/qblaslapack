#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dppsv_(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qppsv(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qppsv_(char *uplo, int *n, int *nrhs, LONG DOUBLE 
#endif

	*ap, LONG DOUBLE *b, int *ldb, int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPPSV computes the solution to a real system of linear equations   
       A * X = B,   
    where A is an N-by-N symmetric positive definite matrix stored in   
    packed format and X and B are N-by-NRHS matrices.   

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

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            See below for further details.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T, in the same storage   
            format as A.   

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

    Further Details   
    ===============   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the symmetric matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = conjg(aji))   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dpptrf_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpptrf(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpptrf_(
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    char *, int *, LONG DOUBLE *, int *), dpptrs_(char 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    char *, int *, LONG DOUBLE *, int *), qpptrs(char 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    char *, int *, LONG DOUBLE *, int *), qpptrs_(char 
#endif

	    *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    int *);


#define AP(I) ap[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < MAX(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPPSV ", &i__1);
	return;
    }

/*     Compute the Cholesky factorization A = U'*U or A = L*L'. */


#ifdef PETSC_PREFIX_SUFFIX
    dpptrf_(uplo, n, &AP(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpptrf(uplo, n, &AP(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpptrf_(uplo, n, &AP(1), info);
#endif

    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dpptrs_(uplo, n, nrhs, &AP(1), &B(1,1), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qpptrs(uplo, n, nrhs, &AP(1), &B(1,1), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qpptrs_(uplo, n, nrhs, &AP(1), &B(1,1), ldb, info);
#endif


    }
    return;

/*     End of DPPSV */

} /* dppsv_ */

