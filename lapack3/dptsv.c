#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dptsv_(int *n, int *nrhs, LONG DOUBLE *d, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qptsv(int *n, int *nrhs, LONG DOUBLE *d, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qptsv_(int *n, int *nrhs, LONG DOUBLE *d, 
#endif

	LONG DOUBLE *e, LONG DOUBLE *b, int *ldb, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPTSV computes the solution to a real system of linear equations   
    A*X = B, where A is an N-by-N symmetric positive definite tridiagonal 
  
    matrix, and X and B are N-by-NRHS matrices.   

    A is factored as A = L*D*L**T, and the factored form of A is then   
    used to solve the system of equations.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal matrix   
            A.  On exit, the n diagonal elements of the diagonal matrix   
            D from the factorization A = L*D*L**T.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix A.  On exit, the (n-1) subdiagonal elements of the   
            unit bidiagonal factor L from the L*D*L**T factorization of   
            A.  (E can also be regarded as the superdiagonal of the unit 
  
            bidiagonal factor U from the U**T*D*U factorization of A.)   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the N-by-NRHS right hand side matrix B.   
            On exit, if INFO = 0, the N-by-NRHS solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the solution has not been   
                  computed.  The factorization has not been completed   
                  unless i = N.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dpttrf_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpttrf(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpttrf_(
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

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < MAX(1,*n)) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPTSV ", &i__1);
	return;
    }

/*     Compute the L*D*L' (or U'*D*U) factorization of A. */


#ifdef PETSC_PREFIX_SUFFIX
    dpttrf_(n, &D(1), &E(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpttrf(n, &D(1), &E(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpttrf_(n, &D(1), &E(1), info);
#endif

    if (*info == 0) {

/*        Solve the system A*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dpttrs_(n, nrhs, &D(1), &E(1), &B(1,1), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qpttrs(n, nrhs, &D(1), &E(1), &B(1,1), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qpttrs_(n, nrhs, &D(1), &E(1), &B(1,1), ldb, info);
#endif

    }
    return;

/*     End of DPTSV */

} /* dptsv_ */

