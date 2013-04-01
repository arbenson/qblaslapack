#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtptrs_(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtptrs(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtptrs_(char *uplo, char *trans, char *diag, int *n, 
#endif

	int *nrhs, LONG DOUBLE *ap, LONG DOUBLE *b, int *ldb, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DTPTRS solves a triangular system of the form   

       A * X = B  or  A**T * X = B,   

    where A is a triangular matrix of order N stored in packed format,   
    and B is an N-by-NRHS matrix.  A check is made to verify that A is   
    nonsingular.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose = Transpose)   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The upper or lower triangular matrix A, packed columnwise in 
  
            a linear array.  The j-th column of A is stored in the array 
  
            AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2*n-j)/2) = A(i,j) for j<=i<=n. 
  

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, if INFO = 0, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the i-th diagonal element of A is zero,   
                  indicating that the matrix is singular and the   
                  solutions X have not been computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static int j;
    extern long int lsame_(char *, char *);
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtpsv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpsv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpsv_(char *, char *, char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static int jc;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int nounit;



#define AP(I) ap[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    nounit = lsame_(diag, "N");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*nrhs < 0) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTPTRS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Check for singularity. */

    if (nounit) {
	if (upper) {
	    jc = 1;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		if (AP(jc + *info - 1) == 0.) {
		    return;
		}
		jc += *info;
/* L10: */
	    }
	} else {
	    jc = 1;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		if (AP(jc) == 0.) {
		    return;
		}
		jc = jc + *n - *info + 1;
/* L20: */
	    }
	}
    }
    *info = 0;

/*     Solve A * x = b  or  A' * x = b. */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
	dtpsv_(uplo, trans, diag, n, &AP(1), &B(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtpsv(uplo, trans, diag, n, &AP(1), &B(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtpsv_(uplo, trans, diag, n, &AP(1), &B(1,j), &c__1);
#endif

/* L30: */
    }

    return;

/*     End of DTPTRS */

} /* dtptrs_ */

