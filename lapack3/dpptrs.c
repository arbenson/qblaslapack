#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpptrs_(char *uplo, int *n, int *nrhs, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpptrs(char *uplo, int *n, int *nrhs, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpptrs_(char *uplo, int *n, int *nrhs, 
#endif

	LONG DOUBLE *ap, LONG DOUBLE *b, int *ldb, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPPTRS solves a system of linear equations A*X = B with a symmetric   
    positive definite matrix A in packed storage using the Cholesky   
    factorization A = U**T*U or A = L*L**T computed by DPPTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The triangular factor U or L from the Cholesky factorization 
  
            A = U**T*U or A = L*L**T, packed columnwise in a linear   
            array.  The j-th column of U or L is stored in the array AP   
            as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static int i;
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

	    LONG DOUBLE *, LONG DOUBLE *, int *), 
	    xerbla_(char *, int *);



#define AP(I) ap[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
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
	xerbla_("DPPTRS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return;
    }

    if (upper) {

/*        Solve A*X = B where A = U'*U. */

	i__1 = *nrhs;
	for (i = 1; i <= *nrhs; ++i) {

/*           Solve U'*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtpsv_("Upper", "Transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtpsv("Upper", "Transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtpsv_("Upper", "Transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif


/*           Solve U*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtpsv_("Upper", "No transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtpsv("Upper", "No transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtpsv_("Upper", "No transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif

/* L10: */
	}
    } else {

/*        Solve A*X = B where A = L*L'. */

	i__1 = *nrhs;
	for (i = 1; i <= *nrhs; ++i) {

/*           Solve L*Y = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtpsv_("Lower", "No transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtpsv("Lower", "No transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtpsv_("Lower", "No transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif


/*           Solve L'*X = Y, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtpsv_("Lower", "Transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtpsv("Lower", "Transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtpsv_("Lower", "Transpose", "Non-unit", n, &AP(1), &B(1,i), &c__1);
#endif

/* L20: */
	}
    }

    return;

/*     End of DPPTRS */

} /* dpptrs_ */

