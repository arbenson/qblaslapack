#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpotf2_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpotf2(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpotf2_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DPOTF2 computes the Cholesky factorization of a real symmetric   
    positive definite matrix A.   

    The factorization has the form   
       A = U' * U ,  if UPLO = 'U', or   
       A = L  * L',  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    This is the unblocked version of the algorithm, calling Level 2 BLAS. 
  

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored.   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n by n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the factor U or L from the Cholesky   
            factorization A = U'*U  or A = L*L'.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, the leading minor of order k is not   
                 positive definite, and the factorization could not be   
                 completed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b10 = -1.;
    static LONG DOUBLE c_b12 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE ddot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qdot(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qdot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    int *);
    static int j;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemv_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static long int upper;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE ajj;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPOTF2", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness
. */

	    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    ajj = A(j,j) - ddot_(&i__2, &A(1,j), &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    ajj = A(j,j) - qdot(&i__2, &A(1,j), &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    ajj = A(j,j) - qdot_(&i__2, &A(1,j), &c__1, 
#endif

		    &A(1,j), &c__1);
	    if (ajj <= 0.) {
		A(j,j) = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    A(j,j) = ajj;

/*           Compute elements J+1:N of row J. */

	    if (j < *n) {
		i__2 = j - 1;
		i__3 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("Transpose", &i__2, &i__3, &c_b10, &A(1,j+1), lda, &A(1,j), &c__1, &c_b12, &A(j,j+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("Transpose", &i__2, &i__3, &c_b10, &A(1,j+1), lda, &A(1,j), &c__1, &c_b12, &A(j,j+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("Transpose", &i__2, &i__3, &c_b10, &A(1,j+1), lda, &A(1,j), &c__1, &c_b12, &A(j,j+1), lda);
#endif

		i__2 = *n - j;
		d__1 = 1. / ajj;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &d__1, &A(j,j+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &d__1, &A(j,j+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &d__1, &A(j,j+1), lda);
#endif

	    }
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness
. */

	    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    ajj = A(j,j) - ddot_(&i__2, &A(j,1), lda, &A(j,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    ajj = A(j,j) - qdot(&i__2, &A(j,1), lda, &A(j,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    ajj = A(j,j) - qdot_(&i__2, &A(j,1), lda, &A(j,1), lda);
#endif

	    if (ajj <= 0.) {
		A(j,j) = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    A(j,j) = ajj;

/*           Compute elements J+1:N of column J. */

	    if (j < *n) {
		i__2 = *n - j;
		i__3 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("No transpose", &i__2, &i__3, &c_b10, &A(j+1,1), lda, &A(j,1), lda, &c_b12, &A(j+1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("No transpose", &i__2, &i__3, &c_b10, &A(j+1,1), lda, &A(j,1), lda, &c_b12, &A(j+1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("No transpose", &i__2, &i__3, &c_b10, &A(j+1,1), lda, &A(j,1), lda, &c_b12, &A(j+1,j), &c__1);
#endif

		i__2 = *n - j;
		d__1 = 1. / ajj;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &d__1, &A(j+1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &d__1, &A(j+1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &d__1, &A(j+1,j), &c__1);
#endif

	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return;

/*     End of DPOTF2 */

} /* dpotf2_ */

