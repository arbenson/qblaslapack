#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsygs2_(int *itype, char *uplo, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsygs2(int *itype, char *uplo, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsygs2_(int *itype, char *uplo, int *n, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *ldb, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DSYGS2 reduces a real symmetric-definite generalized eigenproblem   
    to standard form.   

    If ITYPE = 1, the problem is A*x = lambda*B*x,   
    and A is overwritten by inv(U')*A*inv(U) or inv(L)*A*inv(L')   

    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or   
    B*A*x = lambda*x, and A is overwritten by U*A*U` or L'*A*L.   

    B must have been previously factorized as U'*U or L*L' by DPOTRF.   

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            = 1: compute inv(U')*A*inv(U) or inv(L)*A*inv(L');   
            = 2 or 3: compute U*A*U' or L'*A*L.   

    UPLO    (input) CHARACTER   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored, and how B has been factorized. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n by n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the transformed matrix, stored in the   
            same format as A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            The triangular factor from the Cholesky factorization of B,   
            as returned by DPOTRF.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b6 = -1.;
    static int c__1 = 1;
    static LONG DOUBLE c_b27 = 1.;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsyr2_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr2(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr2_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static int k;

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
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtrmv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmv_(char *, char *, char *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), dtrsv_(char *, char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtrsv(char *, char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtrsv_(char *, char *, char *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static LONG DOUBLE ct;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE akk, bkk;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYGS2", &i__1);
	return;
    }

    if (*itype == 1) {
	if (upper) {

/*           Compute inv(U')*A*inv(U) */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the upper triangle of A(k:n,k:n) */

		akk = A(k,k);
		bkk = B(k,k);
/* Computing 2nd power */
		d__1 = bkk;
		akk /= d__1 * d__1;
		A(k,k) = akk;
		if (k < *n) {
		    i__2 = *n - k;
		    d__1 = 1. / bkk;

#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(&i__2, &d__1, &A(k,k+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(&i__2, &d__1, &A(k,k+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(&i__2, &d__1, &A(k,k+1), lda);
#endif

		    ct = akk * -.5;
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
#endif

		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    dsyr2_(uplo, &i__2, &c_b6, &A(k,k+1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsyr2(uplo, &i__2, &c_b6, &A(k,k+1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsyr2_(uplo, &i__2, &c_b6, &A(k,k+1), lda, 
#endif

			    &B(k,k+1), ldb, &A(k+1,k+1), lda);
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&i__2, &ct, &B(k,k+1), ldb, &A(k,k+1), lda);
#endif

		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrsv_(uplo, "Transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k,k+1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrsv(uplo, "Transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k,k+1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrsv_(uplo, "Transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k,k+1), 
#endif

			    lda);
		}
/* L10: */
	    }
	} else {

/*           Compute inv(L)*A*inv(L') */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the lower triangle of A(k:n,k:n) */

		akk = A(k,k);
		bkk = B(k,k);
/* Computing 2nd power */
		d__1 = bkk;
		akk /= d__1 * d__1;
		A(k,k) = akk;
		if (k < *n) {
		    i__2 = *n - k;
		    d__1 = 1. / bkk;

#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(&i__2, &d__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(&i__2, &d__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(&i__2, &d__1, &A(k+1,k), &c__1);
#endif

		    ct = akk * -.5;
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
#endif

		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    dsyr2_(uplo, &i__2, &c_b6, &A(k+1,k), &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsyr2(uplo, &i__2, &c_b6, &A(k+1,k), &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsyr2_(uplo, &i__2, &c_b6, &A(k+1,k), &c__1, 
#endif

			    &B(k+1,k), &c__1, &A(k+1,k+1), lda);
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&i__2, &ct, &B(k+1,k), &c__1, &A(k+1,k), &c__1);
#endif

		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrsv_(uplo, "No transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k+1,k), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrsv(uplo, "No transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k+1,k), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrsv_(uplo, "No transpose", "Non-unit", &i__2, &B(k+1,k+1), ldb, &A(k+1,k), 
#endif

			    &c__1);
		}
/* L20: */
	    }
	}
    } else {
	if (upper) {

/*           Compute U*A*U' */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the upper triangle of A(1:k,1:k) */

		akk = A(k,k);
		bkk = B(k,k);
		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtrmv_(uplo, "No transpose", "Non-unit", &i__2, &B(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmv(uplo, "No transpose", "Non-unit", &i__2, &B(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmv_(uplo, "No transpose", "Non-unit", &i__2, &B(1,1), 
#endif

			ldb, &A(1,k), &c__1);
		ct = akk * .5;
		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dsyr2_(uplo, &i__2, &c_b27, &A(1,k), &c__1, &B(1,k), &c__1, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsyr2(uplo, &i__2, &c_b27, &A(1,k), &c__1, &B(1,k), &c__1, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsyr2_(uplo, &i__2, &c_b27, &A(1,k), &c__1, &B(1,k), &c__1, &A(1,1), lda);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &ct, &B(1,k), &c__1, &A(1,k), &c__1);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &bkk, &A(1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &bkk, &A(1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &bkk, &A(1,k), &c__1);
#endif

/* Computing 2nd power */
		d__1 = bkk;
		A(k,k) = akk * (d__1 * d__1);
/* L30: */
	    }
	} else {

/*           Compute L'*A*L */

	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {

/*              Update the lower triangle of A(1:k,1:k) */

		akk = A(k,k);
		bkk = B(k,k);
		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtrmv_(uplo, "Transpose", "Non-unit", &i__2, &B(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmv(uplo, "Transpose", "Non-unit", &i__2, &B(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmv_(uplo, "Transpose", "Non-unit", &i__2, &B(1,1), 
#endif

			ldb, &A(k,1), lda);
		ct = akk * .5;
		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dsyr2_(uplo, &i__2, &c_b27, &A(k,1), lda, &B(k,1), ldb, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsyr2(uplo, &i__2, &c_b27, &A(k,1), lda, &B(k,1), ldb, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsyr2_(uplo, &i__2, &c_b27, &A(k,1), lda, &B(k,1), ldb, &A(1,1), lda);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &ct, &B(k,1), ldb, &A(k,1), lda);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &bkk, &A(k,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &bkk, &A(k,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &bkk, &A(k,1), lda);
#endif

/* Computing 2nd power */
		d__1 = bkk;
		A(k,k) = akk * (d__1 * d__1);
/* L40: */
	    }
	}
    }
    return;

/*     End of DSYGS2 */

} /* dsygs2_ */

