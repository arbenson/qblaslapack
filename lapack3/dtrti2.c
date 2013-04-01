#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtrti2_(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtrti2(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtrti2_(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif

	a, int *lda, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DTRTI2 computes the inverse of a real upper or lower triangular   
    matrix.   

    This is the Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the matrix A is upper or lower triangular. 
  
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    DIAG    (input) CHARACTER*1   
            Specifies whether or not the matrix A is unit triangular.   
            = 'N':  Non-unit triangular   
            = 'U':  Unit triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular matrix A.  If UPLO = 'U', the   
            leading n by n upper triangular part of the array A contains 
  
            the upper triangular matrix, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n by n lower triangular part of the array A contains 
  
            the lower triangular matrix, and the strictly upper   
            triangular part of A is not referenced.  If DIAG = 'U', the   
            diagonal elements of A are also not referenced and are   
            assumed to be 1.   

            On exit, the (triangular) inverse of the original matrix, in 
  
            the same storage format.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
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

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), xerbla_(char *, int *);
    static long int nounit;
    static LONG DOUBLE ajj;




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    nounit = lsame_(diag, "N");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTRTI2", &i__1);
	return;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (nounit) {
		A(j,j) = 1. / A(j,j);
		ajj = -A(j,j);
	    } else {
		ajj = -1.;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dtrmv_("Upper", "No transpose", diag, &i__2, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrmv("Upper", "No transpose", diag, &i__2, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrmv_("Upper", "No transpose", diag, &i__2, &A(1,1), lda, &
#endif

		    A(1,j), &c__1);
	    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&i__2, &ajj, &A(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&i__2, &ajj, &A(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&i__2, &ajj, &A(1,j), &c__1);
#endif

/* L10: */
	}
    } else {

/*        Compute inverse of lower triangular matrix. */

	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		A(j,j) = 1. / A(j,j);
		ajj = -A(j,j);
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dtrmv_("Lower", "No transpose", diag, &i__1, &A(j+1,j+1), lda, &A(j+1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtrmv("Lower", "No transpose", diag, &i__1, &A(j+1,j+1), lda, &A(j+1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtrmv_("Lower", "No transpose", diag, &i__1, &A(j+1,j+1), lda, &A(j+1,j), &c__1);
#endif

		i__1 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__1, &ajj, &A(j+1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__1, &ajj, &A(j+1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__1, &ajj, &A(j+1,j), &c__1);
#endif

	    }
/* L20: */
	}
    }

    return;

/*     End of DTRTI2 */

} /* dtrti2_ */

