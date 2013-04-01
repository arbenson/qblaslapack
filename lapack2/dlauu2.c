#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlauu2_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlauu2(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlauu2_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *info)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAUU2 computes the product U * U' or L' * L, where the triangular   
    factor U or L is stored in the upper or lower triangular part of   
    the array A.   

    If UPLO = 'U' or 'u' then the upper triangle of the result is stored, 
  
    overwriting the factor U in A.   
    If UPLO = 'L' or 'l' then the lower triangle of the result is stored, 
  
    overwriting the factor L in A.   

    This is the unblocked form of the algorithm, calling Level 2 BLAS.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the triangular factor stored in the array A 
  
            is upper or lower triangular:   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the triangular factor U or L.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular factor U or L.   
            On exit, if UPLO = 'U', the upper triangle of A is   
            overwritten with the upper triangle of the product U * U';   
            if UPLO = 'L', the lower triangle of A is overwritten with   
            the lower triangle of the product L' * L.   

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
    static LONG DOUBLE c_b7 = 1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
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
    static int i;

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
    static LONG DOUBLE aii;




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
	xerbla_("DLAUU2", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (upper) {

/*        Compute the product U * U'. */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    aii = A(i,i);
	    if (i < *n) {
		i__2 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		A(i,i) = ddot_(&i__2, &A(i,i), lda, &A(i,i), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(i,i) = qdot(&i__2, &A(i,i), lda, &A(i,i), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(i,i) = qdot_(&i__2, &A(i,i), lda, &A(i,i), lda);
#endif

		i__2 = i - 1;
		i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("No transpose", &i__2, &i__3, &c_b7, &A(1,i+1), lda, &A(i,i+1), lda, &aii,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("No transpose", &i__2, &i__3, &c_b7, &A(1,i+1), lda, &A(i,i+1), lda, &aii,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("No transpose", &i__2, &i__3, &c_b7, &A(1,i+1), lda, &A(i,i+1), lda, &aii,
#endif

			 &A(1,i), &c__1);
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i, &aii, &A(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i, &aii, &A(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i, &aii, &A(1,i), &c__1);
#endif

	    }
/* L10: */
	}

    } else {

/*        Compute the product L' * L. */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    aii = A(i,i);
	    if (i < *n) {
		i__2 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		A(i,i) = ddot_(&i__2, &A(i,i), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(i,i) = qdot(&i__2, &A(i,i), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(i,i) = qdot_(&i__2, &A(i,i), &c__1, &
#endif

			A(i,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("Transpose", &i__2, &i__3, &c_b7, &A(i+1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("Transpose", &i__2, &i__3, &c_b7, &A(i+1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("Transpose", &i__2, &i__3, &c_b7, &A(i+1,1), 
#endif

			lda, &A(i+1,i), &c__1, &aii, &A(i,1), lda);
	    } else {

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i, &aii, &A(i,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i, &aii, &A(i,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i, &aii, &A(i,1), lda);
#endif

	    }
/* L20: */
	}
    }

    return;

/*     End of DLAUU2 */

} /* dlauu2_ */

