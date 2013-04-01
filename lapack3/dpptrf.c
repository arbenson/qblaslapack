#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpptrf_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpptrf(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpptrf_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif

	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPPTRF computes the Cholesky factorization of a real symmetric   
    positive definite matrix A stored in packed format.   

    The factorization has the form   
       A = U**T * U,  if UPLO = 'U', or   
       A = L  * L**T,  if UPLO = 'L',   
    where U is an upper triangular matrix and L is lower triangular.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   
            See below for further details.   

            On exit, if INFO = 0, the triangular factor U or L from the   
            Cholesky factorization A = U**T*U or A = L*L**T, in the same 
  
            storage format as A.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the leading minor of order i is not   
                  positive definite, and the factorization could not be   
                  completed.   

    Further Details   
    ======= =======   

    The packed storage scheme is illustrated by the following example   
    when N = 4, UPLO = 'U':   

    Two-dimensional storage of the symmetric matrix A:   

       a11 a12 a13 a14   
           a22 a23 a24   
               a33 a34     (aij = aji)   
                   a44   

    Packed storage of the upper triangle of A:   

    AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b16 = -1.;
    
    /* System generated locals */
    int i__1, i__2;
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

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dspr_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qspr(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qspr_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *);
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
    extern /* Subroutine */ void dtpsv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpsv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpsv_(char *, char *, char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static int jc, jj;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE ajj;



#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPPTRF", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	jj = 0;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    jc = jj + 1;
	    jj += j;

/*           Compute elements 1:J-1 of column J. */

	    if (j > 1) {
		i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtpsv_("Upper", "Transpose", "Non-unit", &i__2, &AP(1), &AP(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtpsv("Upper", "Transpose", "Non-unit", &i__2, &AP(1), &AP(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtpsv_("Upper", "Transpose", "Non-unit", &i__2, &AP(1), &AP(
#endif

			jc), &c__1);
	    }

/*           Compute U(J,J) and test for non-positive-definiteness
. */

	    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    ajj = AP(jj) - ddot_(&i__2, &AP(jc), &c__1, &AP(jc), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    ajj = AP(jj) - qdot(&i__2, &AP(jc), &c__1, &AP(jc), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    ajj = AP(jj) - qdot_(&i__2, &AP(jc), &c__1, &AP(jc), &c__1);
#endif

	    if (ajj <= 0.) {
		AP(jj) = ajj;
		goto L30;
	    }
	    AP(jj) = sqrt(ajj);
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	jj = 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness
. */

	    ajj = AP(jj);
	    if (ajj <= 0.) {
		AP(jj) = ajj;
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    AP(jj) = ajj;

/*           Compute elements J+1:N of column J and update the tra
iling   
             submatrix. */

	    if (j < *n) {
		i__2 = *n - j;
		d__1 = 1. / ajj;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &d__1, &AP(jj + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &d__1, &AP(jj + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &d__1, &AP(jj + 1), &c__1);
#endif

		i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dspr_("Lower", &i__2, &c_b16, &AP(jj + 1), &c__1, &AP(jj + *n 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspr("Lower", &i__2, &c_b16, &AP(jj + 1), &c__1, &AP(jj + *n 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspr_("Lower", &i__2, &c_b16, &AP(jj + 1), &c__1, &AP(jj + *n 
#endif

			- j + 1));
		jj = jj + *n - j + 1;
	    }
/* L20: */
	}
    }
    goto L40;

L30:
    *info = j;

L40:
    return;

/*     End of DPPTRF */

} /* dpptrf_ */

