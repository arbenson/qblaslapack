#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpbtf2_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpbtf2(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpbtf2_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif

	ab, int *ldab, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DPBTF2 computes the Cholesky factorization of a real symmetric   
    positive definite band matrix A.   

    The factorization has the form   
       A = U' * U ,  if UPLO = 'U', or   
       A = L  * L',  if UPLO = 'L',   
    where U is an upper triangular matrix, U' is the transpose of U, and 
  
    L is lower triangular.   

    This is the unblocked version of the algorithm, calling Level 2 BLAS. 
  

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of super-diagonals of the matrix A if UPLO = 'U', 
  
            or the number of sub-diagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the upper or lower triangle of the symmetric band   
            matrix A, stored in the first KD+1 rows of the array.  The   
            j-th column of A is stored in the j-th column of the array AB 
  
            as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  

            On exit, if INFO = 0, the triangular factor U or L from the   
            Cholesky factorization A = U'*U or A = L*L' of the band   
            matrix A, in the same storage format as A.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, the leading minor of order k is not   
                 positive definite, and the factorization could not be   
                 completed.   

    Further Details   
    ===============   

    The band storage scheme is illustrated by the following example, when 
  
    N = 6, KD = 2, and UPLO = 'U':   

    On entry:                       On exit:   

        *    *   a13  a24  a35  a46      *    *   u13  u24  u35  u46   
        *   a12  a23  a34  a45  a56      *   u12  u23  u34  u45  u56   
       a11  a22  a33  a44  a55  a66     u11  u22  u33  u44  u55  u66   

    Similarly, if UPLO = 'L' the format of A is as follows:   

    On entry:                       On exit:   

       a11  a22  a33  a44  a55  a66     l11  l22  l33  l44  l55  l66   
       a21  a32  a43  a54  a65   *      l21  l32  l43  l54  l65   *   
       a31  a42  a53  a64   *    *      l31  l42  l53  l64   *    *   

    Array elements marked * are not used by the routine.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b8 = -1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Builtin functions */
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsyr_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *);
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
    static int kn;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE ajj;
    static int kld;




#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kd < 0) {
	*info = -3;
    } else if (*ldab < *kd + 1) {
	*info = -5;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPBTF2", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/* Computing MAX */
    i__1 = 1, i__2 = *ldab - 1;
    kld = MAX(i__1,i__2);

    if (upper) {

/*        Compute the Cholesky factorization A = U'*U. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute U(J,J) and test for non-positive-definiteness
. */

	    ajj = AB(*kd+1,j);
	    if (ajj <= 0.) {
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    AB(*kd+1,j) = ajj;

/*           Compute elements J+1:J+KN of row J and update the   
             trailing submatrix within the band.   

   Computing MIN */
	    i__2 = *kd, i__3 = *n - j;
	    kn = MIN(i__2,i__3);
	    if (kn > 0) {
		d__1 = 1. / ajj;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&kn, &d__1, &AB(*kd,j+1), &kld);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&kn, &d__1, &AB(*kd,j+1), &kld);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&kn, &d__1, &AB(*kd,j+1), &kld);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dsyr_("Upper", &kn, &c_b8, &AB(*kd,j+1), &kld,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsyr("Upper", &kn, &c_b8, &AB(*kd,j+1), &kld,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsyr_("Upper", &kn, &c_b8, &AB(*kd,j+1), &kld,
#endif

			 &AB(*kd+1,j+1), &kld);
	    }
/* L10: */
	}
    } else {

/*        Compute the Cholesky factorization A = L*L'. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {

/*           Compute L(J,J) and test for non-positive-definiteness
. */

	    ajj = AB(1,j);
	    if (ajj <= 0.) {
		goto L30;
	    }
	    ajj = sqrt(ajj);
	    AB(1,j) = ajj;

/*           Compute elements J+1:J+KN of column J and update the 
  
             trailing submatrix within the band.   

   Computing MIN */
	    i__2 = *kd, i__3 = *n - j;
	    kn = MIN(i__2,i__3);
	    if (kn > 0) {
		d__1 = 1. / ajj;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&kn, &d__1, &AB(2,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&kn, &d__1, &AB(2,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&kn, &d__1, &AB(2,j), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		dsyr_("Lower", &kn, &c_b8, &AB(2,j), &c__1, &AB(1,j+1), &kld);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsyr("Lower", &kn, &c_b8, &AB(2,j), &c__1, &AB(1,j+1), &kld);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsyr_("Lower", &kn, &c_b8, &AB(2,j), &c__1, &AB(1,j+1), &kld);
#endif

	    }
/* L20: */
	}
    }
    return;

L30:
    *info = j;
    return;

/*     End of DPBTF2 */

} /* dpbtf2_ */

