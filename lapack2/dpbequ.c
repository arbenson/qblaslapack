#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpbequ_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpbequ(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpbequ_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif

	ab, int *ldab, LONG DOUBLE *s, LONG DOUBLE *scond, LONG DOUBLE *amax,
	 int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPBEQU computes row and column scalings intended to equilibrate a   
    symmetric positive definite band matrix A and reduce its condition   
    number (with respect to the two-norm).  S contains the scale factors, 
  
    S(i) = 1/sqrt(A(i,i)), chosen so that the scaled matrix B with   
    elements B(i,j) = S(i)*A(i,j)*S(j) has ones on the diagonal.  This   
    choice of S puts the condition number of B within a factor N of the   
    smallest possible condition number over all possible diagonal   
    scalings.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangular of A is stored;   
            = 'L':  Lower triangular of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The upper or lower triangle of the symmetric band matrix A,   
            stored in the first KD+1 rows of the array.  The j-th column 
  
            of A is stored in the j-th column of the array AB as follows: 
  
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  

    LDAB     (input) INTEGER   
            The leading dimension of the array A.  LDAB >= KD+1.   

    S       (output) LONG DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, S contains the scale factors for A.   

    SCOND   (output) LONG DOUBLE PRECISION   
            If INFO = 0, S contains the ratio of the smallest S(i) to   
            the largest S(i).  If SCOND >= 0.1 and AMAX is neither too   
            large nor too small, it is not worth scaling by S.   

    AMAX    (output) LONG DOUBLE PRECISION   
            Absolute value of largest matrix element.  If AMAX is very   
            close to overflow or very close to underflow, the matrix   
            should be scaled.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, the i-th diagonal element is nonpositive. 
  

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    LONG DOUBLE d__1, d__2;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE smin;
    static int i, j;
    extern long int lsame_(char *, char *);
    static long int upper;
    extern /* Subroutine */ void xerbla_(char *, int *);


#define S(I) s[(I)-1]

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
	xerbla_("DPBEQU", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	*scond = 1.;
	*amax = 0.;
	return;
    }

    if (upper) {
	j = *kd + 1;
    } else {
	j = 1;
    }

/*     Initialize SMIN and AMAX. */

    S(1) = AB(j,1);
    smin = S(1);
    *amax = S(1);

/*     Find the minimum and maximum diagonal elements. */

    i__1 = *n;
    for (i = 2; i <= *n; ++i) {
	S(i) = AB(j,i);
/* Computing MIN */
	d__1 = smin, d__2 = S(i);
	smin = MIN(d__1,d__2);
/* Computing MAX */
	d__1 = *amax, d__2 = S(i);
	*amax = MAX(d__1,d__2);
/* L10: */
    }

    if (smin <= 0.) {

/*        Find the first non-positive diagonal element and return. */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (S(i) <= 0.) {
		*info = i;
		return;
	    }
/* L20: */
	}
    } else {

/*        Set the scale factors to the reciprocals   
          of the diagonal elements. */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    S(i) = 1. / sqrt(S(i));
/* L30: */
	}

/*        Compute SCOND = MIN(S(I)) / MAX(S(I)) */

	*scond = sqrt(smin) / sqrt(*amax);
    }
    return;

/*     End of DPBEQU */

} /* dpbequ_ */

