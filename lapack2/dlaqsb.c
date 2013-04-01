#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaqsb_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaqsb(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaqsb_(char *uplo, int *n, int *kd, LONG DOUBLE *
#endif

	ab, int *ldab, LONG DOUBLE *s, LONG DOUBLE *scond, LONG DOUBLE *amax,
	 char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAQSB equilibrates a symmetric band matrix A using the scaling   
    factors in the vector S.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored.   
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

    S       (output) LONG DOUBLE PRECISION array, dimension (N)   
            The scale factors for A.   

    SCOND   (input) LONG DOUBLE PRECISION   
            Ratio of the smallest S(i) to the largest S(i).   

    AMAX    (input) LONG DOUBLE PRECISION   
            Absolute value of largest matrix entry.   

    EQUED   (output) CHARACTER*1   
            Specifies whether or not equilibration was done.   
            = 'N':  No equilibration.   
            = 'Y':  Equilibration was done, i.e., A has been replaced by 
  
                    diag(S) * A * diag(S).   

    Internal Parameters   
    ===================   

    THRESH is a threshold value used to decide if scaling should be done 
  
    based on the ratio of the scaling factors.  If SCOND < THRESH,   
    scaling is done.   

    LARGE and SMALL are threshold values used to decide if scaling should 
  
    be done based on the absolute size of the largest matrix element.   
    If AMAX > LARGE or AMAX < SMALL, scaling is done.   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    /* Local variables */
    static int i, j;
    static LONG DOUBLE large;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE small, cj;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif



#define S(I) s[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    if (*n <= 0) {
	*(unsigned char *)equed = 'N';
	return;
    }

/*     Initialize LARGE and SMALL. */


#ifdef PETSC_PREFIX_SUFFIX
    small = dlamch_("Safe minimum") / dlamch_("Precision");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    small = qlamch("Safe minimum") / dlamch_("Precision");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    small = qlamch_("Safe minimum") / dlamch_("Precision");
#endif

    large = 1. / small;

    if (*scond >= .1 && *amax >= small && *amax <= large) {

/*        No equilibration */

	*(unsigned char *)equed = 'N';
    } else {

/*        Replace A by diag(S) * A * diag(S). */

	if (lsame_(uplo, "U")) {

/*           Upper triangle of A is stored in band format. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
/* Computing MAX */
		i__2 = 1, i__3 = j - *kd;
		i__4 = j;
		for (i = MAX(1,j-*kd); i <= j; ++i) {
		    AB(*kd+1+i-j,j) = cj * S(i) * AB(*kd+1+i-j,j);
/* L10: */
		}
/* L20: */
	    }
	} else {

/*           Lower triangle of A is stored. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
/* Computing MIN */
		i__2 = *n, i__3 = j + *kd;
		i__4 = MIN(i__2,i__3);
		for (i = j; i <= MIN(*n,j+*kd); ++i) {
		    AB(i+1-j,j) = cj * S(i) * AB(i+1-j,j);
/* L30: */
		}
/* L40: */
	    }
	}
	*(unsigned char *)equed = 'Y';
    }

    return;

/*     End of DLAQSB */

} /* dlaqsb_ */

