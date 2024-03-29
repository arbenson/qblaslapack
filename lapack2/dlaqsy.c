#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaqsy_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaqsy(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaqsy_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *s, LONG DOUBLE *scond, LONG DOUBLE *amax, char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAQSY equilibrates a symmetric matrix A using the scaling factors   
    in the vector S.   

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

            On exit, if EQUED = 'Y', the equilibrated matrix:   
            diag(S) * A * diag(S).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(N,1).   

    S       (input) LONG DOUBLE PRECISION array, dimension (N)   
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
    int  i__1, i__2;
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

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

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

/*           Upper triangle of A is stored. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
		i__2 = j;
		for (i = 1; i <= j; ++i) {
		    A(i,j) = cj * S(i) * A(i,j);
/* L10: */
		}
/* L20: */
	    }
	} else {

/*           Lower triangle of A is stored. */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = S(j);
		i__2 = *n;
		for (i = j; i <= *n; ++i) {
		    A(i,j) = cj * S(i) * A(i,j);
/* L30: */
		}
/* L40: */
	    }
	}
	*(unsigned char *)equed = 'Y';
    }

    return;

/*     End of DLAQSY */

} /* dlaqsy_ */

