#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaqge_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaqge(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaqge_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *r, LONG DOUBLE *c, LONG DOUBLE *rowcnd, LONG DOUBLE *
	colcnd, LONG DOUBLE *amax, char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAQGE equilibrates a general M by N matrix A using the row and   
    scaling factors in the vectors R and C.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M by N matrix A.   
            On exit, the equilibrated matrix.  See EQUED for the form of 
  
            the equilibrated matrix.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(M,1).   

    R       (input) LONG DOUBLE PRECISION array, dimension (M)   
            The row scale factors for A.   

    C       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The column scale factors for A.   

    ROWCND  (input) LONG DOUBLE PRECISION   
            Ratio of the smallest R(i) to the largest R(i).   

    COLCND  (input) LONG DOUBLE PRECISION   
            Ratio of the smallest C(i) to the largest C(i).   

    AMAX    (input) LONG DOUBLE PRECISION   
            Absolute value of largest matrix entry.   

    EQUED   (output) CHARACTER*1   
            Specifies the form of equilibration that was done.   
            = 'N':  No equilibration   
            = 'R':  Row equilibration, i.e., A has been premultiplied by 
  
                    diag(R).   
            = 'C':  Column equilibration, i.e., A has been postmultiplied 
  
                    by diag(C).   
            = 'B':  Both row and column equilibration, i.e., A has been   
                    replaced by diag(R) * A * diag(C).   

    Internal Parameters   
    ===================   

    THRESH is a threshold value used to decide if row or column scaling   
    should be done based on the ratio of the row or column scaling   
    factors.  If ROWCND < THRESH, row scaling is done, and if   
    COLCND < THRESH, column scaling is done.   

    LARGE and SMALL are threshold values used to decide if row scaling   
    should be done based on the absolute size of the largest matrix   
    element.  If AMAX > LARGE or AMAX < SMALL, row scaling is done.   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static int i, j;
    static LONG DOUBLE large, small, cj;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif



#define R(I) r[(I)-1]
#define C(I) c[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (*m <= 0 || *n <= 0) {
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

    if (*rowcnd >= .1 && *amax >= small && *amax <= large) {

/*        No row scaling */

	if (*colcnd >= .1) {

/*           No column scaling */

	    *(unsigned char *)equed = 'N';
	} else {

/*           Column scaling */

	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		cj = C(j);
		i__2 = *m;
		for (i = 1; i <= *m; ++i) {
		    A(i,j) = cj * A(i,j);
/* L10: */
		}
/* L20: */
	    }
	    *(unsigned char *)equed = 'C';
	}
    } else if (*colcnd >= .1) {

/*        Row scaling, no column scaling */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		A(i,j) = R(i) * A(i,j);
/* L30: */
	    }
/* L40: */
	}
	*(unsigned char *)equed = 'R';
    } else {

/*        Row and column scaling */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    cj = C(j);
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		A(i,j) = cj * R(i) * A(i,j);
/* L50: */
	    }
/* L60: */
	}
	*(unsigned char *)equed = 'B';
    }

    return;

/*     End of DLAQGE */

} /* dlaqge_ */

