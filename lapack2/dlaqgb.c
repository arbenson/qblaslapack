#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaqgb_(int *m, int *n, int *kl, int *ku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaqgb(int *m, int *n, int *kl, int *ku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaqgb_(int *m, int *n, int *kl, int *ku,
#endif

	 LONG DOUBLE *ab, int *ldab, LONG DOUBLE *r, LONG DOUBLE *c, 
	LONG DOUBLE *rowcnd, LONG DOUBLE *colcnd, LONG DOUBLE *amax, char *equed)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DLAQGB equilibrates a general M by N band matrix A with KL   
    subdiagonals and KU superdiagonals using the row and scaling factors 
  
    in the vectors R and C.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    AB      (input/output) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            On entry, the matrix A in band storage, in rows 1 to KL+KU+1. 
  
            The j-th column of A is stored in the j-th column of the   
            array AB as follows:   
            AB(ku+1+i-j,j) = A(i,j) for MAX(1,j-ku)<=i<=MIN(m,j+kl)   

            On exit, the equilibrated matrix, in the same storage format 
  
            as A.  See EQUED for the form of the equilibrated matrix.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDA >= KL+KU+1.   

    R       (output) LONG DOUBLE PRECISION array, dimension (M)   
            The row scale factors for A.   

    C       (output) LONG DOUBLE PRECISION array, dimension (N)   
            The column scale factors for A.   

    ROWCND  (output) LONG DOUBLE PRECISION   
            Ratio of the smallest R(i) to the largest R(i).   

    COLCND  (output) LONG DOUBLE PRECISION   
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
    int i__1, i__2, i__3, i__4, i__5, i__6;
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

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

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
/* Computing MAX */
		i__2 = 1, i__3 = j - *ku;
/* Computing MIN */
		i__5 = *m, i__6 = j + *kl;
		i__4 = MIN(i__5,i__6);
		for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i) {
		    AB(*ku+1+i-j,j) = cj * AB(*ku+1+i-j,j);
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
/* Computing MAX */
	    i__4 = 1, i__2 = j - *ku;
/* Computing MIN */
	    i__5 = *m, i__6 = j + *kl;
	    i__3 = MIN(i__5,i__6);
	    for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i) {
		AB(*ku+1+i-j,j) = R(i) * AB(*ku+1+i-j,j);
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
/* Computing MAX */
	    i__3 = 1, i__4 = j - *ku;
/* Computing MIN */
	    i__5 = *m, i__6 = j + *kl;
	    i__2 = MIN(i__5,i__6);
	    for (i = MAX(1,j-*ku); i <= MIN(*m,j+*kl); ++i) {
		AB(*ku+1+i-j,j) = cj * R(i) * AB(*ku+1+i-j,j);
/* L50: */
	    }
/* L60: */
	}
	*(unsigned char *)equed = 'B';
    }

    return;

/*     End of DLAQGB */

} /* dlaqgb_ */

