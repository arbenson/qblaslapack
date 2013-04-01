#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbequ_(int *m, int *n, int *kl, int *ku,
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbequ(int *m, int *n, int *kl, int *ku,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbequ_(int *m, int *n, int *kl, int *ku,
#endif

	 LONG DOUBLE *ab, int *ldab, LONG DOUBLE *r, LONG DOUBLE *c, 
	LONG DOUBLE *rowcnd, LONG DOUBLE *colcnd, LONG DOUBLE *amax, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGBEQU computes row and column scalings intended to equilibrate an   
    M-by-N band matrix A and reduce its condition number.  R returns the 
  
    row scale factors and C the column scale factors, chosen to try to   
    make the largest element in each row and column of the matrix B with 
  
    elements B(i,j)=R(i)*A(i,j)*C(j) have absolute value 1.   

    R(i) and C(j) are restricted to be between SMLNUM = smallest safe   
    number and BIGNUM = largest safe number.  Use of these scaling   
    factors is not guaranteed to reduce the condition number of A but   
    works well in practice.   

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
  

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The band matrix A, stored in rows 1 to KL+KU+1.  The j-th   
            column of A is stored in the j-th column of the array AB as   
            follows:   
            AB(ku+1+i-j,j) = A(i,j) for MAX(1,j-ku)<=i<=MIN(m,j+kl).   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KL+KU+1.   

    R       (output) LONG DOUBLE PRECISION array, dimension (M)   
            If INFO = 0, or INFO > M, R contains the row scale factors   
            for A.   

    C       (output) LONG DOUBLE PRECISION array, dimension (N)   
            If INFO = 0, C contains the column scale factors for A.   

    ROWCND  (output) LONG DOUBLE PRECISION   
            If INFO = 0 or INFO > M, ROWCND contains the ratio of the   
            smallest R(i) to the largest R(i).  If ROWCND >= 0.1 and   
            AMAX is neither too large nor too small, it is not worth   
            scaling by R.   

    COLCND  (output) LONG DOUBLE PRECISION   
            If INFO = 0, COLCND contains the ratio of the smallest   
            C(i) to the largest C(i).  If COLCND >= 0.1, it is not   
            worth scaling by C.   

    AMAX    (output) LONG DOUBLE PRECISION   
            Absolute value of largest matrix element.  If AMAX is very   
            close to overflow or very close to underflow, the matrix   
            should be scaled.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, and i is   
                  <= M:  the i-th row of A is exactly zero   
                  >  M:  the (i-M)-th column of A is exactly zero   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int i__1, i__2, i__3, i__4;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static int i, j;
    static LONG DOUBLE rcmin, rcmax;
    static int kd;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE bignum, smlnum;


#define R(I) r[(I)-1]
#define C(I) c[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBEQU", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	*rowcnd = 1.;
	*colcnd = 1.;
	*amax = 0.;
	return;
    }

/*     Get machine constants. */


#ifdef PETSC_PREFIX_SUFFIX
    smlnum = dlamch_("S");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    smlnum = qlamch("S");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    smlnum = qlamch_("S");
#endif

    bignum = 1. / smlnum;

/*     Compute row scale factors. */

    i__1 = *m;
    for (i = 1; i <= *m; ++i) {
	R(i) = 0.;
/* L10: */
    }

/*     Find the maximum element in each row. */

    kd = *ku + 1;
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	i__2 = j - *ku;
/* Computing MIN */
	i__4 = j + *kl;
	i__3 = MIN(i__4,*m);
	for (i = MAX(j-*ku,1); i <= MIN(j+*kl,*m); ++i) {
/* Computing MAX */
	    d__2 = R(i), d__3 = (d__1 = AB(kd+i-j,j), ABS(
		    d__1));
	    R(i) = MAX(d__2,d__3);
/* L20: */
	}
/* L30: */
    }

/*     Find the maximum and minimum scale factors. */

    rcmin = bignum;
    rcmax = 0.;
    i__1 = *m;
    for (i = 1; i <= *m; ++i) {
/* Computing MAX */
	d__1 = rcmax, d__2 = R(i);
	rcmax = MAX(d__1,d__2);
/* Computing MIN */
	d__1 = rcmin, d__2 = R(i);
	rcmin = MIN(d__1,d__2);
/* L40: */
    }
    *amax = rcmax;

    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. 
*/

	i__1 = *m;
	for (i = 1; i <= *m; ++i) {
	    if (R(i) == 0.) {
		*info = i;
		return;
	    }
/* L50: */
	}
    } else {

/*        Invert the scale factors. */

	i__1 = *m;
	for (i = 1; i <= *m; ++i) {
/* Computing MIN   
   Computing MAX */
	    d__2 = R(i);
	    d__1 = MAX(d__2,smlnum);
	    R(i) = 1. / MIN(d__1,bignum);
/* L60: */
	}

/*        Compute ROWCND = MIN(R(I)) / MAX(R(I)) */

	*rowcnd = MAX(rcmin,smlnum) / MIN(rcmax,bignum);
    }

/*     Compute column scale factors */

    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
	C(j) = 0.;
/* L70: */
    }

/*     Find the maximum element in each column,   
       assuming the row scaling computed above. */

    kd = *ku + 1;
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
/* Computing MAX */
	i__3 = j - *ku;
/* Computing MIN */
	i__4 = j + *kl;
	i__2 = MIN(i__4,*m);
	for (i = MAX(j-*ku,1); i <= MIN(j+*kl,*m); ++i) {
/* Computing MAX */
	    d__2 = C(j), d__3 = (d__1 = AB(kd+i-j,j), ABS(
		    d__1)) * R(i);
	    C(j) = MAX(d__2,d__3);
/* L80: */
	}
/* L90: */
    }

/*     Find the maximum and minimum scale factors. */

    rcmin = bignum;
    rcmax = 0.;
    i__1 = *n;
    for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	d__1 = rcmin, d__2 = C(j);
	rcmin = MIN(d__1,d__2);
/* Computing MAX */
	d__1 = rcmax, d__2 = C(j);
	rcmax = MAX(d__1,d__2);
/* L100: */
    }

    if (rcmin == 0.) {

/*        Find the first zero scale factor and return an error code. 
*/

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (C(j) == 0.) {
		*info = *m + j;
		return;
	    }
/* L110: */
	}
    } else {

/*        Invert the scale factors. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN   
   Computing MAX */
	    d__2 = C(j);
	    d__1 = MAX(d__2,smlnum);
	    C(j) = 1. / MIN(d__1,bignum);
/* L120: */
	}

/*        Compute COLCND = MIN(C(J)) / MAX(C(J)) */

	*colcnd = MAX(rcmin,smlnum) / MIN(rcmax,bignum);
    }

    return;

/*     End of DGBEQU */

} /* dgbequ_ */

