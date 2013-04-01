#include <math.h>
#define MIN(a,b)           ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)           ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)             ( ((a)<0.0)   ? -(a) : (a) )



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dorghr_(int *n, int *ilo, int *ihi, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qorghr(int *n, int *ilo, int *ihi, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qorghr_(int *n, int *ilo, int *ihi, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *work, 
	int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGHR generates a real orthogonal matrix Q which is defined as the   
    product of IHI-ILO elementary reflectors of order N, as returned by   
    DGEHRD:   

    Q = H(ilo) H(ilo+1) . . . H(ihi-1).   

    Arguments   
    =========   

    N       (input) INT   
            The order of the matrix Q. N >= 0.   

    ILO     (input) INT   
    IHI     (input) INT   
            ILO and IHI must have the same values as in the previous call 
  
            of DGEHRD. Q is equal to the unit matrix except in the   
            submatrix Q(ilo+1:ihi,ilo+1:ihi).   
            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors, 
  
            as returned by DGEHRD.   
            On exit, the N-by-N orthogonal matrix Q.   

    LDA     (input) INT   
            The leading dimension of the array A. LDA >= MAX(1,N).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEHRD.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INT   
            The dimension of the array WORK. LWORK >= IHI-ILO.   
            For optimum performance LWORK >= (IHI-ILO)*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INT   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static int i, j, iinfo, nh;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dorgqr_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qorgqr(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qorgqr_(
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*ilo < 1 || *ilo > MAX(1,*n)) {
	*info = -2;
    } else if (*ihi < MIN(*ilo,*n) || *ihi > *n) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *ihi - *ilo;
	if (*lwork < MAX(i__1,i__2)) {
	    *info = -8;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGHR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return;
    }

/*     Shift the vectors which define the elementary reflectors one   
       column to the right, and set the first ilo and the last n-ihi   
       rows and columns to those of the unit matrix */

    i__1 = *ilo + 1;
    for (j = *ihi; j >= *ilo+1; --j) {
	i__2 = j - 1;
	for (i = 1; i <= j-1; ++i) {
	    A(i,j) = 0.;
/* L10: */
	}
	i__2 = *ihi;
	for (i = j + 1; i <= *ihi; ++i) {
	    A(i,j) = A(i,j-1);
/* L20: */
	}
	i__2 = *n;
	for (i = *ihi + 1; i <= *n; ++i) {
	    A(i,j) = 0.;
/* L30: */
	}
/* L40: */
    }
    i__1 = *ilo;
    for (j = 1; j <= *ilo; ++j) {
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    A(i,j) = 0.;
/* L50: */
	}
	A(j,j) = 1.;
/* L60: */
    }
    i__1 = *n;
    for (j = *ihi + 1; j <= *n; ++j) {
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    A(i,j) = 0.;
/* L70: */
	}
	A(j,j) = 1.;
/* L80: */
    }

    nh = *ihi - *ilo;
    if (nh > 0) {

/*        Generate Q(ilo+1:ihi,ilo+1:ihi) */


#ifdef PETSC_PREFIX_SUFFIX
	dorgqr_(&nh, &nh, &nh, &A(*ilo+1,*ilo+1), lda, &TAU(*
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorgqr(&nh, &nh, &nh, &A(*ilo+1,*ilo+1), lda, &TAU(*
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorgqr_(&nh, &nh, &nh, &A(*ilo+1,*ilo+1), lda, &TAU(*
#endif

		ilo), &WORK(1), lwork, &iinfo);
    }
    return;

/*     End of DORGHR */

} /* dorghr_ */

