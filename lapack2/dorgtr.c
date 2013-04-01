#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dorgtr_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qorgtr(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qorgtr_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGTR generates a real orthogonal matrix Q which is defined as the   
    product of n-1 elementary reflectors of order N, as returned by   
    DSYTRD:   

    if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),   

    if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangle of A contains elementary reflectors   
                   from DSYTRD;   
            = 'L': Lower triangle of A contains elementary reflectors   
                   from DSYTRD.   

    N       (input) INTEGER   
            The order of the matrix Q. N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the vectors which define the elementary reflectors, 
  
            as returned by DSYTRD.   
            On exit, the N-by-N orthogonal matrix Q.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,N).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DSYTRD.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= MAX(1,N-1).   
            For optimum performance LWORK >= (N-1)*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int i, j;
    extern long int lsame_(char *, char *);
    static int iinfo;
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dorgql_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qorgql(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qorgql_(
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), dorgqr_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qorgqr(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *), qorgqr_(
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, int *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = *n - 1;
	if (*lwork < MAX(i__1,i__2)) {
	    *info = -7;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGTR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return;
    }

    if (upper) {

/*        Q was determined by a call to DSYTRD with UPLO = 'U'   

          Shift the vectors which define the elementary reflectors one
   
          column to the left, and set the last row and column of Q to 
  
          those of the unit matrix */

	i__1 = *n - 1;
	for (j = 1; j <= *n-1; ++j) {
	    i__2 = j - 1;
	    for (i = 1; i <= j-1; ++i) {
		A(i,j) = A(i,j+1);
/* L10: */
	    }
	    A(*n,j) = 0.;
/* L20: */
	}
	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {
	    A(i,*n) = 0.;
/* L30: */
	}
	A(*n,*n) = 1.;

/*        Generate Q(1:n-1,1:n-1) */

	i__1 = *n - 1;
	i__2 = *n - 1;
	i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dorgql_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorgql(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorgql_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), 
#endif

		lwork, &iinfo);

    } else {

/*        Q was determined by a call to DSYTRD with UPLO = 'L'.   

          Shift the vectors which define the elementary reflectors one
   
          column to the right, and set the first row and column of Q t
o   
          those of the unit matrix */

	for (j = *n; j >= 2; --j) {
	    A(1,j) = 0.;
	    i__1 = *n;
	    for (i = j + 1; i <= *n; ++i) {
		A(i,j) = A(i,j-1);
/* L40: */
	    }
/* L50: */
	}
	A(1,1) = 1.;
	i__1 = *n;
	for (i = 2; i <= *n; ++i) {
	    A(i,1) = 0.;
/* L60: */
	}
	if (*n > 1) {

/*           Generate Q(2:n,2:n) */

	    i__1 = *n - 1;
	    i__2 = *n - 1;
	    i__3 = *n - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dorgqr_(&i__1, &i__2, &i__3, &A(2,2), lda, &TAU(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorgqr(&i__1, &i__2, &i__3, &A(2,2), lda, &TAU(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorgqr_(&i__1, &i__2, &i__3, &A(2,2), lda, &TAU(1), 
#endif

		    &WORK(1), lwork, &iinfo);
	}
    }
    return;

/*     End of DORGTR */

} /* dorgtr_ */

