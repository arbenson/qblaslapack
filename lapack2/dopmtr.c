#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dopmtr_(char *side, char *uplo, char *trans, int *m, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qopmtr(char *side, char *uplo, char *trans, int *m, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qopmtr_(char *side, char *uplo, char *trans, int *m, 
#endif

	int *n, LONG DOUBLE *ap, LONG DOUBLE *tau, LONG DOUBLE *c, int *
	ldc, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DOPMTR overwrites the general real M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    where Q is a real orthogonal matrix of order nq, with nq = m if   
    SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of   
    nq-1 elementary reflectors, as returned by DSPTRD using packed   
    storage:   

    if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);   

    if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q**T from the Left;   
            = 'R': apply Q or Q**T from the Right.   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangular packed storage used in previous   
                   call to DSPTRD;   
            = 'L': Lower triangular packed storage used in previous   
                   call to DSPTRD.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q**T.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    AP      (input) LONG DOUBLE PRECISION array, dimension   
                                 (M*(M+1)/2) if SIDE = 'L'   
                                 (N*(N+1)/2) if SIDE = 'R'   
            The vectors which define the elementary reflectors, as   
            returned by DSPTRD.  AP is modified by the routine but   
            restored on exit.   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (M-1) if SIDE = 'L' 
  
                                       or (N-1) if SIDE = 'R'   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DSPTRD.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. 
  

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= MAX(1,M).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension   
                                     (N) if SIDE = 'L'   
                                     (M) if SIDE = 'R'   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static long int left;
    static int i;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarf_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarf(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarf_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *);
    extern long int lsame_(char *, char *);
    static int i1;
    static long int upper;
    static int i2, i3, ic, jc, ii, mi, ni, nq;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int notran, forwrd;
    static LONG DOUBLE aii;



#define AP(I) ap[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");
    upper = lsame_(uplo, "U");

/*     NQ is the order of Q */

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (! left && ! lsame_(side, "R")) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (! notran && ! lsame_(trans, "T")) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*ldc < MAX(1,*m)) {
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DOPMTR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return;
    }

    if (upper) {

/*        Q was determined by a call to DSPTRD with UPLO = 'U' */

	forwrd = (left && notran) || (! left && ! notran);

	if (forwrd) {
	    i1 = 1;
	    i2 = nq - 1;
	    i3 = 1;
	    ii = 2;
	} else {
	    i1 = nq - 1;
	    i2 = 1;
	    i3 = -1;
	    ii = nq * (nq + 1) / 2 - 1;
	}

	if (left) {
	    ni = *n;
	} else {
	    mi = *m;
	}

	i__1 = i2;
	i__2 = i3;
	for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {
	    if (left) {

/*              H(i) is applied to C(1:i,1:n) */

		mi = i;
	    } else {

/*              H(i) is applied to C(1:m,1:i) */

		ni = i;
	    }

/*           Apply H(i) */

	    aii = AP(ii);
	    AP(ii) = 1.;

#ifdef PETSC_PREFIX_SUFFIX
	    dlarf_(side, &mi, &ni, &AP(ii - i + 1), &c__1, &TAU(i), &C(1,1), ldc, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarf(side, &mi, &ni, &AP(ii - i + 1), &c__1, &TAU(i), &C(1,1), ldc, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarf_(side, &mi, &ni, &AP(ii - i + 1), &c__1, &TAU(i), &C(1,1), ldc, &WORK(1));
#endif

	    AP(ii) = aii;

	    if (forwrd) {
		ii = ii + i + 2;
	    } else {
		ii = ii - i - 1;
	    }
/* L10: */
	}
    } else {

/*        Q was determined by a call to DSPTRD with UPLO = 'L'. */

	forwrd = (left && ! notran) || (! left && notran);

	if (forwrd) {
	    i1 = 1;
	    i2 = nq - 1;
	    i3 = 1;
	    ii = 2;
	} else {
	    i1 = nq - 1;
	    i2 = 1;
	    i3 = -1;
	    ii = nq * (nq + 1) / 2 - 1;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	i__2 = i2;
	i__1 = i3;
	for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {
	    aii = AP(ii);
	    AP(ii) = 1.;
	    if (left) {

/*              H(i) is applied to C(i+1:m,1:n) */

		mi = *m - i;
		ic = i + 1;
	    } else {

/*              H(i) is applied to C(1:m,i+1:n) */

		ni = *n - i;
		jc = i + 1;
	    }

/*           Apply H(i) */


#ifdef PETSC_PREFIX_SUFFIX
	    dlarf_(side, &mi, &ni, &AP(ii), &c__1, &TAU(i), &C(ic,jc), ldc, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarf(side, &mi, &ni, &AP(ii), &c__1, &TAU(i), &C(ic,jc), ldc, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarf_(side, &mi, &ni, &AP(ii), &c__1, &TAU(i), &C(ic,jc), ldc, &WORK(1));
#endif

	    AP(ii) = aii;

	    if (forwrd) {
		ii = ii + nq - i + 1;
	    } else {
		ii = ii - nq + i - 2;
	    }
/* L20: */
	}
    }
    return;

/*     End of DOPMTR */

} /* dopmtr_ */

