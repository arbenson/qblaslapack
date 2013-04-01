#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dorm2l_(char *side, char *trans, int *m, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qorm2l(char *side, char *trans, int *m, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qorm2l_(char *side, char *trans, int *m, int *n, 
#endif

	int *k, LONG DOUBLE *a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *
	c, int *ldc, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORM2L overwrites the general real m by n matrix C with   

          Q * C  if SIDE = 'L' and TRANS = 'N', or   

          Q'* C  if SIDE = 'L' and TRANS = 'T', or   

          C * Q  if SIDE = 'R' and TRANS = 'N', or   

          C * Q' if SIDE = 'R' and TRANS = 'T',   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(k) . . . H(2) H(1)   

    as returned by DGEQLF. Q is of order m if SIDE = 'L' and of order n   
    if SIDE = 'R'.   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q' from the Left   
            = 'R': apply Q or Q' from the Right   

    TRANS   (input) CHARACTER*1   
            = 'N': apply Q  (No transpose)   
            = 'T': apply Q' (Transpose)   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA,K)   
            The i-th column must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by 
  
            DGEQLF in the last k columns of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            If SIDE = 'L', LDA >= MAX(1,M);   
            if SIDE = 'R', LDA >= MAX(1,N).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQLF.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the m by n matrix C.   
            On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.   

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= MAX(1,M).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension   
                                     (N) if SIDE = 'L',   
                                     (M) if SIDE = 'R'   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

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
    static int i1, i2, i3, mi, ni, nq;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int notran;
    static LONG DOUBLE aii;



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");

/*     NQ is the order of Q */

    if (left) {
	nq = *m;
    } else {
	nq = *n;
    }
    if (! left && ! lsame_(side, "R")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < MAX(1,nq)) {
	*info = -7;
    } else if (*ldc < MAX(1,*m)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORM2L", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	return;
    }

    if ((left && notran) || (! left && ! notran)) {
	i1 = 1;
	i2 = *k;
	i3 = 1;
    } else {
	i1 = *k;
	i2 = 1;
	i3 = -1;
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

/*           H(i) is applied to C(1:m-k+i,1:n) */

	    mi = *m - *k + i;
	} else {

/*           H(i) is applied to C(1:m,1:n-k+i) */

	    ni = *n - *k + i;
	}

/*        Apply H(i) */

	aii = A(nq-*k+i,i);
	A(nq-*k+i,i) = 1.;

#ifdef PETSC_PREFIX_SUFFIX
	dlarf_(side, &mi, &ni, &A(1,i), &c__1, &TAU(i), &C(1,1), ldc, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarf(side, &mi, &ni, &A(1,i), &c__1, &TAU(i), &C(1,1), ldc, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarf_(side, &mi, &ni, &A(1,i), &c__1, &TAU(i), &C(1,1), ldc, &WORK(1));
#endif

	A(nq-*k+i,i) = aii;
/* L10: */
    }
    return;

/*     End of DORM2L */

} /* dorm2l_ */

