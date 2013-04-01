#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dormtr_(char *side, char *uplo, char *trans, int *m, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qormtr(char *side, char *uplo, char *trans, int *m, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qormtr_(char *side, char *uplo, char *trans, int *m, 
#endif

	int *n, LONG DOUBLE *a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *
	c, int *ldc, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORMTR overwrites the general real M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    where Q is a real orthogonal matrix of order nq, with nq = m if   
    SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of   
    nq-1 elementary reflectors, as returned by DSYTRD:   

    if UPLO = 'U', Q = H(nq-1) . . . H(2) H(1);   

    if UPLO = 'L', Q = H(1) H(2) . . . H(nq-1).   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q**T from the Left;   
            = 'R': apply Q or Q**T from the Right.   

    UPLO    (input) CHARACTER*1   
            = 'U': Upper triangle of A contains elementary reflectors   
                   from DSYTRD;   
            = 'L': Lower triangle of A contains elementary reflectors   
                   from DSYTRD.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q**T.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension   
                                 (LDA,M) if SIDE = 'L'   
                                 (LDA,N) if SIDE = 'R'   
            The vectors which define the elementary reflectors, as   
            returned by DSYTRD.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            LDA >= MAX(1,M) if SIDE = 'L'; LDA >= MAX(1,N) if SIDE = 'R'. 
  

    TAU     (input) LONG DOUBLE PRECISION array, dimension   
                                 (M-1) if SIDE = 'L'   
                                 (N-1) if SIDE = 'R'   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DSYTRD.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q. 
  

    LDC     (input) INTEGER   
            The leading dimension of the array C. LDC >= MAX(1,M).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.   
            If SIDE = 'L', LWORK >= MAX(1,N);   
            if SIDE = 'R', LWORK >= MAX(1,M).   
            For optimum performance LWORK >= N*NB if SIDE = 'L', and   
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal   
            blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static long int left;
    extern long int lsame_(char *, char *);
    static int iinfo, i1;
    static long int upper;
    static int i2, mi, ni, nq, nw;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dormql_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qormql(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qormql_(
#endif

	    char *, char *, int *, int *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *, int *), dormqr_(char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, int *), qormqr(char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, int *), qormqr_(char *, char *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    int *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = lsame_(side, "L");
    upper = lsame_(uplo, "U");

/*     NQ is the order of Q and NW is the minimum dimension of WORK */

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }
    if (! left && ! lsame_(side, "R")) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T")) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*lda < MAX(1,nq)) {
	*info = -7;
    } else if (*ldc < MAX(1,*m)) {
	*info = -10;
    } else if (*lwork < MAX(1,nw)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORMTR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || nq == 1) {
	WORK(1) = 1.;
	return;
    }

    if (left) {
	mi = *m - 1;
	ni = *n;
    } else {
	mi = *m;
	ni = *n - 1;
    }

    if (upper) {

/*        Q was determined by a call to DSYTRD with UPLO = 'U' */

	i__1 = nq - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dormql_(side, trans, &mi, &ni, &i__1, &A(1,2), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qormql(side, trans, &mi, &ni, &i__1, &A(1,2), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qormql_(side, trans, &mi, &ni, &i__1, &A(1,2), lda, &
#endif

		TAU(1), &C(1,1), ldc, &WORK(1), lwork, &iinfo);
    } else {

/*        Q was determined by a call to DSYTRD with UPLO = 'L' */

	if (left) {
	    i1 = 2;
	    i2 = 1;
	} else {
	    i1 = 1;
	    i2 = 2;
	}
	i__1 = nq - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dormqr_(side, trans, &mi, &ni, &i__1, &A(2,1), lda, &TAU(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qormqr(side, trans, &mi, &ni, &i__1, &A(2,1), lda, &TAU(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qormqr_(side, trans, &mi, &ni, &i__1, &A(2,1), lda, &TAU(1), &
#endif

		C(i1,i2), ldc, &WORK(1), lwork, &iinfo);
    }
    return;

/*     End of DORMTR */

} /* dormtr_ */

