#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dormhr_(char *side, char *trans, int *m, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qormhr(char *side, char *trans, int *m, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qormhr_(char *side, char *trans, int *m, int *n, 
#endif

	int *ilo, int *ihi, LONG DOUBLE *a, int *lda, LONG DOUBLE *
	tau, LONG DOUBLE *c, int *ldc, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORMHR overwrites the general real M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    where Q is a real orthogonal matrix of order nq, with nq = m if   
    SIDE = 'L' and nq = n if SIDE = 'R'. Q is defined as the product of   
    IHI-ILO elementary reflectors, as returned by DGEHRD:   

    Q = H(ilo) H(ilo+1) . . . H(ihi-1).   

    Arguments   
    =========   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q or Q**T from the Left;   
            = 'R': apply Q or Q**T from the Right.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q;   
            = 'T':  Transpose, apply Q**T.   

    M       (input) INTEGER   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix C. N >= 0.   

    ILO     (input) INTEGER   
    IHI     (input) INTEGER   
            ILO and IHI must have the same values as in the previous call 
  
            of DGEHRD. Q is equal to the unit matrix except in the   
            submatrix Q(ilo+1:ihi,ilo+1:ihi).   
            If SIDE = 'L', then 1 <= ILO <= IHI <= M, if M > 0, and   
            ILO = 1 and IHI = 0, if M = 0;   
            if SIDE = 'R', then 1 <= ILO <= IHI <= N, if N > 0, and   
            ILO = 1 and IHI = 0, if N = 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension   
                                 (LDA,M) if SIDE = 'L'   
                                 (LDA,N) if SIDE = 'R'   
            The vectors which define the elementary reflectors, as   
            returned by DGEHRD.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.   
            LDA >= MAX(1,M) if SIDE = 'L'; LDA >= MAX(1,N) if SIDE = 'R'. 
  

    TAU     (input) LONG DOUBLE PRECISION array, dimension   
                                 (M-1) if SIDE = 'L'   
                                 (N-1) if SIDE = 'R'   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEHRD.   

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
    static int iinfo, i1, i2, mi, nh, ni, nq, nw;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dormqr_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qormqr(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qormqr_(
#endif

	    char *, char *, int *, int *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, int *);


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = lsame_(side, "L");

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
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*ilo < 1 || *ilo > MAX(1,nq)) {
	*info = -5;
    } else if (*ihi < MIN(*ilo,nq) || *ihi > nq) {
	*info = -6;
    } else if (*lda < MAX(1,nq)) {
	*info = -8;
    } else if (*ldc < MAX(1,*m)) {
	*info = -11;
    } else if (*lwork < MAX(1,nw)) {
	*info = -13;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORMHR", &i__1);
	return;
    }

/*     Quick return if possible */

    nh = *ihi - *ilo;
    if (*m == 0 || *n == 0 || nh == 0) {
	WORK(1) = 1.;
	return;
    }

    if (left) {
	mi = nh;
	ni = *n;
	i1 = *ilo + 1;
	i2 = 1;
    } else {
	mi = *m;
	ni = nh;
	i1 = 1;
	i2 = *ilo + 1;
    }


#ifdef PETSC_PREFIX_SUFFIX
    dormqr_(side, trans, &mi, &ni, &nh, &A(*ilo+1,*ilo), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormqr(side, trans, &mi, &ni, &nh, &A(*ilo+1,*ilo), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormqr_(side, trans, &mi, &ni, &nh, &A(*ilo+1,*ilo), lda, &
#endif

	    TAU(*ilo), &C(i1,i2), ldc, &WORK(1), lwork, &iinfo);
    return;

/*     End of DORMHR */

} /* dormhr_ */

