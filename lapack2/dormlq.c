#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dormlq_(char *side, char *trans, int *m, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qormlq(char *side, char *trans, int *m, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qormlq_(char *side, char *trans, int *m, int *n, 
#endif

	int *k, LONG DOUBLE *a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *
	c, int *ldc, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORMLQ overwrites the general real M-by-N matrix C with   

                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    where Q is a real orthogonal matrix defined as the product of k   
    elementary reflectors   

          Q = H(k) . . . H(2) H(1)   

    as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N   
    if SIDE = 'R'.   

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

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines   
            the matrix Q.   
            If SIDE = 'L', M >= K >= 0;   
            if SIDE = 'R', N >= K >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension   
                                 (LDA,M) if SIDE = 'L',   
                                 (LDA,N) if SIDE = 'R'   
            The i-th row must contain the vector which defines the   
            elementary reflector H(i), for i = 1,2,...,k, as returned by 
  
            DGELQF in the first k rows of its array argument A.   
            A is modified by the routine but restored on exit.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,K).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGELQF.   

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
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static int c__2 = 2;
    static int c__65 = 65;
    
    /* System generated locals */
    int  i__1, i__2, i__4, 
	    i__5;
    char ch__1[3];
    /* Builtin functions   */
    /* Local variables */
    static long int left;
    static int i;
    static LONG DOUBLE t[4160]	/* was [65][64] */;
    extern long int lsame_(char *, char *);
    static int nbmin, iinfo, i1, i2, i3;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorml2_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorml2(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorml2_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *);
    static int ib, ic, jc, nb, mi, ni;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarfb_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfb(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarfb_(char *, char *, char *, char *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static int nq, nw;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarft_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarft(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarft_(char *, char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);
    static long int notran;
    static int ldwork;
    static char transt[1];
    static int iws;



#define T(I) t[(I)]
#define WAS(I) was[(I)]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");

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
    } else if (! notran && ! lsame_(trans, "T")) {
	*info = -2;
    } else if (*m < 0) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*k < 0 || *k > nq) {
	*info = -5;
    } else if (*lda < MAX(1,*k)) {
	*info = -7;
    } else if (*ldc < MAX(1,*m)) {
	*info = -10;
    } else if (*lwork < MAX(1,nw)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORMLQ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0 || *k == 0) {
	WORK(1) = 1.;
	return;
    }

/*     Determine the block size.  NB may be at most NBMAX, where NBMAX   
       is used to define the local array T.   

   Computing MIN   
   Writing concatenation */
    /*i__3[0] = 1, a__1[0] = side;
    i__3[1] = 1, a__1[1] = trans;
    s_cat(ch__1, a__1, i__3, &c__2, 2L); */
    ch__1[0] = *side; ch__1[1] = *trans; ch__1[2] = 0;
    i__1 = 64, i__2 = ilaenv_(&c__1, "DORMLQ", ch__1, m, n, k, &c_n1, 6L, 2L);
    nb = MIN(i__1,i__2);
    nbmin = 2;
    ldwork = nw;
    if (nb > 1 && nb < *k) {
	iws = nw * nb;
	if (*lwork < iws) {
	    nb = *lwork / ldwork;
/* Computing MAX   
   Writing concatenation */
	    /* i__3[0] = 1, a__1[0] = side;
	    i__3[1] = 1, a__1[1] = trans;
	    s_cat(ch__1, a__1, i__3, &c__2, 2L); */
            ch__1[0] = *side; ch__1[1] = *trans; ch__1[2] = 0;
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DORMLQ", ch__1, m, n, k, &c_n1, 
		    6L, 2L);
	    nbmin = MAX(i__1,i__2);
	}
    } else {
	iws = nw;
    }

    if (nb < nbmin || nb >= *k) {

/*        Use unblocked code */


#ifdef PETSC_PREFIX_SUFFIX
	dorml2_(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorml2(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorml2_(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1)
#endif

		, ldc, &WORK(1), &iinfo);
    } else {

/*        Use blocked code */

	if ((left && notran) || (! left && ! notran)) {
	    i1 = 1;
	    i2 = *k;
	    i3 = nb;
	} else {
	    i1 = (*k - 1) / nb * nb + 1;
	    i2 = 1;
	    i3 = -nb;
	}

	if (left) {
	    ni = *n;
	    jc = 1;
	} else {
	    mi = *m;
	    ic = 1;
	}

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}

	i__1 = i2;
	i__2 = i3;
	for (i = i1; i3 < 0 ? i >= i2 : i <= i2; i += i3) {
/* Computing MIN */
	    i__4 = nb, i__5 = *k - i + 1;
	    ib = MIN(i__4,i__5);

/*           Form the triangular factor of the block reflector   
             H = H(i) H(i+1) . . . H(i+ib-1) */

	    i__4 = nq - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlarft_("Forward", "Rowwise", &i__4, &ib, &A(i,i), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarft("Forward", "Rowwise", &i__4, &ib, &A(i,i), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarft_("Forward", "Rowwise", &i__4, &ib, &A(i,i), lda,
#endif

		     &TAU(i), t, &c__65);
	    if (left) {

/*              H or H' is applied to C(i:m,1:n) */

		mi = *m - i + 1;
		ic = i;
	    } else {

/*              H or H' is applied to C(1:m,i:n) */

		ni = *n - i + 1;
		jc = i;
	    }

/*           Apply H or H' */


#ifdef PETSC_PREFIX_SUFFIX
	    dlarfb_(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &A(i,i), lda, t, &c__65, &C(ic,jc), ldc, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfb(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &A(i,i), lda, t, &c__65, &C(ic,jc), ldc, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfb_(side, transt, "Forward", "Rowwise", &mi, &ni, &ib, &A(i,i), lda, t, &c__65, &C(ic,jc), ldc, &
#endif

		    WORK(1), &ldwork);
/* L10: */
	}
    }
    WORK(1) = (LONG DOUBLE) iws;
    return;

/*     End of DORMLQ */

} /* dormlq_ */

