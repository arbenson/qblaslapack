#include <math.h>
#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )



#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dormbr_(char *vect, char *side, char *trans, int *m, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qormbr(char *vect, char *side, char *trans, int *m, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qormbr_(char *vect, char *side, char *trans, int *m, 
#endif

	int *n, int *k, LONG DOUBLE *a, int *lda, LONG DOUBLE *tau, 
	LONG DOUBLE *c, int *ldc, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C   
    with   
                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      Q * C          C * Q   
    TRANS = 'T':      Q**T * C       C * Q**T   

    If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C   
    with   
                    SIDE = 'L'     SIDE = 'R'   
    TRANS = 'N':      P * C          C * P   
    TRANS = 'T':      P**T * C       C * P**T   

    Here Q and P**T are the orthogonal matrices determined by DGEBRD when 
  
    reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and 
  
    P**T are defined as products of elementary reflectors H(i) and G(i)   
    respectively.   

    Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the   
    order of the orthogonal matrix Q or P**T that is applied.   

    If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:   
    if nq >= k, Q = H(1) H(2) . . . H(k);   
    if nq < k, Q = H(1) H(2) . . . H(nq-1).   

    If VECT = 'P', A is assumed to have been a K-by-NQ matrix:   
    if k < nq, P = G(1) G(2) . . . G(k);   
    if k >= nq, P = G(1) G(2) . . . G(nq-1).   

    Arguments   
    =========   

    VECT    (input) CHARACTER*1   
            = 'Q': apply Q or Q**T;   
            = 'P': apply P or P**T.   

    SIDE    (input) CHARACTER*1   
            = 'L': apply Q, Q**T, P or P**T from the Left;   
            = 'R': apply Q, Q**T, P or P**T from the Right.   

    TRANS   (input) CHARACTER*1   
            = 'N':  No transpose, apply Q  or P;   
            = 'T':  Transpose, apply Q**T or P**T.   

    M       (input) INT   
            The number of rows of the matrix C. M >= 0.   

    N       (input) INT   
            The number of columns of the matrix C. N >= 0.   

    K       (input) INT   
            If VECT = 'Q', the number of columns in the original   
            matrix reduced by DGEBRD.   
            If VECT = 'P', the number of rows in the original   
            matrix reduced by DGEBRD.   
            K >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension   
                                  (LDA,MIN(nq,K)) if VECT = 'Q'   
                                  (LDA,nq)        if VECT = 'P'   
            The vectors which define the elementary reflectors H(i) and   
            G(i), whose products determine the matrices Q and P, as   
            returned by DGEBRD.   

    LDA     (input) INT   
            The leading dimension of the array A.   
            If VECT = 'Q', LDA >= MAX(1,nq);   
            if VECT = 'P', LDA >= MAX(1,MIN(nq,K)).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (MIN(nq,K))   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i) or G(i) which determines Q or P, as returned   
            by DGEBRD in the array argument TAUQ or TAUP.   

    C       (input/output) LONG DOUBLE PRECISION array, dimension (LDC,N)   
            On entry, the M-by-N matrix C.   
            On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q   
            or P*C or P**T*C or C*P or C*P**T.   

    LDC     (input) INT   
            The leading dimension of the array C. LDC >= MAX(1,M).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INT   
            The dimension of the array WORK.   
            If SIDE = 'L', LWORK >= MAX(1,N);   
            if SIDE = 'R', LWORK >= MAX(1,M).   
            For optimum performance LWORK >= N*NB if SIDE = 'L', and   
            LWORK >= M*NB if SIDE = 'R', where NB is the optimal   
            blocksize.   

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
    static long int left;
    extern long int lsame_(char *, char *);
    static int iinfo, i1, i2, mi, ni, nq, nw;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dormlq_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qormlq(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qormlq_(
#endif

	    char *, char *, int *, int *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, int *);
    static long int notran;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dormqr_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qormqr(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qormqr_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, int *);
    static long int applyq;
    static char transt[1];


#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define C(I,J) c[(I)-1 + ((J)-1)* ( *ldc)]

    *info = 0;
    applyq = lsame_(vect, "Q");
    left = lsame_(side, "L");
    notran = lsame_(trans, "N");

/*     NQ is the order of Q or P and NW is the minimum dimension of WORK 
*/

    if (left) {
	nq = *m;
	nw = *n;
    } else {
	nq = *n;
	nw = *m;
    }
    if (! applyq && ! lsame_(vect, "P")) {
	*info = -1;
    } else if (! left && ! lsame_(side, "R")) {
	*info = -2;
    } else if (! notran && ! lsame_(trans, "T")) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*k < 0) {
	*info = -6;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = 1, i__2 = MIN(nq,*k);
	if ((applyq && *lda < MAX(1,nq)) || (! applyq && *lda < MAX(i__1,i__2))) {
	    *info = -8;
	} else if (*ldc < MAX(1,*m)) {
	    *info = -11;
	} else if (*lwork < MAX(1,nw)) {
	    *info = -13;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORMBR", &i__1);
	return;
    }

/*     Quick return if possible */

    WORK(1) = 1.;
    if (*m == 0 || *n == 0) {
	return;
    }

    if (applyq) {

/*        Apply Q */

	if (nq >= *k) {

/*           Q was determined by a call to DGEBRD with nq >= k */


#ifdef PETSC_PREFIX_SUFFIX
	    dormqr_(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1), ldc, &WORK(1), lwork, &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormqr(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1), ldc, &WORK(1), lwork, &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormqr_(side, trans, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1), ldc, &WORK(1), lwork, &iinfo);
#endif

	} else if (nq > 1) {

/*           Q was determined by a call to DGEBRD with nq < k */

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dormqr_(side, trans, &mi, &ni, &i__1, &A(2,1), lda, &TAU(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormqr(side, trans, &mi, &ni, &i__1, &A(2,1), lda, &TAU(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormqr_(side, trans, &mi, &ni, &i__1, &A(2,1), lda, &TAU(1)
#endif

		    , &C(i1,i2), ldc, &WORK(1), lwork, &iinfo);
	}
    } else {

/*        Apply P */

	if (notran) {
	    *(unsigned char *)transt = 'T';
	} else {
	    *(unsigned char *)transt = 'N';
	}
	if (nq > *k) {

/*           P was determined by a call to DGEBRD with nq > k */


#ifdef PETSC_PREFIX_SUFFIX
	    dormlq_(side, transt, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1), ldc, &WORK(1), lwork, &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormlq(side, transt, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1), ldc, &WORK(1), lwork, &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormlq_(side, transt, m, n, k, &A(1,1), lda, &TAU(1), &C(1,1), ldc, &WORK(1), lwork, &iinfo);
#endif

	} else if (nq > 1) {

/*           P was determined by a call to DGEBRD with nq <= k */

	    if (left) {
		mi = *m - 1;
		ni = *n;
		i1 = 2;
		i2 = 1;
	    } else {
		mi = *m;
		ni = *n - 1;
		i1 = 1;
		i2 = 2;
	    }
	    i__1 = nq - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dormlq_(side, transt, &mi, &ni, &i__1, &A(1,2), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormlq(side, transt, &mi, &ni, &i__1, &A(1,2), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormlq_(side, transt, &mi, &ni, &i__1, &A(1,2), lda,
#endif

		     &TAU(1), &C(i1,i2), ldc, &WORK(1), lwork, &
		    iinfo);
	}
    }
    return;

/*     End of DORMBR */

} /* dormbr_ */

