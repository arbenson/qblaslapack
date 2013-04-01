#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qggsvp(char *jobu, char *jobv, char *jobq, int *m, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qggsvp_(char *jobu, char *jobv, char *jobq, int *m, 
#endif

	int *p, int *n, LONG DOUBLE *a, int *lda, LONG DOUBLE *b, 
	int *ldb, LONG DOUBLE *tola, LONG DOUBLE *tolb, int *k, int 
	*l, LONG DOUBLE *u, int *ldu, LONG DOUBLE *v, int *ldv, 
	LONG DOUBLE *q, int *ldq, int *iwork, LONG DOUBLE *tau, 
	LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGSVP computes orthogonal matrices U, V and Q such that   

                     N-K-L  K    L   
     U'*A*Q =     K ( 0    A12  A13 )  if M-K-L >= 0;   
                  L ( 0     0   A23 )   
              M-K-L ( 0     0    0  )   

                     N-K-L  K    L   
            =     K ( 0    A12  A13 )  if M-K-L < 0;   
                M-K ( 0     0   A23 )   

                   N-K-L  K    L   
     V'*B*Q =   L ( 0     0   B13 )   
              P-L ( 0     0    0  )   

    where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular   
    upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,   
    otherwise A23 is (M-K)-by-L upper trapezoidal.  K+L = the effective   
    numerical rank of the (M+P)-by-N matrix (A',B')'.  Z' denotes the   
    transpose of Z.   

    This decomposition is the preprocessing step for computing the   
    Generalized Singular Value Decomposition (GSVD), see subroutine   
    DGGSVD.   

    Arguments   
    =========   

    JOBU    (input) CHARACTER*1   
            = 'U':  Orthogonal matrix U is computed;   
            = 'N':  U is not computed.   

    JOBV    (input) CHARACTER*1   
            = 'V':  Orthogonal matrix V is computed;   
            = 'N':  V is not computed.   

    JOBQ    (input) CHARACTER*1   
            = 'Q':  Orthogonal matrix Q is computed;   
            = 'N':  Q is not computed.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B.  P >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrices A and B.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A contains the triangular (or trapezoidal) matrix   
            described in the Purpose section.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,M).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, B contains the triangular matrix described in   
            the Purpose section.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,P).   

    TOLA    (input) LONG DOUBLE PRECISION   
    TOLB    (input) LONG DOUBLE PRECISION   
            TOLA and TOLB are the thresholds to determine the effective   
            numerical rank of matrix B and a subblock of A. Generally,   
            they are set to   
               TOLA = MAX(M,N)*norm(A)*MAZHEPS,   
               TOLB = MAX(P,N)*norm(B)*MAZHEPS.   
            The size of TOLA and TOLB may affect the size of backward   
            errors of the decomposition.   

    K       (output) INTEGER   
    L       (output) INTEGER   
            On exit, K and L specify the dimension of the subblocks   
            described in Purpose.   
            K + L = effective numerical rank of (A',B')'.   

    U       (output) LONG DOUBLE PRECISION array, dimension (LDU,M)   
            If JOBU = 'U', U contains the orthogonal matrix U.   
            If JOBU = 'N', U is not referenced.   

    LDU     (input) INTEGER   
            The leading dimension of the array U. LDU >= MAX(1,M) if   
            JOBU = 'U'; LDU >= 1 otherwise.   

    V       (output) LONG DOUBLE PRECISION array, dimension (LDV,M)   
            If JOBV = 'V', V contains the orthogonal matrix V.   
            If JOBV = 'N', V is not referenced.   

    LDV     (input) INTEGER   
            The leading dimension of the array V. LDV >= MAX(1,P) if   
            JOBV = 'V'; LDV >= 1 otherwise.   

    Q       (output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            If JOBQ = 'Q', Q contains the orthogonal matrix Q.   
            If JOBQ = 'N', Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q. LDQ >= MAX(1,N) if   
            JOBQ = 'Q'; LDQ >= 1 otherwise.   

    IWORK   (workspace) INTEGER array, dimension (N)   

    TAU     (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (MAX(3*N,M,P)) 
  

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   


    Further Details   
    ===============   

    The subroutine uses LAPACK subroutine DGEQPF for the QR factorization 
  
    with column pivoting to detect the effective numerical rank of the   
    a matrix. It may be replaced by a better rank determination strategy. 
  

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b12 = 0.;
    static LONG DOUBLE c_b22 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Local variables */
    static int i, j;
    extern long int lsame_(char *, char *);
    static long int wantq, wantu, wantv;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgeqr2_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeqr2(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeqr2_(int *, int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *), dgerq2_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *), qgerq2(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *), qgerq2_(
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dorg2r_(int *, int *, int *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qorg2r(int *, int *, int *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qorg2r_(int *, int *, int *,
#endif

	     LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dorm2r_(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorm2r(char *, char *, int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorm2r_(char *, char *, int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dormr2_(char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qormr2(char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qormr2_(char *, char *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), dgeqpf_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qgeqpf(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qgeqpf_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_(char *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy(char *, int *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_(char *, int *, int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dlaset_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlaset(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlaset_(char *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), xerbla_(char *, int *), dlapmt_(long int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), xerbla_(char *, int *), qlapmt(long int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), xerbla_(char *, int *), qlapmt_(long int *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, int *);
    static long int forwrd;



#define IWORK(I) iwork[(I)-1]
#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define U(I,J) u[(I)-1 + ((J)-1)* ( *ldu)]
#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    wantu = lsame_(jobu, "U");
    wantv = lsame_(jobv, "V");
    wantq = lsame_(jobq, "Q");
    forwrd = 1;

    *info = 0;
    if (! (wantu || lsame_(jobu, "N"))) {
	*info = -1;
    } else if (! (wantv || lsame_(jobv, "N"))) {
	*info = -2;
    } else if (! (wantq || lsame_(jobq, "N"))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*lda < MAX(1,*m)) {
	*info = -8;
    } else if (*ldb < MAX(1,*p)) {
	*info = -10;
    } else if (*ldu < 1 || (wantu && *ldu < *m)) {
	*info = -16;
    } else if (*ldv < 1 || (wantv && *ldv < *p)) {
	*info = -18;
    } else if (*ldq < 1 || (wantq && *ldq < *n)) {
	*info = -20;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGSVP", &i__1);
	return;
    }

/*     QR with column pivoting of B: B*P = V*( S11 S12 )   
                                             (  0   0  ) */

    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	IWORK(i) = 0;
/* L10: */
    }

#ifdef PETSC_PREFIX_SUFFIX
    dgeqpf_(p, n, &B(1,1), ldb, &IWORK(1), &TAU(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgeqpf(p, n, &B(1,1), ldb, &IWORK(1), &TAU(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgeqpf_(p, n, &B(1,1), ldb, &IWORK(1), &TAU(1), &WORK(1), info);
#endif


/*     Update A := A*P */


#ifdef PETSC_PREFIX_SUFFIX
    dlapmt_(&forwrd, m, n, &A(1,1), lda, &IWORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlapmt(&forwrd, m, n, &A(1,1), lda, &IWORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlapmt_(&forwrd, m, n, &A(1,1), lda, &IWORK(1));
#endif


/*     Determine the effective rank of matrix B. */

    *l = 0;
    i__1 = MIN(*p,*n);
    for (i = 1; i <= MIN(*p,*n); ++i) {
	if ((d__1 = B(i,i), ABS(d__1)) > *tolb) {
	    ++(*l);
	}
/* L20: */
    }

    if (wantv) {

/*        Copy the details of V, and form V. */


#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", p, p, &c_b12, &c_b12, &V(1,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", p, p, &c_b12, &c_b12, &V(1,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", p, p, &c_b12, &c_b12, &V(1,1), ldv);
#endif

	if (*p > 1) {
	    i__1 = *p - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("Lower", &i__1, n, &B(2,1), ldb, &V(2,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("Lower", &i__1, n, &B(2,1), ldb, &V(2,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("Lower", &i__1, n, &B(2,1), ldb, &V(2,1), 
#endif

		    ldv);
	}
	i__1 = MIN(*p,*n);

#ifdef PETSC_PREFIX_SUFFIX
	dorg2r_(p, p, &i__1, &V(1,1), ldv, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorg2r(p, p, &i__1, &V(1,1), ldv, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorg2r_(p, p, &i__1, &V(1,1), ldv, &TAU(1), &WORK(1), info);
#endif

    }

/*     Clean up B */

    i__1 = *l - 1;
    for (j = 1; j <= *l-1; ++j) {
	i__2 = *l;
	for (i = j + 1; i <= *l; ++i) {
	    B(i,j) = 0.;
/* L30: */
	}
/* L40: */
    }
    if (*p > *l) {
	i__1 = *p - *l;

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", &i__1, n, &c_b12, &c_b12, &B(*l+1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", &i__1, n, &c_b12, &c_b12, &B(*l+1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", &i__1, n, &c_b12, &c_b12, &B(*l+1,1), ldb);
#endif

    }

    if (wantq) {

/*        Set Q = I and Update Q := Q*P */


#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b12, &c_b22, &Q(1,1), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b12, &c_b22, &Q(1,1), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b12, &c_b22, &Q(1,1), ldq);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlapmt_(&forwrd, n, n, &Q(1,1), ldq, &IWORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlapmt(&forwrd, n, n, &Q(1,1), ldq, &IWORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlapmt_(&forwrd, n, n, &Q(1,1), ldq, &IWORK(1));
#endif

    }

    if (*p >= *l && *n != *l) {

/*        RQ factorization of (S11 S12): ( S11 S12 ) = ( 0 S12 )*Z */


#ifdef PETSC_PREFIX_SUFFIX
	dgerq2_(l, n, &B(1,1), ldb, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgerq2(l, n, &B(1,1), ldb, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgerq2_(l, n, &B(1,1), ldb, &TAU(1), &WORK(1), info);
#endif


/*        Update A := A*Z' */


#ifdef PETSC_PREFIX_SUFFIX
	dormr2_("Right", "Transpose", m, n, l, &B(1,1), ldb, &TAU(1), &A(1,1), lda, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qormr2("Right", "Transpose", m, n, l, &B(1,1), ldb, &TAU(1), &A(1,1), lda, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qormr2_("Right", "Transpose", m, n, l, &B(1,1), ldb, &TAU(1), &A(1,1), lda, &WORK(1), info);
#endif


	if (wantq) {

/*           Update Q := Q*Z' */


#ifdef PETSC_PREFIX_SUFFIX
	    dormr2_("Right", "Transpose", n, n, l, &B(1,1), ldb, &TAU(1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormr2("Right", "Transpose", n, n, l, &B(1,1), ldb, &TAU(1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormr2_("Right", "Transpose", n, n, l, &B(1,1), ldb, &TAU(1),
#endif

		     &Q(1,1), ldq, &WORK(1), info);
	}

/*        Clean up B */

	i__1 = *n - *l;

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", l, &i__1, &c_b12, &c_b12, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", l, &i__1, &c_b12, &c_b12, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", l, &i__1, &c_b12, &c_b12, &B(1,1), ldb);
#endif

	i__1 = *n;
	for (j = *n - *l + 1; j <= *n; ++j) {
	    i__2 = *l;
	    for (i = j - *n + *l + 1; i <= *l; ++i) {
		B(i,j) = 0.;
/* L50: */
	    }
/* L60: */
	}

    }

/*     Let              N-L     L   
                  A = ( A11    A12 ) M,   

       then the following does the complete QR decomposition of A11:   

                A11 = U*(  0  T12 )*P1'   
                        (  0   0  ) */

    i__1 = *n - *l;
    for (i = 1; i <= *n-*l; ++i) {
	IWORK(i) = 0;
/* L70: */
    }
    i__1 = *n - *l;

#ifdef PETSC_PREFIX_SUFFIX
    dgeqpf_(m, &i__1, &A(1,1), lda, &IWORK(1), &TAU(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgeqpf(m, &i__1, &A(1,1), lda, &IWORK(1), &TAU(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgeqpf_(m, &i__1, &A(1,1), lda, &IWORK(1), &TAU(1), &WORK(1), info);
#endif


/*     Determine the effective rank of A11 */

    *k = 0;
/* Computing MIN */
    i__2 = *m, i__3 = *n - *l;
    i__1 = MIN(i__2,i__3);
    for (i = 1; i <= MIN(*m,*n-*l); ++i) {
	if ((d__1 = A(i,i), ABS(d__1)) > *tola) {
	    ++(*k);
	}
/* L80: */
    }

/*     Update A12 := U'*A12, where A12 = A( 1:M, N-L+1:N )   

   Computing MIN */
    i__2 = *m, i__3 = *n - *l;
    i__1 = MIN(i__2,i__3);

#ifdef PETSC_PREFIX_SUFFIX
    dorm2r_("Left", "Transpose", m, l, &i__1, &A(1,1), lda, &TAU(1), &A(1,*n-*l+1), lda, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qorm2r("Left", "Transpose", m, l, &i__1, &A(1,1), lda, &TAU(1), &A(1,*n-*l+1), lda, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qorm2r_("Left", "Transpose", m, l, &i__1, &A(1,1), lda, &TAU(1), &A(1,*n-*l+1), lda, &WORK(1), info);
#endif


    if (wantu) {

/*        Copy the details of U, and form U */


#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", m, m, &c_b12, &c_b12, &U(1,1), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", m, m, &c_b12, &c_b12, &U(1,1), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", m, m, &c_b12, &c_b12, &U(1,1), ldu);
#endif

	if (*m > 1) {
	    i__1 = *m - 1;
	    i__2 = *n - *l;

#ifdef PETSC_PREFIX_SUFFIX
	    dlacpy_("Lower", &i__1, &i__2, &A(2,1), lda, &U(2,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlacpy("Lower", &i__1, &i__2, &A(2,1), lda, &U(2,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlacpy_("Lower", &i__1, &i__2, &A(2,1), lda, &U(2,1)
#endif

		    , ldu);
	}
/* Computing MIN */
	i__2 = *m, i__3 = *n - *l;
	i__1 = MIN(i__2,i__3);

#ifdef PETSC_PREFIX_SUFFIX
	dorg2r_(m, m, &i__1, &U(1,1), ldu, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorg2r(m, m, &i__1, &U(1,1), ldu, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorg2r_(m, m, &i__1, &U(1,1), ldu, &TAU(1), &WORK(1), info);
#endif

    }

    if (wantq) {

/*        Update Q( 1:N, 1:N-L )  = Q( 1:N, 1:N-L )*P1 */

	i__1 = *n - *l;

#ifdef PETSC_PREFIX_SUFFIX
	dlapmt_(&forwrd, n, &i__1, &Q(1,1), ldq, &IWORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlapmt(&forwrd, n, &i__1, &Q(1,1), ldq, &IWORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlapmt_(&forwrd, n, &i__1, &Q(1,1), ldq, &IWORK(1));
#endif

    }

/*     Clean up A: set the strictly lower triangular part of   
       A(1:K, 1:K) = 0, and A( K+1:M, 1:N-L ) = 0. */

    i__1 = *k - 1;
    for (j = 1; j <= *k-1; ++j) {
	i__2 = *k;
	for (i = j + 1; i <= *k; ++i) {
	    A(i,j) = 0.;
/* L90: */
	}
/* L100: */
    }
    if (*m > *k) {
	i__1 = *m - *k;
	i__2 = *n - *l;

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", &i__1, &i__2, &c_b12, &c_b12, &A(*k+1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", &i__1, &i__2, &c_b12, &c_b12, &A(*k+1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", &i__1, &i__2, &c_b12, &c_b12, &A(*k+1,1), 
#endif

		lda);
    }

    if (*n - *l > *k) {

/*        RQ factorization of ( T11 T12 ) = ( 0 T12 )*Z1 */

	i__1 = *n - *l;

#ifdef PETSC_PREFIX_SUFFIX
	dgerq2_(k, &i__1, &A(1,1), lda, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgerq2(k, &i__1, &A(1,1), lda, &TAU(1), &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgerq2_(k, &i__1, &A(1,1), lda, &TAU(1), &WORK(1), info);
#endif


	if (wantq) {

/*           Update Q( 1:N,1:N-L ) = Q( 1:N,1:N-L )*Z1' */

	    i__1 = *n - *l;

#ifdef PETSC_PREFIX_SUFFIX
	    dormr2_("Right", "Transpose", n, &i__1, k, &A(1,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qormr2("Right", "Transpose", n, &i__1, k, &A(1,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qormr2_("Right", "Transpose", n, &i__1, k, &A(1,1), lda, &
#endif

		    TAU(1), &Q(1,1), ldq, &WORK(1), info);
	}

/*        Clean up A */

	i__1 = *n - *l - *k;

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", k, &i__1, &c_b12, &c_b12, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", k, &i__1, &c_b12, &c_b12, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", k, &i__1, &c_b12, &c_b12, &A(1,1), lda);
#endif

	i__1 = *n - *l;
	for (j = *n - *l - *k + 1; j <= *n-*l; ++j) {
	    i__2 = *k;
	    for (i = j - *n + *l + *k + 1; i <= *k; ++i) {
		A(i,j) = 0.;
/* L110: */
	    }
/* L120: */
	}

    }

    if (*m > *k) {

/*        QR factorization of A( K+1:M,N-L+1:N ) */

	i__1 = *m - *k;

#ifdef PETSC_PREFIX_SUFFIX
	dgeqr2_(&i__1, l, &A(*k+1,*n-*l+1), lda, &TAU(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgeqr2(&i__1, l, &A(*k+1,*n-*l+1), lda, &TAU(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgeqr2_(&i__1, l, &A(*k+1,*n-*l+1), lda, &TAU(1), &
#endif

		WORK(1), info);

	if (wantu) {

/*           Update U(:,K+1:M) := U(:,K+1:M)*U1 */

	    i__1 = *m - *k;
/* Computing MIN */
	    i__3 = *m - *k;
	    i__2 = MIN(i__3,*l);

#ifdef PETSC_PREFIX_SUFFIX
	    dorm2r_("Right", "No transpose", m, &i__1, &i__2, &A(*k+1,*n-*l+1), lda, &TAU(1), &U(1,*k+1), ldu, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorm2r("Right", "No transpose", m, &i__1, &i__2, &A(*k+1,*n-*l+1), lda, &TAU(1), &U(1,*k+1), ldu, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorm2r_("Right", "No transpose", m, &i__1, &i__2, &A(*k+1,*n-*l+1), lda, &TAU(1), &U(1,*k+1), ldu, &WORK(1), info);
#endif

	}

/*        Clean up */

	i__1 = *n;
	for (j = *n - *l + 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = j - *n + *k + *l + 1; i <= *m; ++i) {
		A(i,j) = 0.;
/* L130: */
	    }
/* L140: */
	}

    }

    return;

/*     End of DGGSVP */

} /* dggsvp_ */

