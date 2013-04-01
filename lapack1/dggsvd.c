#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qggsvd(char *jobu, char *jobv, char *jobq, int *m, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qggsvd_(char *jobu, char *jobv, char *jobq, int *m, 
#endif

	int *n, int *p, int *k, int *l, LONG DOUBLE *a, 
	int *lda, LONG DOUBLE *b, int *ldb, LONG DOUBLE *alpha, 
	LONG DOUBLE *beta, LONG DOUBLE *u, int *ldu, LONG DOUBLE *v, int 
	*ldv, LONG DOUBLE *q, int *ldq, LONG DOUBLE *work, int *iwork, 
	int *info)
{
/*  -- LAPACK driver routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGSVD computes the generalized singular value decomposition (GSVD)   
    of an M-by-N real matrix A and P-by-N real matrix B:   

        U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R )   

    where U, V and Q are orthogonal matrices, and Z' is the transpose   
    of Z.  Let K+L = the effective numerical rank of the matrix (A',B')', 
  
    then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and   
    D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the   
    following structures, respectively:   

    If M-K-L >= 0,   

                        K  L   
           D1 =     K ( I  0 )   
                    L ( 0  C )   
                M-K-L ( 0  0 )   

                      K  L   
           D2 =   L ( 0  S )   
                P-L ( 0  0 )   

                    N-K-L  K    L   
      ( 0 R ) = K (  0   R11  R12 )   
                L (  0    0   R22 )   

    where   

      C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),   
      S = diag( BETA(K+1),  ... , BETA(K+L) ),   
      C**2 + S**2 = I.   

      R is stored in A(1:K+L,N-K-L+1:N) on exit.   

    If M-K-L < 0,   

                      K M-K K+L-M   
           D1 =   K ( I  0    0   )   
                M-K ( 0  C    0   )   

                        K M-K K+L-M   
           D2 =   M-K ( 0  S    0  )   
                K+L-M ( 0  0    I  )   
                  P-L ( 0  0    0  )   

                       N-K-L  K   M-K  K+L-M   
      ( 0 R ) =     K ( 0    R11  R12  R13  )   
                  M-K ( 0     0   R22  R23  )   
                K+L-M ( 0     0    0   R33  )   

    where   

      C = diag( ALPHA(K+1), ... , ALPHA(M) ),   
      S = diag( BETA(K+1),  ... , BETA(M) ),   
      C**2 + S**2 = I.   

      (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored   
      ( 0  R22 R23 )   
      in B(M-K+1:L,N+M-K-L+1:N) on exit.   

    The routine computes C, S, R, and optionally the orthogonal   
    transformation matrices U, V and Q.   

    In particular, if B is an N-by-N nonsingular matrix, then the GSVD of 
  
    A and B implicitly gives the SVD of A*inv(B):   
                         A*inv(B) = U*(D1*inv(D2))*V'.   
    If ( A',B')' has orthonormal columns, then the GSVD of A and B is   
    also equal to the CS decomposition of A and B. Furthermore, the GSVD 
  
    can be used to derive the solution of the eigenvalue problem:   
                         A'*A x = lambda* B'*B x.   
    In some literature, the GSVD of A and B is presented in the form   
                     U'*A*X = ( 0 D1 ),   V'*B*X = ( 0 D2 )   
    where U and V are orthogonal and X is nonsingular, D1 and D2 are   
    ``diagonal''.  The former GSVD form can be converted to the latter   
    form by taking the nonsingular matrix X as   

                         X = Q*( I   0    )   
                               ( 0 inv(R) ).   

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

    N       (input) INTEGER   
            The number of columns of the matrices A and B.  N >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B.  P >= 0.   

    K       (output) INTEGER   
    L       (output) INTEGER   
            On exit, K and L specify the dimension of the subblocks   
            described in the Purpose section.   
            K + L = effective numerical rank of (A',B')'.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A contains the triangular matrix R, or part of R.   
            See Purpose for details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,M).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, B contains the triangular matrix R if M-K-L < 0.   
            See Purpose for details.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDA >= MAX(1,P).   

    ALPHA   (output) LONG DOUBLE PRECISION array, dimension (N)   
    BETA    (output) LONG DOUBLE PRECISION array, dimension (N)   
            On exit, ALPHA and BETA contain the generalized singular   
            value pairs of A and B;   
              ALPHA(1:K) = 1,   
              BETA(1:K)  = 0,   
            and if M-K-L >= 0,   
              ALPHA(K+1:K+L) = C,   
              BETA(K+1:K+L)  = S,   
            or if M-K-L < 0,   
              ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0   
              BETA(K+1:M) =S, BETA(M+1:K+L) =1   
            and   
              ALPHA(K+L+1:N) = 0   
              BETA(K+L+1:N)  = 0   

    U       (output) LONG DOUBLE PRECISION array, dimension (LDU,M)   
            If JOBU = 'U', U contains the M-by-M orthogonal matrix U.   
            If JOBU = 'N', U is not referenced.   

    LDU     (input) INTEGER   
            The leading dimension of the array U. LDU >= MAX(1,M) if   
            JOBU = 'U'; LDU >= 1 otherwise.   

    V       (output) LONG DOUBLE PRECISION array, dimension (LDV,P)   
            If JOBV = 'V', V contains the P-by-P orthogonal matrix V.   
            If JOBV = 'N', V is not referenced.   

    LDV     (input) INTEGER   
            The leading dimension of the array V. LDV >= MAX(1,P) if   
            JOBV = 'V'; LDV >= 1 otherwise.   

    Q       (output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q.   
            If JOBQ = 'N', Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q. LDQ >= MAX(1,N) if   
            JOBQ = 'Q'; LDQ >= 1 otherwise.   

    WORK    (workspace) LONG DOUBLE PRECISION array,   
                        dimension (MAX(3*N,M,P)+N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

    INFO    (output)INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = 1, the Jacobi-type procedure failed to   
                  converge.  For further details, see subroutine DTGSJA. 
  

    Internal Parameters   
    ===================   

    TOLA    LONG DOUBLE PRECISION   
    TOLB    LONG DOUBLE PRECISION   
            TOLA and TOLB are the thresholds to determine the effective   
            rank of (A',B')'. Generally, they are set to   
                     TOLA = MAX(M,N)*norm(A)*MAZHEPS,   
                     TOLB = MAX(P,N)*norm(B)*MAZHEPS.   
            The size of TOLA and TOLB may affect the size of backward   
            errors of the decomposition.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static LONG DOUBLE tola, tolb, unfl;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE anorm, bnorm;
    static long int wantq, wantu, wantv;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *), dlange_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *), dlange_(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *), dlange_(char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtgsja_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtgsja(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtgsja_(char *, char *, char *, int *, 
#endif

	    int *, int *, int *, int *, LONG DOUBLE *, int 
	    *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *,
	     int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    int *);
    static int ncycle;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dggsvp_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qggsvp(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qggsvp_(
#endif

	    char *, char *, char *, int *, int *, int *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static LONG DOUBLE ulp;


#define ALPHA(I) alpha[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define U(I,J) u[(I)-1 + ((J)-1)* ( *ldu)]
#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    wantu = lsame_(jobu, "U");
    wantv = lsame_(jobv, "V");
    wantq = lsame_(jobq, "Q");

    *info = 0;
    if (! (wantu || lsame_(jobu, "N"))) {
	*info = -1;
    } else if (! (wantv || lsame_(jobv, "N"))) {
	*info = -2;
    } else if (! (wantq || lsame_(jobq, "N"))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*n < 0) {
	*info = -5;
    } else if (*p < 0) {
	*info = -6;
    } else if (*lda < MAX(1,*m)) {
	*info = -10;
    } else if (*ldb < MAX(1,*p)) {
	*info = -12;
    } else if (*ldu < 1 || (wantu && *ldu < *m)) {
	*info = -16;
    } else if (*ldv < 1 || (wantv && *ldv < *p)) {
	*info = -18;
    } else if (*ldq < 1 || (wantq && *ldq < *n)) {
	*info = -20;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGSVD", &i__1);
	return;
    }

/*     Compute the Frobenius norm of matrices A and B */


#ifdef PETSC_PREFIX_SUFFIX
    anorm = dlange_("1", m, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    anorm = qlange("1", m, n, &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    anorm = qlange_("1", m, n, &A(1,1), lda, &WORK(1));
#endif


#ifdef PETSC_PREFIX_SUFFIX
    bnorm = dlange_("1", p, n, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
    bnorm = qlange("1", p, n, &B(1,1), ldb, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    bnorm = qlange_("1", p, n, &B(1,1), ldb, &WORK(1));
#endif


/*     Get machine precision and set up threshold for determining   
       the effective numerical rank of the matrices A and B. */


#ifdef PETSC_PREFIX_SUFFIX
    ulp = dlamch_("Precision");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    ulp = qlamch("Precision");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    ulp = qlamch_("Precision");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    unfl = dlamch_("Safe Minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    unfl = qlamch("Safe Minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    unfl = qlamch_("Safe Minimum");
#endif

    tola = MAX(*m,*n) * MAX(anorm,unfl) * ulp;
    tolb = MAX(*p,*n) * MAX(bnorm,unfl) * ulp;

/*     Preprocessing */


#ifdef PETSC_PREFIX_SUFFIX
    dggsvp_(jobu, jobv, jobq, m, p, n, &A(1,1), lda, &B(1,1), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qggsvp(jobu, jobv, jobq, m, p, n, &A(1,1), lda, &B(1,1), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qggsvp_(jobu, jobv, jobq, m, p, n, &A(1,1), lda, &B(1,1), ldb, &
#endif

	    tola, &tolb, k, l, &U(1,1), ldu, &V(1,1), ldv, &Q(1,1), ldq, &IWORK(1), &WORK(1), &WORK(*n + 1), info);

/*     Compute the GSVD of two upper "triangular" matrices */


#ifdef PETSC_PREFIX_SUFFIX
    dtgsja_(jobu, jobv, jobq, m, p, n, k, l, &A(1,1), lda, &B(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtgsja(jobu, jobv, jobq, m, p, n, k, l, &A(1,1), lda, &B(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtgsja_(jobu, jobv, jobq, m, p, n, k, l, &A(1,1), lda, &B(1,1), 
#endif

	    ldb, &tola, &tolb, &ALPHA(1), &BETA(1), &U(1,1), ldu, &V(1,1), ldv, &Q(1,1), ldq, &WORK(1), &ncycle, info);

    return;

/*     End of DGGSVD */

} /* dggsvd_ */

