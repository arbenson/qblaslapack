#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtgsja_(char *jobu, char *jobv, char *jobq, int *m, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtgsja(char *jobu, char *jobv, char *jobq, int *m, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtgsja_(char *jobu, char *jobv, char *jobq, int *m, 
#endif

	int *p, int *n, int *k, int *l, LONG DOUBLE *a, 
	int *lda, LONG DOUBLE *b, int *ldb, LONG DOUBLE *tola, 
	LONG DOUBLE *tolb, LONG DOUBLE *alpha, LONG DOUBLE *beta, LONG DOUBLE *u, 
	int *ldu, LONG DOUBLE *v, int *ldv, LONG DOUBLE *q, int *
	ldq, LONG DOUBLE *work, int *ncycle, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DTGSJA computes the generalized singular value decomposition (GSVD)   
    of two real upper triangular (or trapezoidal) matrices A and B.   

    On entry, it is assumed that matrices A and B have the following   
    forms, which may be obtained by the preprocessing subroutine DGGSVP   
    from a general M-by-N matrix A and P-by-N matrix B:   

                 N-K-L  K    L   
       A =    K ( 0    A12  A13 ) if M-K-L >= 0;   
              L ( 0     0   A23 )   
          M-K-L ( 0     0    0  )   

               N-K-L  K    L   
       A =  K ( 0    A12  A13 ) if M-K-L < 0;   
          M-K ( 0     0   A23 )   

               N-K-L  K    L   
       B =  L ( 0     0   B13 )   
          P-L ( 0     0    0  )   

    where the K-by-K matrix A12 and L-by-L matrix B13 are nonsingular   
    upper triangular; A23 is L-by-L upper triangular if M-K-L >= 0,   
    otherwise A23 is (M-K)-by-L upper trapezoidal.   

    On exit,   

                U'*A*Q = D1*( 0 R ),    V'*B*Q = D2*( 0 R ),   

    where U, V and Q are orthogonal matrices, Z' denotes the transpose   
    of Z, R is a nonsingular upper triangular matrix, and D1 and D2 are   
    ``diagonal'' matrices, which are of the following structures:   

    If M-K-L >= 0,   

                        K  L   
           D1 =     K ( I  0 )   
                    L ( 0  C )   
                M-K-L ( 0  0 )   

                      K  L   
           D2 = L   ( 0  S )   
                P-L ( 0  0 )   

                   N-K-L  K    L   
      ( 0 R ) = K (  0   R11  R12 ) K   
                L (  0    0   R22 ) L   

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
        D2 =   M-K ( 0  S    0   )   
             K+L-M ( 0  0    I   )   
               P-L ( 0  0    0   )   

                   N-K-L  K   M-K  K+L-M   
   ( 0 R ) =    K ( 0    R11  R12  R13  )   
              M-K ( 0     0   R22  R23  )   
            K+L-M ( 0     0    0   R33  )   

    where   
    C = diag( ALPHA(K+1), ... , ALPHA(M) ),   
    S = diag( BETA(K+1),  ... , BETA(M) ),   
    C**2 + S**2 = I.   

    R = ( R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N) and R33 is stored 
  
        (  0  R22 R23 )   
    in B(M-K+1:L,N+M-K-L+1:N) on exit.   

    The computation of the orthogonal transformation matrices U, V or Q   
    is optional.  These matrices may either be formed explicitly, or they 
  
    may be postmultiplied into input matrices U1, V1, or Q1.   

    Arguments   
    =========   

    JOBU    (input) CHARACTER*1   
            = 'U':  U must contain an orthogonal matrix U1 on entry, and 
  
                    the product U1*U is returned;   
            = 'I':  U is initialized to the unit matrix, and the   
                    orthogonal matrix U is returned;   
            = 'N':  U is not computed.   

    JOBV    (input) CHARACTER*1   
            = 'V':  V must contain an orthogonal matrix V1 on entry, and 
  
                    the product V1*V is returned;   
            = 'I':  V is initialized to the unit matrix, and the   
                    orthogonal matrix V is returned;   
            = 'N':  V is not computed.   

    JOBQ    (input) CHARACTER*1   
            = 'Q':  Q must contain an orthogonal matrix Q1 on entry, and 
  
                    the product Q1*Q is returned;   
            = 'I':  Q is initialized to the unit matrix, and the   
                    orthogonal matrix Q is returned;   
            = 'N':  Q is not computed.   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B.  P >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrices A and B.  N >= 0.   

    K       (input) INTEGER   
    L       (input) INTEGER   
            K and L specify the subblocks in the input matrices A and B: 
  
            A23 = A(K+1:MIN(K+L,M),N-L+1:N) and B13 = B(1:L,N-L+1:N)   
            of A and B, whose GSVD is going to be computed by DTGSJA.   
            See Further details.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, A(N-K+1:N,1:MIN(K+L,M) ) contains the triangular   
            matrix R or part of R.  See Purpose for details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,M).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, if necessary, B(M-K+1:L,N+M-K-L+1:N) contains   
            a part of R.  See Purpose for details.   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,P).   

    TOLA    (input) LONG DOUBLE PRECISION   
    TOLB    (input) LONG DOUBLE PRECISION   
            TOLA and TOLB are the convergence criteria for the Jacobi-   
            Kogbetliantz iteration procedure. Generally, they are the   
            same as used in the preprocessing step, say   
                TOLA = MAX(M,N)*norm(A)*MAZHEPS,   
                TOLB = MAX(P,N)*norm(B)*MAZHEPS.   

    ALPHA   (output) LONG DOUBLE PRECISION array, dimension (N)   
    BETA    (output) LONG DOUBLE PRECISION array, dimension (N)   
            On exit, ALPHA and BETA contain the generalized singular   
            value pairs of A and B;   
              ALPHA(1:K) = 1,   
              BETA(1:K)  = 0,   
            and if M-K-L >= 0,   
              ALPHA(K+1:K+L) = diag(C),   
              BETA(K+1:K+L)  = diag(S),   
            or if M-K-L < 0,   
              ALPHA(K+1:M)= C, ALPHA(M+1:K+L)= 0   
              BETA(K+1:M) = S, BETA(M+1:K+L) = 1.   
            Furthermore, if K+L < N,   
              ALPHA(K+L+1:N) = 0 and   
              BETA(K+L+1:N)  = 0.   

    U       (input/output) LONG DOUBLE PRECISION array, dimension (LDU,M)   
            On entry, if JOBU = 'U', U must contain a matrix U1 (usually 
  
            the orthogonal matrix returned by DGGSVP).   
            On exit,   
            if JOBU = 'I', U contains the orthogonal matrix U;   
            if JOBU = 'U', U contains the product U1*U.   
            If JOBU = 'N', U is not referenced.   

    LDU     (input) INTEGER   
            The leading dimension of the array U. LDU >= MAX(1,M) if   
            JOBU = 'U'; LDU >= 1 otherwise.   

    V       (input/output) LONG DOUBLE PRECISION array, dimension (LDV,P)   
            On entry, if JOBV = 'V', V must contain a matrix V1 (usually 
  
            the orthogonal matrix returned by DGGSVP).   
            On exit,   
            if JOBV = 'I', V contains the orthogonal matrix V;   
            if JOBV = 'V', V contains the product V1*V.   
            If JOBV = 'N', V is not referenced.   

    LDV     (input) INTEGER   
            The leading dimension of the array V. LDV >= MAX(1,P) if   
            JOBV = 'V'; LDV >= 1 otherwise.   

    Q       (input/output) LONG DOUBLE PRECISION array, dimension (LDQ,N)   
            On entry, if JOBQ = 'Q', Q must contain a matrix Q1 (usually 
  
            the orthogonal matrix returned by DGGSVP).   
            On exit,   
            if JOBQ = 'I', Q contains the orthogonal matrix Q;   
            if JOBQ = 'Q', Q contains the product Q1*Q.   
            If JOBQ = 'N', Q is not referenced.   

    LDQ     (input) INTEGER   
            The leading dimension of the array Q. LDQ >= MAX(1,N) if   
            JOBQ = 'Q'; LDQ >= 1 otherwise.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (2*N)   

    NCYCLE  (output) INTEGER   
            The number of cycles required for convergence.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            = 1:  the procedure does not converge after MAXIT cycles.   

    Internal Parameters   
    ===================   

    MAXIT   INTEGER   
            MAXIT specifies the total loops that the iterative procedure 
  
            may take. If after MAXIT cycles, the routine fails to   
            converge, we return INFO = 1.   

    Further Details   
    ===============   

    DTGSJA essentially uses a variant of Kogbetliantz algorithm to reduce 
  
    MIN(L,M-K)-by-L triangular (or trapezoidal) matrix A23 and L-by-L   
    matrix B13 to the form:   

             U1'*A13*Q1 = C1*R1; V1'*B13*Q1 = S1*R1,   

    where U1, V1 and Q1 are orthogonal matrix, and Z' is the transpose   
    of Z.  C1 and S1 are diagonal matrices satisfying   

                  C1**2 + S1**2 = I,   

    and R1 is an L-by-L nonsingular upper triangular matrix.   

    ===================================================================== 
  



       Decode and test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b13 = 0.;
    static LONG DOUBLE c_b14 = 1.;
    static int c__1 = 1;
    static LONG DOUBLE c_b43 = -1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    LONG DOUBLE d__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void drot_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qrot_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *);
    static int i, j;
    static LONG DOUBLE P_gamma;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, int *);
    static LONG DOUBLE a1;
    static long int initq;
    static LONG DOUBLE a2, a3, b1;
    static long int initu, initv, wantq, upper;
    static LONG DOUBLE b2, b3;
    static long int wantu, wantv;
    static LONG DOUBLE error, ssmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlags2_(long int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlags2(long int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlags2_(long int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *), dlapll_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *), qlapll(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *), qlapll_(int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *);
    static int kcycle;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg(LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlartg_(LONG DOUBLE *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), dlaset_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), qlaset(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *), qlaset_(char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *), xerbla_(char *, int *);
    static LONG DOUBLE csq, csu, csv, snq, rwk, snu, snv;



#define ALPHA(I) alpha[(I)-1]
#define BETA(I) beta[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define U(I,J) u[(I)-1 + ((J)-1)* ( *ldu)]
#define V(I,J) v[(I)-1 + ((J)-1)* ( *ldv)]
#define Q(I,J) q[(I)-1 + ((J)-1)* ( *ldq)]

    initu = lsame_(jobu, "I");
    wantu = initu || lsame_(jobu, "U");

    initv = lsame_(jobv, "I");
    wantv = initv || lsame_(jobv, "V");

    initq = lsame_(jobq, "I");
    wantq = initq || lsame_(jobq, "Q");

    *info = 0;
    if (! (initu || wantu || lsame_(jobu, "N"))) {
	*info = -1;
    } else if (! (initv || wantv || lsame_(jobv, "N"))) {
	*info = -2;
    } else if (! (initq || wantq || lsame_(jobq, "N"))) {
	*info = -3;
    } else if (*m < 0) {
	*info = -4;
    } else if (*p < 0) {
	*info = -5;
    } else if (*n < 0) {
	*info = -6;
    } else if (*lda < MAX(1,*m)) {
	*info = -10;
    } else if (*ldb < MAX(1,*p)) {
	*info = -12;
    } else if (*ldu < 1 || (wantu && *ldu < *m)) {
	*info = -18;
    } else if (*ldv < 1 || (wantv && *ldv < *p)) {
	*info = -20;
    } else if (*ldq < 1 || (wantq && *ldq < *n)) {
	*info = -22;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTGSJA", &i__1);
	return;
    }

/*     Initialize U, V and Q, if necessary */

    if (initu) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", m, m, &c_b13, &c_b14, &U(1,1), ldu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", m, m, &c_b13, &c_b14, &U(1,1), ldu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", m, m, &c_b13, &c_b14, &U(1,1), ldu);
#endif

    }
    if (initv) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", p, p, &c_b13, &c_b14, &V(1,1), ldv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", p, p, &c_b13, &c_b14, &V(1,1), ldv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", p, p, &c_b13, &c_b14, &V(1,1), ldv);
#endif

    }
    if (initq) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b13, &c_b14, &Q(1,1), ldq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b13, &c_b14, &Q(1,1), ldq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b13, &c_b14, &Q(1,1), ldq);
#endif

    }

/*     Loop until convergence */

    upper = 0;
    for (kcycle = 1; kcycle <= 40; ++kcycle) {

	upper = ! upper;

	i__1 = *l - 1;
	for (i = 1; i <= *l-1; ++i) {
	    i__2 = *l;
	    for (j = i + 1; j <= *l; ++j) {

		a1 = 0.;
		a2 = 0.;
		a3 = 0.;
		if (*k + i <= *m) {
		    a1 = A(*k+i,*n-*l+i);
		}
		if (*k + j <= *m) {
		    a3 = A(*k+j,*n-*l+j);
		}

		b1 = B(i,*n-*l+i);
		b3 = B(j,*n-*l+j);

		if (upper) {
		    if (*k + i <= *m) {
			a2 = A(*k+i,*n-*l+j);
		    }
		    b2 = B(i,*n-*l+j);
		} else {
		    if (*k + j <= *m) {
			a2 = A(*k+j,*n-*l+i);
		    }
		    b2 = B(j,*n-*l+i);
		}


#ifdef PETSC_PREFIX_SUFFIX
		dlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlags2(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlags2_(&upper, &a1, &a2, &a3, &b1, &b2, &b3, &csu, &snu, &
#endif

			csv, &snv, &csq, &snq);

/*              Update (K+I)-th and (K+J)-th rows of matrix A:
 U'*A */

		if (*k + j <= *m) {

#ifdef PETSC_PREFIX_SUFFIX
		    drot_(l, &A(*k+j,*n-*l+1), lda, &A(*k+i,*n-*l+1), lda, &csu, &snu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qrot(l, &A(*k+j,*n-*l+1), lda, &A(*k+i,*n-*l+1), lda, &csu, &snu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qrot_(l, &A(*k+j,*n-*l+1), lda, &A(*k+i,*n-*l+1), lda, &csu, &snu);
#endif

		}

/*              Update I-th and J-th rows of matrix B: V'*B */


#ifdef PETSC_PREFIX_SUFFIX
		drot_(l, &B(j,*n-*l+1), ldb, &B(i,*n-*l+1), ldb, &csv, &snv);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(l, &B(j,*n-*l+1), ldb, &B(i,*n-*l+1), ldb, &csv, &snv);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(l, &B(j,*n-*l+1), ldb, &B(i,*n-*l+1), ldb, &csv, &snv);
#endif


/*              Update (N-L+I)-th and (N-L+J)-th columns of ma
trices   
                A and B: A*Q and B*Q   

   Computing MIN */
		i__4 = *k + *l;
		i__3 = MIN(i__4,*m);

#ifdef PETSC_PREFIX_SUFFIX
		drot_(&i__3, &A(1,*n-*l+j), &c__1, &A(1,*n-*l+i), &c__1, &csq, &snq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(&i__3, &A(1,*n-*l+j), &c__1, &A(1,*n-*l+i), &c__1, &csq, &snq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(&i__3, &A(1,*n-*l+j), &c__1, &A(1,*n-*l+i), &c__1, &csq, &snq);
#endif



#ifdef PETSC_PREFIX_SUFFIX
		drot_(l, &B(1,*n-*l+j), &c__1, &B(1,*n-*l+i), &c__1, &csq, &snq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qrot(l, &B(1,*n-*l+j), &c__1, &B(1,*n-*l+i), &c__1, &csq, &snq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qrot_(l, &B(1,*n-*l+j), &c__1, &B(1,*n-*l+i), &c__1, &csq, &snq);
#endif


		if (upper) {
		    if (*k + i <= *m) {
			A(*k+i,*n-*l+j) = 0.;
		    }
		    B(i,*n-*l+j) = 0.;
		} else {
		    if (*k + j <= *m) {
			A(*k+j,*n-*l+i) = 0.;
		    }
		    B(j,*n-*l+i) = 0.;
		}

/*              Update orthogonal matrices U, V, Q, if desired
. */

		if (wantu && *k + j <= *m) {

#ifdef PETSC_PREFIX_SUFFIX
		    drot_(m, &U(1,*k+j), &c__1, &U(1,*k+i), &c__1, &csu, &snu);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qrot(m, &U(1,*k+j), &c__1, &U(1,*k+i), &c__1, &csu, &snu);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qrot_(m, &U(1,*k+j), &c__1, &U(1,*k+i), &c__1, &csu, &snu);
#endif

		}

		if (wantv) {

#ifdef PETSC_PREFIX_SUFFIX
		    drot_(p, &V(1,j), &c__1, &V(1,i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qrot(p, &V(1,j), &c__1, &V(1,i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qrot_(p, &V(1,j), &c__1, &V(1,i), &
#endif

			    c__1, &csv, &snv);
		}

		if (wantq) {

#ifdef PETSC_PREFIX_SUFFIX
		    drot_(n, &Q(1,*n-*l+j), &c__1, &Q(1,*n-*l+i), &c__1, &csq, &snq);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qrot(n, &Q(1,*n-*l+j), &c__1, &Q(1,*n-*l+i), &c__1, &csq, &snq);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qrot_(n, &Q(1,*n-*l+j), &c__1, &Q(1,*n-*l+i), &c__1, &csq, &snq);
#endif

		}

/* L10: */
	    }
/* L20: */
	}

	if (! upper) {

/*           The matrices A13 and B13 were lower triangular at the
 start   
             of the cycle, and are now upper triangular.   

             Convergence test: test the parallelism of the corresp
onding   
             rows of A and B. */

	    error = 0.;
/* Computing MIN */
	    i__2 = *l, i__3 = *m - *k;
	    i__1 = MIN(i__2,i__3);
	    for (i = 1; i <= MIN(*l,*m-*k); ++i) {
		i__2 = *l - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__2, &A(*k+i,*n-*l+i), lda, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__2, &A(*k+i,*n-*l+i), lda, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__2, &A(*k+i,*n-*l+i), lda, &WORK(
#endif

			1), &c__1);
		i__2 = *l - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__2, &B(i,*n-*l+i), ldb, &WORK(*l + 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__2, &B(i,*n-*l+i), ldb, &WORK(*l + 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__2, &B(i,*n-*l+i), ldb, &WORK(*l + 
#endif

			1), &c__1);
		i__2 = *l - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlapll_(&i__2, &WORK(1), &c__1, &WORK(*l + 1), &c__1, &ssmin);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlapll(&i__2, &WORK(1), &c__1, &WORK(*l + 1), &c__1, &ssmin);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlapll_(&i__2, &WORK(1), &c__1, &WORK(*l + 1), &c__1, &ssmin);
#endif

		error = MAX(error,ssmin);
/* L30: */
	    }

	    if (ABS(error) <= (LONG DOUBLE) (*n) * MIN(*tola,*tolb)) {
		goto L50;
	    }
	}

/*        End of cycle loop   

   L40: */
    }

/*     The algorithm has not converged after MAXIT cycles. */

    *info = 1;
    goto L100;

L50:

/*     If ERROR <= N*MIN(TOLA,TOLB), then the algorithm has converged.   
       Compute the generalized singular value pairs (ALPHA, BETA), and   
       set the triangular matrix R to array A. */

    i__1 = *k;
    for (i = 1; i <= *k; ++i) {
	ALPHA(i) = 1.;
	BETA(i) = 0.;
/* L60: */
    }

/* Computing MIN */
    i__2 = *l, i__3 = *m - *k;
    i__1 = MIN(i__2,i__3);
    for (i = 1; i <= MIN(*l,*m-*k); ++i) {

	a1 = A(*k+i,*n-*l+i);
	b1 = B(i,*n-*l+i);

	if (a1 != 0.) {
	    P_gamma = b1 / a1;

/*           change sign if necessary */

	    if (P_gamma < 0.) {
		i__2 = *l - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &c_b43, &B(i,*n-*l+i), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &c_b43, &B(i,*n-*l+i), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &c_b43, &B(i,*n-*l+i), ldb);
#endif

		if (wantv) {

#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(p, &c_b43, &V(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(p, &c_b43, &V(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(p, &c_b43, &V(1,i), &c__1);
#endif

		}
	    }

	    d__1 = ABS(P_gamma);

#ifdef PETSC_PREFIX_SUFFIX
	    dlartg_(&d__1, &c_b14, &BETA(*k + i), &ALPHA(*k + i), &rwk);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlartg(&d__1, &c_b14, &BETA(*k + i), &ALPHA(*k + i), &rwk);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlartg_(&d__1, &c_b14, &BETA(*k + i), &ALPHA(*k + i), &rwk);
#endif


	    if (ALPHA(*k + i) >= BETA(*k + i)) {
		i__2 = *l - i + 1;
		d__1 = 1. / ALPHA(*k + i);

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &d__1, &A(*k+i,*n-*l+i), lda)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &d__1, &A(*k+i,*n-*l+i), lda)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &d__1, &A(*k+i,*n-*l+i), lda)
#endif

			;
	    } else {
		i__2 = *l - i + 1;
		d__1 = 1. / BETA(*k + i);

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &d__1, &B(i,*n-*l+i), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &d__1, &B(i,*n-*l+i), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &d__1, &B(i,*n-*l+i), ldb);
#endif

		i__2 = *l - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
#endif

	    }

	} else {

	    ALPHA(*k + i) = 0.;
	    BETA(*k + i) = 1.;
	    i__2 = *l - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dcopy_(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qcopy(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qcopy_(&i__2, &B(i,*n-*l+i), ldb, &A(*k+i,*n-*l+i), lda);
#endif


	}

/* L70: */
    }

/*     Post-assignment */

    i__1 = *k + *l;
    for (i = *m + 1; i <= *k+*l; ++i) {
	ALPHA(i) = 0.;
	BETA(i) = 1.;
/* L80: */
    }

    if (*k + *l < *n) {
	i__1 = *n;
	for (i = *k + *l + 1; i <= *n; ++i) {
	    ALPHA(i) = 0.;
	    BETA(i) = 0.;
/* L90: */
	}
    }

L100:
    *ncycle = kcycle;
    return;

/*     End of DTGSJA */

} /* dtgsja_ */

