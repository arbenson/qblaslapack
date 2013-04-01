#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dggrqf_(int *m, int *p, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qggrqf(int *m, int *p, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qggrqf_(int *m, int *p, int *n, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *taua, LONG DOUBLE *b, int *ldb, 
	LONG DOUBLE *taub, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGGRQF computes a generalized RQ factorization of an M-by-N matrix A 
  
    and a P-by-N matrix B:   

                A = R*Q,        B = Z*T*Q,   

    where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal   
    matrix, and R and T assume one of the forms:   

    if M <= N,  R = ( 0  R12 ) M,   or if M > N,  R = ( R11 ) M-N,   
                     N-M  M                           ( R21 ) N   
                                                         N   

    where R12 or R21 is upper triangular, and   

    if P >= N,  T = ( T11 ) N  ,   or if P < N,  T = ( T11  T12 ) P,   
                    (  0  ) P-N                         P   N-P   
                       N   

    where T11 is upper triangular.   

    In particular, if B is square and nonsingular, the GRQ factorization 
  
    of A and B implicitly gives the RQ factorization of A*inv(B):   

                 A*inv(B) = (R*inv(T))*Z'   

    where inv(B) denotes the inverse of the matrix B, and Z' denotes the 
  
    transpose of the matrix Z.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    P       (input) INTEGER   
            The number of rows of the matrix B.  P >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrices A and B. N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, if M <= N, the upper triangle of the subarray   
            A(1:M,N-M+1:N) contains the M-by-M upper triangular matrix R; 
  
            if M > N, the elements on and above the (M-N)-th subdiagonal 
  
            contain the M-by-N upper trapezoidal matrix R; the remaining 
  
            elements, with the array TAUA, represent the orthogonal   
            matrix Q as a product of elementary reflectors (see Further   
            Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,M).   

    TAUA    (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Q (see Further Details).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            On entry, the P-by-N matrix B.   
            On exit, the elements on and above the diagonal of the array 
  
            contain the MIN(P,N)-by-N upper trapezoidal matrix T (T is   
            upper triangular if P >= N); the elements below the diagonal, 
  
            with the array TAUB, represent the orthogonal matrix Z as a   
            product of elementary reflectors (see Further Details).   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,P).   

    TAUB    (output) LONG DOUBLE PRECISION array, dimension (MIN(P,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Z (see Further Details).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= MAX(1,N,M,P).   
            For optimum performance LWORK >= MAX(N,M,P)*MAX(NB1,NB2,NB3), 
  
            where NB1 is the optimal blocksize for the RQ factorization   
            of an M-by-N matrix, NB2 is the optimal blocksize for the   
            QR factorization of a P-by-N matrix, and NB3 is the optimal   
            blocksize for a call of DORMRQ.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INF0= -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = MIN(m,n).   

    Each H(i) has the form   

       H(i) = I - taua * v * v'   

    where taua is a real scalar, and v is a real vector with   
    v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in 
  
    A(m-k+i,1:n-k+i-1), and taua in TAUA(i).   
    To form Q explicitly, use LAPACK subroutine DORGRQ.   
    To use Q to update another matrix, use LAPACK subroutine DORMRQ.   

    The matrix Z is represented as a product of elementary reflectors   

       Z = H(1) H(2) . . . H(k), where k = MIN(p,n).   

    Each H(i) has the form   

       H(i) = I - taub * v * v'   

    where taub is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:p) is stored on exit in B(i+1:p,i), 
  
    and taub in TAUB(i).   
    To form Z explicitly, use LAPACK subroutine DORGQR.   
    To use Z to update another matrix, use LAPACK subroutine DORMQR.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int lopt;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgeqrf_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeqrf(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgeqrf_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dgerqf_(int *, int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgerqf(int *, int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgerqf_(int *, int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, int *), xerbla_(char *, int *), dormrq_(char *, char *, int *, int *, int *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, int *), xerbla_(char *, int *), qormrq(char *, char *, int *, int *, int *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, int *), xerbla_(char *, int *), qormrq_(char *, char *, int *, int *, int *,
#endif

	     LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, int *);


#define TAUA(I) taua[(I)-1]
#define TAUB(I) taub[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*p < 0) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*m)) {
	*info = -5;
    } else if (*ldb < MAX(1,*p)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = MAX(1,*m), i__1 = MAX(i__1,*p);
	if (*lwork < MAX(i__1,*n)) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGRQF", &i__1);
	return;
    }

/*     RQ factorization of M-by-N matrix A: A = R*Q */


#ifdef PETSC_PREFIX_SUFFIX
    dgerqf_(m, n, &A(1,1), lda, &TAUA(1), &WORK(1), lwork, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgerqf(m, n, &A(1,1), lda, &TAUA(1), &WORK(1), lwork, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgerqf_(m, n, &A(1,1), lda, &TAUA(1), &WORK(1), lwork, info);
#endif

    lopt = (int) WORK(1);

/*     Update B := B*Q' */

    i__1 = MIN(*m,*n);
/* Computing MAX */
    i__2 = 1, i__3 = *m - *n + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dormrq_("Right", "Transpose", p, n, &i__1, &A(MAX(1,*m-*n+1),1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormrq("Right", "Transpose", p, n, &i__1, &A(MAX(1,*m-*n+1),1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormrq_("Right", "Transpose", p, n, &i__1, &A(MAX(1,*m-*n+1),1), 
#endif

	    lda, &TAUA(1), &B(1,1), ldb, &WORK(1), lwork, info);
/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(1);
    lopt = MAX(i__1,i__2);

/*     QR factorization of P-by-N matrix B: B = Z*T */


#ifdef PETSC_PREFIX_SUFFIX
    dgeqrf_(p, n, &B(1,1), ldb, &TAUB(1), &WORK(1), lwork, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgeqrf(p, n, &B(1,1), ldb, &TAUB(1), &WORK(1), lwork, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgeqrf_(p, n, &B(1,1), ldb, &TAUB(1), &WORK(1), lwork, info);
#endif

/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(1);
    WORK(1) = (LONG DOUBLE) MAX(i__1,i__2);

    return;

/*     End of DGGRQF */

} /* dggrqf_ */

