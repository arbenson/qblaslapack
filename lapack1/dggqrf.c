#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dggqrf_(int *n, int *m, int *p, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qggqrf(int *n, int *m, int *p, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qggqrf_(int *n, int *m, int *p, LONG DOUBLE *
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

    DGGQRF computes a generalized QR factorization of an N-by-M matrix A 
  
    and an N-by-P matrix B:   

                A = Q*R,        B = Q*T*Z,   

    where Q is an N-by-N orthogonal matrix, Z is a P-by-P orthogonal   
    matrix, and R and T assume one of the forms:   

    if N >= M,  R = ( R11 ) M  ,   or if N < M,  R = ( R11  R12 ) N,   
                    (  0  ) N-M                         N   M-N   
                       M   

    where R11 is upper triangular, and   

    if N <= P,  T = ( 0  T12 ) N,   or if N > P,  T = ( T11 ) N-P,   
                     P-N  N                           ( T21 ) P   
                                                         P   

    where T12 or T21 is upper triangular.   

    In particular, if B is square and nonsingular, the GQR factorization 
  
    of A and B implicitly gives the QR factorization of inv(B)*A:   

                 inv(B)*A = Z'*(inv(T)*R)   

    where inv(B) denotes the inverse of the matrix B, and Z' denotes the 
  
    transpose of the matrix Z.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The number of rows of the matrices A and B. N >= 0.   

    M       (input) INTEGER   
            The number of columns of the matrix A.  M >= 0.   

    P       (input) INTEGER   
            The number of columns of the matrix B.  P >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,M)   
            On entry, the N-by-M matrix A.   
            On exit, the elements on and above the diagonal of the array 
  
            contain the MIN(N,M)-by-M upper trapezoidal matrix R (R is   
            upper triangular if N >= M); the elements below the diagonal, 
  
            with the array TAUA, represent the orthogonal matrix Q as a   
            product of MIN(N,M) elementary reflectors (see Further   
            Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A. LDA >= MAX(1,N).   

    TAUA    (output) LONG DOUBLE PRECISION array, dimension (MIN(N,M))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Q (see Further Details).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,P)   
            On entry, the N-by-P matrix B.   
            On exit, if N <= P, the upper triangle of the subarray   
            B(1:N,P-N+1:P) contains the N-by-N upper triangular matrix T; 
  
            if N > P, the elements on and above the (N-P)-th subdiagonal 
  
            contain the N-by-P upper trapezoidal matrix T; the remaining 
  
            elements, with the array TAUB, represent the orthogonal   
            matrix Z as a product of elementary reflectors (see Further   
            Details).   

    LDB     (input) INTEGER   
            The leading dimension of the array B. LDB >= MAX(1,N).   

    TAUB    (output) LONG DOUBLE PRECISION array, dimension (MIN(N,P))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Z (see Further Details).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= MAX(1,N,M,P).   
            For optimum performance LWORK >= MAX(N,M,P)*MAX(NB1,NB2,NB3), 
  
            where NB1 is the optimal blocksize for the QR factorization   
            of an N-by-M matrix, NB2 is the optimal blocksize for the   
            RQ factorization of an N-by-P matrix, and NB3 is the optimal 
  
            blocksize for a call of DORMQR.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = MIN(n,m).   

    Each H(i) has the form   

       H(i) = I - taua * v * v'   

    where taua is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i+1:n,i), 
  
    and taua in TAUA(i).   
    To form Q explicitly, use LAPACK subroutine DORGQR.   
    To use Q to update another matrix, use LAPACK subroutine DORMQR.   

    The matrix Z is represented as a product of elementary reflectors   

       Z = H(1) H(2) . . . H(k), where k = MIN(n,p).   

    Each H(i) has the form   

       H(i) = I - taub * v * v'   

    where taub is a real scalar, and v is a real vector with   
    v(p-k+i+1:p) = 0 and v(p-k+i) = 1; v(1:p-k+i-1) is stored on exit in 
  
    B(n-k+i,1:p-k+i-1), and taub in TAUB(i).   
    To form Z explicitly, use LAPACK subroutine DORGRQ.   
    To use Z to update another matrix, use LAPACK subroutine DORMRQ.   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2;
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
	    *, LONG DOUBLE *, int *, int *), xerbla_(char *, int *), dormqr_(char *, char *, int *, int *, int *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, int *), xerbla_(char *, int *), qormqr(char *, char *, int *, int *, int *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    *, LONG DOUBLE *, int *, int *), xerbla_(char *, int *), qormqr_(char *, char *, int *, int *, int *,
#endif

	     LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, int *);


#define TAUA(I) taua[(I)-1]
#define TAUB(I) taub[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*m < 0) {
	*info = -2;
    } else if (*p < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = MAX(1,*n), i__1 = MAX(i__1,*m);
	if (*lwork < MAX(i__1,*p)) {
	    *info = -11;
	}
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGGQRF", &i__1);
	return;
    }

/*     QR factorization of N-by-M matrix A: A = Q*R */


#ifdef PETSC_PREFIX_SUFFIX
    dgeqrf_(n, m, &A(1,1), lda, &TAUA(1), &WORK(1), lwork, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgeqrf(n, m, &A(1,1), lda, &TAUA(1), &WORK(1), lwork, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgeqrf_(n, m, &A(1,1), lda, &TAUA(1), &WORK(1), lwork, info);
#endif

    lopt = (int) WORK(1);

/*     Update B := Q'*B. */

    i__1 = MIN(*n,*m);

#ifdef PETSC_PREFIX_SUFFIX
    dormqr_("Left", "Transpose", n, p, &i__1, &A(1,1), lda, &TAUA(1), &B(1,1), ldb, &WORK(1), lwork, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qormqr("Left", "Transpose", n, p, &i__1, &A(1,1), lda, &TAUA(1), &B(1,1), ldb, &WORK(1), lwork, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qormqr_("Left", "Transpose", n, p, &i__1, &A(1,1), lda, &TAUA(1), &B(1,1), ldb, &WORK(1), lwork, info);
#endif

/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(1);
    lopt = MAX(i__1,i__2);

/*     RQ factorization of N-by-P matrix B: B = T*Z. */


#ifdef PETSC_PREFIX_SUFFIX
    dgerqf_(n, p, &B(1,1), ldb, &TAUB(1), &WORK(1), lwork, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgerqf(n, p, &B(1,1), ldb, &TAUB(1), &WORK(1), lwork, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgerqf_(n, p, &B(1,1), ldb, &TAUB(1), &WORK(1), lwork, info);
#endif

/* Computing MAX */
    i__1 = lopt, i__2 = (int) WORK(1);
    WORK(1) = (LONG DOUBLE) MAX(i__1,i__2);

    return;

/*     End of DGGQRF */

} /* dggqrf_ */

