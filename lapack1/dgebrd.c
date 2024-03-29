#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgebrd_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgebrd(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgebrd_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *d, LONG DOUBLE *e, LONG DOUBLE *tauq, LONG DOUBLE *taup,
	 LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGEBRD reduces a general real M-by-N matrix A to upper or lower   
    bidiagonal form B by an orthogonal transformation: Q**T * A * P = B. 
  

    If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows in the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns in the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N general matrix to be reduced.   
            On exit,   
            if m >= n, the diagonal and the first superdiagonal are   
              overwritten with the upper bidiagonal matrix B; the   
              elements below the diagonal, with the array TAUQ, represent 
  
              the orthogonal matrix Q as a product of elementary   
              reflectors, and the elements above the first superdiagonal, 
  
              with the array TAUP, represent the orthogonal matrix P as   
              a product of elementary reflectors;   
            if m < n, the diagonal and the first subdiagonal are   
              overwritten with the lower bidiagonal matrix B; the   
              elements below the first subdiagonal, with the array TAUQ, 
  
              represent the orthogonal matrix Q as a product of   
              elementary reflectors, and the elements above the diagonal, 
  
              with the array TAUP, represent the orthogonal matrix P as   
              a product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    D       (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The diagonal elements of the bidiagonal matrix B:   
            D(i) = A(i,i).   

    E       (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N)-1)   
            The off-diagonal elements of the bidiagonal matrix B:   
            if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;   
            if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.   

    TAUQ    (output) LONG DOUBLE PRECISION array dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix Q. See Further Details.   

    TAUP    (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors which   
            represent the orthogonal matrix P. See Further Details.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of the array WORK.  LWORK >= MAX(1,M,N).   
            For optimum performance LWORK >= (M+N)*NB, where NB   
            is the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

    Further Details   
    ===============   

    The matrices Q and P are represented as products of elementary   
    reflectors:   

    If m >= n,   

       Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are real scalars, and v and u are real vectors;   
    v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i); 
  
    u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n); 
  
    tauq is stored in TAUQ(i) and taup in TAUP(i).   

    If m < n,   

       Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)   

    Each H(i) and G(i) has the form:   

       H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'   

    where tauq and taup are real scalars, and v and u are real vectors;   
    v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i); 
  
    u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n); 
  
    tauq is stored in TAUQ(i) and taup in TAUP(i).   

    The contents of A on exit are illustrated by the following examples: 
  

    m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):   

      (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )   
      (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )   
      (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )   
      (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )   
      (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )   
      (  v1  v2  v3  v4  v5 )   

    where d and e denote diagonal and off-diagonal elements of B, vi   
    denotes an element of the vector defining H(i), and ui an element of 
  
    the vector defining G(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static int c__3 = 3;
    static int c__2 = 2;
    static LONG DOUBLE c_b21 = -1.;
    static LONG DOUBLE c_b22 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    /* Local variables */
    static int i, j;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemm_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemm_(char *, char *, int *, int *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int nbmin, iinfo, minmn;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgebd2_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgebd2(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgebd2_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     LONG DOUBLE *, int *);
    static int nb;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlabrd_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabrd(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlabrd_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static int nx;
    static LONG DOUBLE ws;
    extern /* Subroutine */ void xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);
    static int ldwrkx, ldwrky;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAUQ(I) tauq[(I)-1]
#define TAUP(I) taup[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*m)) {
	*info = -4;
    } else /* if(complicated condition) */ {
/* Computing MAX */
	i__1 = MAX(1,*m);
	if (*lwork < MAX(i__1,*n)) {
	    *info = -10;
	}
    }
    if (*info < 0) {
	i__1 = -(*info);
	xerbla_("DGEBRD", &i__1);
	return;
    }

/*     Quick return if possible */

    minmn = MIN(*m,*n);
    if (minmn == 0) {
	WORK(1) = 1.;
	return;
    }

    ws = (LONG DOUBLE) MAX(*m,*n);
    ldwrkx = *m;
    ldwrky = *n;

/*     Set the block size NB and the crossover point NX.   

   Computing MAX */
    i__1 = 1, i__2 = ilaenv_(&c__1, "DGEBRD", " ", m, n, &c_n1, &c_n1, 6L, 1L)
	    ;
    nb = MAX(i__1,i__2);

    if (nb > 1 && nb < minmn) {

/*        Determine when to switch from blocked to unblocked code.   

   Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__3, "DGEBRD", " ", m, n, &c_n1, &c_n1, 
		6L, 1L);
	nx = MAX(i__1,i__2);
	if (nx < minmn) {
	    ws = (LONG DOUBLE) ((*m + *n) * nb);
	    if ((LONG DOUBLE) (*lwork) < ws) {

/*              Not enough work space for the optimal NB, cons
ider using   
                a smaller block size. */

		nbmin = ilaenv_(&c__2, "DGEBRD", " ", m, n, &c_n1, &c_n1, 6L, 
			1L);
		if (*lwork >= (*m + *n) * nbmin) {
		    nb = *lwork / (*m + *n);
		} else {
		    nb = 1;
		    nx = minmn;
		}
	    }
	}
    } else {
	nx = minmn;
    }

    i__1 = minmn - nx;
    i__2 = nb;
    for (i = 1; nb < 0 ? i >= minmn-nx : i <= minmn-nx; i += nb) {

/*        Reduce rows and columns i:i+nb-1 to bidiagonal form and retu
rn   
          the matrices X and Y which are needed to update the unreduce
d   
          part of the matrix */

	i__3 = *m - i + 1;
	i__4 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlabrd_(&i__3, &i__4, &nb, &A(i,i), lda, &D(i), &E(i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlabrd(&i__3, &i__4, &nb, &A(i,i), lda, &D(i), &E(i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlabrd_(&i__3, &i__4, &nb, &A(i,i), lda, &D(i), &E(i), &
#endif

		TAUQ(i), &TAUP(i), &WORK(1), &ldwrkx, &WORK(ldwrkx * nb + 1), 
		&ldwrky);

/*        Update the trailing submatrix A(i+nb:m,i+nb:n), using an upd
ate   
          of the form  A := A - V*Y' - X*U' */

	i__3 = *m - i - nb + 1;
	i__4 = *n - i - nb + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgemm_("No transpose", "Transpose", &i__3, &i__4, &nb, &c_b21, &A(i+nb,i), lda, &WORK(ldwrkx * nb + nb + 1), &ldwrky, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemm("No transpose", "Transpose", &i__3, &i__4, &nb, &c_b21, &A(i+nb,i), lda, &WORK(ldwrkx * nb + nb + 1), &ldwrky, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemm_("No transpose", "Transpose", &i__3, &i__4, &nb, &c_b21, &A(i+nb,i), lda, &WORK(ldwrkx * nb + nb + 1), &ldwrky, &
#endif

		c_b22, &A(i+nb,i+nb), lda);
	i__3 = *m - i - nb + 1;
	i__4 = *n - i - nb + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgemm_("No transpose", "No transpose", &i__3, &i__4, &nb, &c_b21, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgemm("No transpose", "No transpose", &i__3, &i__4, &nb, &c_b21, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgemm_("No transpose", "No transpose", &i__3, &i__4, &nb, &c_b21, &
#endif

		WORK(nb + 1), &ldwrkx, &A(i,i+nb), lda, &c_b22,
		 &A(i+nb,i+nb), lda);

/*        Copy diagonal and off-diagonal elements of B back into A */

	if (*m >= *n) {
	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		A(j,j) = D(j);
		A(j,j+1) = E(j);
/* L10: */
	    }
	} else {
	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		A(j,j) = D(j);
		A(j+1,j) = E(j);
/* L20: */
	    }
	}
/* L30: */
    }

/*     Use unblocked code to reduce the remainder of the matrix */

    i__2 = *m - i + 1;
    i__1 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
    dgebd2_(&i__2, &i__1, &A(i,i), lda, &D(i), &E(i), &TAUQ(i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qgebd2(&i__2, &i__1, &A(i,i), lda, &D(i), &E(i), &TAUQ(i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qgebd2_(&i__2, &i__1, &A(i,i), lda, &D(i), &E(i), &TAUQ(i), &
#endif

	    TAUP(i), &WORK(1), &iinfo);
    WORK(1) = ws;
    return;

/*     End of DGEBRD */

} /* dgebrd_ */

