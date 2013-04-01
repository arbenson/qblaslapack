#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgerqf_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgerqf(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgerqf_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGERQF computes an RQ factorization of a real M-by-N matrix A:   
    A = R * Q.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit,   
            if m <= n, the upper triangle of the subarray   
            A(1:m,n-m+1:n) contains the M-by-M upper triangular matrix R; 
  
            if m >= n, the elements on and above the (m-n)-th subdiagonal 
  
            contain the M-by-N upper trapezoidal matrix R;   
            the remaining elements, with the array TAU, represent the   
            orthogonal matrix Q as a product of MIN(m,n) elementary   
            reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= MAX(1,M).   
            For optimum performance LWORK >= M*NB, where NB is   
            the optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(1) H(2) . . . H(k), where k = MIN(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(n-k+i+1:n) = 0 and v(n-k+i) = 1; v(1:n-k+i-1) is stored on exit in 
  
    A(m-k+i,1:n-k+i-1), and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static int c__3 = 3;
    static int c__2 = 2;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    /* Local variables */
    static int i, k, nbmin, iinfo;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgerq2_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgerq2(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgerq2_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int ib, nb, ki, kk;

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
    static int mu, nu, nx;

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
    static int ldwork, iws;



#define TAU(I) tau[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (*m < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*m)) {
	*info = -4;
    } else if (*lwork < MAX(1,*m)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGERQF", &i__1);
	return;
    }

/*     Quick return if possible */

    k = MIN(*m,*n);
    if (k == 0) {
	WORK(1) = 1.;
	return;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DGERQF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 1;
    iws = *m;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGERQF", " ", m, n, &c_n1, &c_n1, 6L,
		 1L);
	nx = MAX(i__1,i__2);
	if (nx < k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *m;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGERQF", " ", m, n, &c_n1, &
			c_n1, 6L, 1L);
		nbmin = MAX(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially.   
          The last kk rows are handled by the block method. */

	ki = (k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = k, i__2 = ki + nb;
	kk = MIN(i__1,i__2);

	i__1 = k - kk + 1;
	i__2 = -nb;
	for (i = k - kk + ki + 1; -nb < 0 ? i >= k-kk+1 : i <= k-kk+1; i += -nb)
		 {
/* Computing MIN */
	    i__3 = k - i + 1;
	    ib = MIN(i__3,nb);

/*           Compute the RQ factorization of the current block   
             A(m-k+i:m-k+i+ib-1,1:n-k+i+ib-1) */

	    i__3 = *n - k + i + ib - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgerq2_(&ib, &i__3, &A(*m-k+i,1), lda, &TAU(i), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgerq2(&ib, &i__3, &A(*m-k+i,1), lda, &TAU(i), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgerq2_(&ib, &i__3, &A(*m-k+i,1), lda, &TAU(i), &WORK(
#endif

		    1), &iinfo);
	    if (*m - k + i > 1) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i+ib-1) . . . H(i+1) H(i) */

		i__3 = *n - k + i + ib - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarft_("Backward", "Rowwise", &i__3, &ib, &A(*m-k+i,1), lda, &TAU(i), &WORK(1), &ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarft("Backward", "Rowwise", &i__3, &ib, &A(*m-k+i,1), lda, &TAU(i), &WORK(1), &ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarft_("Backward", "Rowwise", &i__3, &ib, &A(*m-k+i,1), lda, &TAU(i), &WORK(1), &ldwork);
#endif


/*              Apply H to A(1:m-k+i-1,1:n-k+i+ib-1) from the 
right */

		i__3 = *m - k + i - 1;
		i__4 = *n - k + i + ib - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarfb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfb("Right", "No transpose", "Backward", "Rowwise", &i__3,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfb_("Right", "No transpose", "Backward", "Rowwise", &i__3,
#endif

			 &i__4, &ib, &A(*m-k+i,1), lda, &WORK(1), &
			ldwork, &A(1,1), lda, &WORK(ib + 1), &ldwork);
	    }
/* L10: */
	}
	mu = *m - k + i + nb - 1;
	nu = *n - k + i + nb - 1;
    } else {
	mu = *m;
	nu = *n;
    }

/*     Use unblocked code to factor the last or only block */

    if (mu > 0 && nu > 0) {

#ifdef PETSC_PREFIX_SUFFIX
	dgerq2_(&mu, &nu, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgerq2(&mu, &nu, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgerq2_(&mu, &nu, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo);
#endif

    }

    WORK(1) = (LONG DOUBLE) iws;
    return;

/*     End of DGERQF */

} /* dgerqf_ */

