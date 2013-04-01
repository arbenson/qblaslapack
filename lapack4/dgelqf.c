
#include <math.h>
#define MIN(a,b)           ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)           ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)             ( ((a)<0.0)   ? -(a) : (a) )


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgelqf_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgelqf(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgelqf_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGELQF computes an LQ factorization of a real M-by-N matrix A:   
    A = L * Q.   

    Arguments   
    =========   

    M       (input) INT   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INT   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the M-by-N matrix A.   
            On exit, the elements on and below the diagonal of the array 
  
            contain the m-by-MIN(m,n) lower trapezoidal matrix L (L is   
            lower triangular if m <= n); the elements above the diagonal, 
  
            with the array TAU, represent the orthogonal matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INT   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INT   
            The dimension of the array WORK.  LWORK >= MAX(1,M).   
            For optimum performance LWORK >= M*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INT   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(k) . . . H(2) H(1), where k = MIN(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n), 
  
    and tau in TAU(i).   

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
    extern /* Subroutine */ void dgelq2_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgelq2(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgelq2_(int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int ib, nb;

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
    static int nx;

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
	xerbla_("DGELQF", &i__1);
	return;
    }

/*     Quick return if possible */

    k = MIN(*m,*n);
    if (k == 0) {
	WORK(1) = 1.;
	return;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DGELQF", " ", m, n, &c_n1, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DGELQF", " ", m, n, &c_n1, &c_n1, 6L,
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
		i__1 = 2, i__2 = ilaenv_(&c__2, "DGELQF", " ", m, n, &c_n1, &
			c_n1, 6L, 1L);
		nbmin = MAX(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < k && nx < k) {

/*        Use blocked code initially */

	i__1 = k - nx;
	i__2 = nb;
	for (i = 1; nb < 0 ? i >= k-nx : i <= k-nx; i += nb) {
/* Computing MIN */
	    i__3 = k - i + 1;
	    ib = MIN(i__3,nb);

/*           Compute the LQ factorization of the current block   
             A(i:i+ib-1,i:n) */

	    i__3 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgelq2_(&ib, &i__3, &A(i,i), lda, &TAU(i), &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgelq2(&ib, &i__3, &A(i,i), lda, &TAU(i), &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgelq2_(&ib, &i__3, &A(i,i), lda, &TAU(i), &WORK(1), &
#endif

		    iinfo);
	    if (i + ib <= *m) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__3 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarft_("Forward", "Rowwise", &i__3, &ib, &A(i,i), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarft("Forward", "Rowwise", &i__3, &ib, &A(i,i), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarft_("Forward", "Rowwise", &i__3, &ib, &A(i,i), 
#endif

			lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H to A(i+ib:m,i:n) from the right */

		i__3 = *m - i - ib + 1;
		i__4 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarfb_("Right", "No transpose", "Forward", "Rowwise", &i__3, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfb("Right", "No transpose", "Forward", "Rowwise", &i__3, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfb_("Right", "No transpose", "Forward", "Rowwise", &i__3, 
#endif

			&i__4, &ib, &A(i,i), lda, &WORK(1), &
			ldwork, &A(i+ib,i), lda, &WORK(ib + 1), &
			ldwork);
	    }
/* L10: */
	}
    } else {
	i = 1;
    }

/*     Use unblocked code to factor the last or only block. */

    if (i <= k) {
	i__2 = *m - i + 1;
	i__1 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dgelq2_(&i__2, &i__1, &A(i,i), lda, &TAU(i), &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgelq2(&i__2, &i__1, &A(i,i), lda, &TAU(i), &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgelq2_(&i__2, &i__1, &A(i,i), lda, &TAU(i), &WORK(1), &
#endif

		iinfo);
    }

    WORK(1) = (LONG DOUBLE) iws;
    return;

/*     End of DGELQF */

} /* dgelqf_ */

