#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dorgqr_(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qorgqr(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qorgqr_(int *m, int *n, int *k, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *lwork, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DORGQR generates an M-by-N real matrix Q with orthonormal columns,   
    which is defined as the first N columns of a product of K elementary 
  
    reflectors of order M   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGEQRF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. M >= N >= 0.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. N >= K >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the i-th column must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGEQRF in the first k columns of its array   
            argument A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= MAX(1,M).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= MAX(1,N).   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument has an illegal value   

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
    int  i__1, i__2, i__3;
    /* Local variables */
    static int i, j, l, nbmin, iinfo;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorg2r_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorg2r(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorg2r_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *);
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
    } else if (*n < 0 || *n > *m) {
	*info = -2;
    } else if (*k < 0 || *k > *n) {
	*info = -3;
    } else if (*lda < MAX(1,*m)) {
	*info = -5;
    } else if (*lwork < MAX(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGQR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	WORK(1) = 1.;
	return;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *n;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DORGQR", " ", m, n, k, &c_n1, 6L, 1L)
		;
	nx = MAX(i__1,i__2);
	if (nx < *k) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  reduc
e NB and   
                determine the minimum value of NB. */

		nb = *lwork / ldwork;
/* Computing MAX */
		i__1 = 2, i__2 = ilaenv_(&c__2, "DORGQR", " ", m, n, k, &c_n1,
			 6L, 1L);
		nbmin = MAX(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the last block.   
          The first kk columns are handled by the block method. */

	ki = (*k - nx - 1) / nb * nb;
/* Computing MIN */
	i__1 = *k, i__2 = ki + nb;
	kk = MIN(i__1,i__2);

/*        Set A(1:kk,kk+1:n) to zero. */

	i__1 = *n;
	for (j = kk + 1; j <= *n; ++j) {
	    i__2 = kk;
	    for (i = 1; i <= kk; ++i) {
		A(i,j) = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the last or only block. */

    if (kk < *n) {
	i__1 = *m - kk;
	i__2 = *n - kk;
	i__3 = *k - kk;

#ifdef PETSC_PREFIX_SUFFIX
	dorg2r_(&i__1, &i__2, &i__3, &A(kk+1,kk+1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qorg2r(&i__1, &i__2, &i__3, &A(kk+1,kk+1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qorg2r_(&i__1, &i__2, &i__3, &A(kk+1,kk+1), lda, &
#endif

		TAU(kk + 1), &WORK(1), &iinfo);
    }

    if (kk > 0) {

/*        Use blocked code */

	i__1 = -nb;
	for (i = ki + 1; -nb < 0 ? i >= 1 : i <= 1; i += -nb) {
/* Computing MIN */
	    i__2 = nb, i__3 = *k - i + 1;
	    ib = MIN(i__2,i__3);
	    if (i + ib <= *n) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i) H(i+1) . . . H(i+ib-1) */

		i__2 = *m - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarft_("Forward", "Columnwise", &i__2, &ib, &A(i,i), lda, &TAU(i), &WORK(1), &ldwork);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarft("Forward", "Columnwise", &i__2, &ib, &A(i,i), lda, &TAU(i), &WORK(1), &ldwork);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarft_("Forward", "Columnwise", &i__2, &ib, &A(i,i), lda, &TAU(i), &WORK(1), &ldwork);
#endif


/*              Apply H to A(i:m,i+ib:n) from the left */

		i__2 = *m - i + 1;
		i__3 = *n - i - ib + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarfb_("Left", "No transpose", "Forward", "Columnwise", &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfb("Left", "No transpose", "Forward", "Columnwise", &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfb_("Left", "No transpose", "Forward", "Columnwise", &
#endif

			i__2, &i__3, &ib, &A(i,i), lda, &WORK(1), &
			ldwork, &A(i,i+ib), lda, &WORK(ib + 1),
			 &ldwork);
	    }

/*           Apply H to rows i:m of current block */

	    i__2 = *m - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dorg2r_(&i__2, &ib, &ib, &A(i,i), lda, &TAU(i), &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorg2r(&i__2, &ib, &ib, &A(i,i), lda, &TAU(i), &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorg2r_(&i__2, &ib, &ib, &A(i,i), lda, &TAU(i), &WORK(
#endif

		    1), &iinfo);

/*           Set rows 1:i-1 of current block to zero */

	    i__2 = i + ib - 1;
	    for (j = i; j <= i+ib-1; ++j) {
		i__3 = i - 1;
		for (l = 1; l <= i-1; ++l) {
		    A(l,j) = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    WORK(1) = (LONG DOUBLE) iws;
    return;

/*     End of DORGQR */

} /* dorgqr_ */

