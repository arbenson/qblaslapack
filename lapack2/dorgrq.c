#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dorgrq_(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qorgrq(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qorgrq_(int *m, int *n, int *k, LONG DOUBLE *
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

    DORGRQ generates an M-by-N real matrix Q with orthonormal rows,   
    which is defined as the last M rows of a product of K elementary   
    reflectors of order N   

          Q  =  H(1) H(2) . . . H(k)   

    as returned by DGERQF.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix Q. M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix Q. N >= M.   

    K       (input) INTEGER   
            The number of elementary reflectors whose product defines the 
  
            matrix Q. M >= K >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the (m-k+i)-th row must contain the vector which   
            defines the elementary reflector H(i), for i = 1,2,...,k, as 
  
            returned by DGERQF in the last k rows of its array argument   
            A.   
            On exit, the M-by-N matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= MAX(1,M).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGERQF.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK. LWORK >= MAX(1,M).   
            For optimum performance LWORK >= M*NB, where NB is the   
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
    int  i__1, i__2, i__3, i__4;
    /* Local variables */
    static int i, j, l, nbmin, iinfo;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dorgr2_(int *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qorgr2(int *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qorgr2_(int *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static int ib, nb, ii, kk;

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
    } else if (*n < *m) {
	*info = -2;
    } else if (*k < 0 || *k > *m) {
	*info = -3;
    } else if (*lda < MAX(1,*m)) {
	*info = -5;
    } else if (*lwork < MAX(1,*m)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGRQ", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m <= 0) {
	WORK(1) = 1.;
	return;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DORGRQ", " ", m, n, k, &c_n1, 6L, 1L);
    nbmin = 2;
    nx = 0;
    iws = *m;
    if (nb > 1 && nb < *k) {

/*        Determine when to cross over from blocked to unblocked code.
   

   Computing MAX */
	i__1 = 0, i__2 = ilaenv_(&c__3, "DORGRQ", " ", m, n, k, &c_n1, 6L, 1L)
		;
	nx = MAX(i__1,i__2);
	if (nx < *k) {

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
		i__1 = 2, i__2 = ilaenv_(&c__2, "DORGRQ", " ", m, n, k, &c_n1,
			 6L, 1L);
		nbmin = MAX(i__1,i__2);
	    }
	}
    }

    if (nb >= nbmin && nb < *k && nx < *k) {

/*        Use blocked code after the first block.   
          The last kk rows are handled by the block method.   

   Computing MIN */
	i__1 = *k, i__2 = (*k - nx + nb - 1) / nb * nb;
	kk = MIN(i__1,i__2);

/*        Set A(1:m-kk,n-kk+1:n) to zero. */

	i__1 = *n;
	for (j = *n - kk + 1; j <= *n; ++j) {
	    i__2 = *m - kk;
	    for (i = 1; i <= *m-kk; ++i) {
		A(i,j) = 0.;
/* L10: */
	    }
/* L20: */
	}
    } else {
	kk = 0;
    }

/*     Use unblocked code for the first or only block. */

    i__1 = *m - kk;
    i__2 = *n - kk;
    i__3 = *k - kk;

#ifdef PETSC_PREFIX_SUFFIX
    dorgr2_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo)
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qorgr2(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qorgr2_(&i__1, &i__2, &i__3, &A(1,1), lda, &TAU(1), &WORK(1), &iinfo)
#endif

	    ;

    if (kk > 0) {

/*        Use blocked code */

	i__1 = *k;
	i__2 = nb;
	for (i = *k - kk + 1; nb < 0 ? i >= *k : i <= *k; i += nb) {
/* Computing MIN */
	    i__3 = nb, i__4 = *k - i + 1;
	    ib = MIN(i__3,i__4);
	    ii = *m - *k + i;
	    if (ii > 1) {

/*              Form the triangular factor of the block reflec
tor   
                H = H(i+ib-1) . . . H(i+1) H(i) */

		i__3 = *n - *k + i + ib - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarft_("Backward", "Rowwise", &i__3, &ib, &A(ii,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarft("Backward", "Rowwise", &i__3, &ib, &A(ii,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarft_("Backward", "Rowwise", &i__3, &ib, &A(ii,1), 
#endif

			lda, &TAU(i), &WORK(1), &ldwork);

/*              Apply H' to A(1:m-k+i-1,1:n-k+i+ib-1) from the
 right */

		i__3 = ii - 1;
		i__4 = *n - *k + i + ib - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarfb_("Right", "Transpose", "Backward", "Rowwise", &i__3, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfb("Right", "Transpose", "Backward", "Rowwise", &i__3, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfb_("Right", "Transpose", "Backward", "Rowwise", &i__3, &
#endif

			i__4, &ib, &A(ii,1), lda, &WORK(1), &ldwork, &
			A(1,1), lda, &WORK(ib + 1), &ldwork);
	    }

/*           Apply H' to columns 1:n-k+i+ib-1 of current block */

	    i__3 = *n - *k + i + ib - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dorgr2_(&ib, &i__3, &ib, &A(ii,1), lda, &TAU(i), &WORK(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qorgr2(&ib, &i__3, &ib, &A(ii,1), lda, &TAU(i), &WORK(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qorgr2_(&ib, &i__3, &ib, &A(ii,1), lda, &TAU(i), &WORK(1), 
#endif

		    &iinfo);

/*           Set columns n-k+i+ib:n of current block to zero */

	    i__3 = *n;
	    for (l = *n - *k + i + ib; l <= *n; ++l) {
		i__4 = ii + ib - 1;
		for (j = ii; j <= ii+ib-1; ++j) {
		    A(j,l) = 0.;
/* L30: */
		}
/* L40: */
	    }
/* L50: */
	}
    }

    WORK(1) = (LONG DOUBLE) iws;
    return;

/*     End of DORGRQ */

} /* dorgrq_ */

