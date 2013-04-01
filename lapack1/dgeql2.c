#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgeql2_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgeql2(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgeql2_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGEQL2 computes a QL factorization of a real m by n matrix A:   
    A = Q * L.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, if m >= n, the lower triangle of the subarray   
            A(m-n+1:m,1:n) contains the n by n lower triangular matrix L; 
  
            if m <= n, the elements on and below the (n-m)-th   
            superdiagonal contain the m by n lower trapezoidal matrix L; 
  
            the remaining elements, with the array TAU, represent the   
            orthogonal matrix Q as a product of elementary reflectors   
            (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    The matrix Q is represented as a product of elementary reflectors   

       Q = H(k) . . . H(2) H(1), where k = MIN(m,n).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(m-k+i+1:m) = 0 and v(m-k+i) = 1; v(1:m-k+i-1) is stored on exit in 
  
    A(1:m-k+i-1,n-k+i), and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static int i, k;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlarf_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarf(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlarf_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *), dlarfg_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *), qlarfg(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *), qlarfg_(int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *), xerbla_(char *, int *);
    static LONG DOUBLE aii;



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
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGEQL2", &i__1);
	return;
    }

    k = MIN(*m,*n);

    for (i = k; i >= 1; --i) {

/*        Generate elementary reflector H(i) to annihilate   
          A(1:m-k+i-1,n-k+i) */

	i__1 = *m - k + i;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&i__1, &A(*m-k+i,*n-k+i), &A(1,*n-k+i), &c__1, &TAU(i));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&i__1, &A(*m-k+i,*n-k+i), &A(1,*n-k+i), &c__1, &TAU(i));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&i__1, &A(*m-k+i,*n-k+i), &A(1,*n-k+i), &c__1, &TAU(i));
#endif


/*        Apply H(i) to A(1:m-k+i,1:n-k+i-1) from the left */

	aii = A(*m-k+i,*n-k+i);
	A(*m-k+i,*n-k+i) = 1.;
	i__1 = *m - k + i;
	i__2 = *n - k + i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlarf_("Left", &i__1, &i__2, &A(1,*n-k+i), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarf("Left", &i__1, &i__2, &A(1,*n-k+i), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarf_("Left", &i__1, &i__2, &A(1,*n-k+i), &c__1, &
#endif

		TAU(i), &A(1,1), lda, &WORK(1));
	A(*m-k+i,*n-k+i) = aii;
/* L10: */
    }
    return;

/*     End of DGEQL2 */

} /* dgeql2_ */

