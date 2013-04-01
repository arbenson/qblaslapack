#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgelq2_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgelq2(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgelq2_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DGELQ2 computes an LQ factorization of a real m by n matrix A:   
    A = L * Q.   

    Arguments   
    =========   

    M       (input) INTEGER   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INTEGER   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix A.   
            On exit, the elements on and below the diagonal of the array 
  
            contain the m by MIN(m,n) lower trapezoidal matrix L (L is   
            lower triangular if m <= n); the elements above the diagonal, 
  
            with the array TAU, represent the orthogonal matrix Q as a   
            product of elementary reflectors (see Further Details).   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (MIN(M,N))   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (M)   

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
    v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n), 
  
    and tau in TAU(i).   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3;
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
	xerbla_("DGELQ2", &i__1);
	return;
    }

    k = MIN(*m,*n);

    i__1 = k;
    for (i = 1; i <= k; ++i) {

/*        Generate elementary reflector H(i) to annihilate A(i,i+1:n) 
*/

	i__2 = *n - i + 1;
/* Computing MIN */
	i__3 = i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dlarfg_(&i__2, &A(i,i), &A(i,MIN(i+1,*n)), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarfg(&i__2, &A(i,i), &A(i,MIN(i+1,*n)), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarfg_(&i__2, &A(i,i), &A(i,MIN(i+1,*n)), lda,
#endif

		 &TAU(i));
	if (i < *m) {

/*           Apply H(i) to A(i+1:m,i:n) from the right */

	    aii = A(i,i);
	    A(i,i) = 1.;
	    i__2 = *m - i;
	    i__3 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlarf_("Right", &i__2, &i__3, &A(i,i), lda, &TAU(i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarf("Right", &i__2, &i__3, &A(i,i), lda, &TAU(i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarf_("Right", &i__2, &i__3, &A(i,i), lda, &TAU(i), &
#endif

		    A(i+1,i), lda, &WORK(1));
	    A(i,i) = aii;
	}
/* L10: */
    }
    return;

/*     End of DGELQ2 */

} /* dgelq2_ */

