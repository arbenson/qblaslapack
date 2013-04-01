#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dorgr2_(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qorgr2(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qorgr2_(int *m, int *n, int *k, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORGR2 generates an m by n real matrix Q with orthonormal rows,   
    which is defined as the last m rows of a product of k elementary   
    reflectors of order n   

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
            On exit, the m by n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= MAX(1,M).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGERQF.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (M)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Local variables */
    static int i, j, l;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *), dlarf_(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qlarf(char *, int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qlarf_(char *, int *, int *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *);
    static int ii;
    extern /* Subroutine */ void xerbla_(char *, int *);


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
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORGR2", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m <= 0) {
	return;
    }

    if (*k < *m) {

/*        Initialise rows 1:m-k to rows of the unit matrix */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m - *k;
	    for (l = 1; l <= *m-*k; ++l) {
		A(l,j) = 0.;
/* L10: */
	    }
	    if (j > *n - *m && j <= *n - *k) {
		A(*m-*n+j,j) = 1.;
	    }
/* L20: */
	}
    }

    i__1 = *k;
    for (i = 1; i <= *k; ++i) {
	ii = *m - *k + i;

/*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the right */

	A(ii,*n-*m+ii) = 1.;
	i__2 = ii - 1;
	i__3 = *n - *m + ii;

#ifdef PETSC_PREFIX_SUFFIX
	dlarf_("Right", &i__2, &i__3, &A(ii,1), lda, &TAU(i), &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlarf("Right", &i__2, &i__3, &A(ii,1), lda, &TAU(i), &A(1,1), lda, &WORK(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlarf_("Right", &i__2, &i__3, &A(ii,1), lda, &TAU(i), &A(1,1), lda, &WORK(1));
#endif

	i__2 = *n - *m + ii - 1;
	d__1 = -TAU(i);

#ifdef PETSC_PREFIX_SUFFIX
	dscal_(&i__2, &d__1, &A(ii,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qscal(&i__2, &d__1, &A(ii,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qscal_(&i__2, &d__1, &A(ii,1), lda);
#endif

	A(ii,*n-*m+ii) = 1. - TAU(i);

/*        Set A(m-k+i,n-k+i+1:n) to zero */

	i__2 = *n;
	for (l = *n - *m + ii + 1; l <= *n; ++l) {
	    A(ii,l) = 0.;
/* L30: */
	}
/* L40: */
    }
    return;

/*     End of DORGR2 */

} /* dorgr2_ */

