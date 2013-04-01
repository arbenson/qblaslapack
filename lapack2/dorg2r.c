#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dorg2r_(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qorg2r(int *m, int *n, int *k, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qorg2r_(int *m, int *n, int *k, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *tau, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       February 29, 1992   


    Purpose   
    =======   

    DORG2R generates an m by n real matrix Q with orthonormal columns,   
    which is defined as the first n columns of a product of k elementary 
  
    reflectors of order m   

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
            On exit, the m-by-n matrix Q.   

    LDA     (input) INTEGER   
            The first dimension of the array A. LDA >= MAX(1,M).   

    TAU     (input) LONG DOUBLE PRECISION array, dimension (K)   
            TAU(i) must contain the scalar factor of the elementary   
            reflector H(i), as returned by DGEQRF.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument has an illegal value   

    ===================================================================== 
  


       Test the input arguments   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2;
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

	    int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *), xerbla_(char *, int *);



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
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DORG2R", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return;
    }

/*     Initialise columns k+1:n to columns of the unit matrix */

    i__1 = *n;
    for (j = *k + 1; j <= *n; ++j) {
	i__2 = *m;
	for (l = 1; l <= *m; ++l) {
	    A(l,j) = 0.;
/* L10: */
	}
	A(j,j) = 1.;
/* L20: */
    }

    for (i = *k; i >= 1; --i) {

/*        Apply H(i) to A(i:m,i:n) from the left */

	if (i < *n) {
	    A(i,i) = 1.;
	    i__1 = *m - i + 1;
	    i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
	    dlarf_("Left", &i__1, &i__2, &A(i,i), &c__1, &TAU(i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarf("Left", &i__1, &i__2, &A(i,i), &c__1, &TAU(i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarf_("Left", &i__1, &i__2, &A(i,i), &c__1, &TAU(i), &
#endif

		    A(i,i+1), lda, &WORK(1));
	}
	if (i < *m) {
	    i__1 = *m - i;
	    d__1 = -TAU(i);

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&i__1, &d__1, &A(i+1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&i__1, &d__1, &A(i+1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&i__1, &d__1, &A(i+1,i), &c__1);
#endif

	}
	A(i,i) = 1. - TAU(i);

/*        Set A(1:i-1,i) to zero */

	i__1 = i - 1;
	for (l = 1; l <= i-1; ++l) {
	    A(l,i) = 0.;
/* L30: */
	}
/* L40: */
    }
    return;

/*     End of DORG2R */

} /* dorg2r_ */

