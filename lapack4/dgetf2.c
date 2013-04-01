
#include <math.h>
#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgetf2_(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgetf2(int *m, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgetf2_(int *m, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *ipiv, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       June 30, 1992   


    Purpose   
    =======   

    DGETF2 computes an LU factorization of a general m-by-n matrix A   
    using partial pivoting with row interchanges.   

    The factorization has the form   
       A = P * L * U   
    where P is a permutation matrix, L is lower triangular with unit   
    diagonal elements (lower trapezoidal if m > n), and U is upper   
    triangular (upper trapezoidal if m < n).   

    This is the right-looking Level 2 BLAS version of the algorithm.   

    Arguments   
    =========   

    M       (input) INT   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INT   
            The number of columns of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the m by n matrix to be factored.   
            On exit, the factors L and U from the factorization   
            A = P*L*U; the unit diagonal elements of L are not stored.   

    LDA     (input) INT   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

    IPIV    (output) INT array, dimension (MIN(M,N))   
            The pivot indices; for 1 <= i <= MIN(M,N), row i of the   
            matrix was interchanged with row IPIV(i).   

    INFO    (output) INT   
            = 0: successful exit   
            < 0: if INFO = -k, the k-th argument had an illegal value   
            > 0: if INFO = k, U(k,k) is exactly zero. The factorization   
                 has been completed, but the factor U is exactly   
                 singular, and division by zero will occur if it is used 
  
                 to solve a system of equations.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b6 = -1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    LONG DOUBLE d__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dger_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qger(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qger_(int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static int j;

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
	    int *), dswap_(int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qswap(int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qswap_(int *, LONG DOUBLE *, int *, LONG DOUBLE 
#endif

	    *, int *);
    static int jp;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    extern /* Subroutine */ void xerbla_(char *, int *);



#define IPIV(I) ipiv[(I)-1]

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
	xerbla_("DGETF2", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*m == 0 || *n == 0) {
	return;
    }

    i__1 = MIN(*m,*n);
    for (j = 1; j <= MIN(*m,*n); ++j) {

/*        Find pivot and test for singularity. */

	i__2 = *m - j + 1;

#ifdef PETSC_PREFIX_SUFFIX
	jp = j - 1 + idamax_(&i__2, &A(j,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	jp = j - 1 + iqamax(&i__2, &A(j,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	jp = j - 1 + iqamax_(&i__2, &A(j,j), &c__1);
#endif

	IPIV(j) = jp;
	if (A(jp,j) != 0.) {

/*           Apply the interchange to columns 1:N. */

	    if (jp != j) {

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(n, &A(j,1), lda, &A(jp,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(n, &A(j,1), lda, &A(jp,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(n, &A(j,1), lda, &A(jp,1), lda);
#endif

	    }

/*           Compute elements J+1:M of J-th column. */

	    if (j < *m) {
		i__2 = *m - j;
		d__1 = 1. / A(j,j);

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &d__1, &A(j+1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &d__1, &A(j+1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &d__1, &A(j+1,j), &c__1);
#endif

	    }

	} else if (*info == 0) {

	    *info = j;
	}

	if (j < MIN(*m,*n)) {

/*           Update trailing submatrix. */

	    i__2 = *m - j;
	    i__3 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
	    dger_(&i__2, &i__3, &c_b6, &A(j+1,j), &c__1, &A(j,j+1), lda, &A(j+1,j+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qger(&i__2, &i__3, &c_b6, &A(j+1,j), &c__1, &A(j,j+1), lda, &A(j+1,j+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qger_(&i__2, &i__3, &c_b6, &A(j+1,j), &c__1, &A(j,j+1), lda, &A(j+1,j+1), lda);
#endif

	}
/* L10: */
    }
    return;

/*     End of DGETF2 */

} /* dgetf2_ */

