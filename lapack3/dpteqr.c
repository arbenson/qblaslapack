#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpteqr_(char *compz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpteqr(char *compz, int *n, LONG DOUBLE *d, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpteqr_(char *compz, int *n, LONG DOUBLE *d, 
#endif

	LONG DOUBLE *e, LONG DOUBLE *z, int *ldz, LONG DOUBLE *work, int 
	*info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DPTEQR computes all eigenvalues and, optionally, eigenvectors of a   
    symmetric positive definite tridiagonal matrix by first factoring the 
  
    matrix using DPTTRF, and then calling DBDSQR to compute the singular 
  
    values of the bidiagonal factor.   

    This routine computes the eigenvalues of the positive definite   
    tridiagonal matrix to high relative accuracy.  This means that if the 
  
    eigenvalues range over many orders of magnitude in size, then the   
    small eigenvalues and corresponding eigenvectors will be computed   
    more accurately than, for example, with the standard QR method.   

    The eigenvectors of a full or band symmetric positive definite matrix 
  
    can also be found if DSYTRD, DSPTRD, or DSBTRD has been used to   
    reduce this matrix to tridiagonal form. (The reduction to tridiagonal 
  
    form, however, may preclude the possibility of obtaining high   
    relative accuracy in the small eigenvalues of the original matrix, if 
  
    these eigenvalues range over many orders of magnitude.)   

    Arguments   
    =========   

    COMPZ   (input) CHARACTER*1   
            = 'N':  Compute eigenvalues only.   
            = 'V':  Compute eigenvectors of original symmetric   
                    matrix also.  Array Z contains the orthogonal   
                    matrix used to reduce the original matrix to   
                    tridiagonal form.   
            = 'I':  Compute eigenvectors of tridiagonal matrix also.   

    N       (input) INTEGER   
            The order of the matrix.  N >= 0.   

    D       (input/output) LONG DOUBLE PRECISION array, dimension (N)   
            On entry, the n diagonal elements of the tridiagonal   
            matrix.   
            On normal exit, D contains the eigenvalues, in descending   
            order.   

    E       (input/output) LONG DOUBLE PRECISION array, dimension (N-1)   
            On entry, the (n-1) subdiagonal elements of the tridiagonal   
            matrix.   
            On exit, E has been destroyed.   

    Z       (input/output) LONG DOUBLE PRECISION array, dimension (LDZ, N)   
            On entry, if COMPZ = 'V', the orthogonal matrix used in the   
            reduction to tridiagonal form.   
            On exit, if COMPZ = 'V', the orthonormal eigenvectors of the 
  
            original symmetric matrix;   
            if COMPZ = 'I', the orthonormal eigenvectors of the   
            tridiagonal matrix.   
            If INFO > 0 on exit, Z contains the eigenvectors associated   
            with only the stored eigenvalues.   
            If  COMPZ = 'N', then Z is not referenced.   

    LDZ     (input) INTEGER   
            The leading dimension of the array Z.  LDZ >= 1, and if   
            COMPZ = 'V' or 'I', LDZ >= MAX(1,N).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LWORK)   
            If  COMPZ = 'N', then LWORK = 2*N   
            If  COMPZ = 'V' or 'I', then LWORK = MAX(1,4*N-4)   

    INFO    (output) INTEGER   
            = 0:  successful exit.   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   
            > 0:  if INFO = i, and i is:   
                  <= N  the Cholesky factorization of the matrix could   
                        not be performed because the i-th principal minor 
  
                        was not positive definite.   
                  > N   the SVD algorithm failed to converge;   
                        if INFO = N+i, i off-diagonal elements of the   
                        bidiagonal factor did not converge to zero.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b7 = 0.;
    static LONG DOUBLE c_b8 = 1.;
    static int c__0 = 0;
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE c[1]	/* was [1][1] */;
    static int i;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE vt[1]	/* was [1][1] */;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlaset_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlaset_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    xerbla_(char *, int *), dbdsqr_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    xerbla_(char *, int *), qbdsqr(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    xerbla_(char *, int *), qbdsqr_(char *, int *, 
#endif

	    int *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *);
    static int icompz;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dpttrf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qpttrf(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qpttrf_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *);
    static int nru;



#define C(I) c[(I)]
#define WAS(I) was[(I)]
#define VT(I) vt[(I)]
#define WAS(I) was[(I)]
#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define WORK(I) work[(I)-1]

#define Z(I,J) z[(I)-1 + ((J)-1)* ( *ldz)]

    *info = 0;

    if (lsame_(compz, "N")) {
	icompz = 0;
    } else if (lsame_(compz, "V")) {
	icompz = 1;
    } else if (lsame_(compz, "I")) {
	icompz = 2;
    } else {
	icompz = -1;
    }
    if (icompz < 0) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*ldz < 1 || (icompz > 0 && *ldz < MAX(1,*n))) {
	*info = -6;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPTEQR", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

    if (*n == 1) {
	if (icompz > 0) {
	    Z(1,1) = 1.;
	}
	return;
    }
    if (icompz == 2) {

#ifdef PETSC_PREFIX_SUFFIX
	dlaset_("Full", n, n, &c_b7, &c_b8, &Z(1,1), ldz);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaset("Full", n, n, &c_b7, &c_b8, &Z(1,1), ldz);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaset_("Full", n, n, &c_b7, &c_b8, &Z(1,1), ldz);
#endif

    }

/*     Call DPTTRF to factor the matrix. */


#ifdef PETSC_PREFIX_SUFFIX
    dpttrf_(n, &D(1), &E(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qpttrf(n, &D(1), &E(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qpttrf_(n, &D(1), &E(1), info);
#endif

    if (*info != 0) {
	return;
    }
    i__1 = *n;
    for (i = 1; i <= *n; ++i) {
	D(i) = sqrt(D(i));
/* L10: */
    }
    i__1 = *n - 1;
    for (i = 1; i <= *n-1; ++i) {
	E(i) *= D(i);
/* L20: */
    }

/*     Call DBDSQR to compute the singular values/vectors of the   
       bidiagonal factor. */

    if (icompz > 0) {
	nru = *n;
    } else {
	nru = 0;
    }

#ifdef PETSC_PREFIX_SUFFIX
    dbdsqr_("Lower", n, &c__0, &nru, &c__0, &D(1), &E(1), vt, &c__1, &Z(1,1), ldz, c, &c__1, &WORK(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qbdsqr("Lower", n, &c__0, &nru, &c__0, &D(1), &E(1), vt, &c__1, &Z(1,1), ldz, c, &c__1, &WORK(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qbdsqr_("Lower", n, &c__0, &nru, &c__0, &D(1), &E(1), vt, &c__1, &Z(1,1), ldz, c, &c__1, &WORK(1), info);
#endif


/*     Square the singular values. */

    if (*info == 0) {
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    D(i) *= D(i);
/* L30: */
	}
    } else {
	*info = *n + *info;
    }

    return;

/*     End of DPTEQR */

} /* dpteqr_ */

