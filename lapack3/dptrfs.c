#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dptrfs_(int *n, int *nrhs, LONG DOUBLE *d, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qptrfs(int *n, int *nrhs, LONG DOUBLE *d, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qptrfs_(int *n, int *nrhs, LONG DOUBLE *d, 
#endif

	LONG DOUBLE *e, LONG DOUBLE *df, LONG DOUBLE *ef, LONG DOUBLE *b, int 
	*ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *ferr, LONG DOUBLE *berr,
	 LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DPTRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is symmetric positive definite 
  
    and tridiagonal, and provides error bounds and backward error   
    estimates for the solution.   

    Arguments   
    =========   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the tridiagonal matrix A.   

    E       (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) subdiagonal elements of the tridiagonal matrix A.   

    DF      (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the diagonal matrix D from the   
            factorization computed by DPTTRF.   

    EF      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) subdiagonal elements of the unit bidiagonal factor 
  
            L from the factorization computed by DPTTRF.   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (input/output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by DPTTRS.   
            On exit, the improved solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= MAX(1,N).   

    FERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).   

    BERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (2*N)   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Internal Parameters   
    ===================   

    ITMAX is the maximum number of steps of iterative refinement.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b11 = 1.;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static LONG DOUBLE safe1, safe2;
    static int i, j;
    static LONG DOUBLE s;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static int count;
    static LONG DOUBLE bi;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif

    static LONG DOUBLE cx, dx, ex;
    static int ix;

#ifdef PETSC_PREFIX_SUFFIX
    extern int idamax_(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern int iqamax(int *, LONG DOUBLE *, int *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern int iqamax_(int *, LONG DOUBLE *, int *);
#endif

    static int nz;
    static LONG DOUBLE safmin;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE lstres;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dpttrs_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qpttrs(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qpttrs_(int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, int *);
    static LONG DOUBLE eps;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define DF(I) df[(I)-1]
#define EF(I) ef[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    if (*n < 0) {
	*info = -1;
    } else if (*nrhs < 0) {
	*info = -2;
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    } else if (*ldx < MAX(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPTRFS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {
	    FERR(j) = 0.;
	    BERR(j) = 0.;
/* L10: */
	}
	return;
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

    nz = 4;

#ifdef PETSC_PREFIX_SUFFIX
    eps = dlamch_("Epsilon");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    eps = qlamch("Epsilon");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    eps = qlamch_("Epsilon");
#endif


#ifdef PETSC_PREFIX_SUFFIX
    safmin = dlamch_("Safe minimum");
#endif
#ifdef Q_C_PREFIX_SUFFIX
    safmin = qlamch("Safe minimum");
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    safmin = qlamch_("Safe minimum");
#endif

    safe1 = nz * safmin;
    safe2 = safe1 / eps;

/*     Do for each right hand side */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

	count = 1;
	lstres = 3.;
L20:

/*        Loop until stopping criterion is satisfied.   

          Compute residual R = B - A * X.  Also compute   
          ABS(A)*ABS(x) + ABS(b) for use in the backward error bound. 
*/

	if (*n == 1) {
	    bi = B(1,j);
	    dx = D(1) * X(1,j);
	    WORK(*n + 1) = bi - dx;
	    WORK(1) = ABS(bi) + ABS(dx);
	} else {
	    bi = B(1,j);
	    dx = D(1) * X(1,j);
	    ex = E(1) * X(2,j);
	    WORK(*n + 1) = bi - dx - ex;
	    WORK(1) = ABS(bi) + ABS(dx) + ABS(ex);
	    i__2 = *n - 1;
	    for (i = 2; i <= *n-1; ++i) {
		bi = B(i,j);
		cx = E(i - 1) * X(i-1,j);
		dx = D(i) * X(i,j);
		ex = E(i) * X(i+1,j);
		WORK(*n + i) = bi - cx - dx - ex;
		WORK(i) = ABS(bi) + ABS(cx) + ABS(dx) + ABS(ex);
/* L30: */
	    }
	    bi = B(*n,j);
	    cx = E(*n - 1) * X(*n-1,j);
	    dx = D(*n) * X(*n,j);
	    WORK(*n + *n) = bi - cx - dx;
	    WORK(*n) = ABS(bi) + ABS(cx) + ABS(dx);
	}

/*        Compute componentwise relative backward error from formula 
  

          MAX(i) ( ABS(R(i)) / ( ABS(A)*ABS(X) + ABS(B) )(i) )   

          where ABS(Z) is the componentwise absolute value of the matr
ix   
          or vector Z.  If the i-th component of the denominator is le
ss   
          than SAFE2, then SAFE1 is added to the i-th components of th
e   
          numerator and denominator before dividing. */

	s = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
/* Computing MAX */
		d__2 = s, d__3 = (d__1 = WORK(*n + i), ABS(d__1)) / WORK(i);
		s = MAX(d__2,d__3);
	    } else {
/* Computing MAX */
		d__2 = s, d__3 = ((d__1 = WORK(*n + i), ABS(d__1)) + safe1) / 
			(WORK(i) + safe1);
		s = MAX(d__2,d__3);
	    }
/* L40: */
	}
	BERR(j) = s;

/*        Test stopping criterion. Continue iterating if   
             1) The residual BERR(J) is larger than machine epsilon, a
nd   
             2) BERR(J) decreased by at least a factor of 2 during the
   
                last iteration, and   
             3) At most ITMAX iterations tried. */

	if (BERR(j) > eps && BERR(j) * 2. <= lstres && count <= 5) {

/*           Update solution and try again. */


#ifdef PETSC_PREFIX_SUFFIX
	    dpttrs_(n, &c__1, &DF(1), &EF(1), &WORK(*n + 1), n, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qpttrs(n, &c__1, &DF(1), &EF(1), &WORK(*n + 1), n, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qpttrs_(n, &c__1, &DF(1), &EF(1), &WORK(*n + 1), n, info);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    daxpy_(n, &c_b11, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qaxpy(n, &c_b11, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qaxpy_(n, &c_b11, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif

		    ;
	    lstres = BERR(j);
	    ++count;
	    goto L20;
	}

/*        Bound error from formula   

          norm(X - XTRUE) / norm(X) .le. FERR =   
          norm( ABS(inv(A))*   
             ( ABS(R) + NZ*EPS*( ABS(A)*ABS(X)+ABS(B) ))) / norm(X)   

          where   
            norm(Z) is the magnitude of the largest component of Z   
            inv(A) is the inverse of A   
            ABS(Z) is the componentwise absolute value of the matrix o
r   
               vector Z   
            NZ is the maximum number of nonzeros in any row of A, plus
 1   
            EPS is machine epsilon   

          The i-th component of ABS(R)+NZ*EPS*(ABS(A)*ABS(X)+ABS(B)) 
  
          is incremented by SAFE1 if the i-th component of   
          ABS(A)*ABS(X) + ABS(B) is less than SAFE2. */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
		WORK(i) = (d__1 = WORK(*n + i), ABS(d__1)) + nz * eps * WORK(
			i);
	    } else {
		WORK(i) = (d__1 = WORK(*n + i), ABS(d__1)) + nz * eps * WORK(
			i) + safe1;
	    }
/* L50: */
	}

#ifdef PETSC_PREFIX_SUFFIX
	ix = idamax_(n, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	ix = iqamax(n, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	ix = iqamax_(n, &WORK(1), &c__1);
#endif

	FERR(j) = WORK(ix);

/*        Estimate the norm of inv(A).   

          Solve M(A) * x = e, where M(A) = (m(i,j)) is given by   

             m(i,j) =  ABS(A(i,j)), i = j,   
             m(i,j) = -ABS(A(i,j)), i .ne. j,   

          and e = [ 1, 1, ..., 1 ]'.  Note M(A) = M(L)*D*M(L)'.   

          Solve M(L) * x = e. */

	WORK(1) = 1.;
	i__2 = *n;
	for (i = 2; i <= *n; ++i) {
	    WORK(i) = WORK(i - 1) * (d__1 = EF(i - 1), ABS(d__1)) + 1.;
/* L60: */
	}

/*        Solve D * M(L)' * x = b. */

	WORK(*n) /= DF(*n);
	for (i = *n - 1; i >= 1; --i) {
	    WORK(i) = WORK(i) / DF(i) + WORK(i + 1) * (d__1 = EF(i), ABS(d__1)
		    );
/* L70: */
	}

/*        Compute norm(inv(A)) = MAX(x(i)), 1<=i<=n. */


#ifdef PETSC_PREFIX_SUFFIX
	ix = idamax_(n, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	ix = iqamax(n, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	ix = iqamax_(n, &WORK(1), &c__1);
#endif

	FERR(j) *= (d__1 = WORK(ix), ABS(d__1));

/*        Normalize error. */

	lstres = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__2 = lstres, d__3 = (d__1 = X(i,j), ABS(d__1));
	    lstres = MAX(d__2,d__3);
/* L80: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L90: */
    }

    return;

/*     End of DPTRFS */

} /* dptrfs_ */

