#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgtrfs_(char *trans, int *n, int *nrhs, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgtrfs(char *trans, int *n, int *nrhs, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgtrfs_(char *trans, int *n, int *nrhs, 
#endif

	LONG DOUBLE *dl, LONG DOUBLE *d, LONG DOUBLE *du, LONG DOUBLE *dlf, 
	LONG DOUBLE *df, LONG DOUBLE *duf, LONG DOUBLE *du2, int *ipiv, 
	LONG DOUBLE *b, int *ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *
	ferr, LONG DOUBLE *berr, LONG DOUBLE *work, int *iwork, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGTRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is tridiagonal, and provides   
    error bounds and backward error estimates for the solution.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B     (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    DL      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) subdiagonal elements of A.   

    D       (input) LONG DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of A.   

    DU      (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) superdiagonal elements of A.   

    DLF     (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) multipliers that define the matrix L from the   
            LU factorization of A as computed by DGTTRF.   

    DF      (input) LONG DOUBLE PRECISION array, dimension (N)   
            The n diagonal elements of the upper triangular matrix U from 
  
            the LU factorization of A.   

    DUF     (input) LONG DOUBLE PRECISION array, dimension (N-1)   
            The (n-1) elements of the first superdiagonal of U.   

    DU2     (input) LONG DOUBLE PRECISION array, dimension (N-2)   
            The (n-2) elements of the second superdiagonal of U.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= n, row i of the matrix was   
            interchanged with row IPIV(i).  IPIV(i) will always be either 
  
            i or i+1; IPIV(i) = i indicates a row interchange was not   
            required.   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (input/output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by DGTTRS.   
            On exit, the improved solution matrix X.   

    LDX     (input) INTEGER   
            The leading dimension of the array X.  LDX >= MAX(1,N).   

    FERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The estimated forward error bound for each solution vector   
            X(j) (the j-th column of the solution matrix X).   
            If XTRUE is the true solution corresponding to X(j), FERR(j) 
  
            is an estimated upper bound for the magnitude of the largest 
  
            element in (X(j) - XTRUE) divided by the magnitude of the   
            largest element in X(j).  The estimate is as reliable as   
            the estimate for RCOND, and is almost always a slight   
            overestimate of the true error.   

    BERR    (output) LONG DOUBLE PRECISION array, dimension (NRHS)   
            The componentwise relative backward error of each solution   
            vector X(j) (i.e., the smallest relative change in   
            any element of A or B that makes X(j) an exact solution).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (3*N)   

    IWORK   (workspace) INTEGER array, dimension (N)   

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
    static LONG DOUBLE c_b18 = -1.;
    static LONG DOUBLE c_b19 = 1.;
    
    /* System generated locals */
    int  i__1, i__2;
    LONG DOUBLE d__1, d__2, d__3, d__4;
    /* Local variables */
    static int kase;
    static LONG DOUBLE safe1, safe2;
    static int i, j;
    static LONG DOUBLE s;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dcopy_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qcopy_(int *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), daxpy_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qaxpy(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qaxpy_(int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static int count;

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE dlamch_(char *);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch(char *);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qlamch_(char *);
#endif


#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlacon_(int *, LONG DOUBLE *, LONG DOUBLE *,
#endif

	     int *, LONG DOUBLE *, int *);
    static int nz;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlagtm_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlagtm(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlagtm_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *);
    static LONG DOUBLE safmin;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int notran;
    static char transn[1];

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgttrs_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgttrs(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgttrs_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *,
	     LONG DOUBLE *, int *, int *);
    static char transt[1];
    static LONG DOUBLE lstres, eps;



#define DL(I) dl[(I)-1]
#define D(I) d[(I)-1]
#define DU(I) du[(I)-1]
#define DLF(I) dlf[(I)-1]
#define DF(I) df[(I)-1]
#define DUF(I) duf[(I)-1]
#define DU2(I) du2[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*ldb < MAX(1,*n)) {
	*info = -13;
    } else if (*ldx < MAX(1,*n)) {
	*info = -15;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGTRFS", &i__1);
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

    if (notran) {
	*(unsigned char *)transn = 'N';
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transn = 'T';
	*(unsigned char *)transt = 'N';
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

          Compute residual R = B - op(A) * X,   
          where op(A) = A, A**T, or A**H, depending on TRANS. */


#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dlagtm_(trans, n, &c__1, &c_b18, &DL(1), &D(1), &DU(1), &X(1,j), ldx, &c_b19, &WORK(*n + 1), n);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlagtm(trans, n, &c__1, &c_b18, &DL(1), &D(1), &DU(1), &X(1,j), ldx, &c_b19, &WORK(*n + 1), n);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlagtm_(trans, n, &c__1, &c_b18, &DL(1), &D(1), &DU(1), &X(1,j), ldx, &c_b19, &WORK(*n + 1), n);
#endif


/*        Compute ABS(op(A))*ABS(x) + ABS(b) for use in the backward 
  
          error bound. */

	if (notran) {
	    if (*n == 1) {
		WORK(1) = (d__1 = B(1,j), ABS(d__1)) + (d__2 = D(1)
			 * X(1,j), ABS(d__2));
	    } else {
		WORK(1) = (d__1 = B(1,j), ABS(d__1)) + (d__2 = D(1)
			 * X(1,j), ABS(d__2)) + (d__3 = DU(1) * X(2,j), ABS(d__3));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    WORK(i) = (d__1 = B(i,j), ABS(d__1)) + (d__2 = 
			    DL(i - 1) * X(i-1,j), ABS(d__2)) + (
			    d__3 = D(i) * X(i,j), ABS(d__3)) + (
			    d__4 = DU(i) * X(i+1,j), ABS(d__4));
/* L30: */
		}
		WORK(*n) = (d__1 = B(*n,j), ABS(d__1)) + (d__2 = 
			DL(*n - 1) * X(*n-1,j), ABS(d__2)) + (
			d__3 = D(*n) * X(*n,j), ABS(d__3));
	    }
	} else {
	    if (*n == 1) {
		WORK(1) = (d__1 = B(1,j), ABS(d__1)) + (d__2 = D(1)
			 * X(1,j), ABS(d__2));
	    } else {
		WORK(1) = (d__1 = B(1,j), ABS(d__1)) + (d__2 = D(1)
			 * X(1,j), ABS(d__2)) + (d__3 = DL(1) * X(2,j), ABS(d__3));
		i__2 = *n - 1;
		for (i = 2; i <= *n-1; ++i) {
		    WORK(i) = (d__1 = B(i,j), ABS(d__1)) + (d__2 = 
			    DU(i - 1) * X(i-1,j), ABS(d__2)) + (
			    d__3 = D(i) * X(i,j), ABS(d__3)) + (
			    d__4 = DL(i) * X(i+1,j), ABS(d__4));
/* L40: */
		}
		WORK(*n) = (d__1 = B(*n,j), ABS(d__1)) + (d__2 = 
			DU(*n - 1) * X(*n-1,j), ABS(d__2)) + (
			d__3 = D(*n) * X(*n,j), ABS(d__3));
	    }
	}

/*        Compute componentwise relative backward error from formula 
  

          MAX(i) ( ABS(R(i)) / ( ABS(op(A))*ABS(X) + ABS(B) )(i) )   

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
/* L50: */
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
	    dgttrs_(trans, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgttrs(trans, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgttrs_(trans, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &IPIV(
#endif

		    1), &WORK(*n + 1), n, info);

#ifdef PETSC_PREFIX_SUFFIX
	    daxpy_(n, &c_b19, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qaxpy(n, &c_b19, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qaxpy_(n, &c_b19, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif

		    ;
	    lstres = BERR(j);
	    ++count;
	    goto L20;
	}

/*        Bound error from formula   

          norm(X - XTRUE) / norm(X) .le. FERR =   
          norm( ABS(inv(op(A)))*   
             ( ABS(R) + NZ*EPS*( ABS(op(A))*ABS(X)+ABS(B) ))) / norm(X
)   

          where   
            norm(Z) is the magnitude of the largest component of Z   
            inv(op(A)) is the inverse of op(A)   
            ABS(Z) is the componentwise absolute value of the matrix o
r   
               vector Z   
            NZ is the maximum number of nonzeros in any row of A, plus
 1   
            EPS is machine epsilon   

          The i-th component of ABS(R)+NZ*EPS*(ABS(op(A))*ABS(X)+ABS(B
))   
          is incremented by SAFE1 if the i-th component of   
          ABS(op(A))*ABS(X) + ABS(B) is less than SAFE2.   

          Use DLACON to estimate the infinity-norm of the matrix   
             inv(op(A)) * diag(W),   
          where W = ABS(R) + NZ*EPS*( ABS(op(A))*ABS(X)+ABS(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
		WORK(i) = (d__1 = WORK(*n + i), ABS(d__1)) + nz * eps * WORK(
			i);
	    } else {
		WORK(i) = (d__1 = WORK(*n + i), ABS(d__1)) + nz * eps * WORK(
			i) + safe1;
	    }
/* L60: */
	}

	kase = 0;
L70:

#ifdef PETSC_PREFIX_SUFFIX
	dlacon_(n, &WORK((*n << 1) + 1), &WORK(*n + 1), &IWORK(1), &FERR(j), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlacon(n, &WORK((*n << 1) + 1), &WORK(*n + 1), &IWORK(1), &FERR(j), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlacon_(n, &WORK((*n << 1) + 1), &WORK(*n + 1), &IWORK(1), &FERR(j), &
#endif

		kase);
	if (kase != 0) {
	    if (kase == 1) {

/*              Multiply by diag(W)*inv(op(A)**T). */


#ifdef PETSC_PREFIX_SUFFIX
		dgttrs_(transt, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgttrs(transt, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgttrs_(transt, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
#endif

			IPIV(1), &WORK(*n + 1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) = WORK(i) * WORK(*n + i);
/* L80: */
		}
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) = WORK(i) * WORK(*n + i);
/* L90: */
		}

#ifdef PETSC_PREFIX_SUFFIX
		dgttrs_(transn, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgttrs(transn, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgttrs_(transn, n, &c__1, &DLF(1), &DF(1), &DUF(1), &DU2(1), &
#endif

			IPIV(1), &WORK(*n + 1), n, info);
	    }
	    goto L70;
	}

/*        Normalize error. */

	lstres = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__2 = lstres, d__3 = (d__1 = X(i,j), ABS(d__1));
	    lstres = MAX(d__2,d__3);
/* L100: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L110: */
    }

    return;

/*     End of DGTRFS */

} /* dgtrfs_ */

