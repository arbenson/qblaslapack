#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpbrfs_(char *uplo, int *n, int *kd, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpbrfs(char *uplo, int *n, int *kd, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpbrfs_(char *uplo, int *n, int *kd, int *
#endif

	nrhs, LONG DOUBLE *ab, int *ldab, LONG DOUBLE *afb, int *ldafb, 
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

    DPBRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is symmetric positive definite 
  
    and banded, and provides error bounds and backward error estimates   
    for the solution.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The upper or lower triangle of the symmetric band matrix A,   
            stored in the first KD+1 rows of the array.  The j-th column 
  
            of A is stored in the j-th column of the array AB as follows: 
  
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    AFB     (input) LONG DOUBLE PRECISION array, dimension (LDAFB,N)   
            The triangular factor U or L from the Cholesky factorization 
  
            A = U**T*U or A = L*L**T of the band matrix A as computed by 
  
            DPBTRF, in the same storage format as A (see AB).   

    LDAFB   (input) INTEGER   
            The leading dimension of the array AFB.  LDAFB >= KD+1.   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (input/output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by DPBTRS.   
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
    static LONG DOUBLE c_b12 = -1.;
    static LONG DOUBLE c_b14 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4, i__5;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static int kase;
    static LONG DOUBLE safe1, safe2;
    static int i, j, k, l;
    static LONG DOUBLE s;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsbmv_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsbmv(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsbmv_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), dcopy_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qcopy_(int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), daxpy_(int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy(int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qaxpy_(int 
#endif

	    *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *)
	    ;
    static int count;
    static long int upper;

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
    static LONG DOUBLE xk;
    static int nz;
    static LONG DOUBLE safmin;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dpbtrs_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpbtrs(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qpbtrs_(
#endif

	    char *, int *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *, int *);
    static LONG DOUBLE lstres, eps;



#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define AFB(I,J) afb[(I)-1 + ((J)-1)* ( *ldafb)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kd < 0) {
	*info = -3;
    } else if (*nrhs < 0) {
	*info = -4;
    } else if (*ldab < *kd + 1) {
	*info = -6;
    } else if (*ldafb < *kd + 1) {
	*info = -8;
    } else if (*ldb < MAX(1,*n)) {
	*info = -10;
    } else if (*ldx < MAX(1,*n)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPBRFS", &i__1);
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

/*     NZ = maximum number of nonzero elements in each row of A, plus 1   

   Computing MIN */
    i__1 = *n + 1, i__2 = (*kd << 1) + 2;
    nz = MIN(i__1,i__2);

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

          Compute residual R = B - A * X */


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
	dsbmv_(uplo, n, kd, &c_b12, &AB(1,1), ldab, &X(1,j), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsbmv(uplo, n, kd, &c_b12, &AB(1,1), ldab, &X(1,j), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsbmv_(uplo, n, kd, &c_b12, &AB(1,1), ldab, &X(1,j), 
#endif

		&c__1, &c_b14, &WORK(*n + 1), &c__1);

/*        Compute componentwise relative backward error from formula 
  

          MAX(i) ( ABS(R(i)) / ( ABS(A)*ABS(X) + ABS(B) )(i) )   

          where ABS(Z) is the componentwise absolute value of the matr
ix   
          or vector Z.  If the i-th component of the denominator is le
ss   
          than SAFE2, then SAFE1 is added to the i-th components of th
e   
          numerator and denominator before dividing. */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    WORK(i) = (d__1 = B(i,j), ABS(d__1));
/* L30: */
	}

/*        Compute ABS(A)*ABS(X) + ABS(B). */

	if (upper) {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		s = 0.;
		xk = (d__1 = X(k,j), ABS(d__1));
		l = *kd + 1 - k;
/* Computing MAX */
		i__3 = 1, i__4 = k - *kd;
		i__5 = k - 1;
		for (i = MAX(1,k-*kd); i <= k-1; ++i) {
		    WORK(i) += (d__1 = AB(l+i,k), ABS(d__1)) * 
			    xk;
		    s += (d__1 = AB(l+i,k), ABS(d__1)) * (d__2 =
			     X(i,j), ABS(d__2));
/* L40: */
		}
		WORK(k) = WORK(k) + (d__1 = AB(*kd+1,k), ABS(
			d__1)) * xk + s;
/* L50: */
	    }
	} else {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		s = 0.;
		xk = (d__1 = X(k,j), ABS(d__1));
		WORK(k) += (d__1 = AB(1,k), ABS(d__1)) * xk;
		l = 1 - k;
/* Computing MIN */
		i__3 = *n, i__4 = k + *kd;
		i__5 = MIN(i__3,i__4);
		for (i = k + 1; i <= MIN(*n,k+*kd); ++i) {
		    WORK(i) += (d__1 = AB(l+i,k), ABS(d__1)) * 
			    xk;
		    s += (d__1 = AB(l+i,k), ABS(d__1)) * (d__2 =
			     X(i,j), ABS(d__2));
/* L60: */
		}
		WORK(k) += s;
/* L70: */
	    }
	}
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
/* L80: */
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
	    dpbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n + 1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qpbtrs(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n + 1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qpbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n + 1)
#endif

		    , n, info);

#ifdef PETSC_PREFIX_SUFFIX
	    daxpy_(n, &c_b14, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qaxpy(n, &c_b14, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qaxpy_(n, &c_b14, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
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
          ABS(A)*ABS(X) + ABS(B) is less than SAFE2.   

          Use DLACON to estimate the infinity-norm of the matrix   
             inv(A) * diag(W),   
          where W = ABS(R) + NZ*EPS*( ABS(A)*ABS(X)+ABS(B) ))) */

	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (WORK(i) > safe2) {
		WORK(i) = (d__1 = WORK(*n + i), ABS(d__1)) + nz * eps * WORK(
			i);
	    } else {
		WORK(i) = (d__1 = WORK(*n + i), ABS(d__1)) + nz * eps * WORK(
			i) + safe1;
	    }
/* L90: */
	}

	kase = 0;
L100:

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

/*              Multiply by diag(W)*inv(A'). */


#ifdef PETSC_PREFIX_SUFFIX
		dpbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qpbtrs(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qpbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
#endif

			+ 1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) *= WORK(i);
/* L110: */
		}
	    } else if (kase == 2) {

/*              Multiply by inv(A)*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) *= WORK(i);
/* L120: */
		}

#ifdef PETSC_PREFIX_SUFFIX
		dpbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qpbtrs(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qpbtrs_(uplo, n, kd, &c__1, &AFB(1,1), ldafb, &WORK(*n 
#endif

			+ 1), n, info);
	    }
	    goto L100;
	}

/*        Normalize error. */

	lstres = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__2 = lstres, d__3 = (d__1 = X(i,j), ABS(d__1));
	    lstres = MAX(d__2,d__3);
/* L130: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L140: */
    }

    return;

/*     End of DPBRFS */

} /* dpbrfs_ */

