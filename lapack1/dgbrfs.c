#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbrfs_(char *trans, int *n, int *kl, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbrfs(char *trans, int *n, int *kl, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbrfs_(char *trans, int *n, int *kl, int *
#endif

	ku, int *nrhs, LONG DOUBLE *ab, int *ldab, LONG DOUBLE *afb, 
	int *ldafb, int *ipiv, LONG DOUBLE *b, int *ldb, 
	LONG DOUBLE *x, int *ldx, LONG DOUBLE *ferr, LONG DOUBLE *berr, 
	LONG DOUBLE *work, int *iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DGBRFS improves the computed solution to a system of linear   
    equations when the coefficient matrix is banded, and provides   
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

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The original band matrix A, stored in rows 1 to KL+KU+1.   
            The j-th column of A is stored in the j-th column of the   
            array AB as follows:   
            AB(ku+1+i-j,j) = A(i,j) for MAX(1,j-ku)<=i<=MIN(n,j+kl).   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KL+KU+1.   

    AFB     (input) LONG DOUBLE PRECISION array, dimension (LDAFB,N)   
            Details of the LU factorization of the band matrix A, as   
            computed by DGBTRF.  U is stored as an upper triangular band 
  
            matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and   
            the multipliers used during the factorization are stored in   
            rows KL+KU+2 to 2*KL+KU+1.   

    LDAFB   (input) INTEGER   
            The leading dimension of the array AFB.  LDAFB >= 2*KL*KU+1. 
  

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices from DGBTRF; for 1<=i<=N, row i of the   
            matrix was interchanged with row IPIV(i).   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (input/output) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            On entry, the solution matrix X, as computed by DGBTRS.   
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
    static LONG DOUBLE c_b15 = -1.;
    static LONG DOUBLE c_b17 = 1.;
    
    /* System generated locals */
    int   i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static int kase;
    static LONG DOUBLE safe1, safe2;
    static int i, j, k;
    static LONG DOUBLE s;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgbmv_(char *, int *, int *, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbmv(char *, int *, int *, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgbmv_(char *, int *, int *, int *
#endif

	    , int *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, LONG DOUBLE *, int *);
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
    static int count, kk;

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
    extern /* Subroutine */ void xerbla_(char *, int *), dgbtrs_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qgbtrs(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qgbtrs_(
#endif

	    char *, int *, int *, int *, int *, LONG DOUBLE *, 
	    int *, int *, LONG DOUBLE *, int *, int *);
    static long int notran;
    static char transt[1];
    static LONG DOUBLE lstres, eps;



#define IPIV(I) ipiv[(I)-1]
#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define AFB(I,J) afb[(I)-1 + ((J)-1)* ( *ldafb)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*kl < 0) {
	*info = -3;
    } else if (*ku < 0) {
	*info = -4;
    } else if (*nrhs < 0) {
	*info = -5;
    } else if (*ldab < *kl + *ku + 1) {
	*info = -7;
    } else if (*ldafb < (*kl << 1) + *ku + 1) {
	*info = -9;
    } else if (*ldb < MAX(1,*n)) {
	*info = -12;
    } else if (*ldx < MAX(1,*n)) {
	*info = -14;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBRFS", &i__1);
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
	*(unsigned char *)transt = 'T';
    } else {
	*(unsigned char *)transt = 'N';
    }

/*     NZ = maximum number of nonzero elements in each row of A, plus 1   

   Computing MIN */
    i__1 = *kl + *ku + 2, i__2 = *n + 1;
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
	dgbmv_(trans, n, n, kl, ku, &c_b15, &AB(1,1), ldab, &X(1,j), &c__1, &c_b17, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qgbmv(trans, n, n, kl, ku, &c_b15, &AB(1,1), ldab, &X(1,j), &c__1, &c_b17, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qgbmv_(trans, n, n, kl, ku, &c_b15, &AB(1,1), ldab, &X(1,j), &c__1, &c_b17, &WORK(*n + 1), &c__1);
#endif


/*        Compute componentwise relative backward error from formula 
  

          MAX(i) ( ABS(R(i)) / ( ABS(op(A))*ABS(X) + ABS(B) )(i) )   

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

/*        Compute ABS(op(A))*ABS(X) + ABS(B). */

	if (notran) {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		kk = *ku + 1 - k;
		xk = (d__1 = X(k,j), ABS(d__1));
/* Computing MAX */
		i__3 = 1, i__4 = k - *ku;
/* Computing MIN */
		i__6 = *n, i__7 = k + *kl;
		i__5 = MIN(i__6,i__7);
		for (i = MAX(1,k-*ku); i <= MIN(*n,k+*kl); ++i) {
		    WORK(i) += (d__1 = AB(kk+i,k), ABS(d__1)) * 
			    xk;
/* L40: */
		}
/* L50: */
	    }
	} else {
	    i__2 = *n;
	    for (k = 1; k <= *n; ++k) {
		s = 0.;
		kk = *ku + 1 - k;
/* Computing MAX */
		i__5 = 1, i__3 = k - *ku;
/* Computing MIN */
		i__6 = *n, i__7 = k + *kl;
		i__4 = MIN(i__6,i__7);
		for (i = MAX(1,k-*ku); i <= MIN(*n,k+*kl); ++i) {
		    s += (d__1 = AB(kk+i,k), ABS(d__1)) * (d__2 
			    = X(i,j), ABS(d__2));
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
	    dgbtrs_(trans, n, kl, ku, &c__1, &AFB(1,1), ldafb, &IPIV(1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgbtrs(trans, n, kl, ku, &c__1, &AFB(1,1), ldafb, &IPIV(1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgbtrs_(trans, n, kl, ku, &c__1, &AFB(1,1), ldafb, &IPIV(1)
#endif

		    , &WORK(*n + 1), n, info);

#ifdef PETSC_PREFIX_SUFFIX
	    daxpy_(n, &c_b17, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qaxpy(n, &c_b17, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qaxpy_(n, &c_b17, &WORK(*n + 1), &c__1, &X(1,j), &c__1)
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

/*              Multiply by diag(W)*inv(op(A)**T). */


#ifdef PETSC_PREFIX_SUFFIX
		dgbtrs_(transt, n, kl, ku, &c__1, &AFB(1,1), ldafb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgbtrs(transt, n, kl, ku, &c__1, &AFB(1,1), ldafb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgbtrs_(transt, n, kl, ku, &c__1, &AFB(1,1), ldafb, &
#endif

			IPIV(1), &WORK(*n + 1), n, info);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) *= WORK(i);
/* L110: */
		}
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) *= WORK(i);
/* L120: */
		}

#ifdef PETSC_PREFIX_SUFFIX
		dgbtrs_(trans, n, kl, ku, &c__1, &AFB(1,1), ldafb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgbtrs(trans, n, kl, ku, &c__1, &AFB(1,1), ldafb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgbtrs_(trans, n, kl, ku, &c__1, &AFB(1,1), ldafb, &
#endif

			IPIV(1), &WORK(*n + 1), n, info);
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

/*     End of DGBRFS */

} /* dgbrfs_ */

