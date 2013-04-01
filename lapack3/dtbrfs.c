#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtbrfs_(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtbrfs(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtbrfs_(char *uplo, char *trans, char *diag, int *n, 
#endif

	int *kd, int *nrhs, LONG DOUBLE *ab, int *ldab, LONG DOUBLE 
	*b, int *ldb, LONG DOUBLE *x, int *ldx, LONG DOUBLE *ferr, 
	LONG DOUBLE *berr, LONG DOUBLE *work, int *iwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DTBRFS provides error bounds and backward error estimates for the   
    solution to a system of linear equations with a triangular band   
    coefficient matrix.   

    The solution matrix X must be computed by DTBTRS or some other   
    means before entering this routine.  DTBRFS does not do iterative   
    refinement because doing so cannot improve the backward error.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A**T * X = B  (Transpose)   
            = 'C':  A**H * X = B  (Conjugate transpose = Transpose)   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals or subdiagonals of the   
            triangular band matrix A.  KD >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrices B and X.  NRHS >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The upper or lower triangular band matrix A, stored in the   
            first kd+1 rows of the array. The j-th column of A is stored 
  
            in the j-th column of the array AB as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  
            If DIAG = 'U', the diagonal elements of A are not referenced 
  
            and are assumed to be 1.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            The right hand side matrix B.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    X       (input) LONG DOUBLE PRECISION array, dimension (LDX,NRHS)   
            The solution matrix X.   

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

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b19 = -1.;
    
    /* System generated locals */
    int  i__1, 
	    i__2, i__3, i__4, i__5;
    LONG DOUBLE d__1, d__2, d__3;
    /* Local variables */
    static int kase;
    static LONG DOUBLE safe1, safe2;
    static int i, j, k;
    static LONG DOUBLE s;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtbmv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtbmv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtbmv_(char *, char *, char *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), dcopy_(int *, LONG DOUBLE *, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qcopy(int *, LONG DOUBLE *, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), qcopy_(int *, LONG DOUBLE *, int *
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    , LONG DOUBLE *, int *), dtbsv_(char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    , LONG DOUBLE *, int *), qtbsv(char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    , LONG DOUBLE *, int *), qtbsv_(char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    int *), daxpy_(int *, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *), qaxpy(int *, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *), qaxpy_(int *, LONG DOUBLE *
#endif

	    , LONG DOUBLE *, int *, LONG DOUBLE *, int *);
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
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int notran;
    static char transt[1];
    static long int nounit;
    static LONG DOUBLE lstres, eps;



#define FERR(I) ferr[(I)-1]
#define BERR(I) berr[(I)-1]
#define WORK(I) work[(I)-1]
#define IWORK(I) iwork[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]
#define X(I,J) x[(I)-1 + ((J)-1)* ( *ldx)]

    *info = 0;
    upper = lsame_(uplo, "U");
    notran = lsame_(trans, "N");
    nounit = lsame_(diag, "N");

    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, 
	    "C")) {
	*info = -2;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -3;
    } else if (*n < 0) {
	*info = -4;
    } else if (*kd < 0) {
	*info = -5;
    } else if (*nrhs < 0) {
	*info = -6;
    } else if (*ldab < *kd + 1) {
	*info = -8;
    } else if (*ldb < MAX(1,*n)) {
	*info = -10;
    } else if (*ldx < MAX(1,*n)) {
	*info = -12;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTBRFS", &i__1);
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

/*     NZ = maximum number of nonzero elements in each row of A, plus 1 */

    nz = *kd + 2;

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

/*        Compute residual R = B - op(A) * X,   
          where op(A) = A or A', depending on TRANS. */


#ifdef PETSC_PREFIX_SUFFIX
	dcopy_(n, &X(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qcopy(n, &X(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qcopy_(n, &X(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
	dtbmv_(uplo, trans, diag, n, kd, &AB(1,1), ldab, &WORK(*n + 1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtbmv(uplo, trans, diag, n, kd, &AB(1,1), ldab, &WORK(*n + 1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtbmv_(uplo, trans, diag, n, kd, &AB(1,1), ldab, &WORK(*n + 1), 
#endif

		&c__1);

#ifdef PETSC_PREFIX_SUFFIX
	daxpy_(n, &c_b19, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qaxpy(n, &c_b19, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qaxpy_(n, &c_b19, &B(1,j), &c__1, &WORK(*n + 1), &c__1);
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
/* L20: */
	}

	if (notran) {

/*           Compute ABS(A)*ABS(X) + ABS(B). */

	    if (upper) {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			xk = (d__1 = X(k,j), ABS(d__1));
/* Computing MAX */
			i__3 = 1, i__4 = k - *kd;
			i__5 = k;
			for (i = MAX(1,k-*kd); i <= k; ++i) {
			    WORK(i) += (d__1 = AB(*kd+1+i-k,k), ABS(d__1)) * xk;
/* L30: */
			}
/* L40: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			xk = (d__1 = X(k,j), ABS(d__1));
/* Computing MAX */
			i__5 = 1, i__3 = k - *kd;
			i__4 = k - 1;
			for (i = MAX(1,k-*kd); i <= k-1; ++i) {
			    WORK(i) += (d__1 = AB(*kd+1+i-k,k), ABS(d__1)) * xk;
/* L50: */
			}
			WORK(k) += xk;
/* L60: */
		    }
		}
	    } else {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			xk = (d__1 = X(k,j), ABS(d__1));
/* Computing MIN */
			i__5 = *n, i__3 = k + *kd;
			i__4 = MIN(i__5,i__3);
			for (i = k; i <= MIN(*n,k+*kd); ++i) {
			    WORK(i) += (d__1 = AB(i+1-k,k), 
				    ABS(d__1)) * xk;
/* L70: */
			}
/* L80: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			xk = (d__1 = X(k,j), ABS(d__1));
/* Computing MIN */
			i__5 = *n, i__3 = k + *kd;
			i__4 = MIN(i__5,i__3);
			for (i = k + 1; i <= MIN(*n,k+*kd); ++i) {
			    WORK(i) += (d__1 = AB(i+1-k,k), 
				    ABS(d__1)) * xk;
/* L90: */
			}
			WORK(k) += xk;
/* L100: */
		    }
		}
	    }
	} else {

/*           Compute ABS(A')*ABS(X) + ABS(B). */

	    if (upper) {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.;
/* Computing MAX */
			i__4 = 1, i__5 = k - *kd;
			i__3 = k;
			for (i = MAX(1,k-*kd); i <= k; ++i) {
			    s += (d__1 = AB(*kd+1+i-k,k), 
				    ABS(d__1)) * (d__2 = X(i,j), 
				    ABS(d__2));
/* L110: */
			}
			WORK(k) += s;
/* L120: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = (d__1 = X(k,j), ABS(d__1));
/* Computing MAX */
			i__3 = 1, i__4 = k - *kd;
			i__5 = k - 1;
			for (i = MAX(1,k-*kd); i <= k-1; ++i) {
			    s += (d__1 = AB(*kd+1+i-k,k), 
				    ABS(d__1)) * (d__2 = X(i,j), 
				    ABS(d__2));
/* L130: */
			}
			WORK(k) += s;
/* L140: */
		    }
		}
	    } else {
		if (nounit) {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = 0.;
/* Computing MIN */
			i__3 = *n, i__4 = k + *kd;
			i__5 = MIN(i__3,i__4);
			for (i = k; i <= MIN(*n,k+*kd); ++i) {
			    s += (d__1 = AB(i+1-k,k), ABS(
				    d__1)) * (d__2 = X(i,j), ABS(
				    d__2));
/* L150: */
			}
			WORK(k) += s;
/* L160: */
		    }
		} else {
		    i__2 = *n;
		    for (k = 1; k <= *n; ++k) {
			s = (d__1 = X(k,j), ABS(d__1));
/* Computing MIN */
			i__3 = *n, i__4 = k + *kd;
			i__5 = MIN(i__3,i__4);
			for (i = k + 1; i <= MIN(*n,k+*kd); ++i) {
			    s += (d__1 = AB(i+1-k,k), ABS(
				    d__1)) * (d__2 = X(i,j), ABS(
				    d__2));
/* L170: */
			}
			WORK(k) += s;
/* L180: */
		    }
		}
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
/* L190: */
	}
	BERR(j) = s;

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
/* L200: */
	}

	kase = 0;
L210:

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

/*              Multiply by diag(W)*inv(op(A)'). */


#ifdef PETSC_PREFIX_SUFFIX
		dtbsv_(uplo, transt, diag, n, kd, &AB(1,1), ldab, &WORK(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtbsv(uplo, transt, diag, n, kd, &AB(1,1), ldab, &WORK(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtbsv_(uplo, transt, diag, n, kd, &AB(1,1), ldab, &WORK(
#endif

			*n + 1), &c__1);
		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) = WORK(i) * WORK(*n + i);
/* L220: */
		}
	    } else {

/*              Multiply by inv(op(A))*diag(W). */

		i__2 = *n;
		for (i = 1; i <= *n; ++i) {
		    WORK(*n + i) = WORK(i) * WORK(*n + i);
/* L230: */
		}

#ifdef PETSC_PREFIX_SUFFIX
		dtbsv_(uplo, trans, diag, n, kd, &AB(1,1), ldab, &WORK(*
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtbsv(uplo, trans, diag, n, kd, &AB(1,1), ldab, &WORK(*
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtbsv_(uplo, trans, diag, n, kd, &AB(1,1), ldab, &WORK(*
#endif

			n + 1), &c__1);
	    }
	    goto L210;
	}

/*        Normalize error. */

	lstres = 0.;
	i__2 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__2 = lstres, d__3 = (d__1 = X(i,j), ABS(d__1));
	    lstres = MAX(d__2,d__3);
/* L240: */
	}
	if (lstres != 0.) {
	    FERR(j) /= lstres;
	}

/* L250: */
    }

    return;

/*     End of DTBRFS */

} /* dtbrfs_ */

