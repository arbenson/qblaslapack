#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsptri_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsptri(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsptri_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif

	ipiv, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSPTRI computes the inverse of a real symmetric indefinite matrix   
    A in packed storage using the factorization A = U*D*U**T or   
    A = L*D*L**T computed by DSPTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the details of the factorization are stored 
  
            as an upper or lower triangular matrix.   
            = 'U':  Upper triangular, form is A = U*D*U**T;   
            = 'L':  Lower triangular, form is A = L*D*L**T.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the block diagonal matrix D and the multipliers   
            used to obtain the factor U or L as computed by DSPTRF,   
            stored as a packed triangular matrix.   

            On exit, if INFO = 0, the (symmetric) inverse of the original 
  
            matrix, stored as a packed triangular matrix. The j-th column 
  
            of inv(A) is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = inv(A)(i,j) for 1<=i<=j;   
            if UPLO = 'L',   
               AP(i + (j-1)*(2n-j)/2) = inv(A)(i,j) for j<=i<=n.   

    IPIV    (input) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D   
            as determined by DSPTRF.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (N)   

    INFO    (output) INTEGER   
            = 0: successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   
            > 0: if INFO = i, D(i,i) = 0; the matrix is singular and its 
  
                 inverse could not be computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b11 = -1.;
    static LONG DOUBLE c_b13 = 0.;
    
    /* System generated locals */
    int i__1;
    LONG DOUBLE d__1;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern LONG DOUBLE ddot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern LONG DOUBLE qdot(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern LONG DOUBLE qdot_(int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
#endif

	    int *);
    static LONG DOUBLE temp, akkp1, d;
    static int j, k;
    static LONG DOUBLE t;
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
	    LONG DOUBLE *, int *), dswap_(int *, LONG DOUBLE *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qswap(int *, LONG DOUBLE *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qswap_(int *, LONG DOUBLE *, int 
#endif

	    *, LONG DOUBLE *, int *);
    static int kstep;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dspmv_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qspmv(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qspmv_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *,
	     int *);
    static long int upper;
    static LONG DOUBLE ak;
    static int kc, kp, kx;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static int kcnext, kpc, npp;
    static LONG DOUBLE akp1;



#define WORK(I) work[(I)-1]
#define IPIV(I) ipiv[(I)-1]
#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPTRI", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Check that the diagonal matrix D is nonsingular. */

    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

	kp = *n * (*n + 1) / 2;
	for (*info = *n; *info >= 1; --(*info)) {
	    if (IPIV(*info) > 0 && AP(kp) == 0.) {
		return;
	    }
	    kp -= *info;
/* L10: */
	}
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

	kp = 1;
	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (IPIV(*info) > 0 && AP(kp) == 0.) {
		return;
	    }
	    kp = kp + *n - *info + 1;
/* L20: */
	}
    }
    *info = 0;

    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
	kc = 1;
L30:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L50;
	}

	kcnext = kc + k;
	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    AP(kc + k - 1) = 1. / AP(kc + k - 1);

/*           Compute column K of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
#endif

		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif

			AP(kc), &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kc + k - 1) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kc + k - 1) -= qdot(&i__1, &WORK(1), &c__1, &AP(kc), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kc + k - 1) -= qdot_(&i__1, &WORK(1), &c__1, &AP(kc), &
#endif

			c__1);
	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    t = (d__1 = AP(kcnext + k - 1), ABS(d__1));
	    ak = AP(kc + k - 1) / t;
	    akp1 = AP(kcnext + k) / t;
	    akkp1 = AP(kcnext + k - 1) / t;
	    d = t * (ak * akp1 - 1.);
	    AP(kc + k - 1) = akp1 / d;
	    AP(kcnext + k) = ak / d;
	    AP(kcnext + k - 1) = -akkp1 / d;

/*           Compute columns K and K+1 of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &AP(kc), &c__1, &WORK(1), &c__1);
#endif

		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif

			AP(kc), &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kc + k - 1) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kc + k - 1) -= qdot(&i__1, &WORK(1), &c__1, &AP(kc), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kc + k - 1) -= qdot_(&i__1, &WORK(1), &c__1, &AP(kc), &
#endif

			c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kcnext + k - 1) -= ddot_(&i__1, &AP(kc), &c__1, &AP(kcnext)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kcnext + k - 1) -= qdot(&i__1, &AP(kc), &c__1, &AP(kcnext)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kcnext + k - 1) -= qdot_(&i__1, &AP(kc), &c__1, &AP(kcnext)
#endif

			, &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &AP(kcnext), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &AP(kcnext), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &AP(kcnext), &c__1, &WORK(1), &c__1);
#endif

		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__1, &c_b11, &AP(1), &WORK(1), &c__1, &c_b13, &
#endif

			AP(kcnext), &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kcnext + k) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kcnext), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kcnext + k) -= qdot(&i__1, &WORK(1), &c__1, &AP(kcnext), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kcnext + k) -= qdot_(&i__1, &WORK(1), &c__1, &AP(kcnext), &
#endif

			c__1);
	    }
	    kstep = 2;
	    kcnext = kcnext + k + 1;
	}

	kp = (i__1 = IPIV(k), ABS(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading 
  
             submatrix A(1:k+1,1:k+1) */

	    kpc = (kp - 1) * kp / 2 + 1;
	    i__1 = kp - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dswap_(&i__1, &AP(kc), &c__1, &AP(kpc), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qswap(&i__1, &AP(kc), &c__1, &AP(kpc), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qswap_(&i__1, &AP(kc), &c__1, &AP(kpc), &c__1);
#endif

	    kx = kpc + kp - 1;
	    i__1 = k - 1;
	    for (j = kp + 1; j <= k-1; ++j) {
		kx = kx + j - 1;
		temp = AP(kc + j - 1);
		AP(kc + j - 1) = AP(kx);
		AP(kx) = temp;
/* L40: */
	    }
	    temp = AP(kc + k - 1);
	    AP(kc + k - 1) = AP(kpc + kp - 1);
	    AP(kpc + kp - 1) = temp;
	    if (kstep == 2) {
		temp = AP(kc + k + k - 1);
		AP(kc + k + k - 1) = AP(kc + k + kp - 1);
		AP(kc + k + kp - 1) = temp;
	    }
	}

	k += kstep;
	kc = kcnext;
	goto L30;
L50:

	;
    } else {

/*        Compute inv(A) from the factorization A = L*D*L'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	npp = *n * (*n + 1) / 2;
	k = *n;
	kc = npp;
L60:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L80;
	}

	kcnext = kc - (*n - k + 2);
	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    AP(kc) = 1. / AP(kc);

/*           Compute column K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
#endif

		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__1, &c_b11, &AP(kc + *n - k + 1), &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__1, &c_b11, &AP(kc + *n - k + 1), &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__1, &c_b11, &AP(kc + *n - k + 1), &WORK(1), &
#endif

			c__1, &c_b13, &AP(kc + 1), &c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kc) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kc) -= qdot(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kc) -= qdot_(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
#endif

	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    t = (d__1 = AP(kcnext + 1), ABS(d__1));
	    ak = AP(kcnext) / t;
	    akp1 = AP(kc) / t;
	    akkp1 = AP(kcnext + 1) / t;
	    d = t * (ak * akp1 - 1.);
	    AP(kcnext) = akp1 / d;
	    AP(kc) = ak / d;
	    AP(kcnext + 1) = -akkp1 / d;

/*           Compute columns K-1 and K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &AP(kc + 1), &c__1, &WORK(1), &c__1);
#endif

		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
#endif

			&c__1, &c_b13, &AP(kc + 1), &c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kc) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kc) -= qdot(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kc) -= qdot_(&i__1, &WORK(1), &c__1, &AP(kc + 1), &c__1);
#endif

		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kcnext + 1) -= ddot_(&i__1, &AP(kc + 1), &c__1, &AP(kcnext 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kcnext + 1) -= qdot(&i__1, &AP(kc + 1), &c__1, &AP(kcnext 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kcnext + 1) -= qdot_(&i__1, &AP(kc + 1), &c__1, &AP(kcnext 
#endif

			+ 2), &c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &AP(kcnext + 2), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &AP(kcnext + 2), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &AP(kcnext + 2), &c__1, &WORK(1), &c__1);
#endif

		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__1, &c_b11, &AP(kc + (*n - k + 1)), &WORK(1), 
#endif

			&c__1, &c_b13, &AP(kcnext + 2), &c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		AP(kcnext) -= ddot_(&i__1, &WORK(1), &c__1, &AP(kcnext + 2), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(kcnext) -= qdot(&i__1, &WORK(1), &c__1, &AP(kcnext + 2), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(kcnext) -= qdot_(&i__1, &WORK(1), &c__1, &AP(kcnext + 2), &
#endif

			c__1);
	    }
	    kstep = 2;
	    kcnext -= *n - k + 3;
	}

	kp = (i__1 = IPIV(k), ABS(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing
   
             submatrix A(k-1:n,k-1:n) */

	    kpc = npp - (*n - kp + 1) * (*n - kp + 2) / 2 + 1;
	    if (kp < *n) {
		i__1 = *n - kp;

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(&i__1, &AP(kc + kp - k + 1), &c__1, &AP(kpc + 1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(&i__1, &AP(kc + kp - k + 1), &c__1, &AP(kpc + 1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(&i__1, &AP(kc + kp - k + 1), &c__1, &AP(kpc + 1), &
#endif

			c__1);
	    }
	    kx = kc + kp - k;
	    i__1 = kp - 1;
	    for (j = k + 1; j <= kp-1; ++j) {
		kx = kx + *n - j + 1;
		temp = AP(kc + j - k);
		AP(kc + j - k) = AP(kx);
		AP(kx) = temp;
/* L70: */
	    }
	    temp = AP(kc);
	    AP(kc) = AP(kpc);
	    AP(kpc) = temp;
	    if (kstep == 2) {
		temp = AP(kc - *n + k - 1);
		AP(kc - *n + k - 1) = AP(kc - *n + kp - 1);
		AP(kc - *n + kp - 1) = temp;
	    }
	}

	k -= kstep;
	kc = kcnext;
	goto L60;
L80:
	;
    }

    return;

/*     End of DSPTRI */

} /* dsptri_ */

