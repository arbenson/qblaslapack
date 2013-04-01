#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsytri_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsytri(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsytri_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *ipiv, LONG DOUBLE *work, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSYTRI computes the inverse of a real symmetric indefinite matrix   
    A using the factorization A = U*D*U**T or A = L*D*L**T computed by   
    DSYTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the details of the factorization are stored 
  
            as an upper or lower triangular matrix.   
            = 'U':  Upper triangular, form is A = U*D*U**T;   
            = 'L':  Lower triangular, form is A = L*D*L**T.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the block diagonal matrix D and the multipliers   
            used to obtain the factor U or L as computed by DSYTRF.   

            On exit, if INFO = 0, the (symmetric) inverse of the original 
  
            matrix.  If UPLO = 'U', the upper triangular part of the   
            inverse is formed and the part of A below the diagonal is not 
  
            referenced; if UPLO = 'L' the lower triangular part of the   
            inverse is formed and the part of A above the diagonal is   
            not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    IPIV    (input) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D   
            as determined by DSYTRF.   

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
    int  i__1;
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
    static int k;
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
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsymv_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsymv(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsymv_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    LONG DOUBLE *, int *);
    static LONG DOUBLE ak;
    static int kp;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE akp1;



#define IPIV(I) ipiv[(I)-1]
#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYTRI", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Check that the diagonal matrix D is nonsingular. */

    if (upper) {

/*        Upper triangular storage: examine D from bottom to top */

	for (*info = *n; *info >= 1; --(*info)) {
	    if (IPIV(*info) > 0 && A(*info,*info) == 0.) {
		return;
	    }
/* L10: */
	}
    } else {

/*        Lower triangular storage: examine D from top to bottom. */

	i__1 = *n;
	for (*info = 1; *info <= i__1; ++(*info)) {
	    if (IPIV(*info) > 0 && A(*info,*info) == 0.) {
		return;
	    }
/* L20: */
	}
    }
    *info = 0;

    if (upper) {

/*        Compute inv(A) from the factorization A = U*D*U'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = 1;
L30:

/*        If K > N, exit from loop. */

	if (k > *n) {
	    goto L40;
	}

	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    A(k,k) = 1. / A(k,k);

/*           Compute column K of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
#endif

		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif

			c__1, &c_b13, &A(1,k), &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		A(k,k) -= ddot_(&i__1, &WORK(1), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k,k) -= qdot(&i__1, &WORK(1), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k,k) -= qdot_(&i__1, &WORK(1), &c__1, &A(1,k), &c__1);
#endif

	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    t = (d__1 = A(k,k+1), ABS(d__1));
	    ak = A(k,k) / t;
	    akp1 = A(k+1,k+1) / t;
	    akkp1 = A(k,k+1) / t;
	    d = t * (ak * akp1 - 1.);
	    A(k,k) = akp1 / d;
	    A(k+1,k+1) = ak / d;
	    A(k,k+1) = -akkp1 / d;

/*           Compute columns K and K+1 of the inverse. */

	    if (k > 1) {
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &A(1,k), &c__1, &WORK(1), &c__1);
#endif

		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif

			c__1, &c_b13, &A(1,k), &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		A(k,k) -= ddot_(&i__1, &WORK(1), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k,k) -= qdot(&i__1, &WORK(1), &c__1, &A(1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k,k) -= qdot_(&i__1, &WORK(1), &c__1, &A(1,k), &c__1);
#endif

		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		A(k,k+1) -= ddot_(&i__1, &A(1,k), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k,k+1) -= qdot(&i__1, &A(1,k), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k,k+1) -= qdot_(&i__1, &A(1,k), &
#endif

			c__1, &A(1,k+1), &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &A(1,k+1), &c__1, &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &A(1,k+1), &c__1, &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &A(1,k+1), &c__1, &WORK(1), &
#endif

			c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i__1, &c_b11, &A(1,1), lda, &WORK(1), &
#endif

			c__1, &c_b13, &A(1,k+1), &c__1);
		i__1 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		A(k+1,k+1) -= ddot_(&i__1, &WORK(1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k+1,k+1) -= qdot(&i__1, &WORK(1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k+1,k+1) -= qdot_(&i__1, &WORK(1), &c__1, &
#endif

			A(1,k+1), &c__1);
	    }
	    kstep = 2;
	}

	kp = (i__1 = IPIV(k), ABS(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the leading 
  
             submatrix A(1:k+1,1:k+1) */

	    i__1 = kp - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dswap_(&i__1, &A(1,k), &c__1, &A(1,kp), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qswap(&i__1, &A(1,k), &c__1, &A(1,kp), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qswap_(&i__1, &A(1,k), &c__1, &A(1,kp), &
#endif

		    c__1);
	    i__1 = k - kp - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dswap_(&i__1, &A(kp+1,k), &c__1, &A(kp,kp+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qswap(&i__1, &A(kp+1,k), &c__1, &A(kp,kp+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qswap_(&i__1, &A(kp+1,k), &c__1, &A(kp,kp+1), lda);
#endif

	    temp = A(k,k);
	    A(k,k) = A(kp,kp);
	    A(kp,kp) = temp;
	    if (kstep == 2) {
		temp = A(k,k+1);
		A(k,k+1) = A(kp,k+1);
		A(kp,k+1) = temp;
	    }
	}

	k += kstep;
	goto L30;
L40:

	;
    } else {

/*        Compute inv(A) from the factorization A = L*D*L'.   

          K is the main loop index, increasing from 1 to N in steps of
   
          1 or 2, depending on the size of the diagonal blocks. */

	k = *n;
L50:

/*        If K < 1, exit from loop. */

	if (k < 1) {
	    goto L60;
	}

	if (IPIV(k) > 0) {

/*           1 x 1 diagonal block   

             Invert the diagonal block. */

	    A(k,k) = 1. / A(k,k);

/*           Compute column K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
#endif

		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif

			 &WORK(1), &c__1, &c_b13, &A(k+1,k), &
			c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		A(k,k) -= ddot_(&i__1, &WORK(1), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k,k) -= qdot(&i__1, &WORK(1), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k,k) -= qdot_(&i__1, &WORK(1), &c__1, &A(k+1,k), &c__1);
#endif

	    }
	    kstep = 1;
	} else {

/*           2 x 2 diagonal block   

             Invert the diagonal block. */

	    t = (d__1 = A(k,k-1), ABS(d__1));
	    ak = A(k-1,k-1) / t;
	    akp1 = A(k,k) / t;
	    akkp1 = A(k,k-1) / t;
	    d = t * (ak * akp1 - 1.);
	    A(k-1,k-1) = akp1 / d;
	    A(k,k) = ak / d;
	    A(k,k-1) = -akkp1 / d;

/*           Compute columns K-1 and K of the inverse. */

	    if (k < *n) {
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &A(k+1,k), &c__1, &WORK(1), &c__1);
#endif

		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif

			 &WORK(1), &c__1, &c_b13, &A(k+1,k), &
			c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		A(k,k) -= ddot_(&i__1, &WORK(1), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k,k) -= qdot(&i__1, &WORK(1), &c__1, &A(k+1,k), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k,k) -= qdot_(&i__1, &WORK(1), &c__1, &A(k+1,k), &c__1);
#endif

		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		A(k,k-1) -= ddot_(&i__1, &A(k+1,k)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k,k-1) -= qdot(&i__1, &A(k+1,k)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k,k-1) -= qdot_(&i__1, &A(k+1,k)
#endif

			, &c__1, &A(k+1,k-1), &c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dcopy_(&i__1, &A(k+1,k-1), &c__1, &WORK(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qcopy(&i__1, &A(k+1,k-1), &c__1, &WORK(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qcopy_(&i__1, &A(k+1,k-1), &c__1, &WORK(1), &
#endif

			c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i__1, &c_b11, &A(k+1,k+1), lda,
#endif

			 &WORK(1), &c__1, &c_b13, &A(k+1,k-1)
			, &c__1);
		i__1 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		A(k-1,k-1) -= ddot_(&i__1, &WORK(1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		A(k-1,k-1) -= qdot(&i__1, &WORK(1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		A(k-1,k-1) -= qdot_(&i__1, &WORK(1), &c__1, &
#endif

			A(k+1,k-1), &c__1);
	    }
	    kstep = 2;
	}

	kp = (i__1 = IPIV(k), ABS(i__1));
	if (kp != k) {

/*           Interchange rows and columns K and KP in the trailing
   
             submatrix A(k-1:n,k-1:n) */

	    if (kp < *n) {
		i__1 = *n - kp;

#ifdef PETSC_PREFIX_SUFFIX
		dswap_(&i__1, &A(kp+1,k), &c__1, &A(kp+1,kp), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qswap(&i__1, &A(kp+1,k), &c__1, &A(kp+1,kp), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qswap_(&i__1, &A(kp+1,k), &c__1, &A(kp+1,kp), &c__1);
#endif

	    }
	    i__1 = kp - k - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dswap_(&i__1, &A(k+1,k), &c__1, &A(kp,k+1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qswap(&i__1, &A(k+1,k), &c__1, &A(kp,k+1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qswap_(&i__1, &A(k+1,k), &c__1, &A(kp,k+1), lda);
#endif

	    temp = A(k,k);
	    A(k,k) = A(kp,kp);
	    A(kp,kp) = temp;
	    if (kstep == 2) {
		temp = A(k,k-1);
		A(k,k-1) = A(kp,k-1);
		A(kp,k-1) = temp;
	    }
	}

	k -= kstep;
	goto L50;
L60:
	;
    }

    return;

/*     End of DSYTRI */

} /* dsytri_ */

