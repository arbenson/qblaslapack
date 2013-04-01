#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dspgst_(int *itype, char *uplo, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qspgst(int *itype, char *uplo, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qspgst_(int *itype, char *uplo, int *n, 
#endif

	LONG DOUBLE *ap, LONG DOUBLE *bp, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DSPGST reduces a real symmetric-definite generalized eigenproblem   
    to standard form, using packed storage.   

    If ITYPE = 1, the problem is A*x = lambda*B*x,   
    and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)   

    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or   
    B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.   

    B must have been previously factorized as U**T*U or L*L**T by DPPTRF. 
  

    Arguments   
    =========   

    ITYPE   (input) INTEGER   
            = 1: compute inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T);   
            = 2 or 3: compute U*A*U**T or L**T*A*L.   

    UPLO    (input) CHARACTER   
            = 'U':  Upper triangle of A is stored and B is factored as   
                    U**T*U;   
            = 'L':  Lower triangle of A is stored and B is factored as   
                    L*L**T.   

    N       (input) INTEGER   
            The order of the matrices A and B.  N >= 0.   

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the upper or lower triangle of the symmetric matrix 
  
            A, packed columnwise in a linear array.  The j-th column of A 
  
            is stored in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   

            On exit, if INFO = 0, the transformed matrix, stored in the   
            same format as A.   

    BP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The triangular factor from the Cholesky factorization of B,   
            stored in the same format as A, as returned by DPPTRF.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b9 = -1.;
    static LONG DOUBLE c_b11 = 1.;
    
    /* System generated locals */
    int i__1, i__2;
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

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dspr2_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qspr2(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qspr2_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *);
    static int j, k;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qscal_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *);
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), dspmv_(char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qspmv(char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, int *), qspmv_(char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *,
	     LONG DOUBLE *, int *);
    static long int upper;
    static int P_j1, k1;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtpmv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpmv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpmv_(char *, char *, char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dtpsv_(char *, char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtpsv(char *, char *, char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtpsv_(char *, char *, char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *);
    static int jj, kk;
    static LONG DOUBLE ct;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static LONG DOUBLE ajj;
    static int P_j1P_j1;
    static LONG DOUBLE akk;
    static int k1k1;
    static LONG DOUBLE bjj, bkk;



#define BP(I) bp[(I)-1]
#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSPGST", &i__1);
	return;
    }

    if (*itype == 1) {
	if (upper) {

/*           Compute inv(U')*A*inv(U)   

             P_J1 and JJ are the indices of A(1,j) and A(j,j) */

	    jj = 0;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		P_j1 = jj + 1;
		jj += j;

/*              Compute the j-th column of the upper triangle 
of A */

		bjj = BP(jj);

#ifdef PETSC_PREFIX_SUFFIX
		dtpsv_(uplo, "Transpose", "Nonunit", &j, &BP(1), &AP(P_j1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtpsv(uplo, "Transpose", "Nonunit", &j, &BP(1), &AP(P_j1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtpsv_(uplo, "Transpose", "Nonunit", &j, &BP(1), &AP(P_j1), &
#endif

			c__1);
		i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__2, &c_b9, &AP(1), &BP(P_j1), &c__1, &c_b11, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__2, &c_b9, &AP(1), &BP(P_j1), &c__1, &c_b11, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__2, &c_b9, &AP(1), &BP(P_j1), &c__1, &c_b11, &
#endif

			AP(P_j1), &c__1);
		i__2 = j - 1;
		d__1 = 1. / bjj;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &d__1, &AP(P_j1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &d__1, &AP(P_j1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &d__1, &AP(P_j1), &c__1);
#endif

		i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		AP(jj) = (AP(jj) - ddot_(&i__2, &AP(P_j1), &c__1, &BP(P_j1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(jj) = (AP(jj) - qdot(&i__2, &AP(P_j1), &c__1, &BP(P_j1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(jj) = (AP(jj) - qdot_(&i__2, &AP(P_j1), &c__1, &BP(P_j1), &
#endif

			c__1)) / bjj;
/* L10: */
	    }
	} else {

/*           Compute inv(L)*A*inv(L')   

             KK and K1K1 are the indices of A(k,k) and A(k+1,k+1) 
*/

	    kk = 1;
	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {
		k1k1 = kk + *n - k + 1;

/*              Update the lower triangle of A(k:n,k:n) */

		akk = AP(kk);
		bkk = BP(kk);
/* Computing 2nd power */
		d__1 = bkk;
		akk /= d__1 * d__1;
		AP(kk) = akk;
		if (k < *n) {
		    i__2 = *n - k;
		    d__1 = 1. / bkk;

#ifdef PETSC_PREFIX_SUFFIX
		    dscal_(&i__2, &d__1, &AP(kk + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qscal(&i__2, &d__1, &AP(kk + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qscal_(&i__2, &d__1, &AP(kk + 1), &c__1);
#endif

		    ct = akk * -.5;
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&i__2, &ct, &BP(kk + 1), &c__1, &AP(kk + 1), &c__1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&i__2, &ct, &BP(kk + 1), &c__1, &AP(kk + 1), &c__1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&i__2, &ct, &BP(kk + 1), &c__1, &AP(kk + 1), &c__1)
#endif

			    ;
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    dspr2_(uplo, &i__2, &c_b9, &AP(kk + 1), &c__1, &BP(kk + 1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qspr2(uplo, &i__2, &c_b9, &AP(kk + 1), &c__1, &BP(kk + 1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qspr2_(uplo, &i__2, &c_b9, &AP(kk + 1), &c__1, &BP(kk + 1)
#endif

			    , &c__1, &AP(k1k1));
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    daxpy_(&i__2, &ct, &BP(kk + 1), &c__1, &AP(kk + 1), &c__1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qaxpy(&i__2, &ct, &BP(kk + 1), &c__1, &AP(kk + 1), &c__1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qaxpy_(&i__2, &ct, &BP(kk + 1), &c__1, &AP(kk + 1), &c__1)
#endif

			    ;
		    i__2 = *n - k;

#ifdef PETSC_PREFIX_SUFFIX
		    dtpsv_(uplo, "No transpose", "Non-unit", &i__2, &BP(k1k1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtpsv(uplo, "No transpose", "Non-unit", &i__2, &BP(k1k1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtpsv_(uplo, "No transpose", "Non-unit", &i__2, &BP(k1k1),
#endif

			     &AP(kk + 1), &c__1);
		}
		kk = k1k1;
/* L20: */
	    }
	}
    } else {
	if (upper) {

/*           Compute U*A*U'   

             K1 and KK are the indices of A(1,k) and A(k,k) */

	    kk = 0;
	    i__1 = *n;
	    for (k = 1; k <= *n; ++k) {
		k1 = kk + 1;
		kk += k;

/*              Update the upper triangle of A(1:k,1:k) */

		akk = AP(kk);
		bkk = BP(kk);
		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtpmv_(uplo, "No transpose", "Non-unit", &i__2, &BP(1), &AP(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtpmv(uplo, "No transpose", "Non-unit", &i__2, &BP(1), &AP(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtpmv_(uplo, "No transpose", "Non-unit", &i__2, &BP(1), &AP(
#endif

			k1), &c__1);
		ct = akk * .5;
		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &ct, &BP(k1), &c__1, &AP(k1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &ct, &BP(k1), &c__1, &AP(k1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &ct, &BP(k1), &c__1, &AP(k1), &c__1);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dspr2_(uplo, &i__2, &c_b11, &AP(k1), &c__1, &BP(k1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspr2(uplo, &i__2, &c_b11, &AP(k1), &c__1, &BP(k1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspr2_(uplo, &i__2, &c_b11, &AP(k1), &c__1, &BP(k1), &c__1, &
#endif

			AP(1));
		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &ct, &BP(k1), &c__1, &AP(k1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &ct, &BP(k1), &c__1, &AP(k1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &ct, &BP(k1), &c__1, &AP(k1), &c__1);
#endif

		i__2 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &bkk, &AP(k1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &bkk, &AP(k1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &bkk, &AP(k1), &c__1);
#endif

/* Computing 2nd power */
		d__1 = bkk;
		AP(kk) = akk * (d__1 * d__1);
/* L30: */
	    }
	} else {

/*           Compute L'*A*L   

             JJ and P_J1P_J1 are the indices of A(j,j) and A(j+1,j+1) 
*/

	    jj = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		P_j1P_j1 = jj + *n - j + 1;

/*              Compute the j-th column of the lower triangle 
of A */

		ajj = AP(jj);
		bjj = BP(jj);
		i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		AP(jj) = ajj * bjj + ddot_(&i__2, &AP(jj + 1), &c__1, &BP(jj 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		AP(jj) = ajj * bjj + qdot(&i__2, &AP(jj + 1), &c__1, &BP(jj 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		AP(jj) = ajj * bjj + qdot_(&i__2, &AP(jj + 1), &c__1, &BP(jj 
#endif

			+ 1), &c__1);
		i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &bjj, &AP(jj + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &bjj, &AP(jj + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &bjj, &AP(jj + 1), &c__1);
#endif

		i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dspmv_(uplo, &i__2, &c_b11, &AP(P_j1P_j1), &BP(jj + 1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspmv(uplo, &i__2, &c_b11, &AP(P_j1P_j1), &BP(jj + 1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspmv_(uplo, &i__2, &c_b11, &AP(P_j1P_j1), &BP(jj + 1), &c__1, &
#endif

			c_b11, &AP(jj + 1), &c__1);
		i__2 = *n - j + 1;

#ifdef PETSC_PREFIX_SUFFIX
		dtpmv_(uplo, "Transpose", "Non-unit", &i__2, &BP(jj), &AP(jj),
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtpmv(uplo, "Transpose", "Non-unit", &i__2, &BP(jj), &AP(jj),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtpmv_(uplo, "Transpose", "Non-unit", &i__2, &BP(jj), &AP(jj),
#endif

			 &c__1);
		jj = P_j1P_j1;
/* L40: */
	    }
	}
    }
    return;

/*     End of DSPGST */

} /* dspgst_ */

