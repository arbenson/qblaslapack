#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpptri_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpptri(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpptri_(char *uplo, int *n, LONG DOUBLE *ap, int *
#endif

	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPPTRI computes the inverse of a real symmetric positive definite   
    matrix A using the Cholesky factorization A = U**T*U or A = L*L**T   
    computed by DPPTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangular factor is stored in AP;   
            = 'L':  Lower triangular factor is stored in AP.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the triangular factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T, packed columnwise as 
  
            a linear array.  The j-th column of U or L is stored in the   
            array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = U(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = L(i,j) for j<=i<=n.   

            On exit, the upper or lower triangle of the (symmetric)   
            inverse of A, overwriting the input factor U or L.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the (i,i) element of the factor U or L is 
  
                  zero, and the inverse could not be computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b8 = 1.;
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1, i__2;
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
    extern /* Subroutine */ void dspr_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qspr(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qspr_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *);
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

	    int *);
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtpmv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpmv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtpmv_(char *, char *, char *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static long int upper;
    static int jc, jj;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dtptri_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qtptri(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qtptri_(
#endif

	    char *, char *, int *, LONG DOUBLE *, int *);
    static LONG DOUBLE ajj;
    static int jjn;



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
	xerbla_("DPPTRI", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Invert the triangular Cholesky factor U or L. */


#ifdef PETSC_PREFIX_SUFFIX
    dtptri_(uplo, "Non-unit", n, &AP(1), info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtptri(uplo, "Non-unit", n, &AP(1), info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtptri_(uplo, "Non-unit", n, &AP(1), info);
#endif

    if (*info > 0) {
	return;
    }

    if (upper) {

/*        Compute the product inv(U) * inv(U)'. */

	jj = 0;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    jc = jj + 1;
	    jj += j;
	    if (j > 1) {
		i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dspr_("Upper", &i__2, &c_b8, &AP(jc), &c__1, &AP(1));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qspr("Upper", &i__2, &c_b8, &AP(jc), &c__1, &AP(1));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qspr_("Upper", &i__2, &c_b8, &AP(jc), &c__1, &AP(1));
#endif

	    }
	    ajj = AP(jj);

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&j, &ajj, &AP(jc), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&j, &ajj, &AP(jc), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&j, &ajj, &AP(jc), &c__1);
#endif

/* L10: */
	}

    } else {

/*        Compute the product inv(L)' * inv(L). */

	jj = 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    jjn = jj + *n - j + 1;
	    i__2 = *n - j + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    AP(jj) = ddot_(&i__2, &AP(jj), &c__1, &AP(jj), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    AP(jj) = qdot(&i__2, &AP(jj), &c__1, &AP(jj), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    AP(jj) = qdot_(&i__2, &AP(jj), &c__1, &AP(jj), &c__1);
#endif

	    if (j < *n) {
		i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dtpmv_("Lower", "Transpose", "Non-unit", &i__2, &AP(jjn), &AP(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtpmv("Lower", "Transpose", "Non-unit", &i__2, &AP(jjn), &AP(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtpmv_("Lower", "Transpose", "Non-unit", &i__2, &AP(jjn), &AP(
#endif

			jj + 1), &c__1);
	    }
	    jj = jjn;
/* L20: */
	}
    }

    return;

/*     End of DPPTRI */

} /* dpptri_ */

