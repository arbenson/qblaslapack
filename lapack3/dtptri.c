#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtptri_(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtptri(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtptri_(char *uplo, char *diag, int *n, LONG DOUBLE *
#endif

	ap, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DTPTRI computes the inverse of a real upper or lower triangular   
    matrix A stored in packed format.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    DIAG    (input) CHARACTER*1   
            = 'N':  A is non-unit triangular;   
            = 'U':  A is unit triangular.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    AP      (input/output) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2) 
  
            On entry, the upper or lower triangular matrix A, stored   
            columnwise in a linear array.  The j-th column of A is stored 
  
            in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*((2*n-j)/2) = A(i,j) for j<=i<=n. 
  
            See below for further details.   
            On exit, the (triangular) inverse of the original matrix, in 
  
            the same packed storage format.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, A(i,i) is exactly zero.  The triangular   
                  matrix is singular and its inverse can not be computed. 
  

    Further Details   
    ===============   

    A triangular matrix A can be transferred to packed storage using one 
  
    of the following program segments:   

    UPLO = 'U':                      UPLO = 'L':   

          JC = 1                           JC = 1   
          DO 2 J = 1, N                    DO 2 J = 1, N   
             DO 1 I = 1, J                    DO 1 I = J, N   
                AP(JC+I-1) = A(I,J)              AP(JC+I-J) = A(I,J)   
        1    CONTINUE                    1    CONTINUE   
             JC = JC + J                      JC = JC + N - J + 1   
        2 CONTINUE                       2 CONTINUE   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1, i__2;
    /* Local variables */
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
    extern /* Subroutine */ void xerbla_(char *, int *);
    static int jclast;
    static long int nounit;
    static LONG DOUBLE ajj;



#define AP(I) ap[(I)-1]


    *info = 0;
    upper = lsame_(uplo, "U");
    nounit = lsame_(diag, "N");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! nounit && ! lsame_(diag, "U")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTPTRI", &i__1);
	return;
    }

/*     Check for singularity if non-unit. */

    if (nounit) {
	if (upper) {
	    jj = 0;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		jj += *info;
		if (AP(jj) == 0.) {
		    return;
		}
/* L10: */
	    }
	} else {
	    jj = 1;
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		if (AP(jj) == 0.) {
		    return;
		}
		jj = jj + *n - *info + 1;
/* L20: */
	    }
	}
	*info = 0;
    }

    if (upper) {

/*        Compute inverse of upper triangular matrix. */

	jc = 1;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    if (nounit) {
		AP(jc + j - 1) = 1. / AP(jc + j - 1);
		ajj = -AP(jc + j - 1);
	    } else {
		ajj = -1.;
	    }

/*           Compute elements 1:j-1 of j-th column. */

	    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dtpmv_("Upper", "No transpose", diag, &i__2, &AP(1), &AP(jc), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtpmv("Upper", "No transpose", diag, &i__2, &AP(1), &AP(jc), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtpmv_("Upper", "No transpose", diag, &i__2, &AP(1), &AP(jc), &
#endif

		    c__1);
	    i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dscal_(&i__2, &ajj, &AP(jc), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qscal(&i__2, &ajj, &AP(jc), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qscal_(&i__2, &ajj, &AP(jc), &c__1);
#endif

	    jc += j;
/* L30: */
	}

    } else {

/*        Compute inverse of lower triangular matrix. */

	jc = *n * (*n + 1) / 2;
	for (j = *n; j >= 1; --j) {
	    if (nounit) {
		AP(jc) = 1. / AP(jc);
		ajj = -AP(jc);
	    } else {
		ajj = -1.;
	    }
	    if (j < *n) {

/*              Compute elements j+1:n of j-th column. */

		i__1 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dtpmv_("Lower", "No transpose", diag, &i__1, &AP(jclast), &AP(
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qtpmv("Lower", "No transpose", diag, &i__1, &AP(jclast), &AP(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qtpmv_("Lower", "No transpose", diag, &i__1, &AP(jclast), &AP(
#endif

			jc + 1), &c__1);
		i__1 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__1, &ajj, &AP(jc + 1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__1, &ajj, &AP(jc + 1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__1, &ajj, &AP(jc + 1), &c__1);
#endif

	    }
	    jclast = jc;
	    jc = jc - *n + j - 2;
/* L40: */
	}
    }

    return;

/*     End of DTPTRI */

} /* dtptri_ */

