#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dtbtrs_(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qtbtrs(char *uplo, char *trans, char *diag, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qtbtrs_(char *uplo, char *trans, char *diag, int *n, 
#endif

	int *kd, int *nrhs, LONG DOUBLE *ab, int *ldab, LONG DOUBLE 
	*b, int *ldb, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DTBTRS solves a triangular system of the form   

       A * X = B  or  A**T * X = B,   

    where A is a triangular band matrix of order N, and B is an   
    N-by NRHS matrix.  A check is made to verify that A is nonsingular.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  A is upper triangular;   
            = 'L':  A is lower triangular.   

    TRANS   (input) CHARACTER*1   
            Specifies the form the system of equations:   
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
            of the matrix B.  NRHS >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The upper or lower triangular band matrix A, stored in the   
            first kd+1 rows of AB.  The j-th column of A is stored   
            in the j-th column of the array AB as follows:   
            if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=MIN(n,j+kd). 
  
            If DIAG = 'U', the diagonal elements of A are not referenced 
  
            and are assumed to be 1.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, if INFO = 0, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the i-th diagonal element of A is zero,   
                  indicating that the matrix is singular and the   
                  solutions X have not been computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1;
    /* Local variables */
    static int j;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtbsv_(char *, char *, char *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtbsv(char *, char *, char *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtbsv_(char *, char *, char *, int *, 
#endif

	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, int *);
    static long int upper;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int nounit;




#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    nounit = lsame_(diag, "N");
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
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
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DTBTRS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Check for singularity. */

    if (nounit) {
	if (upper) {
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		if (AB(*kd+1,*info) == 0.) {
		    return;
		}
/* L10: */
	    }
	} else {
	    i__1 = *n;
	    for (*info = 1; *info <= i__1; ++(*info)) {
		if (AB(1,*info) == 0.) {
		    return;
		}
/* L20: */
	    }
	}
    }
    *info = 0;

/*     Solve A * X = B  or  A' * X = B. */

    i__1 = *nrhs;
    for (j = 1; j <= *nrhs; ++j) {

#ifdef PETSC_PREFIX_SUFFIX
	dtbsv_(uplo, trans, diag, n, kd, &AB(1,1), ldab, &B(1,j), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtbsv(uplo, trans, diag, n, kd, &AB(1,1), ldab, &B(1,j), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtbsv_(uplo, trans, diag, n, kd, &AB(1,1), ldab, &B(1,j), &c__1);
#endif

/* L30: */
    }

    return;

/*     End of DTBTRS */

} /* dtbtrs_ */

