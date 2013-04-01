#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpbtrs_(char *uplo, int *n, int *kd, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpbtrs(char *uplo, int *n, int *kd, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpbtrs_(char *uplo, int *n, int *kd, int *
#endif

	nrhs, LONG DOUBLE *ab, int *ldab, LONG DOUBLE *b, int *ldb, 
	int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DPBTRS solves a system of linear equations A*X = B with a symmetric   
    positive definite band matrix A using the Cholesky factorization   
    A = U**T*U or A = L*L**T computed by DPBTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangular factor stored in AB;   
            = 'L':  Lower triangular factor stored in AB.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KD      (input) INTEGER   
            The number of superdiagonals of the matrix A if UPLO = 'U',   
            or the number of subdiagonals if UPLO = 'L'.  KD >= 0.   

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The triangular factor U or L from the Cholesky factorization 
  
            A = U**T*U or A = L*L**T of the band matrix A, stored in the 
  
            first KD+1 rows of the array.  The j-th column of U or L is   
            stored in the j-th column of the array AB as follows:   
            if UPLO ='U', AB(kd+1+i-j,j) = U(i,j) for MAX(1,j-kd)<=i<=j; 
  
            if UPLO ='L', AB(1+i-j,j)    = L(i,j) for j<=i<=MIN(n,j+kd). 
  

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= KD+1.   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

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




#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

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
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPBTRS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return;
    }

    if (upper) {

/*        Solve A*X = B where A = U'*U. */

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve U'*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtbsv_("Upper", "Transpose", "Non-unit", n, kd, &AB(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtbsv("Upper", "Transpose", "Non-unit", n, kd, &AB(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtbsv_("Upper", "Transpose", "Non-unit", n, kd, &AB(1,1), 
#endif

		    ldab, &B(1,j), &c__1);

/*           Solve U*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtbsv_("Upper", "No transpose", "Non-unit", n, kd, &AB(1,1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtbsv("Upper", "No transpose", "Non-unit", n, kd, &AB(1,1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtbsv_("Upper", "No transpose", "Non-unit", n, kd, &AB(1,1),
#endif

		     ldab, &B(1,j), &c__1);
/* L10: */
	}
    } else {

/*        Solve A*X = B where A = L*L'. */

	i__1 = *nrhs;
	for (j = 1; j <= *nrhs; ++j) {

/*           Solve L*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtbsv_("Lower", "No transpose", "Non-unit", n, kd, &AB(1,1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtbsv("Lower", "No transpose", "Non-unit", n, kd, &AB(1,1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtbsv_("Lower", "No transpose", "Non-unit", n, kd, &AB(1,1),
#endif

		     ldab, &B(1,j), &c__1);

/*           Solve L'*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	    dtbsv_("Lower", "Transpose", "Non-unit", n, kd, &AB(1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtbsv("Lower", "Transpose", "Non-unit", n, kd, &AB(1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtbsv_("Lower", "Transpose", "Non-unit", n, kd, &AB(1,1), 
#endif

		    ldab, &B(1,j), &c__1);
/* L20: */
	}
    }

    return;

/*     End of DPBTRS */

} /* dpbtrs_ */

