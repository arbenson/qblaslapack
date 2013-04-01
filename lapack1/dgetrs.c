#include <math.h>

#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgetrs_(char *trans, int *n, int *nrhs, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgetrs(char *trans, int *n, int *nrhs, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgetrs_(char *trans, int *n, int *nrhs, 
#endif

	LONG DOUBLE *a, int *lda, int *ipiv, LONG DOUBLE *b, int *
	ldb, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGETRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general N-by-N matrix A using the LU factorization computed   
    by DGETRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations:   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INT   
            The order of the matrix A.  N >= 0.   

    NRHS    (input) INT   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            The factors L and U from the factorization A = P*L*U   
            as computed by DGETRF.   

    LDA     (input) INT   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    IPIV    (input) INT array, dimension (N)   
            The pivot indices from DGETRF; for 1<=i<=N, row i of the   
            matrix was interchanged with row IPIV(i).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INT   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INT   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static LONG DOUBLE c_b12 = 1.;
    static int c_n1 = -1;
    
    /* System generated locals */
    int i__1;
    /* Local variables */
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtrsm_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrsm(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrsm_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *), xerbla_(

#ifdef PETSC_PREFIX_SUFFIX
	    char *, int *), dlaswp_(int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    char *, int *), qlaswp(int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    char *, int *), qlaswp_(int *, LONG DOUBLE *, 
#endif

	    int *, int *, int *, int *, int *);
    static long int notran;



#define IPIV(I) ipiv[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    notran = lsame_(trans, "N");
    if (! notran && ! lsame_(trans, "T") && ! lsame_(trans, "C")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*nrhs < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -8;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGETRS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return;
    }

    if (notran) {

/*        Solve A * X = B.   

          Apply row interchanges to the right hand sides. */


#ifdef PETSC_PREFIX_SUFFIX
	dlaswp_(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaswp(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaswp_(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c__1);
#endif


/*        Solve L*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dtrsm_("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrsm("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrsm_("Left", "Lower", "No transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif


/*        Solve U*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dtrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrsm("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrsm_("Left", "Upper", "No transpose", "Non-unit", n, nrhs, &c_b12, &
#endif

		A(1,1), lda, &B(1,1), ldb);
    } else {

/*        Solve A' * X = B.   

          Solve U'*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dtrsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrsm("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrsm_("Left", "Upper", "Transpose", "Non-unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif


/*        Solve L'*X = B, overwriting B with X. */


#ifdef PETSC_PREFIX_SUFFIX
	dtrsm_("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qtrsm("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qtrsm_("Left", "Lower", "Transpose", "Unit", n, nrhs, &c_b12, &A(1,1), lda, &B(1,1), ldb);
#endif


/*        Apply row interchanges to the solution vectors. */


#ifdef PETSC_PREFIX_SUFFIX
	dlaswp_(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c_n1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlaswp(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c_n1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlaswp_(nrhs, &B(1,1), ldb, &c__1, n, &IPIV(1), &c_n1);
#endif

    }

    return;

/*     End of DGETRS */

} /* dgetrs_ */

