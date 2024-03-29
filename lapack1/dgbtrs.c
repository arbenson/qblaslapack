#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dgbtrs_(char *trans, int *n, int *kl, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qgbtrs(char *trans, int *n, int *kl, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qgbtrs_(char *trans, int *n, int *kl, int *
#endif

	ku, int *nrhs, LONG DOUBLE *ab, int *ldab, int *ipiv, 
	LONG DOUBLE *b, int *ldb, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DGBTRS solves a system of linear equations   
       A * X = B  or  A' * X = B   
    with a general band matrix A using the LU factorization computed   
    by DGBTRF.   

    Arguments   
    =========   

    TRANS   (input) CHARACTER*1   
            Specifies the form of the system of equations.   
            = 'N':  A * X = B  (No transpose)   
            = 'T':  A'* X = B  (Transpose)   
            = 'C':  A'* X = B  (Conjugate transpose = Transpose)   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    KL      (input) INTEGER   
            The number of subdiagonals within the band of A.  KL >= 0.   

    KU      (input) INTEGER   
            The number of superdiagonals within the band of A.  KU >= 0. 
  

    NRHS    (input) INTEGER   
            The number of right hand sides, i.e., the number of columns   
            of the matrix B.  NRHS >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            Details of the LU factorization of the band matrix A, as   
            computed by DGBTRF.  U is stored as an upper triangular band 
  
            matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and   
            the multipliers used during the factorization are stored in   
            rows KL+KU+2 to 2*KL+KU+1.   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.   

    IPIV    (input) INTEGER array, dimension (N)   
            The pivot indices; for 1 <= i <= N, row i of the matrix was   
            interchanged with row IPIV(i).   

    B       (input/output) LONG DOUBLE PRECISION array, dimension (LDB,NRHS)   
            On entry, the right hand side matrix B.   
            On exit, the solution matrix X.   

    LDB     (input) INTEGER   
            The leading dimension of the array B.  LDB >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0: if INFO = -i, the i-th argument had an illegal value   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b7 = -1.;
    static int c__1 = 1;
    static LONG DOUBLE c_b23 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dger_(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qger(int *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qger_(int *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static int i, j, l;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dgemv_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qgemv_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), dswap_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qswap(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qswap_(int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), dtbsv_(char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtbsv(char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, int *), qtbsv_(char *, 
#endif

	    char *, char *, int *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, int *);
    static long int lnoti;
    static int kd, lm;
    extern /* Subroutine */ void xerbla_(char *, int *);
    static long int notran;



#define IPIV(I) ipiv[(I)-1]

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

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
    } else if (*ldab < (*kl << 1) + *ku + 1) {
	*info = -7;
    } else if (*ldb < MAX(1,*n)) {
	*info = -10;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DGBTRS", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0 || *nrhs == 0) {
	return;
    }

    kd = *ku + *kl + 1;
    lnoti = *kl > 0;

    if (notran) {

/*        Solve  A*X = B.   

          Solve L*X = B, overwriting B with X.   

          L is represented as a product of permutations and unit lower
   
          triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
   
          where each transformation L(i) is a rank-one modification of
   
          the identity matrix. */

	if (lnoti) {
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
/* Computing MIN */
		i__2 = *kl, i__3 = *n - j;
		lm = MIN(i__2,i__3);
		l = IPIV(j);
		if (l != j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dswap_(nrhs, &B(l,1), ldb, &B(j,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qswap(nrhs, &B(l,1), ldb, &B(j,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qswap_(nrhs, &B(l,1), ldb, &B(j,1), ldb);
#endif

		}

#ifdef PETSC_PREFIX_SUFFIX
		dger_(&lm, nrhs, &c_b7, &AB(kd+1,j), &c__1, &B(j,1), ldb, &B(j+1,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qger(&lm, nrhs, &c_b7, &AB(kd+1,j), &c__1, &B(j,1), ldb, &B(j+1,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qger_(&lm, nrhs, &c_b7, &AB(kd+1,j), &c__1, &B(j,1), ldb, &B(j+1,1), ldb);
#endif

/* L10: */
	    }
	}

	i__1 = *nrhs;
	for (i = 1; i <= *nrhs; ++i) {

/*           Solve U*X = B, overwriting B with X. */

	    i__2 = *kl + *ku;

#ifdef PETSC_PREFIX_SUFFIX
	    dtbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &AB(1,1), ldab, &B(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtbsv("Upper", "No transpose", "Non-unit", n, &i__2, &AB(1,1), ldab, &B(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtbsv_("Upper", "No transpose", "Non-unit", n, &i__2, &AB(1,1), ldab, &B(1,i), &c__1);
#endif

/* L20: */
	}

    } else {

/*        Solve A'*X = B. */

	i__1 = *nrhs;
	for (i = 1; i <= *nrhs; ++i) {

/*           Solve U'*X = B, overwriting B with X. */

	    i__2 = *kl + *ku;

#ifdef PETSC_PREFIX_SUFFIX
	    dtbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &AB(1,1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtbsv("Upper", "Transpose", "Non-unit", n, &i__2, &AB(1,1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtbsv_("Upper", "Transpose", "Non-unit", n, &i__2, &AB(1,1),
#endif

		     ldab, &B(1,i), &c__1);
/* L30: */
	}

/*        Solve L'*X = B, overwriting B with X. */

	if (lnoti) {
	    for (j = *n - 1; j >= 1; --j) {
/* Computing MIN */
		i__1 = *kl, i__2 = *n - j;
		lm = MIN(i__1,i__2);

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("Transpose", &lm, nrhs, &c_b7, &B(j+1,1), ldb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("Transpose", &lm, nrhs, &c_b7, &B(j+1,1), ldb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("Transpose", &lm, nrhs, &c_b7, &B(j+1,1), ldb,
#endif

			 &AB(kd+1,j), &c__1, &c_b23, &B(j,1), ldb);
		l = IPIV(j);
		if (l != j) {

#ifdef PETSC_PREFIX_SUFFIX
		    dswap_(nrhs, &B(l,1), ldb, &B(j,1), ldb);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qswap(nrhs, &B(l,1), ldb, &B(j,1), ldb);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qswap_(nrhs, &B(l,1), ldb, &B(j,1), ldb);
#endif

		}
/* L40: */
	    }
	}
    }
    return;

/*     End of DGBTRS */

} /* dgbtrs_ */

