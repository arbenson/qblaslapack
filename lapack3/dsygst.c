#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsygst_(int *itype, char *uplo, int *n, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsygst(int *itype, char *uplo, int *n, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsygst_(int *itype, char *uplo, int *n, 
#endif

	LONG DOUBLE *a, int *lda, LONG DOUBLE *b, int *ldb, int *
	info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYGST reduces a real symmetric-definite generalized eigenproblem   
    to standard form.   

    If ITYPE = 1, the problem is A*x = lambda*B*x,   
    and A is overwritten by inv(U**T)*A*inv(U) or inv(L)*A*inv(L**T)   

    If ITYPE = 2 or 3, the problem is A*B*x = lambda*x or   
    B*A*x = lambda*x, and A is overwritten by U*A*U**T or L**T*A*L.   

    B must have been previously factorized as U**T*U or L*L**T by DPOTRF. 
  

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

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, if INFO = 0, the transformed matrix, stored in the   
            same format as A.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    B       (input) LONG DOUBLE PRECISION array, dimension (LDB,N)   
            The triangular factor from the Cholesky factorization of B,   
            as returned by DPOTRF.   

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
    static int c_n1 = -1;
    static LONG DOUBLE c_b14 = 1.;
    static LONG DOUBLE c_b16 = -.5;
    static LONG DOUBLE c_b19 = -1.;
    static LONG DOUBLE c_b52 = .5;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int k;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dtrmm_(char *, char *, char *, char *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmm(char *, char *, char *, char *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qtrmm_(char *, char *, char *, char *, 
#endif

	    int *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dsymm_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qsymm(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qsymm_(
#endif

	    char *, char *, int *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *, LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, 
	    int *);
    static long int upper;

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

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dsygs2_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qsygs2(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qsygs2_(
#endif

	    int *, char *, int *, LONG DOUBLE *, int *, LONG DOUBLE 
	    *, int *, int *);
    static int kb;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsyr2k_(char *, char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr2k(char *, char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr2k_(char *, char *, int *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, 
	    LONG DOUBLE *, LONG DOUBLE *, int *);
    static int nb;
    extern /* Subroutine */ void xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);




#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define B(I,J) b[(I)-1 + ((J)-1)* ( *ldb)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (*itype < 1 || *itype > 3) {
	*info = -1;
    } else if (! upper && ! lsame_(uplo, "L")) {
	*info = -2;
    } else if (*n < 0) {
	*info = -3;
    } else if (*lda < MAX(1,*n)) {
	*info = -5;
    } else if (*ldb < MAX(1,*n)) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYGST", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Determine the block size for this environment. */

    nb = ilaenv_(&c__1, "DSYGST", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);

    if (nb <= 1 || nb >= *n) {

/*        Use unblocked code */


#ifdef PETSC_PREFIX_SUFFIX
	dsygs2_(itype, uplo, n, &A(1,1), lda, &B(1,1), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsygs2(itype, uplo, n, &A(1,1), lda, &B(1,1), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsygs2_(itype, uplo, n, &A(1,1), lda, &B(1,1), ldb, info);
#endif

    } else {

/*        Use blocked code */

	if (*itype == 1) {
	    if (upper) {

/*              Compute inv(U')*A*inv(U) */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = MIN(i__3,nb);

/*                 Update the upper triangle of A(k:n,k:n)
 */


#ifdef PETSC_PREFIX_SUFFIX
		    dsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsygs2(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif

		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Left", uplo, "Transpose", "Non-unit", &kb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Left", uplo, "Transpose", "Non-unit", &kb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Left", uplo, "Transpose", "Non-unit", &kb, &
#endif

				i__3, &c_b14, &B(k,k), ldb, &A(k,k+kb), lda);
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsymm_("Left", uplo, &kb, &i__3, &c_b16, &A(k,k), lda, &B(k,k+kb), ldb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsymm("Left", uplo, &kb, &i__3, &c_b16, &A(k,k), lda, &B(k,k+kb), ldb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsymm_("Left", uplo, &kb, &i__3, &c_b16, &A(k,k), lda, &B(k,k+kb), ldb, 
#endif

				&c_b14, &A(k,k+kb), lda);
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsyr2k_(uplo, "Transpose", &i__3, &kb, &c_b19, &A(k,k+kb), lda, &B(k,k+kb), ldb, &c_b14, &A(k+kb,k+kb), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsyr2k(uplo, "Transpose", &i__3, &kb, &c_b19, &A(k,k+kb), lda, &B(k,k+kb), ldb, &c_b14, &A(k+kb,k+kb), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsyr2k_(uplo, "Transpose", &i__3, &kb, &c_b19, &A(k,k+kb), lda, &B(k,k+kb), ldb, &c_b14, &A(k+kb,k+kb), lda);
#endif

			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsymm_("Left", uplo, &kb, &i__3, &c_b16, &A(k,k), lda, &B(k,k+kb), ldb, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsymm("Left", uplo, &kb, &i__3, &c_b16, &A(k,k), lda, &B(k,k+kb), ldb, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsymm_("Left", uplo, &kb, &i__3, &c_b16, &A(k,k), lda, &B(k,k+kb), ldb, 
#endif

				&c_b14, &A(k,k+kb), lda);
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Right", uplo, "No transpose", "Non-unit", &kb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Right", uplo, "No transpose", "Non-unit", &kb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Right", uplo, "No transpose", "Non-unit", &kb,
#endif

				 &i__3, &c_b14, &B(k+kb,k+kb)
				, ldb, &A(k,k+kb), lda);
		    }
/* L10: */
		}
	    } else {

/*              Compute inv(L)*A*inv(L') */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = MIN(i__3,nb);

/*                 Update the lower triangle of A(k:n,k:n)
 */


#ifdef PETSC_PREFIX_SUFFIX
		    dsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsygs2(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif

		    if (k + kb <= *n) {
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Right", uplo, "Transpose", "Non-unit", &i__3, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Right", uplo, "Transpose", "Non-unit", &i__3, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Right", uplo, "Transpose", "Non-unit", &i__3, 
#endif

				&kb, &c_b14, &B(k,k), ldb, &A(k+kb,k), lda);
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsymm_("Right", uplo, &i__3, &kb, &c_b16, &A(k,k), lda, &B(k+kb,k), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsymm("Right", uplo, &i__3, &kb, &c_b16, &A(k,k), lda, &B(k+kb,k), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsymm_("Right", uplo, &i__3, &kb, &c_b16, &A(k,k), lda, &B(k+kb,k), ldb, &
#endif

				c_b14, &A(k+kb,k), lda);
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsyr2k_(uplo, "No transpose", &i__3, &kb, &c_b19, &A(k+kb,k), lda, &B(k+kb,k), ldb, &c_b14, &A(k+kb,k+kb), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsyr2k(uplo, "No transpose", &i__3, &kb, &c_b19, &A(k+kb,k), lda, &B(k+kb,k), ldb, &c_b14, &A(k+kb,k+kb), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsyr2k_(uplo, "No transpose", &i__3, &kb, &c_b19, &A(k+kb,k), lda, &B(k+kb,k), ldb, &c_b14, &A(k+kb,k+kb), lda);
#endif

			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dsymm_("Right", uplo, &i__3, &kb, &c_b16, &A(k,k), lda, &B(k+kb,k), ldb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qsymm("Right", uplo, &i__3, &kb, &c_b16, &A(k,k), lda, &B(k+kb,k), ldb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qsymm_("Right", uplo, &i__3, &kb, &c_b16, &A(k,k), lda, &B(k+kb,k), ldb, &
#endif

				c_b14, &A(k+kb,k), lda);
			i__3 = *n - k - kb + 1;

#ifdef PETSC_PREFIX_SUFFIX
			dtrsm_("Left", uplo, "No transpose", "Non-unit", &
#endif
#ifdef Q_C_PREFIX_SUFFIX
			qtrsm("Left", uplo, "No transpose", "Non-unit", &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
			qtrsm_("Left", uplo, "No transpose", "Non-unit", &
#endif

				i__3, &kb, &c_b14, &B(k+kb,k+kb), ldb, &A(k+kb,k), lda);
		    }
/* L20: */
		}
	    }
	} else {
	    if (upper) {

/*              Compute U*A*U' */

		i__1 = *n;
		i__2 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = MIN(i__3,nb);

/*                 Update the upper triangle of A(1:k+kb-1
,1:k+kb-1) */

		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrmm_("Left", uplo, "No transpose", "Non-unit", &i__3, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrmm("Left", uplo, "No transpose", "Non-unit", &i__3, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrmm_("Left", uplo, "No transpose", "Non-unit", &i__3, &
#endif

			    kb, &c_b14, &B(1,1), ldb, &A(1,k),
			     lda);
		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dsymm_("Right", uplo, &i__3, &kb, &c_b52, &A(k,k), lda, &B(1,k), ldb, &c_b14, &A(1,k), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsymm("Right", uplo, &i__3, &kb, &c_b52, &A(k,k), lda, &B(1,k), ldb, &c_b14, &A(1,k), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsymm_("Right", uplo, &i__3, &kb, &c_b52, &A(k,k), lda, &B(1,k), ldb, &c_b14, &A(1,k), lda);
#endif

		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dsyr2k_(uplo, "No transpose", &i__3, &kb, &c_b14, &A(1,k), lda, &B(1,k), ldb, &c_b14,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsyr2k(uplo, "No transpose", &i__3, &kb, &c_b14, &A(1,k), lda, &B(1,k), ldb, &c_b14,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsyr2k_(uplo, "No transpose", &i__3, &kb, &c_b14, &A(1,k), lda, &B(1,k), ldb, &c_b14,
#endif

			     &A(1,1), lda);
		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dsymm_("Right", uplo, &i__3, &kb, &c_b52, &A(k,k), lda, &B(1,k), ldb, &c_b14, &A(1,k), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsymm("Right", uplo, &i__3, &kb, &c_b52, &A(k,k), lda, &B(1,k), ldb, &c_b14, &A(1,k), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsymm_("Right", uplo, &i__3, &kb, &c_b52, &A(k,k), lda, &B(1,k), ldb, &c_b14, &A(1,k), lda);
#endif

		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrmm_("Right", uplo, "Transpose", "Non-unit", &i__3, &kb,
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrmm("Right", uplo, "Transpose", "Non-unit", &i__3, &kb,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrmm_("Right", uplo, "Transpose", "Non-unit", &i__3, &kb,
#endif

			     &c_b14, &B(k,k), ldb, &A(1,k), lda);

#ifdef PETSC_PREFIX_SUFFIX
		    dsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsygs2(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif

/* L30: */
		}
	    } else {

/*              Compute L'*A*L */

		i__2 = *n;
		i__1 = nb;
		for (k = 1; nb < 0 ? k >= *n : k <= *n; k += nb) {
/* Computing MIN */
		    i__3 = *n - k + 1;
		    kb = MIN(i__3,nb);

/*                 Update the lower triangle of A(1:k+kb-1
,1:k+kb-1) */

		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrmm_("Right", uplo, "No transpose", "Non-unit", &kb, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrmm("Right", uplo, "No transpose", "Non-unit", &kb, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrmm_("Right", uplo, "No transpose", "Non-unit", &kb, &
#endif

			    i__3, &c_b14, &B(1,1), ldb, &A(k,1), 
			    lda);
		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dsymm_("Left", uplo, &kb, &i__3, &c_b52, &A(k,k), lda, &B(k,1), ldb, &c_b14, &A(k,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsymm("Left", uplo, &kb, &i__3, &c_b52, &A(k,k), lda, &B(k,1), ldb, &c_b14, &A(k,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsymm_("Left", uplo, &kb, &i__3, &c_b52, &A(k,k), lda, &B(k,1), ldb, &c_b14, &A(k,1), lda);
#endif

		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dsyr2k_(uplo, "Transpose", &i__3, &kb, &c_b14, &A(k,1), lda, &B(k,1), ldb, &c_b14, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsyr2k(uplo, "Transpose", &i__3, &kb, &c_b14, &A(k,1), lda, &B(k,1), ldb, &c_b14, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsyr2k_(uplo, "Transpose", &i__3, &kb, &c_b14, &A(k,1), lda, &B(k,1), ldb, &c_b14, &A(1,1), lda);
#endif

		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dsymm_("Left", uplo, &kb, &i__3, &c_b52, &A(k,k), lda, &B(k,1), ldb, &c_b14, &A(k,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsymm("Left", uplo, &kb, &i__3, &c_b52, &A(k,k), lda, &B(k,1), ldb, &c_b14, &A(k,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsymm_("Left", uplo, &kb, &i__3, &c_b52, &A(k,k), lda, &B(k,1), ldb, &c_b14, &A(k,1), lda);
#endif

		    i__3 = k - 1;

#ifdef PETSC_PREFIX_SUFFIX
		    dtrmm_("Left", uplo, "Transpose", "Non-unit", &kb, &i__3, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qtrmm("Left", uplo, "Transpose", "Non-unit", &kb, &i__3, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qtrmm_("Left", uplo, "Transpose", "Non-unit", &kb, &i__3, 
#endif

			    &c_b14, &B(k,k), ldb, &A(k,1), 
			    lda);

#ifdef PETSC_PREFIX_SUFFIX
		    dsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qsygs2(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qsygs2_(itype, uplo, &kb, &A(k,k), lda, &B(k,k), ldb, info);
#endif

/* L40: */
		}
	    }
	}
    }
    return;

/*     End of DSYGST */

} /* dsygst_ */

