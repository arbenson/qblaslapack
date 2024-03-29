#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dpotri_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qpotri(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qpotri_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       March 31, 1993   


    Purpose   
    =======   

    DPOTRI computes the inverse of a real symmetric positive definite   
    matrix A using the Cholesky factorization A = U**T*U or A = L*L**T   
    computed by DPOTRF.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the triangular factor U or L from the Cholesky   
            factorization A = U**T*U or A = L*L**T, as computed by   
            DPOTRF.   
            On exit, the upper or lower triangle of the (symmetric)   
            inverse of A, overwriting the input factor U or L.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, the (i,i) element of the factor U or L is 
  
                  zero, and the inverse could not be computed.   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1;
    /* Local variables */
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), dlauum_(
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qlauum(
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void xerbla_(char *, int *), qlauum_(
#endif

	    char *, int *, LONG DOUBLE *, int *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dtrtri_(char *, char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qtrtri(char *, char *, int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qtrtri_(char *, char *, int *, LONG DOUBLE *, int *, 
#endif

	    int *);



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DPOTRI", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	return;
    }

/*     Invert the triangular Cholesky factor U or L. */


#ifdef PETSC_PREFIX_SUFFIX
    dtrtri_(uplo, "Non-unit", n, &A(1,1), lda, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qtrtri(uplo, "Non-unit", n, &A(1,1), lda, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qtrtri_(uplo, "Non-unit", n, &A(1,1), lda, info);
#endif

    if (*info > 0) {
	return;
    }

/*     Form inv(U)*inv(U)' or inv(L)'*inv(L). */


#ifdef PETSC_PREFIX_SUFFIX
    dlauum_(uplo, n, &A(1,1), lda, info);
#endif
#ifdef Q_C_PREFIX_SUFFIX
    qlauum(uplo, n, &A(1,1), lda, info);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    qlauum_(uplo, n, &A(1,1), lda, info);
#endif


    return;

/*     End of DPOTRI */

} /* dpotri_ */

