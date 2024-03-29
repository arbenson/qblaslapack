#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsytrf_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsytrf(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsytrf_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, int *ipiv, LONG DOUBLE *work, int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYTRF computes the factorization of a real symmetric matrix A using 
  
    the Bunch-Kaufman diagonal pivoting method.  The form of the   
    factorization is   

       A = U*D*U**T  or  A = L*D*L**T   

    where U (or L) is a product of permutation and unit upper (lower)   
    triangular matrices, and D is symmetric and block diagonal with   
    1-by-1 and 2-by-2 diagonal blocks.   

    This is the blocked version of the algorithm, calling Level 3 BLAS.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            = 'U':  Upper triangle of A is stored;   
            = 'L':  Lower triangle of A is stored.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            N-by-N upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading N-by-N lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   

            On exit, the block diagonal matrix D and the multipliers used 
  
            to obtain the factor U or L (see below for further details). 
  

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    IPIV    (output) INTEGER array, dimension (N)   
            Details of the interchanges and the block structure of D.   
            If IPIV(k) > 0, then rows and columns k and IPIV(k) were   
            interchanged and D(k,k) is a 1-by-1 diagonal block.   
            If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and   
            columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k) 
  
            is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =   
            IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were   
            interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The length of WORK.  LWORK >=1.  For best performance   
            LWORK >= N*NB, where NB is the block size returned by ILAENV. 
  

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   
            > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization 
  
                  has been completed, but the block diagonal matrix D is 
  
                  exactly singular, and division by zero will occur if it 
  
                  is used to solve a system of equations.   

    Further Details   
    ===============   

    If UPLO = 'U', then A = U*D*U', where   
       U = P(n)*U(n)* ... *P(k)U(k)* ...,   
    i.e., U is a product of terms P(k)*U(k), where k decreases from n to 
  
    1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1   
    and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as   
    defined by IPIV(k), and U(k) is a unit upper triangular matrix, such 
  
    that if the diagonal block D(k) is of order s (s = 1 or 2), then   

               (   I    v    0   )   k-s   
       U(k) =  (   0    I    0   )   s   
               (   0    0    I   )   n-k   
                  k-s   s   n-k   

    If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).   
    If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k), 
  
    and A(k,k), and v overwrites A(1:k-2,k-1:k).   

    If UPLO = 'L', then A = L*D*L', where   
       L = P(1)*L(1)* ... *P(k)*L(k)* ...,   
    i.e., L is a product of terms P(k)*L(k), where k increases from 1 to 
  
    n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1   
    and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as   
    defined by IPIV(k), and L(k) is a unit lower triangular matrix, such 
  
    that if the diagonal block D(k) is of order s (s = 1 or 2), then   

               (   I    0     0   )  k-1   
       L(k) =  (   0    I     0   )  s   
               (   0    v     I   )  n-k-s+1   
                  k-1   s  n-k-s+1   

    If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).   
    If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),   
    and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).   

    ===================================================================== 
  


       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static int c__2 = 2;
    
    /* System generated locals */
    int  i__1, i__2;
    /* Local variables */
    static int j, k;
    extern long int lsame_(char *, char *);
    static int nbmin, iinfo;
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsytf2_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsytf2(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsytf2_(char *, int *, LONG DOUBLE *, 
#endif

	    int *, int *, int *);
    static int kb, nb;
    extern /* Subroutine */ void xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlasyf_(char *, int *, int *, int 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasyf(char *, int *, int *, int 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlasyf_(char *, int *, int *, int 
#endif

	    *, LONG DOUBLE *, int *, int *, LONG DOUBLE *, int *, 
	    int *);
    static int ldwork, iws;



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
    } else if (*lwork < 1) {
	*info = -7;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYTRF", &i__1);
	return;
    }

/*     Determine the block size */

    nb = ilaenv_(&c__1, "DSYTRF", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);
    nbmin = 2;
    ldwork = *n;
    if (nb > 1 && nb < *n) {
	iws = ldwork * nb;
	if (*lwork < iws) {
/* Computing MAX */
	    i__1 = *lwork / ldwork;
	    nb = MAX(i__1,1);
/* Computing MAX */
	    i__1 = 2, i__2 = ilaenv_(&c__2, "DSYTRF", uplo, n, &c_n1, &c_n1, &
		    c_n1, 6L, 1L);
	    nbmin = MAX(i__1,i__2);
	}
    } else {
	iws = 1;
    }
    if (nb < nbmin) {
	nb = *n;
    }

    if (upper) {

/*        Factorize A as U*D*U' using the upper triangle of A   

          K is the main loop index, decreasing from N to 1 in steps of
   
          KB, where KB is the number of columns factorized by DLASYF; 
  
          KB is either NB or NB-1, or K for the last block */

	k = *n;
L10:

/*        If K < 1, exit from loop */

	if (k < 1) {
	    goto L40;
	}

	if (k > nb) {

/*           Factorize columns k-kb+1:k of A and use blocked code 
to   
             update columns 1:k-kb */


#ifdef PETSC_PREFIX_SUFFIX
	    dlasyf_(uplo, &k, &nb, &kb, &A(1,1), lda, &IPIV(1), &WORK(1),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasyf(uplo, &k, &nb, &kb, &A(1,1), lda, &IPIV(1), &WORK(1),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasyf_(uplo, &k, &nb, &kb, &A(1,1), lda, &IPIV(1), &WORK(1),
#endif

		     &ldwork, &iinfo);
	} else {

/*           Use unblocked code to factorize columns 1:k of A */


#ifdef PETSC_PREFIX_SUFFIX
	    dsytf2_(uplo, &k, &A(1,1), lda, &IPIV(1), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsytf2(uplo, &k, &A(1,1), lda, &IPIV(1), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsytf2_(uplo, &k, &A(1,1), lda, &IPIV(1), &iinfo);
#endif

	    kb = k;
	}

/*        Set INFO on the first occurrence of a zero pivot */

	if (*info == 0 && iinfo > 0) {
	    *info = iinfo;
	}

/*        Decrease K and return to the start of the main loop */

	k -= kb;
	goto L10;

    } else {

/*        Factorize A as L*D*L' using the lower triangle of A   

          K is the main loop index, increasing from 1 to N in steps of
   
          KB, where KB is the number of columns factorized by DLASYF; 
  
          KB is either NB or NB-1, or N-K+1 for the last block */

	k = 1;
L20:

/*        If K > N, exit from loop */

	if (k > *n) {
	    goto L40;
	}

	if (k <= *n - nb) {

/*           Factorize columns k:k+kb-1 of A and use blocked code 
to   
             update columns k+kb:n */

	    i__1 = *n - k + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlasyf_(uplo, &i__1, &nb, &kb, &A(k,k), lda, &IPIV(k), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlasyf(uplo, &i__1, &nb, &kb, &A(k,k), lda, &IPIV(k), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlasyf_(uplo, &i__1, &nb, &kb, &A(k,k), lda, &IPIV(k), 
#endif

		    &WORK(1), &ldwork, &iinfo);
	} else {

/*           Use unblocked code to factorize columns k:n of A */

	    i__1 = *n - k + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dsytf2_(uplo, &i__1, &A(k,k), lda, &IPIV(k), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsytf2(uplo, &i__1, &A(k,k), lda, &IPIV(k), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsytf2_(uplo, &i__1, &A(k,k), lda, &IPIV(k), &iinfo);
#endif

	    kb = *n - k + 1;
	}

/*        Set INFO on the first occurrence of a zero pivot */

	if (*info == 0 && iinfo > 0) {
	    *info = iinfo + k - 1;
	}

/*        Adjust IPIV */

	i__1 = k + kb - 1;
	for (j = k; j <= k+kb-1; ++j) {
	    if (IPIV(j) > 0) {
		IPIV(j) = IPIV(j) + k - 1;
	    } else {
		IPIV(j) = IPIV(j) - k + 1;
	    }
/* L30: */
	}

/*        Increase K and return to the start of the main loop */

	k += kb;
	goto L20;

    }

L40:
    WORK(1) = (LONG DOUBLE) iws;
    return;

/*     End of DSYTRF */

} /* dsytrf_ */

