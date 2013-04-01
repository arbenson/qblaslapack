#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsytrd_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsytrd(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsytrd_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *d, LONG DOUBLE *e, LONG DOUBLE *tau, LONG DOUBLE *work, 
	int *lwork, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       September 30, 1994   


    Purpose   
    =======   

    DSYTRD reduces a real symmetric matrix A to real symmetric   
    tridiagonal form T by an orthogonal similarity transformation:   
    Q**T * A * Q = T.   

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
            On exit, if UPLO = 'U', the diagonal and first superdiagonal 
  
            of A are overwritten by the corresponding elements of the   
            tridiagonal matrix T, and the elements above the first   
            superdiagonal, with the array TAU, represent the orthogonal   
            matrix Q as a product of elementary reflectors; if UPLO   
            = 'L', the diagonal and first subdiagonal of A are over-   
            written by the corresponding elements of the tridiagonal   
            matrix T, and the elements below the first subdiagonal, with 
  
            the array TAU, represent the orthogonal matrix Q as a product 
  
            of elementary reflectors. See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(1,N).   

    D       (output) LONG DOUBLE PRECISION array, dimension (N)   
            The diagonal elements of the tridiagonal matrix T:   
            D(i) = A(i,i).   

    E       (output) LONG DOUBLE PRECISION array, dimension (N-1)   
            The off-diagonal elements of the tridiagonal matrix T:   
            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'. 
  

    TAU     (output) LONG DOUBLE PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors (see Further 
  
            Details).   

    WORK    (workspace/output) LONG DOUBLE PRECISION array, dimension (LWORK) 
  
            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.   

    LWORK   (input) INTEGER   
            The dimension of the array WORK.  LWORK >= 1.   
            For optimum performance LWORK >= N*NB, where NB is the   
            optimal blocksize.   

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n-1) . . . H(2) H(1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in   
    A(1:i-1,i+1), and tau in TAU(i).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(n-1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i), 
  
    and tau in TAU(i).   

    The contents of A on exit are illustrated by the following examples   
    with n = 5:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  d   e   v2  v3  v4 )              (  d                  )   
      (      d   e   v3  v4 )              (  e   d              )   
      (          d   e   v4 )              (  v1  e   d          )   
      (              d   e  )              (  v1  v2  e   d      )   
      (                  d  )              (  v1  v2  v3  e   d  )   

    where d and e denote diagonal and off-diagonal elements of T, and vi 
  
    denotes an element of the vector defining H(i).   

    ===================================================================== 
  


       Test the input parameters   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    static int c_n1 = -1;
    static int c__3 = 3;
    static int c__2 = 2;
    static LONG DOUBLE c_b22 = -1.;
    static LONG DOUBLE c_b23 = 1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int i, j;
    extern long int lsame_(char *, char *);
    static int nbmin, iinfo;
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsytd2_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsytd2(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsytd2_(char *, int *, LONG DOUBLE *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), dsyr2k_(char *, char *, int *, int *, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qsyr2k(char *, char *, int *, int *, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *, int *), qsyr2k_(char *, char *, int *, int *, LONG DOUBLE 
#endif

	    *, LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *,
	     LONG DOUBLE *, int *);
    static int nb, kk, nx;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlatrd_(char *, int *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatrd(char *, int *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlatrd_(char *, int *, int *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, LONG DOUBLE *,
	     int *), xerbla_(char *, int *);
    extern int ilaenv_(int *, char *, char *, int *, int *, 
	    int *, int *, long int, long int);
    static int ldwork, iws;



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]
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
	*info = -9;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYTRD", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n == 0) {
	WORK(1) = 1.;
	return;
    }

/*     Determine the block size. */

    nb = ilaenv_(&c__1, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1, 6L, 1L);
    nx = *n;
    iws = 1;
    if (nb > 1 && nb < *n) {

/*        Determine when to cross over from blocked to unblocked code 
  
          (last block is always handled by unblocked code).   

   Computing MAX */
	i__1 = nb, i__2 = ilaenv_(&c__3, "DSYTRD", uplo, n, &c_n1, &c_n1, &
		c_n1, 6L, 1L);
	nx = MAX(i__1,i__2);
	if (nx < *n) {

/*           Determine if workspace is large enough for blocked co
de. */

	    ldwork = *n;
	    iws = ldwork * nb;
	    if (*lwork < iws) {

/*              Not enough workspace to use optimal NB:  deter
mine the   
                minimum value of NB, and reduce NB or force us
e of   
                unblocked code by setting NX = N.   

   Computing MAX */
		i__1 = *lwork / ldwork;
		nb = MAX(i__1,1);
		nbmin = ilaenv_(&c__2, "DSYTRD", uplo, n, &c_n1, &c_n1, &c_n1,
			 6L, 1L);
		if (nb < nbmin) {
		    nx = *n;
		}
	    }
	} else {
	    nx = *n;
	}
    } else {
	nb = 1;
    }

    if (upper) {

/*        Reduce the upper triangle of A.   
          Columns 1:kk are handled by the unblocked method. */

	kk = *n - (*n - nx + nb - 1) / nb * nb;
	i__1 = kk + 1;
	i__2 = -nb;
	for (i = *n - nb + 1; -nb < 0 ? i >= kk+1 : i <= kk+1; i += -nb) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form 
the   
             matrix W which is needed to update the unreduced part
 of   
             the matrix */

	    i__3 = i + nb - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlatrd_(uplo, &i__3, &nb, &A(1,1), lda, &E(1), &TAU(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatrd(uplo, &i__3, &nb, &A(1,1), lda, &E(1), &TAU(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatrd_(uplo, &i__3, &nb, &A(1,1), lda, &E(1), &TAU(1), &
#endif

		    WORK(1), &ldwork);

/*           Update the unreduced submatrix A(1:i-1,1:i-1), using 
an   
             update of the form:  A := A - V*W' - W*V' */

	    i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &A(1,i), lda, &WORK(1), &ldwork, &c_b23, &A(1,1), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsyr2k(uplo, "No transpose", &i__3, &nb, &c_b22, &A(1,i), lda, &WORK(1), &ldwork, &c_b23, &A(1,1), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &A(1,i), lda, &WORK(1), &ldwork, &c_b23, &A(1,1), lda);
#endif


/*           Copy superdiagonal elements back into A, and diagonal
   
             elements into D */

	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		A(j-1,j) = E(j - 1);
		D(j) = A(j,j);
/* L10: */
	    }
/* L20: */
	}

/*        Use unblocked code to reduce the last or only block */


#ifdef PETSC_PREFIX_SUFFIX
	dsytd2_(uplo, &kk, &A(1,1), lda, &D(1), &E(1), &TAU(1), &iinfo);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsytd2(uplo, &kk, &A(1,1), lda, &D(1), &E(1), &TAU(1), &iinfo);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsytd2_(uplo, &kk, &A(1,1), lda, &D(1), &E(1), &TAU(1), &iinfo);
#endif

    } else {

/*        Reduce the lower triangle of A */

	i__2 = *n - nx;
	i__1 = nb;
	for (i = 1; nb < 0 ? i >= *n-nx : i <= *n-nx; i += nb) {

/*           Reduce columns i:i+nb-1 to tridiagonal form and form 
the   
             matrix W which is needed to update the unreduced part
 of   
             the matrix */

	    i__3 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dlatrd_(uplo, &i__3, &nb, &A(i,i), lda, &E(i), &TAU(i),
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlatrd(uplo, &i__3, &nb, &A(i,i), lda, &E(i), &TAU(i),
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlatrd_(uplo, &i__3, &nb, &A(i,i), lda, &E(i), &TAU(i),
#endif

		     &WORK(1), &ldwork);

/*           Update the unreduced submatrix A(i+ib:n,i+ib:n), usin
g   
             an update of the form:  A := A - V*W' - W*V' */

	    i__3 = *n - i - nb + 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &A(i+nb,i), lda, &WORK(nb + 1), &ldwork, &c_b23, &A(i+nb,i+nb), lda);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsyr2k(uplo, "No transpose", &i__3, &nb, &c_b22, &A(i+nb,i), lda, &WORK(nb + 1), &ldwork, &c_b23, &A(i+nb,i+nb), lda);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsyr2k_(uplo, "No transpose", &i__3, &nb, &c_b22, &A(i+nb,i), lda, &WORK(nb + 1), &ldwork, &c_b23, &A(i+nb,i+nb), lda);
#endif


/*           Copy subdiagonal elements back into A, and diagonal 
  
             elements into D */

	    i__3 = i + nb - 1;
	    for (j = i; j <= i+nb-1; ++j) {
		A(j+1,j) = E(j);
		D(j) = A(j,j);
/* L30: */
	    }
/* L40: */
	}

/*        Use unblocked code to reduce the last or only block */

	i__1 = *n - i + 1;

#ifdef PETSC_PREFIX_SUFFIX
	dsytd2_(uplo, &i__1, &A(i,i), lda, &D(i), &E(i), &TAU(i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qsytd2(uplo, &i__1, &A(i,i), lda, &D(i), &E(i), &TAU(i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qsytd2_(uplo, &i__1, &A(i,i), lda, &D(i), &E(i), &TAU(i), &
#endif

		iinfo);
    }

    WORK(1) = (LONG DOUBLE) iws;
    return;

/*     End of DSYTRD */

} /* dsytrd_ */

