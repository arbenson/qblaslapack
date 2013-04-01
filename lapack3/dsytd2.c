#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dsytd2_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qsytd2(char *uplo, int *n, LONG DOUBLE *a, int *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qsytd2_(char *uplo, int *n, LONG DOUBLE *a, int *
#endif

	lda, LONG DOUBLE *d, LONG DOUBLE *e, LONG DOUBLE *tau, int *info)
{
/*  -- LAPACK routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DSYTD2 reduces a real symmetric matrix A to symmetric tridiagonal   
    form T by an orthogonal similarity transformation: Q' * A * Q = T.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U':  Upper triangular   
            = 'L':  Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower 
  
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

    INFO    (output) INTEGER   
            = 0:  successful exit   
            < 0:  if INFO = -i, the i-th argument had an illegal value.   

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
    static LONG DOUBLE c_b8 = 0.;
    static LONG DOUBLE c_b14 = -1.;
    
    /* System generated locals */
    int  i__1, i__2, i__3;
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
    static LONG DOUBLE taui;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsyr2_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr2(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsyr2_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 
	    int *);
    static int i;
    static LONG DOUBLE alpha;
    extern long int lsame_(char *, char *);

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void daxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qaxpy_(int *, LONG DOUBLE *, LONG DOUBLE *, 
#endif

	    int *, LONG DOUBLE *, int *);
    static long int upper;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dsymv_(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qsymv(char *, int *, LONG DOUBLE *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qsymv_(char *, int *, LONG DOUBLE *, 
#endif

	    LONG DOUBLE *, int *, LONG DOUBLE *, int *, LONG DOUBLE *, 

#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), dlarfg_(int *, LONG DOUBLE *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlarfg(int *, LONG DOUBLE *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *), qlarfg_(int *, LONG DOUBLE *,
#endif

	     LONG DOUBLE *, int *, LONG DOUBLE *), xerbla_(char *, int *
	    );



#define D(I) d[(I)-1]
#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    *info = 0;
    upper = lsame_(uplo, "U");
    if (! upper && ! lsame_(uplo, "L")) {
	*info = -1;
    } else if (*n < 0) {
	*info = -2;
    } else if (*lda < MAX(1,*n)) {
	*info = -4;
    }
    if (*info != 0) {
	i__1 = -(*info);
	xerbla_("DSYTD2", &i__1);
	return;
    }

/*     Quick return if possible */

    if (*n <= 0) {
	return;
    }

    if (upper) {

/*        Reduce the upper triangle of A */

	for (i = *n - 1; i >= 1; --i) {

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(1:i-1,i+1) */


#ifdef PETSC_PREFIX_SUFFIX
	    dlarfg_(&i, &A(i,i+1), &A(1,i+1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfg(&i, &A(i,i+1), &A(1,i+1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfg_(&i, &A(i,i+1), &A(1,i+1), &
#endif

		    c__1, &taui);
	    E(i) = A(i,i+1);

	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(1:i,1:i) */

		A(i,i+1) = 1.;

/*              Compute  x := tau * A * v  storing x in TAU(1:
i) */


#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i, &taui, &A(1,1), lda, &A(1,i+1), &c__1, &c_b8, &TAU(1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i, &taui, &A(1,1), lda, &A(1,i+1), &c__1, &c_b8, &TAU(1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i, &taui, &A(1,1), lda, &A(1,i+1), &c__1, &c_b8, &TAU(1), &c__1);
#endif


/*              Compute  w := x - 1/2 * tau * (x'*v) * v */


#ifdef PETSC_PREFIX_SUFFIX
		alpha = taui * -.5 * ddot_(&i, &TAU(1), &c__1, &A(1,i+1), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		alpha = taui * -.5 * qdot(&i, &TAU(1), &c__1, &A(1,i+1), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		alpha = taui * -.5 * qdot_(&i, &TAU(1), &c__1, &A(1,i+1), &c__1);
#endif


#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i, &alpha, &A(1,i+1), &c__1, &TAU(1), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i, &alpha, &A(1,i+1), &c__1, &TAU(1), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i, &alpha, &A(1,i+1), &c__1, &TAU(1), &
#endif

			c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */


#ifdef PETSC_PREFIX_SUFFIX
		dsyr2_(uplo, &i, &c_b14, &A(1,i+1), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsyr2(uplo, &i, &c_b14, &A(1,i+1), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsyr2_(uplo, &i, &c_b14, &A(1,i+1), &c__1, &
#endif

			TAU(1), &c__1, &A(1,1), lda);

		A(i,i+1) = E(i);
	    }
	    D(i + 1) = A(i+1,i+1);
	    TAU(i) = taui;
/* L10: */
	}
	D(1) = A(1,1);
    } else {

/*        Reduce the lower triangle of A */

	i__1 = *n - 1;
	for (i = 1; i <= *n-1; ++i) {

/*           Generate elementary reflector H(i) = I - tau * v * v'
   
             to annihilate A(i+2:n,i) */

	    i__2 = *n - i;
/* Computing MIN */
	    i__3 = i + 2;

#ifdef PETSC_PREFIX_SUFFIX
	    dlarfg_(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), &c__1, &taui);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlarfg(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), &c__1, &taui);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlarfg_(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), &c__1, &taui);
#endif

	    E(i) = A(i+1,i);

	    if (taui != 0.) {

/*              Apply H(i) from both sides to A(i+1:n,i+1:n) 
*/

		A(i+1,i) = 1.;

/*              Compute  x := tau * A * v  storing y in TAU(i:
n-1) */

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_(uplo, &i__2, &taui, &A(i+1,i+1), lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv(uplo, &i__2, &taui, &A(i+1,i+1), lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_(uplo, &i__2, &taui, &A(i+1,i+1), lda, 
#endif

			&A(i+1,i), &c__1, &c_b8, &TAU(i), &c__1);

/*              Compute  w := x - 1/2 * tau * (x'*v) * v */

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		alpha = taui * -.5 * ddot_(&i__2, &TAU(i), &c__1, &A(i+1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		alpha = taui * -.5 * qdot(&i__2, &TAU(i), &c__1, &A(i+1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		alpha = taui * -.5 * qdot_(&i__2, &TAU(i), &c__1, &A(i+1,i), &c__1);
#endif

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &TAU(i), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &alpha, &A(i+1,i), &c__1, &TAU(i), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &TAU(i), 
#endif

			&c__1);

/*              Apply the transformation as a rank-2 update: 
  
                   A := A - v * w' - w * v' */

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dsyr2_(uplo, &i__2, &c_b14, &A(i+1,i), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsyr2(uplo, &i__2, &c_b14, &A(i+1,i), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsyr2_(uplo, &i__2, &c_b14, &A(i+1,i), &c__1, &
#endif

			TAU(i), &c__1, &A(i+1,i+1), lda);

		A(i+1,i) = E(i);
	    }
	    D(i) = A(i,i);
	    TAU(i) = taui;
/* L20: */
	}
	D(*n) = A(*n,*n);
    }

    return;

/*     End of DSYTD2 */

} /* dsytd2_ */

