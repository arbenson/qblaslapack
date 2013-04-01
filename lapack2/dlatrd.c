#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlatrd_(char *uplo, int *n, int *nb, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlatrd(char *uplo, int *n, int *nb, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlatrd_(char *uplo, int *n, int *nb, LONG DOUBLE *
#endif

	a, int *lda, LONG DOUBLE *e, LONG DOUBLE *tau, LONG DOUBLE *w, 
	int *ldw)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLATRD reduces NB rows and columns of a real symmetric matrix A to   
    symmetric tridiagonal form by an orthogonal similarity   
    transformation Q' * A * Q, and returns the matrices V and W which are 
  
    needed to apply the transformation to the unreduced part of A.   

    If UPLO = 'U', DLATRD reduces the last NB rows and columns of a   
    matrix, of which the upper triangle is supplied;   
    if UPLO = 'L', DLATRD reduces the first NB rows and columns of a   
    matrix, of which the lower triangle is supplied.   

    This is an auxiliary routine called by DSYTRD.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is stored:   
            = 'U': Upper triangular   
            = 'L': Lower triangular   

    N       (input) INTEGER   
            The order of the matrix A.   

    NB      (input) INTEGER   
            The number of rows and columns to be reduced.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On entry, the symmetric matrix A.  If UPLO = 'U', the leading 
  
            n-by-n upper triangular part of A contains the upper   
            triangular part of the matrix A, and the strictly lower   
            triangular part of A is not referenced.  If UPLO = 'L', the   
            leading n-by-n lower triangular part of A contains the lower 
  
            triangular part of the matrix A, and the strictly upper   
            triangular part of A is not referenced.   
            On exit:   
            if UPLO = 'U', the last NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements above the diagonal 
  
              with the array TAU, represent the orthogonal matrix Q as a 
  
              product of elementary reflectors;   
            if UPLO = 'L', the first NB columns have been reduced to   
              tridiagonal form, with the diagonal elements overwriting   
              the diagonal elements of A; the elements below the diagonal 
  
              with the array TAU, represent the  orthogonal matrix Q as a 
  
              product of elementary reflectors.   
            See Further Details.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= (1,N).   

    E       (output) LONG DOUBLE PRECISION array, dimension (N-1)   
            If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal   
            elements of the last NB columns of the reduced matrix;   
            if UPLO = 'L', E(1:nb) contains the subdiagonal elements of   
            the first NB columns of the reduced matrix.   

    TAU     (output) LONG DOUBLE PRECISION array, dimension (N-1)   
            The scalar factors of the elementary reflectors, stored in   
            TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'. 
  
            See Further Details.   

    W       (output) LONG DOUBLE PRECISION array, dimension (LDW,NB)   
            The n-by-nb matrix W required to update the unreduced part   
            of A.   

    LDW     (input) INTEGER   
            The leading dimension of the array W. LDW >= MAX(1,N).   

    Further Details   
    ===============   

    If UPLO = 'U', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(n) H(n-1) . . . H(n-nb+1).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i), 
  
    and tau in TAU(i-1).   

    If UPLO = 'L', the matrix Q is represented as a product of elementary 
  
    reflectors   

       Q = H(1) H(2) . . . H(nb).   

    Each H(i) has the form   

       H(i) = I - tau * v * v'   

    where tau is a real scalar, and v is a real vector with   
    v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i), 
  
    and tau in TAU(i).   

    The elements of the vectors v together form the n-by-nb matrix V   
    which is needed, with W, to apply the transformation to the unreduced 
  
    part of the matrix, using a symmetric rank-2k update of the form:   
    A := A - V*W' - W*V'.   

    The contents of A on exit are illustrated by the following examples   
    with n = 5 and nb = 2:   

    if UPLO = 'U':                       if UPLO = 'L':   

      (  a   a   a   v4  v5 )              (  d                  )   
      (      a   a   v4  v5 )              (  1   d              )   
      (          a   1   v5 )              (  v1  1   a          )   
      (              d   1  )              (  v1  v2  a   a      )   
      (                  d  )              (  v1  v2  a   a   a  )   

    where d denotes a diagonal element of the reduced matrix, a denotes   
    an element of the original matrix that is unchanged, and vi denotes   
    an element of the vector defining H(i).   

    ===================================================================== 
  


       Quick return if possible   

    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static LONG DOUBLE c_b5 = -1.;
    static LONG DOUBLE c_b6 = 1.;
    static int c__1 = 1;
    static LONG DOUBLE c_b16 = 0.;
    
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
    static int i;
    static LONG DOUBLE alpha;

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
	    LONG DOUBLE *, LONG DOUBLE *, int *), daxpy_(int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qaxpy(int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, LONG DOUBLE *, int *), qaxpy_(int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *, int *, LONG DOUBLE *, int *), 

#ifdef PETSC_PREFIX_SUFFIX
	    dsymv_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qsymv(char *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qsymv_(char *, int *, LONG DOUBLE *, LONG DOUBLE *, int *, 
#endif


#ifdef PETSC_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), dlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *, int *,
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), qlarfg(int *, LONG DOUBLE *, LONG DOUBLE *, int *,
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    LONG DOUBLE *, int *, LONG DOUBLE *, LONG DOUBLE *, int *), qlarfg_(int *, LONG DOUBLE *, LONG DOUBLE *, int *,
#endif

	     LONG DOUBLE *);
    static int iw;



#define E(I) e[(I)-1]
#define TAU(I) tau[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]
#define W(I,J) w[(I)-1 + ((J)-1)* ( *ldw)]

    if (*n <= 0) {
	return;
    }

    if (lsame_(uplo, "U")) {

/*        Reduce last NB columns of upper triangle */

	i__1 = *n - *nb + 1;
	for (i = *n; i >= *n-*nb+1; --i) {
	    iw = i - *n + *nb;
	    if (i < *n) {

/*              Update A(1:i,i) */

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("No transpose", &i, &i__2, &c_b5, &A(1,i+1), lda, &W(i,iw+1), ldw, &c_b6, &A(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("No transpose", &i, &i__2, &c_b5, &A(1,i+1), lda, &W(i,iw+1), ldw, &c_b6, &A(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("No transpose", &i, &i__2, &c_b5, &A(1,i+1), lda, &W(i,iw+1), ldw, &c_b6, &A(1,i), &c__1);
#endif

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("No transpose", &i, &i__2, &c_b5, &W(1,iw+1), ldw, &A(i,i+1), lda, &c_b6, &A(1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("No transpose", &i, &i__2, &c_b5, &W(1,iw+1), ldw, &A(i,i+1), lda, &c_b6, &A(1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("No transpose", &i, &i__2, &c_b5, &W(1,iw+1), ldw, &A(i,i+1), lda, &c_b6, &A(1,i), &c__1);
#endif

	    }
	    if (i > 1) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(1:i-2,i) */

		i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlarfg_(&i__2, &A(i-1,i), &A(1,i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfg(&i__2, &A(i-1,i), &A(1,i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfg_(&i__2, &A(i-1,i), &A(1,i), &
#endif

			c__1, &TAU(i - 1));
		E(i - 1) = A(i-1,i);
		A(i-1,i) = 1.;

/*              Compute W(1:i-1,i) */

		i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_("Upper", &i__2, &c_b6, &A(1,1), lda, &A(1,i), &c__1, &c_b16, &W(1,iw), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv("Upper", &i__2, &c_b6, &A(1,1), lda, &A(1,i), &c__1, &c_b16, &W(1,iw), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_("Upper", &i__2, &c_b6, &A(1,1), lda, &A(1,i), &c__1, &c_b16, &W(1,iw), &
#endif

			c__1);
		if (i < *n) {
		    i__2 = i - 1;
		    i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemv_("Transpose", &i__2, &i__3, &c_b6, &W(1,iw+1), ldw, &A(1,i), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemv("Transpose", &i__2, &i__3, &c_b6, &W(1,iw+1), ldw, &A(1,i), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemv_("Transpose", &i__2, &i__3, &c_b6, &W(1,iw+1), ldw, &A(1,i), &c__1, &
#endif

			    c_b16, &W(i+1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemv_("No transpose", &i__2, &i__3, &c_b5, &A(1,i+1), lda, &W(i+1,iw), &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemv("No transpose", &i__2, &i__3, &c_b5, &A(1,i+1), lda, &W(i+1,iw), &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemv_("No transpose", &i__2, &i__3, &c_b5, &A(1,i+1), lda, &W(i+1,iw), &c__1, 
#endif

			    &c_b6, &W(1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemv_("Transpose", &i__2, &i__3, &c_b6, &A(1,i+1), lda, &A(1,i), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemv("Transpose", &i__2, &i__3, &c_b6, &A(1,i+1), lda, &A(1,i), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemv_("Transpose", &i__2, &i__3, &c_b6, &A(1,i+1), lda, &A(1,i), &c__1, &
#endif

			    c_b16, &W(i+1,iw), &c__1);
		    i__2 = i - 1;
		    i__3 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		    dgemv_("No transpose", &i__2, &i__3, &c_b5, &W(1,iw+1), ldw, &W(i+1,iw), &c__1, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qgemv("No transpose", &i__2, &i__3, &c_b5, &W(1,iw+1), ldw, &W(i+1,iw), &c__1, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qgemv_("No transpose", &i__2, &i__3, &c_b5, &W(1,iw+1), ldw, &W(i+1,iw), &c__1, 
#endif

			    &c_b6, &W(1,iw), &c__1);
		}
		i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &TAU(i - 1), &W(1,iw), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &TAU(i - 1), &W(1,iw), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &TAU(i - 1), &W(1,iw), &c__1);
#endif

		i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		alpha = TAU(i - 1) * -.5 * ddot_(&i__2, &W(1,iw), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		alpha = TAU(i - 1) * -.5 * qdot(&i__2, &W(1,iw), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		alpha = TAU(i - 1) * -.5 * qdot_(&i__2, &W(1,iw), &
#endif

			c__1, &A(1,i), &c__1);
		i__2 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &alpha, &A(1,i), &c__1, &W(1,iw), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &alpha, &A(1,i), &c__1, &W(1,iw), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &alpha, &A(1,i), &c__1, &W(1,iw), &c__1);
#endif

	    }

/* L10: */
	}
    } else {

/*        Reduce first NB columns of lower triangle */

	i__1 = *nb;
	for (i = 1; i <= *nb; ++i) {

/*           Update A(i:n,i) */

	    i__2 = *n - i + 1;
	    i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("No transpose", &i__2, &i__3, &c_b5, &A(i,1), lda, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("No transpose", &i__2, &i__3, &c_b5, &A(i,1), lda, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("No transpose", &i__2, &i__3, &c_b5, &A(i,1), lda, &
#endif

		    W(i,1), ldw, &c_b6, &A(i,i), &c__1)
		    ;
	    i__2 = *n - i + 1;
	    i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
	    dgemv_("No transpose", &i__2, &i__3, &c_b5, &W(i,1), ldw, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qgemv("No transpose", &i__2, &i__3, &c_b5, &W(i,1), ldw, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qgemv_("No transpose", &i__2, &i__3, &c_b5, &W(i,1), ldw, &
#endif

		    A(i,1), lda, &c_b6, &A(i,i), &c__1)
		    ;
	    if (i < *n) {

/*              Generate elementary reflector H(i) to annihila
te   
                A(i+2:n,i) */

		i__2 = *n - i;
/* Computing MIN */
		i__3 = i + 2;

#ifdef PETSC_PREFIX_SUFFIX
		dlarfg_(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), &c__1, &TAU(i));
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlarfg(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), &c__1, &TAU(i));
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlarfg_(&i__2, &A(i+1,i), &A(MIN(i+2,*n),i), &c__1, &TAU(i));
#endif

		E(i) = A(i+1,i);
		A(i+1,i) = 1.;

/*              Compute W(i+1:n,i) */

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dsymv_("Lower", &i__2, &c_b6, &A(i+1,i+1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qsymv("Lower", &i__2, &c_b6, &A(i+1,i+1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qsymv_("Lower", &i__2, &c_b6, &A(i+1,i+1), 
#endif

			lda, &A(i+1,i), &c__1, &c_b16, &W(i+1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("Transpose", &i__2, &i__3, &c_b6, &W(i+1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("Transpose", &i__2, &i__3, &c_b6, &W(i+1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("Transpose", &i__2, &i__3, &c_b6, &W(i+1,1), 
#endif

			ldw, &A(i+1,i), &c__1, &c_b16, &W(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &A(i+1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("No transpose", &i__2, &i__3, &c_b5, &A(i+1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("No transpose", &i__2, &i__3, &c_b5, &A(i+1,1)
#endif

			, lda, &W(1,i), &c__1, &c_b6, &W(i+1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("Transpose", &i__2, &i__3, &c_b6, &A(i+1,1), 
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("Transpose", &i__2, &i__3, &c_b6, &A(i+1,1), 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("Transpose", &i__2, &i__3, &c_b6, &A(i+1,1), 
#endif

			lda, &A(i+1,i), &c__1, &c_b16, &W(1,i), &c__1);
		i__2 = *n - i;
		i__3 = i - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dgemv_("No transpose", &i__2, &i__3, &c_b5, &W(i+1,1)
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qgemv("No transpose", &i__2, &i__3, &c_b5, &W(i+1,1)
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qgemv_("No transpose", &i__2, &i__3, &c_b5, &W(i+1,1)
#endif

			, ldw, &W(1,i), &c__1, &c_b6, &W(i+1,i), &c__1);
		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		dscal_(&i__2, &TAU(i), &W(i+1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qscal(&i__2, &TAU(i), &W(i+1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qscal_(&i__2, &TAU(i), &W(i+1,i), &c__1);
#endif

		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		alpha = TAU(i) * -.5 * ddot_(&i__2, &W(i+1,i), &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		alpha = TAU(i) * -.5 * qdot(&i__2, &W(i+1,i), &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		alpha = TAU(i) * -.5 * qdot_(&i__2, &W(i+1,i), &
#endif

			c__1, &A(i+1,i), &c__1);
		i__2 = *n - i;

#ifdef PETSC_PREFIX_SUFFIX
		daxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &W(i+1,i), &c__1);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qaxpy(&i__2, &alpha, &A(i+1,i), &c__1, &W(i+1,i), &c__1);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qaxpy_(&i__2, &alpha, &A(i+1,i), &c__1, &W(i+1,i), &c__1);
#endif

	    }

/* L20: */
	}
    }

    return;

/*     End of DLATRD */

} /* dlatrd_ */

