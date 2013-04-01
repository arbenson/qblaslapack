
#include <math.h>
#define MIN(a,b)           ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)           ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)             ( ((a)<0.0)   ? -(a) : (a) )


#ifdef PETSC_PREFIX_SUFFIX
/* Subroutine */ void dlaset_(char *uplo, int *m, int *n, LONG DOUBLE *
#endif
#ifdef Q_C_PREFIX_SUFFIX
/* Subroutine */ void qlaset(char *uplo, int *m, int *n, LONG DOUBLE *
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
/* Subroutine */ void qlaset_(char *uplo, int *m, int *n, LONG DOUBLE *
#endif

	alpha, LONG DOUBLE *beta, LONG DOUBLE *a, int *lda)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLASET initializes an m-by-n matrix A to BETA on the diagonal and   
    ALPHA on the offdiagonals.   

    Arguments   
    =========   

    UPLO    (input) CHARACTER*1   
            Specifies the part of the matrix A to be set.   
            = 'U':      Upper triangular part is set; the strictly lower 
  
                        triangular part of A is not changed.   
            = 'L':      Lower triangular part is set; the strictly upper 
  
                        triangular part of A is not changed.   
            Otherwise:  All of the matrix A is set.   

    M       (input) INT   
            The number of rows of the matrix A.  M >= 0.   

    N       (input) INT   
            The number of columns of the matrix A.  N >= 0.   

    ALPHA   (input) LONG DOUBLE PRECISION   
            The constant to which the offdiagonal elements are to be set. 
  

    BETA    (input) LONG DOUBLE PRECISION   
            The constant to which the diagonal elements are to be set.   

    A       (input/output) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            On exit, the leading m-by-n submatrix of A is set as follows: 
  

            if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,   
            if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,   
            otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,   

            and, for all UPLO, A(i,i) = BETA, 1<=i<=MIN(m,n).   

    LDA     (input) INT   
            The leading dimension of the array A.  LDA >= MAX(1,M).   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* System generated locals */
    int  i__1, i__2, i__3;
    /* Local variables */
    static int i, j;
    extern long int lsame_(char *, char *);



#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (lsame_(uplo, "U")) {

/*        Set the strictly upper triangular or trapezoidal part of the
   
          array to ALPHA. */

	i__1 = *n;
	for (j = 2; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = j - 1;
	    i__2 = MIN(i__3,*m);
	    for (i = 1; i <= MIN(j-1,*m); ++i) {
		A(i,j) = *alpha;
/* L10: */
	    }
/* L20: */
	}

    } else if (lsame_(uplo, "L")) {

/*        Set the strictly lower triangular or trapezoidal part of the
   
          array to ALPHA. */

	i__1 = MIN(*m,*n);
	for (j = 1; j <= MIN(*m,*n); ++j) {
	    i__2 = *m;
	    for (i = j + 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L30: */
	    }
/* L40: */
	}

    } else {

/*        Set the leading m-by-n submatrix to ALPHA. */

	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    i__2 = *m;
	    for (i = 1; i <= *m; ++i) {
		A(i,j) = *alpha;
/* L50: */
	    }
/* L60: */
	}
    }

/*     Set the first MIN(M,N) diagonal elements to BETA. */

    i__1 = MIN(*m,*n);
    for (i = 1; i <= MIN(*m,*n); ++i) {
	A(i,i) = *beta;
/* L70: */
    }

    return;

/*     End of DLASET */

} /* dlaset_ */

