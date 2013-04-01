#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlanhs_(char *norm, int *n, LONG DOUBLE *a, int *lda, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlanhs(char *norm, int *n, LONG DOUBLE *a, int *lda, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlanhs_(char *norm, int *n, LONG DOUBLE *a, int *lda, 
#endif

	LONG DOUBLE *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANHS  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    Hessenberg matrix A.   

    Description   
    ===========   

    DLANHS returns the value   

       DLANHS = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'   
                (   
                ( norm1(A),         NORM = '1', 'O' or 'o'   
                (   
                ( normI(A),         NORM = 'I' or 'i'   
                (   
                ( normF(A),         NORM = 'F', 'f', 'E' or 'e'   

    where  norm1  denotes the  one norm of a matrix (maximum column sum), 
  
    normI  denotes the  infinity norm  of a matrix  (maximum row sum) and 
  
    normF  denotes the  Frobenius norm of a matrix (square root of sum of 
  
    squares).  Note that  MAX(ABS(A(i,j)))  is not a  matrix norm.   

    Arguments   
    =========   

    NORM    (input) CHARACTER*1   
            Specifies the value to be returned in DLANHS as described   
            above.   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANHS is   
            set to zero.   

    A       (input) LONG DOUBLE PRECISION array, dimension (LDA,N)   
            The n by n upper Hessenberg matrix A; the part of A below the 
  
            first sub-diagonal is not referenced.   

    LDA     (input) INTEGER   
            The leading dimension of the array A.  LDA >= MAX(N,1).   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I'; otherwise, WORK is not   
            referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int  i__1, i__2, i__3, i__4;
    LONG DOUBLE ret_val, d__1, d__2, d__3;
    /* Builtin functions */

    /* Local variables */
    static int i, j;
    static LONG DOUBLE scale;
    extern long int lsame_(char *, char *);
    static LONG DOUBLE value;

#ifdef PETSC_PREFIX_SUFFIX
    extern /* Subroutine */ void dlassq_(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
    extern /* Subroutine */ void qlassq(int *, LONG DOUBLE *, int *, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
    extern /* Subroutine */ void qlassq_(int *, LONG DOUBLE *, int *, 
#endif

	    LONG DOUBLE *, LONG DOUBLE *);
    static LONG DOUBLE sum;



#define WORK(I) work[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find MAX(ABS(A(i,j))). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = MIN(i__3,i__4);
	    for (i = 1; i <= MIN(*n,j+1); ++i) {
/* Computing MAX */
		d__2 = value, d__3 = (d__1 = A(i,j), ABS(d__1));
		value = MAX(d__2,d__3);
/* L10: */
	    }
/* L20: */
	}
    } else if (lsame_(norm, "O") || *(unsigned char *)norm == '1') {

/*        Find norm1(A). */

	value = 0.;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
	    sum = 0.;
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = MIN(i__3,i__4);
	    for (i = 1; i <= MIN(*n,j+1); ++i) {
		sum += (d__1 = A(i,j), ABS(d__1));
/* L30: */
	    }
	    value = MAX(value,sum);
/* L40: */
	}
    } else if (lsame_(norm, "I")) {

/*        Find normI(A). */

	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    WORK(i) = 0.;
/* L50: */
	}
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = MIN(i__3,i__4);
	    for (i = 1; i <= MIN(*n,j+1); ++i) {
		WORK(i) += (d__1 = A(i,j), ABS(d__1));
/* L60: */
	    }
/* L70: */
	}
	value = 0.;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
/* Computing MAX */
	    d__1 = value, d__2 = WORK(i);
	    value = MAX(d__1,d__2);
/* L80: */
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	i__1 = *n;
	for (j = 1; j <= *n; ++j) {
/* Computing MIN */
	    i__3 = *n, i__4 = j + 1;
	    i__2 = MIN(i__3,i__4);

#ifdef PETSC_PREFIX_SUFFIX
	    dlassq_(&i__2, &A(1,j), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	    qlassq(&i__2, &A(1,j), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	    qlassq_(&i__2, &A(1,j), &c__1, &scale, &sum);
#endif

/* L90: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANHS */

} /* dlanhs_ */

