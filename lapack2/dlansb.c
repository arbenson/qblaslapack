#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlansb_(char *norm, char *uplo, int *n, int *k, LONG DOUBLE 
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlansb(char *norm, char *uplo, int *n, int *k, LONG DOUBLE 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlansb_(char *norm, char *uplo, int *n, int *k, LONG DOUBLE 
#endif

	*ab, int *ldab, LONG DOUBLE *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANSB  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the element of  largest absolute value  of an 
  
    n by n symmetric band matrix A,  with k super-diagonals.   

    Description   
    ===========   

    DLANSB returns the value   

       DLANSB = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in DLANSB as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            band matrix A is supplied.   
            = 'U':  Upper triangular part is supplied   
            = 'L':  Lower triangular part is supplied   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANSB is   
            set to zero.   

    K       (input) INTEGER   
            The number of super-diagonals or sub-diagonals of the   
            band matrix A.  K >= 0.   

    AB      (input) LONG DOUBLE PRECISION array, dimension (LDAB,N)   
            The upper or lower triangle of the symmetric band matrix A,   
            stored in the first K+1 rows of AB.  The j-th column of A is 
  
            stored in the j-th column of the array AB as follows:   
            if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for MAX(1,j-k)<=i<=j;   
            if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=MIN(n,j+k).   

    LDAB    (input) INTEGER   
            The leading dimension of the array AB.  LDAB >= K+1.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,   
            WORK is not referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1, i__2, i__3, i__4;
    LONG DOUBLE ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE absa;
    static int i, j, l;
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

#define AB(I,J) ab[(I)-1 + ((J)-1)* ( *ldab)]

    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find MAX(ABS(A(i,j))). */

	value = 0.;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
/* Computing MAX */
		i__2 = *k + 2 - j;
		i__3 = *k + 1;
		for (i = MAX(*k+2-j,1); i <= *k+1; ++i) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = AB(i,j), ABS(
			    d__1));
		    value = MAX(d__2,d__3);
/* L10: */
		}
/* L20: */
	    }
	} else {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
/* Computing MIN */
		i__2 = *n + 1 - j, i__4 = *k + 1;
		i__3 = MIN(i__2,i__4);
		for (i = 1; i <= MIN(*n+1-j,*k+1); ++i) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = AB(i,j), ABS(
			    d__1));
		    value = MAX(d__2,d__3);
/* L30: */
		}
/* L40: */
	    }
	}
    } else if (lsame_(norm, "I") || lsame_(norm, "O") || *(
	    unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

	value = 0.;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		sum = 0.;
		l = *k + 1 - j;
/* Computing MAX */
		i__3 = 1, i__2 = j - *k;
		i__4 = j - 1;
		for (i = MAX(1,j-*k); i <= j-1; ++i) {
		    absa = (d__1 = AB(l+i,j), ABS(d__1));
		    sum += absa;
		    WORK(i) += absa;
/* L50: */
		}
		WORK(j) = sum + (d__1 = AB(*k+1,j), ABS(d__1));
/* L60: */
	    }
	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
/* Computing MAX */
		d__1 = value, d__2 = WORK(i);
		value = MAX(d__1,d__2);
/* L70: */
	    }
	} else {
	    i__1 = *n;
	    for (i = 1; i <= *n; ++i) {
		WORK(i) = 0.;
/* L80: */
	    }
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		sum = WORK(j) + (d__1 = AB(1,j), ABS(d__1));
		l = 1 - j;
/* Computing MIN */
		i__3 = *n, i__2 = j + *k;
		i__4 = MIN(i__3,i__2);
		for (i = j + 1; i <= MIN(*n,j+*k); ++i) {
		    absa = (d__1 = AB(l+i,j), ABS(d__1));
		    sum += absa;
		    WORK(i) += absa;
/* L90: */
		}
		value = MAX(value,sum);
/* L100: */
	    }
	}
    } else if (lsame_(norm, "F") || lsame_(norm, "E")) {

/*        Find normF(A). */

	scale = 0.;
	sum = 1.;
	if (*k > 0) {
	    if (lsame_(uplo, "U")) {
		i__1 = *n;
		for (j = 2; j <= *n; ++j) {
/* Computing MIN */
		    i__3 = j - 1;
		    i__4 = MIN(i__3,*k);
/* Computing MAX */
		    i__2 = *k + 2 - j;

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&i__4, &AB(MAX(*k+2-j,1),j), &c__1, &
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&i__4, &AB(MAX(*k+2-j,1),j), &c__1, &
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&i__4, &AB(MAX(*k+2-j,1),j), &c__1, &
#endif

			    scale, &sum);
/* L110: */
		}
		l = *k + 1;
	    } else {
		i__1 = *n - 1;
		for (j = 1; j <= *n-1; ++j) {
/* Computing MIN */
		    i__3 = *n - j;
		    i__4 = MIN(i__3,*k);

#ifdef PETSC_PREFIX_SUFFIX
		    dlassq_(&i__4, &AB(2,j), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		    qlassq(&i__4, &AB(2,j), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		    qlassq_(&i__4, &AB(2,j), &c__1, &scale, &sum);
#endif

/* L120: */
		}
		l = 1;
	    }
	    sum *= 2;
	} else {
	    l = 1;
	}

#ifdef PETSC_PREFIX_SUFFIX
	dlassq_(n, &AB(l,1), ldab, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
	qlassq(n, &AB(l,1), ldab, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
	qlassq_(n, &AB(l,1), ldab, &scale, &sum);
#endif

	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANSB */

} /* dlansb_ */

