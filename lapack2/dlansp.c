#define MIN(a,b)      ( ((a)<(b)) ? (a) : (b) )
#define MAX(a,b)      ( ((a)<(b)) ? (b) : (a) )
#define ABS(a)        ( ((a)<0.0) ? -(a) : (a))
#include <math.h>


#ifdef PETSC_PREFIX_SUFFIX
LONG DOUBLE dlansp_(char *norm, char *uplo, int *n, LONG DOUBLE *ap, 
#endif
#ifdef Q_C_PREFIX_SUFFIX
LONG DOUBLE qlansp(char *norm, char *uplo, int *n, LONG DOUBLE *ap, 
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
LONG DOUBLE qlansp_(char *norm, char *uplo, int *n, LONG DOUBLE *ap, 
#endif

	LONG DOUBLE *work)
{
/*  -- LAPACK auxiliary routine (version 2.0) --   
       Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,   
       Courant Institute, Argonne National Lab, and Rice University   
       October 31, 1992   


    Purpose   
    =======   

    DLANSP  returns the value of the one norm,  or the Frobenius norm, or 
  
    the  infinity norm,  or the  element of  largest absolute value  of a 
  
    real symmetric matrix A,  supplied in packed form.   

    Description   
    ===========   

    DLANSP returns the value   

       DLANSP = ( MAX(ABS(A(i,j))), NORM = 'M' or 'm'   
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
            Specifies the value to be returned in DLANSP as described   
            above.   

    UPLO    (input) CHARACTER*1   
            Specifies whether the upper or lower triangular part of the   
            symmetric matrix A is supplied.   
            = 'U':  Upper triangular part of A is supplied   
            = 'L':  Lower triangular part of A is supplied   

    N       (input) INTEGER   
            The order of the matrix A.  N >= 0.  When N = 0, DLANSP is   
            set to zero.   

    AP      (input) LONG DOUBLE PRECISION array, dimension (N*(N+1)/2)   
            The upper or lower triangle of the symmetric matrix A, packed 
  
            columnwise in a linear array.  The j-th column of A is stored 
  
            in the array AP as follows:   
            if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;   
            if UPLO = 'L', AP(i + (j-1)*(2n-j)/2) = A(i,j) for j<=i<=n.   

    WORK    (workspace) LONG DOUBLE PRECISION array, dimension (LWORK),   
            where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,   
            WORK is not referenced.   

   ===================================================================== 
  


    
   Parameter adjustments   
       Function Body */
    /* Table of constant values */
    static int c__1 = 1;
    
    /* System generated locals */
    int i__1, i__2;
    LONG DOUBLE ret_val, d__1, d__2, d__3;
    /* Builtin functions */
    /* Local variables */
    static LONG DOUBLE absa;
    static int i, j, k;
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
#define AP(I) ap[(I)-1]


    if (*n == 0) {
	value = 0.;
    } else if (lsame_(norm, "M")) {

/*        Find MAX(ABS(A(i,j))). */

	value = 0.;
	if (lsame_(uplo, "U")) {
	    k = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = k + j - 1;
		for (i = k; i <= k+j-1; ++i) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = AP(i), ABS(d__1));
		    value = MAX(d__2,d__3);
/* L10: */
		}
		k += j;
/* L20: */
	    }
	} else {
	    k = 1;
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		i__2 = k + *n - j;
		for (i = k; i <= k+*n-j; ++i) {
/* Computing MAX */
		    d__2 = value, d__3 = (d__1 = AP(i), ABS(d__1));
		    value = MAX(d__2,d__3);
/* L30: */
		}
		k = k + *n - j + 1;
/* L40: */
	    }
	}
    } else if (lsame_(norm, "I") || lsame_(norm, "O") || *(
	    unsigned char *)norm == '1') {

/*        Find normI(A) ( = norm1(A), since A is symmetric). */

	value = 0.;
	k = 1;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 1; j <= *n; ++j) {
		sum = 0.;
		i__2 = j - 1;
		for (i = 1; i <= j-1; ++i) {
		    absa = (d__1 = AP(k), ABS(d__1));
		    sum += absa;
		    WORK(i) += absa;
		    ++k;
/* L50: */
		}
		WORK(j) = sum + (d__1 = AP(k), ABS(d__1));
		++k;
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
		sum = WORK(j) + (d__1 = AP(k), ABS(d__1));
		++k;
		i__2 = *n;
		for (i = j + 1; i <= *n; ++i) {
		    absa = (d__1 = AP(k), ABS(d__1));
		    sum += absa;
		    WORK(i) += absa;
		    ++k;
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
	k = 2;
	if (lsame_(uplo, "U")) {
	    i__1 = *n;
	    for (j = 2; j <= *n; ++j) {
		i__2 = j - 1;

#ifdef PETSC_PREFIX_SUFFIX
		dlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlassq(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif

		k += j;
/* L110: */
	    }
	} else {
	    i__1 = *n - 1;
	    for (j = 1; j <= *n-1; ++j) {
		i__2 = *n - j;

#ifdef PETSC_PREFIX_SUFFIX
		dlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_C_PREFIX_SUFFIX
		qlassq(&i__2, &AP(k), &c__1, &scale, &sum);
#endif
#ifdef Q_NORMAL_PREFIX_SUFFIX
		qlassq_(&i__2, &AP(k), &c__1, &scale, &sum);
#endif

		k = k + *n - j + 1;
/* L120: */
	    }
	}
	sum *= 2;
	k = 1;
	i__1 = *n;
	for (i = 1; i <= *n; ++i) {
	    if (AP(k) != 0.) {
		absa = (d__1 = AP(k), ABS(d__1));
		if (scale < absa) {
/* Computing 2nd power */
		    d__1 = scale / absa;
		    sum = sum * (d__1 * d__1) + 1.;
		    scale = absa;
		} else {
/* Computing 2nd power */
		    d__1 = absa / scale;
		    sum += d__1 * d__1;
		}
	    }
	    if (lsame_(uplo, "U")) {
		k = k + i + 1;
	    } else {
		k = k + *n - i + 1;
	    }
/* L130: */
	}
	value = scale * sqrt(sum);
    }

    ret_val = value;
    return ret_val;

/*     End of DLANSP */

} /* dlansp_ */

